# SPDX-License-Identifier: GPL-2.0-only
# Copyright (C) 2025 Neulite Core Team <neulite-core@numericalbrain.org>

import ast
from collections import Counter
from decimal import Decimal, ROUND_HALF_UP
import glob
import h5py
import json
import logging
import numpy as np
import os
import pandas as pd
import pickle
import re
import shutil
from sklearn.decomposition import PCA
import sys
import time

from bmtk.builder import NetworkBuilder
from bmtk.builder.network_adaptors import DenseNetwork
from bmtk.builder.network_adaptors.edges_collator import EdgesCollator
from bmtk.builder.index_builders import create_index_in_memory
from bmtk.builder.builder_utils import comm, mpi_rank, mpi_size, barrier
from bmtk.builder.edges_sorter import sort_edges
from bmtk.builder.builder_utils import add_hdf5_attrs
from bmtk.utils.sonata.config import SonataConfig

logger = logging.getLogger(__name__)

class NeuliteNetwork(DenseNetwork):
    # Connection CSV header definition
    CONNECTION_CSV_HEADER = ["#pre nid", "post nid", "post cid", "weight",
                             "tau_decay", "tau_rise", "erev", "delay", "e/i"]

    def __init__(self, name, convert_morphologies=True, convert_ion_channels=True, neulite_only=True,
                 simulation_config=None, generate_config_h=True, **network_props):
        super(NeuliteNetwork, self).__init__(name, **network_props or {})

        # Detect base_dir from the calling script
        self.my_base_dir = ""
        if len(sys.argv) > 0 and os.path.exists(sys.argv[0]):
            try:
                with open(sys.argv[0], 'r', encoding='utf-8') as f:
                    for line in f:
                        match = re.search(r'base_dir\s*=\s*[\'"]([^\'\"]+)[\'"]', line)
                        if match:
                            self.my_base_dir = match.group(1)
                            logger.debug(f'Detected base_dir from script: {self.my_base_dir}')
                            break
            except Exception as e:
                logger.warning(f"Failed to detect base_dir from script: {e}")

        self._nodes = []
        self._DenseNetwork__edges_tables = []
        self._target_networks = {}
        self.morphologies_dir = ""
        self.biophysical_neuron_models_dir = ""
        self.neulite_dir = f"{self.my_base_dir}_nl" if self.my_base_dir else "neulite"
        self.ion_dir = "data"
        self.swc_dir = "data"
        self.convert_morphologies_flag = convert_morphologies
        self.convert_ion_channels_flag = convert_ion_channels
        self.neulite_only = neulite_only

        # Generate config.h if simulation_config is provided
        if simulation_config is not None and generate_config_h:
            if mpi_rank == 0:
                logger.info("Generating config.h from simulation configuration")
                self.generate_config_h(simulation_config)

        try:
            config = SonataConfig.from_json(self._get_config_path())
            self.morphologies_dir = config["components"]["morphologies_dir"]
            self.biophysical_neuron_models_dir = config["components"]["biophysical_neuron_models_dir"]
            if mpi_rank == 0:
                if self.convert_morphologies_flag and os.path.exists(self.morphologies_dir):
                    logger.info("Converting morphology SWC files for Neulite")
                    self._convert_morphologies()
                elif self.convert_morphologies_flag:
                    logger.warning(f"Morphologies directory does not exists. [{self.morphologies_dir}]")
                if self.convert_ion_channels_flag and self.biophysical_neuron_models_dir:
                    logger.info("Converting ion channel JSON files to CSV format")
                    self.convert_ion_channels()

                # Generate config.h from config.json if simulation_config was not provided
                if simulation_config is None and generate_config_h:
                    logger.info("Generating config.h from config.json")
                    self.generate_config_h(config)
            barrier()
        except Exception as e:
            logger.warning(f"Initialization warning: {e}")
    
    def add_nodes(self, N=1, **properties):
        """Override add_nodes to optionally filter non-biophysical nodes"""
        if self.neulite_only:
            # Neulite-only mode: only accept biophysical nodes
            if properties.get('model_type') == 'biophysical':
                super(NeuliteNetwork, self).add_nodes(N=N, **properties)
            else:
                logger.info(f"[Neulite-only mode] Skipping non-biophysical nodes: model_type={properties.get('model_type')}, N={N}")
        else:
            # Normal mode: accept all nodes
            super(NeuliteNetwork, self).add_nodes(N=N, **properties)
    
    def add_edges(self, source=None, target=None, connection_rule=None, connection_params=None, 
                  iterator='one_to_one', edge_type_properties=None, **properties):
        """Add edges between source and target nodes"""
        super(NeuliteNetwork, self).add_edges(source=source, target=target, 
                                            connection_rule=connection_rule,
                                            connection_params=connection_params,
                                            iterator=iterator,
                                            edge_type_properties=edge_type_properties,
                                            **properties)

    def set_config_path(self, config_path):
        """Set the path for config.json.

        :param config_path: the path for config.json
        """
        config = SonataConfig.from_json(config_path)
        self.morphologies_dir = config["components"]["morphologies_dir"]
        self.biophysical_neuron_models_dir = config["components"]["biophysical_neuron_models_dir"]
        if mpi_rank == 0:
            if self.convert_morphologies_flag and os.path.exists(self.morphologies_dir):
                logger.info("Converting morphology SWC files for Neulite")
                self._convert_morphologies()
            elif self.convert_morphologies_flag:
                logger.error(f"Morphologies directory does not exists. [{self.morphologies_dir}]")
                raise Exception(f"Morphologies directory does not exists. [{self.morphologies_dir}]")
            if self.convert_ion_channels_flag and self.biophysical_neuron_models_dir:
                logger.info("Converting ion channel JSON files to CSV format")
                self.convert_ion_channels()
        barrier()

    def set_dir(self, ion_dir, swc_dir):
        """Set the directory path of ion and swc files for Neulite.

        :param ion_dir: the directory path of ion csv files for Neulite
        :param swc_dir: the directory path of swc files for Neulite
        """
        self.ion_dir = ion_dir
        self.swc_dir = swc_dir

    def _get_config_path(self, config_file="config.json"):
        """Get the full path to config file based on base_dir

        :param config_file: config file name (default: config.json)
        :return: Full path to config file
        """
        if self.my_base_dir:
            return f"{self.my_base_dir}/{config_file}"
        return config_file

    def _convert_morphologies(self):
        # Create output directory if it doesn't exist
        output_dir = os.path.join(self.neulite_dir, self.swc_dir)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        files = glob.glob(os.path.join(self.morphologies_dir, "*.swc"))
        for file in files:
            logger.debug(f"_convert_morphologies: {file}")
            swc = pd.read_csv(file, names=['id', 'type', 'x', 'y', 'z', 'r', 'parent'], comment="#", delimiter=" ", dtype={'id': int, 'parent': int})
            swc = self._replace_zero_based_index(swc)
            swc = self._fix_axon_perisomatic_directed_df(swc)
            swc = self._sort_by_DFS(swc)
            swc = swc.rename(columns={'id': '#id'})
            output_file = os.path.join(output_dir, os.path.basename(file))
            swc.to_csv(output_file, index=False, sep=" ", float_format='%.6f')
            logger.debug(f"Saved converted morphology to: {output_file}")
        
        logger.info(f"Morphology conversion completed. Converted {len(files)} SWC files to {output_dir}")

    def _replace_zero_based_index(self, swc_df):
        if swc_df['id'][0] == 1:
            swc_df['id'] = swc_df['id'] - 1
            swc_df['parent'] = swc_df['parent'].mask(swc_df['parent'] > 0, swc_df['parent']-1)
        return swc_df

    def _sort_ids_by_DFS(self, graph, current, visited, sorted_nodes):

        visited.add(current)
        sorted_nodes.append(current)

        for child in graph.get(current, []):
            if child not in visited:
                # recursive
                self._sort_ids_by_DFS(graph, child, visited, sorted_nodes)

        return sorted_nodes

    def _sort_by_DFS(self, swc_df):
        # create graph
        # key:parent, value:list(id)
        graph = swc_df[swc_df["parent"] != -1].groupby("parent")["id"].apply(list).to_dict()

        # get root node
        root = swc_df[swc_df['parent'] == -1]['id'].values[0]

        # sort ids by DFS
        sorted_ids = self._sort_ids_by_DFS(graph, root, set(), [])

        # create id mapping
        id_mapping = {old_id: new_id for new_id, old_id in enumerate(sorted_ids, 0)}

        # create new DataFrame
        new_df = swc_df.copy()
        new_df['new_id'] = new_df['id'].map(id_mapping)
        new_df['new_parent'] = new_df['parent'].map(lambda x: id_mapping.get(x, -1))

        # sort by new id
        new_df = new_df.sort_values('new_id')

        # create DataFrame
        final_df = pd.DataFrame({
            'id': new_df['new_id'],
            'type': new_df['type'],
            'x': new_df['x'],
            'y': new_df['y'],
            'z': new_df['z'],
            'r': new_df['r'],
            'parent': new_df['new_parent']
        })

        return final_df

    def _get_axon_direction_df(self, df):
        # Extract soma coordinates
        soma_df = df[df['type'] == 1]
        soma_end = soma_df.iloc[-1][['x', 'y', 'z']].values
        soma_mid = soma_df.iloc[len(soma_df)//2][['x', 'y', 'z']].values

        # Extract axon coordinates
        axon_df = df[df['type'] == 2]
        axon_p3d = axon_df[['x', 'y', 'z']].values

        # Combine soma and axon coordinates
        p3d = np.vstack((soma_mid.reshape(1, -1), axon_p3d))

        # Compute PCA
        pca = PCA(n_components=3)
        pca.fit(p3d)
        unit_v = pca.components_[0]

        # Normalize the vector
        unit_v = unit_v / np.linalg.norm(unit_v)

        # Determine direction
        axon_end = axon_p3d[-1] - soma_mid
        if np.dot(unit_v, axon_end) < 0:
            unit_v *= -1

        # Calculate new axon segment coordinates
        axon_seg_coor = np.zeros((2, 3))
        axon_seg_coor[0] = soma_end + (unit_v * 30.)
        axon_seg_coor[1] = soma_end + (unit_v * 60.)

        return axon_seg_coor

    def _fix_axon_perisomatic_directed_df(self, df):
        # Check if axon exists
        if 2 not in df['type'].values:
            raise Exception('There is no axonal reconstruction in swc file.')

        # Get new axon coordinates
        axon_coords = self._get_axon_direction_df(df)

        # Remove existing axon
        df = df[df['type'] != 2].reset_index(drop=True)

        # Add new axon segments
        max_id = df['id'].max()
        new_axon_df = pd.DataFrame({
            'id': range(max_id + 1, max_id + 3),
            'type': [2, 2],
            'x': axon_coords[:, 0],
            'y': axon_coords[:, 1],
            'z': axon_coords[:, 2],
            'r': [0.5, 0.5],
            'parent': [df[df['type'] == 1]['id'].iloc[-1], max_id + 1]
        })

        # Concatenate original df (without axon) and new axon df
        result_df = pd.concat([df, new_axon_df], ignore_index=True)

        return result_df

    def save_nodes(self, nodes_file_name=None, node_types_file_name=None, output_dir='.', force_overwrite=True, mode='w', compression='gzip'):
        """Save the instantiated nodes in SONATA format files and neulite format files.

        :param nodes_file_name: file-name of hdf5 nodes file. By default will use <network.name>_nodes.h5.
        :param node_types_file_name: file-name of the csv node-types file. By default will use <network.name>_node_types.csv
        :param output_dir: Directory where network files will be generated. Default, current working directory.
        :param force_overwrite: Overwrites existing network files.
        :param compression: Compression format.
        """
        logger.info("save_nodes start")
        
        population_file = os.path.join(self.neulite_dir, f'{self.name}_population.csv')
        if not force_overwrite and os.path.exists(population_file):
            raise Exception('File {} already exists. Please delete existing file, use a different name, or use force_overwrite.'.format(population_file))
        nf_dir = os.path.dirname(population_file)

        if mpi_rank == 0:
            # Create neulite directory if it doesn't exist
            if not os.path.exists(self.neulite_dir):
                os.makedirs(self.neulite_dir)
            # Create subdirectory if it doesn't exist    
            data_dir = os.path.join(self.neulite_dir, self.swc_dir)
            if not os.path.exists(data_dir):
                os.makedirs(data_dir)

        super(NeuliteNetwork, self).save_nodes(nodes_file_name, node_types_file_name, output_dir, force_overwrite, mode, compression)

        barrier()

        self._save_population(population_file)

        # Create empty connection CSV file
        if mpi_rank == 0:
            self._create_empty_connection_files()

        barrier()
        logger.info("save_nodes end")

    def _save_population(self, population_file):
        logger.info("_save_population start")
        if mpi_rank == 0:
            logger.debug('Saving {} node-types to {}.'.format(self.name, population_file))

            output_data = pd.DataFrame(columns=["#n_cell", "n_comp", "name", "swc_file", "ion_file"])
            population_list = []

            # make node_type_id list
            node_type_id_list = [node.node_type_id for node in self.nodes()]
            # count make_type_id
            node_type_counter = Counter(node_type_id_list)
            for model_id in node_type_counter:
                cnt = node_type_counter[model_id]
                node_types_property = self._node_types_properties[model_id]
                
                # Skip non-biophysical models
                if node_types_property.get("model_type") != "biophysical":
                    logger.debug(f"Skipping non-biophysical model: {model_id=}, model_type={node_types_property.get('model_type')}")
                    continue

                # One pop_name may use different morphologies.
                pop_name = node_types_property.get("pop_name", f"default_{model_id}")
                model_name = f'{pop_name}_{str(model_id)}'

                # ion file path
                ion_file = node_types_property["dynamics_params"]
                ion_file = f"{ion_file[:-5]}.csv"
                ion_file = os.path.join(self.ion_dir, ion_file)

                n_comp, swc_file = self._get_swc(model_id)
                population = {"#n_cell": cnt, "n_comp": n_comp, "name": model_name, "swc_file": swc_file, "ion_file": ion_file}
                population_list.append(population)

            output_data = pd.concat([output_data, pd.DataFrame(population_list)])
            output_data.to_csv(population_file, index=False)
            logger.info(f"Created Neulite population file: {population_file}")

        barrier()
        logger.info("_save_population end")

    def _create_empty_connection_files(self):
        """Create empty connection CSV file with header only.

        Creates empty connection file for self-connection ({name}_{name}).
        If save_edges() is called later, this file will be overwritten.
        """
        connection_file = os.path.join(self.neulite_dir,
                                        f"{self.name}_{self.name}_connection.csv")

        if not os.path.exists(connection_file):
            logger.info(f"Creating empty connection file: {connection_file}")
            empty_df = pd.DataFrame(columns=self.CONNECTION_CSV_HEADER)
            empty_df.to_csv(connection_file, index=False)
            logger.info(f"Created empty connection file: {connection_file}")
        else:
            logger.debug(f"Connection file already exists, skipping: {connection_file}")

    def _get_swc(self, model_id):
        # set default value
        n_comp = 0
        swc_file = "nan"

        node_types_property = self._node_types_properties[model_id]
        if "morphology" not in node_types_property:
            logger.warning(f"morphology does not exists. [{model_id=}]")
            return n_comp, swc_file

        model_morph = node_types_property["morphology"]
        if model_morph == None:
            logger.warning(f"morphology is None. [{model_id=}]")
            return n_comp, swc_file

        # swc file path
        swc_file = os.path.join(self.morphologies_dir, model_morph)

        if not os.path.exists(swc_file):
            logger.warning(f"swc file does not exists. [{swc_file=}]")
            return n_comp, swc_file

        # use converted swc file path
        converted_swc_file = os.path.join(self.neulite_dir, self.swc_dir, model_morph)

        if os.path.exists(converted_swc_file):
            # if converted swc file exists, use it
            data_morph = pd.read_csv(converted_swc_file, delimiter=" ", comment="#", header=None)
            n_comp = data_morph.shape[0]
        else:
            # if converted swc file does not exist, use original swc file (fallback)
            if not os.path.exists(swc_file):
                logger.warning(f"swc file does not exists. [{swc_file=}]")
                return n_comp, swc_file
            data_morph = pd.read_csv(swc_file, delimiter=" ", comment="#", header=None)
            n_comp = data_morph.shape[0]

        swc_file_for_out = os.path.join(self.swc_dir, model_morph)
        return n_comp, swc_file_for_out

    def _save_edges(self, edges_file_name, src_network, trg_network, pop_name=None, sort_by='target_node_id',
                    index_by=('target_node_id', 'source_node_id'), compression='gzip', mode='w'):
        barrier()

        if compression == 'none':
            compression = None  # legit option for h5py for no compression

        if mpi_rank == 0:
            logger.debug('Saving {} --> {} edges to {}.'.format(src_network, trg_network, edges_file_name))
            logger.info("_save_edges start")
            logger.debug(f"{mpi_size=}")

        filtered_edge_types = [
            # Some edges may not match the source/target population
            et for et in self._DenseNetwork__edges_tables
            if et.source_network == src_network and et.target_network == trg_network
        ]

        merged_edges = EdgesCollator(filtered_edge_types, network_name=self.name)
        merged_edges.process()
        n_total_conns = merged_edges.n_total_edges
        barrier()

        if n_total_conns == 0:
            if mpi_rank == 0:
                logger.warning('Was not able to generate any edges using the "connection_rule". Not saving.')
            return

        # Try to sort before writing file, If edges are split across ranks/files for MPI/size issues then we need to
        # write to disk first then sort the hdf5 file
        sort_on_disk = False
        edges_file_name_final = edges_file_name
        edges_file_dirname = os.path.dirname(edges_file_name)
        edges_file_basename = os.path.basename(edges_file_name)
        if sort_by:
            if merged_edges.can_sort:
                merged_edges.sort(sort_by=sort_by)
            else:
                sort_on_disk = True
                edges_file_name_final = edges_file_name

                edges_file_name = os.path.join(edges_file_dirname, '.unsorted.{}'.format(edges_file_basename))
                if mpi_rank == 0:
                    logger.debug('Unable to sort edges in memory, will temporarly save to {}'.format(edges_file_name) +
                                 ' before sorting hdf5 file.')
        barrier()

        start_time = time.time()
        if mpi_rank == 0:
            logger.debug('Saving {} edges to disk'.format(n_total_conns))

            # neulite start
            start_time = time.time()
            logger.info("Starting the processing...")

            logger.info("Loading edge types data...")
            matching_et = [c.edge_type_properties for c in self._connection_maps
                           if c.source_network_name == src_network and c.target_network_name == trg_network]
            edge_types_list = []
            cols = ["edge_type_id", "syn_weight", "target_sections", "delay", "tau2", "tau1", "erev"]
            for edge_type in matching_et:
                # get basic parameters
                base_params = [edge_type.get(cname, 'NULL') if edge_type.get(cname, 'NULL') is not None else 'NULL' for cname in cols[:-3]]
                
                # extract tau1, tau2, erev from dynamics_params
                dynamics_params = edge_type.get('dynamics_params')
                tau1, tau2, erev = self._extract_synapse_params(dynamics_params)
                
                # add all parameters to list
                base_params.extend([tau1, tau2, erev])
                edge_types_list.append(base_params)
            edge_types_data = pd.DataFrame(edge_types_list, columns=cols)
            edge_types_data["delay"] = edge_types_data["delay"].map(lambda x: Decimal(x).quantize(Decimal('1'), ROUND_HALF_UP))
            logger.debug(f"{len(edge_types_data)=}")

            node_id_table = np.zeros(self._nnodes)  # todo: set dtypes
            node_type_id_table = np.zeros(self._nnodes)

            for i, node in enumerate(self.nodes()):
                node_id_table[i] = node.node_id  # node_id
                node_type_id_table[i] = node.node_type_id  # node_type_id

            h5_node_type = pd.DataFrame(node_type_id_table, columns=["node_type_id"])
            node_ids = pd.DataFrame(node_id_table, columns=["node_id"])
            h5_node_type = pd.concat([node_ids, h5_node_type], axis=1)

            # for MPI
            h5_node_type = h5_node_type.to_dict()

            # neulite end

            pop_name = '{}_to_{}'.format(src_network, trg_network) if pop_name is None else pop_name
            with h5py.File(edges_file_name, 'w') as hf:
                # Initialize the hdf5 groups and datasets
                add_hdf5_attrs(hf)
                pop_grp = hf.create_group('/edges/{}'.format(pop_name))

                pop_grp.create_dataset('source_node_id', (n_total_conns,), dtype='uint64', compression=compression)
                pop_grp['source_node_id'].attrs['node_population'] = src_network
                pop_grp.create_dataset('target_node_id', (n_total_conns,), dtype='uint64', compression=compression)
                pop_grp['target_node_id'].attrs['node_population'] = trg_network
                pop_grp.create_dataset('edge_group_id', (n_total_conns,), dtype='uint16', compression=compression)
                pop_grp.create_dataset('edge_group_index', (n_total_conns,), dtype='uint32', compression=compression)
                pop_grp.create_dataset('edge_type_id', (n_total_conns,), dtype='uint32', compression=compression)

                for group_id in merged_edges.group_ids:
                    # different model-groups will have different datasets/properties depending on what edge information
                    # is being saved for each edges
                    model_grp = pop_grp.create_group(str(group_id))
                    for prop_mdata in merged_edges.get_group_metadata(group_id):
                        model_grp.create_dataset(prop_mdata['name'], shape=prop_mdata['dim'], dtype=prop_mdata['type'], compression=compression)

                # neulite start
                connections_list = []
                # neulite end

                # Uses the collated edges (eg combined edges across all edge-types) to actually write the data to hdf5,
                # potentially in multiple chunks. For small networks doing it this way isn't very effiecent, however
                # this has the benefits:
                #  * For very large networks it won't always be possible to store all the data in memory.
                #  * When using MPI/multi-node the chunks can represent data from different ranks.
                for chunk_id, idx_beg, idx_end in merged_edges.itr_chunks():
                    pop_grp['source_node_id'][idx_beg:idx_end] = merged_edges.get_source_node_ids(chunk_id)
                    pop_grp['target_node_id'][idx_beg:idx_end] = merged_edges.get_target_node_ids(chunk_id)
                    pop_grp['edge_type_id'][idx_beg:idx_end] = merged_edges.get_edge_type_ids(chunk_id)
                    pop_grp['edge_group_id'][idx_beg:idx_end] = merged_edges.get_edge_group_ids(chunk_id)
                    pop_grp['edge_group_index'][idx_beg:idx_end] = merged_edges.get_edge_group_indices(chunk_id)

                    for group_id, prop_name, grp_idx_beg, grp_idx_end in merged_edges.get_group_data(chunk_id):
                        prop_array = merged_edges.get_group_property(prop_name, group_id, chunk_id)
                        pop_grp[str(group_id)][prop_name][grp_idx_beg:grp_idx_end] = prop_array

                    # neulite start
                    logger.info("Loading edge data...")
                    edge_type_id = pd.DataFrame(pop_grp['edge_type_id'][idx_beg:idx_end], columns=["edge_type_id"])
                    source_node_id = pd.DataFrame(pop_grp['source_node_id'][idx_beg:idx_end], columns=["pre nid"])
                    target_node_id = pd.DataFrame(pop_grp['target_node_id'][idx_beg:idx_end], columns=["post nid"])
                    nsyns_data = pd.DataFrame(pop_grp["0/nsyns"][idx_beg:idx_end], columns=["nsyns"], dtype="object")

                    # Merging edge types data on edge_type_id
                    logger.info("Merging edge types data...")
                    connections = pd.merge(edge_type_id, edge_types_data, on='edge_type_id', how='left')
                    logger.debug(f"afrer marge: {chunk_id=}, {len(connections)=}")
                    connections = pd.concat([connections, nsyns_data], axis=1)
                    logger.debug(f"after concat: {chunk_id=}, {len(connections)=}")
                    connections.insert(0, "pre nid", source_node_id)
                    connections.insert(1, "post nid", target_node_id)
                    logger.debug(f"after insert: {chunk_id=}, {len(connections)=}")
                    connections_list.append(connections)

                logger.debug(f"before concat: {len(connections_list)=}")
                merged_connections = pd.concat(connections_list)
                del connections_list
                logger.debug(f"after concat: {len(merged_connections)=}")
                merged_connections["target_sections"] = merged_connections["target_sections"].apply(self._str_to_list)
                merged_connections["target_sections"] = merged_connections["target_sections"].apply(self._replace_strings)
                merged_connections["post cid"] = -1
                merged_connections["ei"] = ""
                merged_connections["nsyns"] = merged_connections["nsyns"].apply(self._to_list)
                merged_connections = merged_connections.explode("nsyns")
                logger.info(f"Finished merging edge types data.")

                # set chunk_size
                chunk_size = (len(merged_connections) // mpi_size ) + 1
                chunks = [merged_connections.iloc[i*chunk_size:(i+1)*chunk_size] for i in range(mpi_size)]
                logger.debug(f"{mpi_size=}, {chunk_size=}, {len(chunks)=}")
                # for MPI
                serialized_chunks = [self._serialize_data(chunk) for chunk in chunks]
                logger.info("Processing chunks...")
                # neulite end

            if sort_on_disk:
                logger.debug('Sorting {} by {} to {}'.format(edges_file_name, sort_by, edges_file_name_final))
                sort_edges(
                    input_edges_path=edges_file_name,
                    output_edges_path=edges_file_name_final,
                    edges_population='/edges/{}'.format(pop_name),
                    sort_by=sort_by,
                    compression=compression,
                    # sort_on_disk=True,
                )
                try:
                    logger.debug('Deleting intermediate edges file {}.'.format(edges_file_name))
                    os.remove(edges_file_name)
                except OSError as e:
                    logger.warning('Unable to remove intermediate edges file {}.'.format(edges_file_name))

            if index_by:
                index_by = index_by if isinstance(index_by, (list, tuple)) else [index_by]
                for index_type in index_by:
                    logger.debug('Creating index {}'.format(index_type))
                    create_index_in_memory(
                        edges_file=edges_file_name_final,
                        edges_population='/edges/{}'.format(pop_name),
                        index_type=index_type,
                        compression=compression
                    )

        else:
            h5_node_type = None
            serialized_chunks = None

        barrier()
        del merged_edges

        # neulite start
        # Use "send" to prevent overflow.
        if mpi_size == 1:
            # Single process mode: no MPI communication needed
            (serialized_chunk,) = serialized_chunks
            chunk = self._deserialize_data(serialized_chunk)
        elif mpi_rank == 0:
            # Multi-process mode: rank 0 distributes chunks
            for i, attr in enumerate(serialized_chunks):
                if i == 0:
                    chunk = self._deserialize_data(attr)
                else:
                    if comm is not None:
                        comm.send(attr, dest=i)
        else:
            # Multi-process mode: other ranks receive chunks
            if comm is not None:
                serialized_chunk = comm.recv(source=0)
                chunk = self._deserialize_data(serialized_chunk)
                logger.debug(f"chunk is received:rank{mpi_rank}")
        # neulite end

        if mpi_size == 1:
            # Single process mode: h5_node_type is already available, just convert to DataFrame
            h5_node_type = pd.DataFrame.from_dict(h5_node_type)
        else:
            # Multi-process mode: use MPI broadcast
            try:
                if comm is not None:
                    h5_node_type = comm.bcast(h5_node_type, root=0)
                # dict to DataFrame
                h5_node_type = pd.DataFrame.from_dict(h5_node_type)
            except Exception as e:
                logger.error(f"Error during bcast: {e}")
                if comm is not None:
                    comm.Abort()
                raise

        barrier()

        # neulite start
        file_cache = {}
        # Calculate path to converted (preprocessed) SWC files
        converted_swc_path = os.path.join(self.neulite_dir, self.swc_dir)
        processed_chunk = self._process_chunk(chunk, file_cache, h5_node_type, self._node_types_properties, self.morphologies_dir, converted_morphologies_path=converted_swc_path)
        end_time = time.time()
        logger.info(f"Processing completed in {end_time - start_time} seconds: from rank{mpi_rank}")

        barrier()

        # Collect processed chunks from all ranks
        if mpi_size == 1:
            # Single process mode: no MPI communication needed
            processed_list = [processed_chunk]
        elif mpi_rank == 0:
            # Multi-process mode: rank 0 collects chunks from other ranks
            processed_list = [processed_chunk]
            for i in range(1, mpi_size):
                # receive from rank i
                if comm is not None:
                    serialized_chunk = comm.recv(source=i)
                    processed_chunk = self._deserialize_data(serialized_chunk)
                    logger.debug(f"processed_chunk is received:rank{mpi_rank}")
                    processed_list.append(processed_chunk)
        else:
            # Multi-process mode: other ranks send to rank 0
            processed_list = None
            # send to rank 0
            if comm is not None:
                comm.send(self._serialize_data(processed_chunk), dest=0)
        # neulite end

        barrier()

        if mpi_rank == 0:
            # Filter out None values from chunks
            valid_chunks = [chunk for chunk in processed_list if chunk is not None]
            
            if not valid_chunks:
                logger.warning(f"No valid edges to save for {src_network} -> {trg_network}")
                logger.info("_save_edges end")
                return
                
            output_data = pd.concat(valid_chunks)
            del processed_list
            if sort_by:
                output_data = output_data.sort_values(["#pre nid", "post nid", "post cid"])
            output_path = os.path.join(self.neulite_dir, f"{src_network}_{trg_network}_connection.csv")
            output_data.to_csv(output_path, mode='w', header=True, index=False)
            logger.info(f"Created Neulite connection file: {output_path}")
            logger.info("_save_edges end")
            logger.debug('Saving completed.')

    def _extract_synapse_params(self, dynamics_params, default_tau1=0.1, default_tau2=1.7, default_erev=0.0):
        """Extract tau1, tau2, and erev from dynamics_params
        
        :param dynamics_params: JSON file path or dictionary containing synapse parameters
        :param default_tau1: Default tau1 value if not found
        :param default_tau2: Default tau2 value if not found
        :param default_erev: Default erev value if not found
        :return: tuple (tau1, tau2, erev)
        """
        if not dynamics_params:
            return default_tau1, default_tau2, default_erev
        
        try:
            if isinstance(dynamics_params, str):
                # if JSON file
                params_path = dynamics_params
                if not os.path.isabs(params_path):
                    # if relative path, combine with synaptic_models_dir
                    config = SonataConfig.from_json(self._get_config_path())
                    synaptic_models_dir = config["components"]["synaptic_models_dir"]
                    params_path = os.path.join(synaptic_models_dir, dynamics_params)
                
                with open(params_path, 'r') as f:
                    params = json.load(f)
            else:
                # if dictionary
                params = dynamics_params
            
            tau1 = params.get('tau1', default_tau1)
            tau2 = params.get('tau2', default_tau2)
            erev = params.get('erev', default_erev)
            return tau1, tau2, erev
            
        except Exception as e:
            logger.error(f"Failed to extract synapse parameters from {dynamics_params}: {e}")
            raise

    @staticmethod
    def _replace_strings(lst):
        if isinstance(lst, list):
            mapping = {"somatic": 1, "axon": 2, "basal": 3, "apical": 4}
            return [mapping.get(item, item) for item in lst]
        return lst

    @staticmethod
    def _str_to_list(text):
        if text == "NULL":
            return None
        if isinstance(text, str):
            return ast.literal_eval(text)
        return text

    @staticmethod
    def _get_cid(row, file_cache, h5_node_type, node_types, morphologies_path, converted_morphologies_path=None):
        try:
            post_node_id = row["post nid"]
            post_node_type = h5_node_type.loc[post_node_id, "node_type_id"]
            sections = row["target_sections"]
            if post_node_type not in node_types or node_types[post_node_type].get("morphology") == None:
                return -1

            morph_path = node_types[post_node_type].get("morphology")

            if morph_path in file_cache:
                morph_data = file_cache[morph_path]
            else:
                # Try to use converted (preprocessed) SWC file first
                if converted_morphologies_path:
                    converted_morph_path_full = os.path.join(converted_morphologies_path, morph_path)
                    if os.path.exists(converted_morph_path_full):
                        morph_path_full = converted_morph_path_full
                        logger.debug(f"Using converted SWC file: {morph_path_full}")
                    else:
                        # Fallback: use original file
                        morph_path_full = os.path.join(morphologies_path, morph_path)
                        logger.debug(f"Converted SWC not found, using original: {morph_path_full}")
                else:
                    # If converted_morphologies_path is not provided, use original file
                    morph_path_full = os.path.join(morphologies_path, morph_path)
                    logger.debug(f"Using original SWC file: {morph_path_full}")

                try:
                    morph_data = pd.read_csv(morph_path_full, delimiter=" ", comment="#", usecols=["id", "type"], names=["id", "type", "x", "y", "z", "r", "p"])
                    file_cache[morph_path] = morph_data
                except Exception as e:
                    logger.error(f"Error reading file {morph_path_full}: {e}")
                    raise

            # random sampling
            ret = morph_data[morph_data["type"].isin(sections)].sample()["id"].iloc[0]
            return ret
        except Exception as e:
            logger.error(f"Error _get_cid: {e}")
            logger.error(row)
            raise

    @staticmethod
    def _get_ei(row, h5_node_type, node_types):
        try:
            pre_node = row["pre nid"]
            pre_node_type = h5_node_type.loc[pre_node, "node_type_id"]
            # Check if 'ei' key exists, default to 'e' if not found
            ei_value = node_types[pre_node_type].get("ei", "e")
            if ei_value not in ["e", "i"]:
                logger.warning(f"Unexpected ei value '{ei_value}' for node_type {pre_node_type}, defaulting to 'e'")
                return "e"
            return ei_value
        except Exception as e:
            logger.error(f"Error in _get_ei: {e}, pre_node={row.get('pre nid', 'unknown')}")
            return "e"  # Default to excitatory

    @staticmethod
    def _is_biophysical_node(node_id, h5_node_type, node_types_data):
        """Check if a node is biophysical type"""
        try:
            node_type_id = h5_node_type.loc[node_id, "node_type_id"]
            return node_types_data.get(node_type_id, {}).get("model_type") == "biophysical"
        except:
            return False

    @staticmethod
    def _process_chunk(chunk, file_cache, h5_node_type, node_types_data, morphologies_path, converted_morphologies_path=None):
        try:
            logger.debug(f"Start processing chunk...")
            chunk = chunk.copy()
            logger.debug(f"End coppy chunk...")
            chunk.loc[:, "post cid"] = chunk.apply(NeuliteNetwork._get_cid, axis=1, file_cache=file_cache, h5_node_type=h5_node_type, node_types=node_types_data, morphologies_path=morphologies_path, converted_morphologies_path=converted_morphologies_path)
            logger.debug(f"End _get_cid chunk...")
            chunk.loc[:, "ei"] = chunk.apply(NeuliteNetwork._get_ei, axis=1, h5_node_type=h5_node_type, node_types=node_types_data)
            logger.debug(f"End _get_ei chunk...")
            
            # Filter edges to keep only biophysical->biophysical connections
            # Check source (pre) nodes
            pre_biophysical_mask = chunk["pre nid"].apply(
                lambda nid: NeuliteNetwork._is_biophysical_node(nid, h5_node_type, node_types_data)
            )
            
            # Check target (post) nodes (post_cid == -1 means non-biophysical)
            post_biophysical_mask = chunk["post cid"] != -1
            
            # Keep only edges where both source and target are biophysical
            biophysical_edges_mask = pre_biophysical_mask & post_biophysical_mask
            
            # Log filtering statistics
            total_edges = len(chunk)
            biophysical_edges = biophysical_edges_mask.sum()
            filtered_count = total_edges - biophysical_edges
            
            if filtered_count > 0:
                logger.info(f"Filtering out {filtered_count} edges involving non-biophysical nodes")
                logger.debug(f"  - Non-biophysical source: {(~pre_biophysical_mask).sum()}")
                logger.debug(f"  - Non-biophysical target: {(~post_biophysical_mask).sum()}")
                chunk = chunk[biophysical_edges_mask]
            
            if len(chunk) == 0:
                logger.info("All edges were filtered out")
                return None
                
            chunk.drop(["nsyns", "target_sections"], axis=1, inplace=True)
            chunk = chunk.rename(columns={"pre nid": "#pre nid", "syn_weight": "weight", "ei":"e/i", "tau2": "tau_decay", "tau1": "tau_rise"})
            chunk = chunk[NeuliteNetwork.CONNECTION_CSV_HEADER]
            logger.debug(f"Finished processing chunk.")
            return chunk
        except Exception as e:
            logger.error(f"Error processing chunk : {e}")
            return None

    @staticmethod
    def _on_task_complete(result):
        if result is not None:
            logger.info("Task complete")
        else:
            logger.error("Task failed")

    @staticmethod
    def _to_list(n):
            return list(range(1, n+1))
        
    @staticmethod
    def _serialize_data(data):
        return pickle.dumps(data)

    @staticmethod
    def _deserialize_data(data):
        return pickle.loads(data)
    
    
    def convert_ion_channels(self, input_dir=None, output_dir=None):
        """Convert ion channel JSON files to CSV format for Neulite.
        
        :param input_dir: Directory containing *_fit.json files. If None, uses biophysical_neuron_models_dir from config
        :param output_dir: Directory where CSV files will be saved. If None, uses ion_dir
        """
        if input_dir is None:
            input_dir = self.biophysical_neuron_models_dir
        if output_dir is None:
            output_dir = self.ion_dir
            
        files = glob.glob(os.path.join(input_dir, "*_fit.json"))
        channel_list = ['Cm', 'Ra', 'leak', 'e_pas', 'gamma', 'decay', 'NaV', 'NaTs', 'NaTa', 'Nap',
                        'Kv2like', 'Kv3_1', 'K_P', 'K_T', 'Kd', 'Im', 'Im_v2', 'Ih', 'SK', 'Ca_HVA', 'Ca_LVA']

        # Create output directory if it doesn't exist
        actual_output_dir = os.path.join(self.neulite_dir, output_dir)
        if not os.path.exists(actual_output_dir):
            os.makedirs(actual_output_dir)

        for file in files:
            csv_data = pd.DataFrame(0.0, index=range(4), columns=channel_list, dtype=np.float64)

            with open(file, 'r') as json_open:
                json_load = json.load(json_open)

            passive = pd.DataFrame(json_load['passive'])
            ion_data = pd.DataFrame(json_load['genome'])

            # Insert values of ra (genome first)
            ra_genome = ion_data[ion_data['name'] == 'ra']
            if not ra_genome.empty:
                for data in ra_genome.iterrows():
                    if data[1]['section'] == 'soma':
                        csv_data.loc[0, 'Ra'] = float(data[1]['value'])
                    elif data[1]['section'] == 'axon':
                        csv_data.loc[1, 'Ra'] = float(data[1]['value'])
                    elif data[1]['section'] == 'dend':
                        csv_data.loc[2, 'Ra'] = float(data[1]['value'])
                        csv_data.loc[3, 'Ra'] = float(data[1]['value'])
                    else:
                        csv_data.loc[3, 'Ra'] = float(data[1]['value'])
            else:
                csv_data['Ra'] = float(passive['ra'].iloc[-1])

            # Insert values of e_pas (genome first)
            e_pas_genome = ion_data[ion_data['name'] == 'e_pas']
            if not e_pas_genome.empty:
                for data in e_pas_genome.iterrows():
                    if data[1]['section'] == 'soma':
                        csv_data.loc[0, 'e_pas'] = float(data[1]['value'])
                    elif data[1]['section'] == 'axon':
                        csv_data.loc[1, 'e_pas'] = float(data[1]['value'])
                    elif data[1]['section'] == 'dend':
                        csv_data.loc[2, 'e_pas'] = float(data[1]['value'])
                        csv_data.loc[3, 'e_pas'] = float(data[1]['value'])
                    else:
                        csv_data.loc[3, 'e_pas'] = float(data[1]['value'])
            else:
                csv_data['e_pas'] = float(passive['e_pas'].iloc[-1])

            # Insert values of cm (genome first)
            cm_genome = ion_data[ion_data['name'] == 'cm']
            if not cm_genome.empty:
                for data in cm_genome.iterrows():
                    if data[1]['section'] == 'soma':
                        csv_data.loc[0, 'Cm'] = float(data[1]['value'])
                    elif data[1]['section'] == 'axon':
                        csv_data.loc[1, 'Cm'] = float(data[1]['value'])
                    elif data[1]['section'] == 'dend':
                        csv_data.loc[2, 'Cm'] = float(data[1]['value'])
                        csv_data.loc[3, 'Cm'] = float(data[1]['value'])
                    else:
                        csv_data.loc[3, 'Cm'] = float(data[1]['value'])
            else:
                for data in passive.loc[0, 'cm']:
                    if data['section'] == 'soma':
                        csv_data.loc[0, 'Cm'] = float(data['cm'])
                    elif data['section'] == 'axon':
                        csv_data.loc[1, 'Cm'] = float(data['cm'])
                    elif data['section'] == 'dend':
                        csv_data.loc[2, 'Cm'] = float(data['cm'])
                        csv_data.loc[3, 'Cm'] = float(data['cm'])
                    else:
                        csv_data.loc[3, 'Cm'] = float(data['cm'])

            # Insert values of leak
            g_pas = ion_data[ion_data['name'] == 'g_pas']
            for data in g_pas.iterrows():
                if data[1]['section'] == 'soma':
                    csv_data.loc[0, 'leak'] = float(data[1]['value'])
                elif data[1]['section'] == 'axon':
                    csv_data.loc[1, 'leak'] = float(data[1]['value'])
                elif data[1]['section'] == 'dend':
                    csv_data.loc[2, 'leak'] = float(data[1]['value'])
                    csv_data.loc[3, 'leak'] = float(data[1]['value'])
                else:
                    csv_data.loc[3, 'leak'] = float(data[1]['value'])

            # Insert values of gamma
            gamma = ion_data[ion_data['name'] == 'gamma_CaDynamics']
            csv_data['gamma'] = 0.05  # Initialize to 0.05
            for data in gamma.iterrows():
                if data[1]['section'] == 'soma':
                    csv_data.loc[0, 'gamma'] = float(data[1]['value'])
                elif data[1]['section'] == 'axon':
                    csv_data.loc[1, 'gamma'] = float(data[1]['value'])
                elif data[1]['section'] == 'dend':
                    csv_data.loc[2, 'gamma'] = float(data[1]['value'])
                    csv_data.loc[3, 'gamma'] = float(data[1]['value'])
                else:
                    csv_data.loc[3, 'gamma'] = float(data[1]['value'])

            # Insert values of decay
            decay = ion_data[ion_data['name'] == 'decay_CaDynamics']
            csv_data['decay'] = 80.0  # Initialize to 80.0
            for data in decay.iterrows():
                if data[1]['section'] == 'soma':
                    csv_data.loc[0, 'decay'] = float(data[1]['value'])
                elif data[1]['section'] == 'axon':
                    csv_data.loc[1, 'decay'] = float(data[1]['value'])
                elif data[1]['section'] == 'dend':
                    csv_data.loc[2, 'decay'] = float(data[1]['value'])
                    csv_data.loc[3, 'decay'] = float(data[1]['value'])
                else:
                    csv_data.loc[3, 'decay'] = float(data[1]['value'])

            # Insert values of ion channels
            for channel in channel_list[6:]:
                value = ion_data[ion_data['mechanism'] == channel]['value']
                csv_data.loc[0, channel] = 0.0 if len(value) == 0 else float(value.iloc[-1])

            csv_data.index = csv_data.index + 1
            output_file = os.path.join(actual_output_dir, os.path.splitext(os.path.basename(file))[0] + ".csv")
            csv_data.to_csv(output_file, header=False)
            logger.debug(f"Created ion channel CSV: {output_file}")
        
        logger.info(f"Ion channel conversion completed. Created {len(files)} CSV files in {actual_output_dir}")

    def get_simulation_params_from_config(self, config):
        """Extract simulation parameters from SonataConfig object

        :param config: SonataConfig object
        :return: dict of simulation parameters
        """
        # Extract run parameters with defaults
        dt = config.run.get('dt', 0.1) if hasattr(config, 'run') else 0.1
        tstop = config.run.get('tstop', 1000.0) if hasattr(config, 'run') else 1000.0
        spike_threshold = config.run.get('spike_threshold', -15.0) if hasattr(config, 'run') else -15.0

        # Check inputs section for current clamp
        current_injection_params = {
            'amp': 0.1,        # Default values
            'delay': 500.0,
            'duration': 500.0
        }

        if hasattr(config, 'inputs'):
            # Look for current_clamp type input
            for input_name, input_params in config.inputs.items():
                if input_params.get('input_type') == 'current_clamp' or \
                   input_params.get('module') == 'IClamp':
                    # Found current_clamp
                    current_injection_params = {
                        'amp': input_params.get('amp', 0.1),
                        'delay': input_params.get('delay', 500.0),
                        'duration': input_params.get('duration', 500.0)
                    }
                    break  # Use first current_clamp found

        return {
            'dt': dt,
            'tstop': tstop,
            'spike_threshold': spike_threshold,
            'current_injection': current_injection_params
        }

    def generate_config_h(self, config, output_path=None):
        """Generate config.h file (v15 compatible)

        :param config: SonataConfig object
        :param output_path: Output file path (if None, uses neulite_dir/kernel/config.h)
        """
        if mpi_rank != 0:
            return

        # Set default output path if not provided
        if output_path is None:
            output_path = os.path.join(self.neulite_dir, "kernel", "config.h")

        # Extract parameters from SonataConfig
        params = self.get_simulation_params_from_config(config)
        logger.info("Loaded parameters from SonataConfig object")

        lines = [
            "// Automatically generated by bionet_lite",
            "",
            "#pragma once",
            "",
            "#undef DEBUG",
            "",
            "// Simulation parameters",
            f"#define TSTOP ( {params['tstop']} )",
            f"#define DT ( {params['dt']} )",
            f"#define INV_DT ( ( int ) ( 1.0 / ( DT ) ) )",
            "",
            "// Neuron parameters",
            f"#define SPIKE_THRESHOLD ( {params['spike_threshold']} )",
            "#define ALLACTIVE ( 0 ) // Set to 1 for allactive models",
            "",
            "// Current injection parameters",
            f"#define I_AMP ( {params['current_injection']['amp']} )",
            f"#define I_DELAY ( {params['current_injection']['delay']} )",
            f"#define I_DURATION ( {params['current_injection']['duration']} )",
        ]

        # Create directory if needed
        os.makedirs(os.path.dirname(output_path) if os.path.dirname(output_path) else '.', exist_ok=True)

        # Write to file
        with open(output_path, 'w') as f:
            f.write('\n'.join(lines))
            f.write('\n')  # Add final newline

        logger.info(f"Generated config.h: {output_path}")


class NeuliteBuilder(NetworkBuilder):

    def __init__(self, name, adaptor_cls=NeuliteNetwork, convert_morphologies=True, convert_ion_channels=True, neulite_only=True,
                 simulation_config=None, generate_config_h=True, **network_props):
        self.adaptor = adaptor_cls(name, convert_morphologies=convert_morphologies, convert_ion_channels=convert_ion_channels,
                                   neulite_only=neulite_only, simulation_config=simulation_config,
                                   generate_config_h=generate_config_h, **network_props)

