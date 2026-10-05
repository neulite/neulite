# SPDX-License-Identifier: GPL-2.0-only
# Copyright (C) 2025 Neulite Core Team <neulite-core@numericalbrain.org>

import ast
from collections import Counter
from decimal import Decimal, ROUND_HALF_UP
import glob
import h5py
import hashlib
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
import traceback

from bmtk.builder import NetworkBuilder
from bmtk.builder.network_adaptors import DenseNetwork
from bmtk.builder.network_adaptors.edges_collator import EdgesCollator
from bmtk.builder.index_builders import create_index_in_memory
from bmtk.builder.builder_utils import comm, mpi_rank, mpi_size, barrier, build_time_uuid
from bmtk.builder.edges_sorter import sort_edges
from bmtk.builder.builder_utils import add_hdf5_attrs
from bmtk.utils.sonata.config import SonataConfig

logger = logging.getLogger(__name__)

# Tutorial-specific simulation parameters
# These are applied automatically based on base_dir pattern matching
# to maintain "import-only change" design philosophy

# Default parameters (used when no pattern matches, and for Tutorial2,3)
BIONETLITE_DEFAULTS = {
    'tstop': 1000.0,
    'dt': 0.1,
    'current_clamp': {
        'amp': 0.1,
        'delay': 500.0,
        'duration': 500.0
    }
}

TUTORIAL_PARAMS = {
    'sim_ch01': {  # Tutorial1: Single cell with current injection
        'tstop': 2000.0,
        'dt': 0.1,
        'report_vars': ['v', 'cai'],
        'current_clamp': {
            'amp': 0.12,
            'delay': 500.0,
            'duration': 1000.0
        }
    },
    'sim_ch02': BIONETLITE_DEFAULTS,  # Tutorial2: No current injection
    'sim_ch03': BIONETLITE_DEFAULTS,  # Tutorial3: No current injection
    'sim_ch04': {  # Tutorial4: E/I network
        'tstop': 3000.0,
        'dt': 0.1,
        'report_vars': ['v'],
        'current_clamp': {
            'amp': 0.1,
            'delay': 500.0,
            'duration': 2000.0
        }
    },
}


_WEIGHT_FUNCTIONS = {}
_DIRECTION_RULE_NAMES = {"DirectionRule_EE", "DirectionRule_others"}


def register_weight_function(name):
    """Register a chunk-level weight function under a SONATA weight_function name."""
    def decorator(func):
        _WEIGHT_FUNCTIONS[name] = func
        return func
    return decorator


def normalize_weight_function_name(name):
    if name is None:
        return "wmax"
    if isinstance(name, str):
        stripped = name.strip()
        if stripped in ("", "NULL", "None", "nan"):
            return "wmax"
        return stripped
    try:
        if pd.isna(name):
            return "wmax"
    except (TypeError, ValueError):
        pass
    return str(name)


def known_weight_functions():
    return tuple(sorted(_WEIGHT_FUNCTIONS))


def _coerce_node_props(node_props=None, node_tuning=None):
    if node_props is None:
        if node_tuning is None:
            return None
        return pd.DataFrame.from_dict(
            {node_id: {"tuning_angle": value} for node_id, value in node_tuning.items()},
            orient="index",
        )
    if isinstance(node_props, pd.DataFrame):
        if "node_id" in node_props.columns:
            props = node_props.set_index("node_id", drop=False)
        else:
            props = node_props
    else:
        props = pd.DataFrame.from_dict(node_props, orient="index")

    if node_tuning is not None:
        tuning = pd.Series(node_tuning)
        tuning_by_node_id = pd.Series(props.index.map(tuning), index=props.index)
        if "tuning_angle" in props.columns:
            props = props.copy()
            props["tuning_angle"] = props["tuning_angle"].fillna(tuning_by_node_id)
        else:
            props = props.copy()
            props["tuning_angle"] = tuning_by_node_id
    return props


def _node_property(rows, node_props, node_column, property_name):
    if node_props is None or property_name not in node_props.columns:
        return pd.Series(np.nan, index=rows.index)
    return rows[node_column].map(node_props[property_name])


def _numeric_edge_property(rows, property_name, default=None):
    if property_name not in rows.columns:
        return pd.Series(np.nan, index=rows.index)

    values = pd.to_numeric(
        rows[property_name].replace("NULL", np.nan),
        errors="coerce",
    )
    if default is not None:
        values = values.fillna(default)
    return values.astype(float)


def _raise_if_missing(function_name, values_by_name):
    missing = pd.DataFrame(values_by_name).isna()
    if not missing.any().any():
        return

    missing_columns = sorted(missing.columns[missing.any()].tolist())
    bad_indices = missing.index[missing.any(axis=1)].tolist()[:5]
    raise ValueError(
        f"{function_name} requires non-missing {', '.join(missing_columns)} "
        f"for row index(es): {bad_indices}"
    )


def _delta_tuning_180(src_tuning, trg_tuning):
    return np.abs(np.abs((np.abs(trg_tuning - src_tuning) % 360.0) - 180.0) - 180.0)


@register_weight_function("wmax")
def wmax(rows, node_props=None, node_tuning=None):
    syn_weight = _numeric_edge_property(rows, "syn_weight")
    _raise_if_missing("wmax", {"syn_weight": syn_weight})
    return syn_weight


register_weight_function("default_weight_fnc")(wmax)
register_weight_function("set_syn_weight")(wmax)


@register_weight_function("gaussianLL")
def gaussian_ll(rows, node_props=None, node_tuning=None):
    node_props = _coerce_node_props(node_props, node_tuning)

    src_tuning = _node_property(rows, node_props, "pre nid", "tuning_angle")
    trg_tuning = _node_property(rows, node_props, "post nid", "tuning_angle")
    sigma = _numeric_edge_property(rows, "weight_sigma")
    syn_weight = _numeric_edge_property(rows, "syn_weight")

    _raise_if_missing(
        "gaussianLL",
        {
            "source tuning_angle": src_tuning,
            "target tuning_angle": trg_tuning,
            "weight_sigma": sigma,
            "syn_weight": syn_weight,
        },
    )

    src_t = src_tuning.astype(float).values
    trg_t = trg_tuning.astype(float).values
    delta_tuning = np.abs(np.abs(np.abs(180.0 - np.abs(trg_t - src_t) % 360.0) - 90.0) - 90.0)
    weights = syn_weight.values * np.exp(-(delta_tuning / sigma.values) ** 2)
    return pd.Series(weights, index=rows.index)


@register_weight_function("DirectionRule_others")
def direction_rule_others(rows, node_props=None, node_tuning=None):
    node_props = _coerce_node_props(node_props, node_tuning)
    src_tuning = _node_property(rows, node_props, "pre nid", "tuning_angle")
    trg_tuning = _node_property(rows, node_props, "post nid", "tuning_angle")
    sigma = _numeric_edge_property(rows, "weight_sigma")
    syn_weight = _numeric_edge_property(rows, "syn_weight")

    _raise_if_missing(
        "DirectionRule_others",
        {
            "source tuning_angle": src_tuning,
            "target tuning_angle": trg_tuning,
            "weight_sigma": sigma,
            "syn_weight": syn_weight,
        },
    )

    delta_tuning = _delta_tuning_180(src_tuning.astype(float).values, trg_tuning.astype(float).values)
    weights = np.exp(-((delta_tuning / sigma.values) ** 2)) * syn_weight.values
    return pd.Series(weights, index=rows.index)


@register_weight_function("DirectionRule_EE")
def direction_rule_ee(rows, node_props=None, node_tuning=None):
    node_props = _coerce_node_props(node_props, node_tuning)
    src_tuning = _node_property(rows, node_props, "pre nid", "tuning_angle")
    trg_tuning = _node_property(rows, node_props, "post nid", "tuning_angle")
    x_src = _node_property(rows, node_props, "pre nid", "x")
    z_src = _node_property(rows, node_props, "pre nid", "z")
    x_trg = _node_property(rows, node_props, "post nid", "x")
    z_trg = _node_property(rows, node_props, "post nid", "z")
    sigma = _numeric_edge_property(rows, "weight_sigma")
    syn_weight = _numeric_edge_property(rows, "syn_weight")

    _raise_if_missing(
        "DirectionRule_EE",
        {
            "source tuning_angle": src_tuning,
            "target tuning_angle": trg_tuning,
            "source x": x_src,
            "source z": z_src,
            "target x": x_trg,
            "target z": z_trg,
            "weight_sigma": sigma,
            "syn_weight": syn_weight,
        },
    )

    src_t = src_tuning.astype(float).values
    trg_t = trg_tuning.astype(float).values
    delta_tuning = _delta_tuning_180(src_t, trg_t)
    w_multiplier = np.exp(-((delta_tuning / sigma.values) ** 2))

    delta_x = (x_trg.astype(float).values - x_src.astype(float).values) * 0.07
    delta_z = (z_trg.astype(float).values - z_src.astype(float).values) * 0.04
    theta_pref = trg_t * (np.pi / 180.0)
    xz = delta_x * np.cos(theta_pref) + delta_z * np.sin(theta_pref)
    phase_scale_ratio = np.exp(-(xz ** 2 / 2.0))

    theta_tar_scale = np.abs(np.abs(np.abs(180.0 - np.abs(trg_t) % 360.0) - 90.0) - 90.0)
    phase_scale_ratio = phase_scale_ratio * (5.5 / 4.0 - 11.0 / 1680.0 * theta_tar_scale)

    weights = w_multiplier * phase_scale_ratio * syn_weight.values
    return pd.Series(weights, index=rows.index)


def apply_weight_functions(chunk, node_tuning=None, node_props=None, strict=False, enable_direction_rules=False):
    """Apply registered weight functions to a connection chunk.

    In strict mode, unknown names raise ValueError. With strict=False they keep
    the historical bionetlite behavior and pass through syn_weight.

    DirectionRule_* names are registered but disabled by default to preserve the
    existing no-direction Neulite workflow. Set enable_direction_rules=True to
    match BMTK's DirectionRule_* runtime weight functions at build time.
    """
    base_weight = _numeric_edge_property(chunk, "syn_weight")
    _raise_if_missing("weight function dispatch", {"syn_weight": base_weight})
    chunk.loc[:, "weight"] = base_weight
    if "weight_function" not in chunk.columns:
        return chunk

    node_props = _coerce_node_props(node_props, node_tuning)
    function_names = chunk["weight_function"].map(normalize_weight_function_name)
    if not enable_direction_rules:
        function_names = function_names.map(lambda name: "wmax" if name in _DIRECTION_RULE_NAMES else name)

    unknown_names = sorted({name for name in function_names.unique() if name not in _WEIGHT_FUNCTIONS})
    if strict and unknown_names:
        raise ValueError(f"Unknown weight_function(s): {', '.join(unknown_names)}")

    for name, indices in function_names.groupby(function_names).groups.items():
        func = _WEIGHT_FUNCTIONS.get(name)
        if func is None:
            continue
        chunk.loc[indices, "weight"] = func(
            chunk.loc[indices],
            node_props=node_props,
            node_tuning=node_tuning,
        )

    return chunk


class NeuliteNetwork(DenseNetwork):
    # Connection CSV header definition
    CONNECTION_CSV_HEADER = ["#pre nid", "post nid", "post cid", "weight",
                             "tau_decay", "tau_rise", "erev", "delay", "e/i"]
    EDGE_TYPE_COLUMNS = ["edge_type_id", "syn_weight", "target_sections", "delay",
                         "tau2", "tau1", "erev", "weight_function", "weight_sigma"]
    NODE_PROPERTY_COLUMNS = ["node_id", "node_type_id", "tuning_angle", "x", "z"]
    SWC_ID_MAP_DIR = "swc_id_map"
    SWC_ID_MAP_COLUMNS = [
        "morphology",
        "bmtk_swc_id",
        "zero_based_id",
        "converted_cid",
        "type",
        "x",
        "y",
        "z",
        "r",
        "parent",
        "status",
    ]

    def __init__(self, name, convert_morphologies=True, convert_ion_channels=True,
                 simulation_config=None, generate_config_h=True,
                 enable_direction_rules=False,
                 preselect_synapse_locations=False,
                 preselect_synapse_seed=None,
                 write_precise_target_debug_properties=True,
                 **network_props):
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
        self.enable_direction_rules = enable_direction_rules
        self.preselect_synapse_locations = preselect_synapse_locations
        self.preselect_synapse_seed = preselect_synapse_seed
        self.write_precise_target_debug_properties = write_precise_target_debug_properties
        self.skipped_node_groups = []
        self._edge_type_order_by_id = {}
        self._next_edge_type_order = 0

        if simulation_config is not None and generate_config_h:
            def generate_from_simulation_config():
                logger.info("Generating config.h from simulation configuration")
                self.generate_config_h(simulation_config)
            self._run_rank0_task(
                generate_from_simulation_config,
                "config.h generation from simulation configuration",
            )

        config_path = self._get_config_path()

        # Auto-copy BMTK example files if config.json does not exist and base_dir is detected.
        # This enables BMTK tutorials to work without running build_env_bionet first.
        if self.my_base_dir and not os.path.exists(config_path):
            def setup_bmtk_environment():
                logger.info(f"config.json not found. Auto-setting up environment in {self.my_base_dir}")
                from bmtk.utils.sim_setup import build_env_bionet

                # Get tutorial-specific parameters based on base_dir pattern
                tutorial_params = self._get_tutorial_params(self.my_base_dir)

                build_env_kwargs = {
                    'base_dir': self.my_base_dir,
                    'include_examples': True,
                    'config_file': 'config.json',
                    'tstop': tutorial_params['tstop'],
                    'dt': tutorial_params['dt'],
                    'current_clamp': tutorial_params['current_clamp'],
                }
                if 'report_vars' in tutorial_params:
                    build_env_kwargs['report_vars'] = tutorial_params['report_vars']
                build_env_bionet(**build_env_kwargs)

            self._run_rank0_task(setup_bmtk_environment, "BMTK environment setup")

        self._initialize_from_config_path(config_path, simulation_config, generate_config_h)

    def build(self, force=False):
        if mpi_size > 1:
            build_time_uuid()
        if self.preselect_synapse_locations and mpi_size > 1:
            from bmtk.builder.bionet import get_swc as _get_swc
            barrier()
        return super(NeuliteNetwork, self).build(force=force)

    @staticmethod
    def _rank0_failure_message(context, failure):
        exception_type, message, formatted_traceback = failure
        return f"{context} failed on rank 0 ({exception_type}: {message})\n{formatted_traceback}"

    def _run_rank0_task(self, task, context):
        failure = None
        original_exception = None

        if mpi_rank == 0:
            try:
                task()
            except Exception as exc:
                original_exception = exc
                failure = (type(exc).__name__, str(exc), traceback.format_exc())

        if mpi_size > 1:
            if comm is None:
                if original_exception is not None:
                    raise original_exception
                raise RuntimeError(f"MPI communicator is not available for {context}")
            failure = comm.bcast(failure, root=0)

        if failure is not None:
            if original_exception is not None:
                raise original_exception
            raise RuntimeError(self._rank0_failure_message(context, failure))

        barrier()

    def _initialize_from_config_path(self, config_path, simulation_config, generate_config_h):
        if not os.path.exists(config_path):
            logger.warning(f"Initialization warning: config file not found: {config_path}")
            return

        config = SonataConfig.from_json(config_path)
        self.morphologies_dir = config["components"]["morphologies_dir"]
        self.biophysical_neuron_models_dir = config["components"]["biophysical_neuron_models_dir"]
        self._run_rank0_task(
            lambda: self._initialize_neulite_assets(config, simulation_config, generate_config_h),
            "Neulite asset initialization",
        )

    def _initialize_neulite_assets(self, config, simulation_config, generate_config_h):
        if self.convert_morphologies_flag:
            if not self.morphologies_dir or not os.path.exists(self.morphologies_dir):
                raise FileNotFoundError(f"Morphologies directory does not exist. [{self.morphologies_dir}]")
            logger.info("Converting morphology SWC files for Neulite")
            self._convert_morphologies()

        if self.convert_ion_channels_flag and self.biophysical_neuron_models_dir:
            logger.info("Converting ion channel JSON files to CSV format")
            self.convert_ion_channels()

        if simulation_config is None and generate_config_h:
            logger.info("Generating config.h from config.json")
            self.generate_config_h(config)
    
    def add_nodes(self, N=1, **properties):
        """Add nodes to the network. Only biophysical nodes are supported."""
        model_type = properties.get('model_type')
        if model_type != 'biophysical':
            skip_record = {
                "reason": "non-biophysical",
                "N": N,
                "model_type": model_type,
                "pop_name": properties.get("pop_name"),
                "node_type_id": properties.get("node_type_id"),
            }
            self.skipped_node_groups.append(skip_record)
            logger.info(
                "Skipping non-biophysical nodes: "
                f"model_type={model_type}, N={N}, "
                f"pop_name={skip_record['pop_name']}, node_type_id={skip_record['node_type_id']}"
            )
            return
        super(NeuliteNetwork, self).add_nodes(N=N, **properties)

    def add_edges(self, source=None, target=None, connection_rule=None, connection_params=None,
                  iterator='one_to_one', edge_type_properties=None, **properties):
        """Add edges between source and target nodes"""
        connection = super(NeuliteNetwork, self).add_edges(source=source, target=target,
                                                           connection_rule=connection_rule,
                                                           connection_params=connection_params,
                                                           iterator=iterator,
                                                           edge_type_properties=edge_type_properties,
                                                           **properties)
        edge_type_id = connection.edge_type_properties.get("edge_type_id")
        if edge_type_id is not None:
            self._edge_type_order_by_id[int(edge_type_id)] = self._next_edge_type_order
            self._next_edge_type_order += 1
        if self.preselect_synapse_locations:
            self._add_precise_targeting_properties(connection, properties)
        return connection

    def _add_precise_targeting_properties(self, connection, edge_properties):
        if connection.source_network_name != self.name or connection.target_network_name != self.name:
            raise ValueError(
                "Automatic precise targeting currently supports self-network biophysical edges only; "
                f"got {connection.source_network_name} -> {connection.target_network_name}"
            )

        if self.write_precise_target_debug_properties:
            names = ["sec_id", "sec_x", "bmtk_swc_id", "bmtk_swc_dist", "neulite_post_cid"]
            dtypes = [np.uint32, np.float64, np.int64, np.float64, np.int64]
        else:
            names = ["sec_id", "sec_x", "neulite_post_cid"]
            dtypes = [np.uint32, np.float64, np.int64]

        if self._is_missing_value(edge_properties.get("target_sections")):
            rule = self._missing_precise_targeting_rule(
                "target_sections is required when preselect_synapse_locations=True"
            )
            connection.add_properties(names, rule=rule, dtypes=dtypes)
            return

        target_sections = self._normalize_precise_target_sections(edge_properties.get("target_sections"))
        distance_range = self._normalize_precise_distance_range(edge_properties.get("distance_range"))
        edge_type_id = int(connection.edge_type_properties["edge_type_id"])
        rule = self._precise_targeting_rule(edge_type_id, target_sections, distance_range)
        connection.add_properties(names, rule=rule, dtypes=dtypes)

    @staticmethod
    def _missing_precise_targeting_rule(message):
        def rule(source, target):
            raise ValueError(message)

        return rule

    def _precise_targeting_rule(self, edge_type_id, target_sections, distance_range):
        synapse_ordinals = Counter()
        swc_cache = {}
        mapping_cache = {}

        def rule(source, target):
            source_node_id = self._node_id_from_builder_node(source)
            target_node_id = self._node_id_from_builder_node(target)
            ordinal_key = (source_node_id, target_node_id)
            synapse_ordinal = synapse_ordinals[ordinal_key]
            synapse_ordinals[ordinal_key] += 1

            placement = self._preselect_synapse_location(
                source,
                target,
                target_sections=target_sections,
                distance_range=distance_range,
                edge_type_id=edge_type_id,
                synapse_ordinal=synapse_ordinal,
                swc_cache=swc_cache,
                mapping_cache=mapping_cache,
            )
            if self.write_precise_target_debug_properties:
                return placement
            return placement[0], placement[1], placement[4]

        return rule

    @staticmethod
    def _node_id_from_builder_node(node):
        if hasattr(node, "node_id"):
            return int(node.node_id)
        return int(node["node_id"])

    @staticmethod
    def _literal_or_value(value):
        if isinstance(value, str):
            try:
                return ast.literal_eval(value)
            except (SyntaxError, ValueError):
                return value
        return value

    @staticmethod
    def _normalize_precise_target_sections(target_sections):
        if NeuliteNetwork._is_missing_value(target_sections):
            raise ValueError("target_sections is required when preselect_synapse_locations=True")

        target_sections = NeuliteNetwork._literal_or_value(target_sections)
        if isinstance(target_sections, str):
            target_sections = [target_sections]
        if isinstance(target_sections, np.ndarray):
            target_sections = target_sections.tolist()
        if not isinstance(target_sections, (list, tuple)) or len(target_sections) == 0:
            raise ValueError("target_sections must be a non-empty list when preselect_synapse_locations=True")

        section_names = []
        section_name_by_swc_type = {
            1: "somatic",
            2: "axon",
            3: "basal",
            4: "apical",
        }
        for section in target_sections:
            if isinstance(section, (int, np.integer)):
                section = section_name_by_swc_type.get(int(section), section)
            elif isinstance(section, str) and section.strip() in ("1", "2", "3", "4"):
                section = section_name_by_swc_type[int(section.strip())]
            section = str(section)
            if section in ("axon", "axonal", "2"):
                raise ValueError("axon precise targeting is not supported")
            section_names.append(section)

        return tuple(section_names)

    @staticmethod
    def _normalize_precise_distance_range(distance_range):
        if NeuliteNetwork._is_missing_value(distance_range):
            logger.info(
                "distance_range not specified for precise targeting; using all-distance default [0.0, 1.0e20]"
            )
            return (0.0, 1.0e20)

        distance_range = NeuliteNetwork._literal_or_value(distance_range)
        if isinstance(distance_range, np.ndarray):
            distance_range = distance_range.tolist()
        if not isinstance(distance_range, (list, tuple)) or len(distance_range) != 2:
            raise ValueError("distance_range must be a two-value sequence when preselect_synapse_locations=True")

        dmin, dmax = float(distance_range[0]), float(distance_range[1])
        if dmin < 0.0 or dmax <= dmin:
            raise ValueError(
                "distance_range must satisfy 0.0 <= min < max when preselect_synapse_locations=True"
            )
        return (dmin, dmax)

    def _stable_preselect_seed(self, edge_type_id, source_node_id, target_node_id, synapse_ordinal):
        payload = "|".join(
            str(value)
            for value in (
                self.preselect_synapse_seed,
                edge_type_id,
                source_node_id,
                target_node_id,
                synapse_ordinal,
            )
        ).encode("utf-8")
        return int.from_bytes(hashlib.sha256(payload).digest()[:8], "little") % (2 ** 32)

    def _bmtk_swc_for_target(self, target, swc_cache):
        target_node_id = self._node_id_from_builder_node(target)
        if target_node_id not in swc_cache:
            from bmtk.builder.bionet import get_swc
            swc_cache[target_node_id] = get_swc(
                target,
                morphology_dir=self.morphologies_dir,
                use_cache=False,
            )
        return swc_cache[target_node_id]

    @staticmethod
    def _validate_candidate_probabilities(candidate_segment_indices, candidate_probabilities):
        candidate_probabilities = np.asarray(candidate_probabilities, dtype=np.float64)
        if len(candidate_segment_indices) != len(candidate_probabilities):
            raise ValueError("BMTK precise targeting returned mismatched candidate segments and probabilities")
        probability_sum = float(np.sum(candidate_probabilities))
        if probability_sum <= 0.0 or not np.isfinite(probability_sum) or not np.all(np.isfinite(candidate_probabilities)):
            raise ValueError("BMTK precise targeting returned invalid candidate probabilities")
        return candidate_probabilities / probability_sum

    def _preselect_synapse_location(self, source, target, target_sections, distance_range,
                                    edge_type_id, synapse_ordinal, swc_cache=None,
                                    mapping_cache=None):
        if swc_cache is None:
            swc_cache = {}
        if mapping_cache is None:
            mapping_cache = {}

        source_node_id = self._node_id_from_builder_node(source)
        target_node_id = self._node_id_from_builder_node(target)
        swc = self._bmtk_swc_for_target(target, swc_cache)
        candidate_segment_indices, candidate_probabilities = swc.find_sections(target_sections, distance_range)
        if len(candidate_segment_indices) == 0:
            raise ValueError(
                "No candidate target sections found for precise targeting: "
                f"target_node_id={target_node_id}, sections={target_sections}, distance_range={distance_range}"
            )
        candidate_probabilities = self._validate_candidate_probabilities(
            candidate_segment_indices,
            candidate_probabilities,
        )

        seed = self._stable_preselect_seed(edge_type_id, source_node_id, target_node_id, synapse_ordinal)
        rng = np.random.default_rng(seed)
        chosen_segment_index = int(rng.choice(candidate_segment_indices, p=candidate_probabilities))
        sec_id = int(swc.seg_props.sec_id[chosen_segment_index])
        sec_x = float(swc.seg_props.x[chosen_segment_index])
        bmtk_swc_id, bmtk_swc_dist = swc.get_swc_id(sec_id, sec_x)
        neulite_post_cid = self._lookup_neulite_cid_from_bmtk_swc_id(
            target.get("morphology"),
            int(bmtk_swc_id),
            mapping_cache=mapping_cache,
        )
        return (
            np.uint32(sec_id),
            float(sec_x),
            int(bmtk_swc_id),
            float(bmtk_swc_dist),
            int(neulite_post_cid),
        )

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

    def _get_tutorial_params(self, base_dir):
        """Get simulation parameters based on base_dir pattern.

        Automatically detects tutorial type from base_dir and returns
        appropriate simulation parameters. This maintains the "import-only
        change" design philosophy by not requiring users to specify parameters.

        :param base_dir: Base directory path (e.g., 'sim_ch01', './sim_ch04')
        :return: dict with tstop, dt, current_clamp parameters
        """
        if not base_dir:
            logger.info("No base_dir specified. Using bionetlite defaults.")
            return BIONETLITE_DEFAULTS.copy()

        # Normalize path for matching
        normalized = os.path.basename(os.path.normpath(base_dir))

        # Also check full path for patterns
        full_path = os.path.abspath(base_dir) if base_dir else ''

        for pattern, params in TUTORIAL_PARAMS.items():
            if pattern in normalized or pattern in full_path:
                logger.info(
                    f"Detected tutorial pattern '{pattern}' in base_dir '{base_dir}'. "
                    f"Using parameters: tstop={params['tstop']}, dt={params['dt']}, "
                    f"i_amp={params['current_clamp']['amp']}"
                )
                return params.copy()

        logger.info(
            f"No tutorial pattern matched for base_dir '{base_dir}'. "
            f"Using bionetlite defaults: tstop={BIONETLITE_DEFAULTS['tstop']}, "
            f"dt={BIONETLITE_DEFAULTS['dt']}"
        )
        return BIONETLITE_DEFAULTS.copy()

    @staticmethod
    def _complete_network_entries(scanned_entries, required_fields):
        return [
            entry for entry in scanned_entries.values()
            if all(entry.get(field) for field in required_fields)
        ]

    @staticmethod
    def _merge_network_entries(existing_entries, scanned_entries, primary_field):
        if not existing_entries:
            return list(scanned_entries), bool(scanned_entries)

        existing_primary_files = {
            entry.get(primary_field)
            for entry in existing_entries
            if isinstance(entry, dict) and entry.get(primary_field)
        }
        scanned_primary_files = {
            entry[primary_field]
            for entry in scanned_entries
            if entry.get(primary_field)
        }

        if not existing_primary_files.intersection(scanned_primary_files):
            return existing_entries, False

        merged_entries = list(existing_entries)
        updated = False
        for entry in scanned_entries:
            primary_file = entry.get(primary_field)
            if primary_file and primary_file not in existing_primary_files:
                merged_entries.append(entry)
                existing_primary_files.add(primary_file)
                updated = True

        return merged_entries, updated

    def update_circuit_config(self):
        """Update circuit_config.json networks section after network files are created.

        This resolves the timing issue where circuit_config.json is created before
        network files exist, resulting in an empty networks section.

        The method scans the network directory for SONATA files and updates the
        networks section of circuit_config.json with the found files.
        """
        if mpi_rank != 0:
            return

        circuit_config_path = self._get_config_path("circuit_config.json")
        if not os.path.exists(circuit_config_path):
            logger.debug(f"circuit_config.json not found, skipping update: {circuit_config_path}")
            return

        # Read existing circuit_config.json
        with open(circuit_config_path, 'r') as f:
            circuit_config = json.load(f)

        # Get network directory path
        network_dir = os.path.join(self.my_base_dir, 'network') if self.my_base_dir else 'network'

        # Scan network directory for SONATA files
        net_nodes = {}
        net_edges = {}

        if os.path.exists(network_dir):
            for fname in sorted(os.listdir(network_dir)):
                file_path = os.path.join(network_dir, fname)
                if not os.path.isfile(file_path) or fname.startswith('.'):
                    continue

                if fname.endswith('_nodes.h5'):
                    net_name = fname[:fname.find('_nodes')]
                    net_nodes.setdefault(net_name, {})['nodes_file'] = f'$NETWORK_DIR/{fname}'
                elif fname.endswith('_node_types.csv'):
                    net_name = fname[:fname.find('_node_types')]
                    net_nodes.setdefault(net_name, {})['node_types_file'] = f'$NETWORK_DIR/{fname}'
                elif fname.endswith('_edges.h5'):
                    net_name = fname[:fname.find('_edges')]
                    net_edges.setdefault(net_name, {})['edges_file'] = f'$NETWORK_DIR/{fname}'
                elif fname.endswith('_edge_types.csv'):
                    net_name = fname[:fname.find('_edge_types')]
                    net_edges.setdefault(net_name, {})['edge_types_file'] = f'$NETWORK_DIR/{fname}'

        # Check if networks section needs update
        networks = circuit_config.setdefault('networks', {})
        existing_nodes = networks.get('nodes', [])
        existing_edges = networks.get('edges', [])
        scanned_nodes = self._complete_network_entries(
            net_nodes,
            required_fields=("nodes_file", "node_types_file"),
        )
        scanned_edges = self._complete_network_entries(
            net_edges,
            required_fields=("edges_file", "edge_types_file"),
        )

        merged_nodes, updated_nodes = self._merge_network_entries(
            existing_nodes,
            scanned_nodes,
            primary_field="nodes_file",
        )
        merged_edges, updated_edges = self._merge_network_entries(
            existing_edges,
            scanned_edges,
            primary_field="edges_file",
        )

        if not updated_nodes and not updated_edges:
            logger.debug("No missing network entries to update in circuit_config.json")
            return

        networks['nodes'] = merged_nodes
        networks['edges'] = merged_edges
        networks.setdefault('gap_juncs', [])

        # Save updated circuit_config.json
        with open(circuit_config_path, 'w') as f:
            json.dump(circuit_config, f, indent=2)

        logger.info(
            "Updated circuit_config.json with "
            f"{len(networks.get('nodes', []))} node population(s) and "
            f"{len(networks.get('edges', []))} edge population(s)"
        )

    def _convert_morphologies(self):
        # Create output directory if it doesn't exist
        output_dir = os.path.join(self.neulite_dir, self.swc_dir)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        files = glob.glob(os.path.join(self.morphologies_dir, "*.swc"))
        for file in files:
            logger.debug(f"_convert_morphologies: {file}")
            swc = pd.read_csv(file, names=['id', 'type', 'x', 'y', 'z', 'r', 'parent'], comment="#", delimiter=" ", dtype={'id': int, 'parent': int})
            swc["bmtk_swc_id"] = swc["id"]
            swc = self._replace_zero_based_index(swc)
            swc["zero_based_id"] = swc["id"]
            swc["status"] = "mapped"
            removed_axon_map = self._removed_axon_mapping_rows(swc, os.path.basename(file))
            swc = self._fix_axon_perisomatic_directed_df(swc)
            swc = self._sort_by_DFS(swc)
            self._write_swc_id_map(os.path.basename(file), swc, removed_axon_map)
            swc = swc[["id", "type", "x", "y", "z", "r", "parent"]].rename(columns={'id': '#id'})
            output_file = os.path.join(output_dir, os.path.basename(file))
            swc.to_csv(output_file, index=False, sep=" ", float_format='%.6f')
            logger.debug(f"Saved converted morphology to: {output_file}")
        
        logger.info(f"Morphology conversion completed. Converted {len(files)} SWC files to {output_dir}")

    def _swc_id_map_path(self, morphology):
        return os.path.join(
            self.neulite_dir,
            self.SWC_ID_MAP_DIR,
            f"{os.path.basename(morphology)}.csv",
        )

    def _write_swc_id_map(self, morphology, converted_swc, removed_axon_map):
        map_dir = os.path.join(self.neulite_dir, self.SWC_ID_MAP_DIR)
        if not os.path.exists(map_dir):
            os.makedirs(map_dir)

        converted_map = converted_swc.copy()
        converted_map["morphology"] = morphology
        converted_map["converted_cid"] = converted_map["id"]

        mapping = pd.concat([converted_map, removed_axon_map], ignore_index=True, sort=False)
        for column in self.SWC_ID_MAP_COLUMNS:
            if column not in mapping.columns:
                mapping[column] = pd.NA

        mapping = mapping[self.SWC_ID_MAP_COLUMNS]
        mapping.to_csv(self._swc_id_map_path(morphology), index=False)

    def _lookup_neulite_cid_from_bmtk_swc_id(self, morphology, bmtk_swc_id, mapping_cache=None):
        if mapping_cache is None:
            mapping_cache = {}
        if self._is_missing_value(morphology):
            raise ValueError("morphology is required for precise targeting SWC id lookup")

        morphology_basename = os.path.basename(morphology)
        if morphology_basename not in mapping_cache:
            map_path = self._swc_id_map_path(morphology_basename)
            if not os.path.exists(map_path):
                raise FileNotFoundError(f"SWC id map not found for morphology {morphology_basename}: {map_path}")
            mapping_cache[morphology_basename] = pd.read_csv(map_path)

        mapping = mapping_cache[morphology_basename]
        matches = mapping[mapping["bmtk_swc_id"] == bmtk_swc_id]
        if matches.empty:
            raise ValueError(
                f"bmtk_swc_id={bmtk_swc_id} is not present in SWC id map for {morphology_basename}"
            )

        mapped = matches[matches["status"] == "mapped"]
        if mapped.empty:
            statuses = sorted(str(status) for status in matches["status"].dropna().unique())
            raise ValueError(
                f"bmtk_swc_id={bmtk_swc_id} is not mapped for {morphology_basename}; "
                f"status={statuses}"
            )

        converted_cid = mapped.iloc[0]["converted_cid"]
        if self._is_missing_value(converted_cid):
            raise ValueError(
                f"bmtk_swc_id={bmtk_swc_id} has no converted_cid for {morphology_basename}"
            )
        return int(converted_cid)

    @staticmethod
    def _removed_axon_mapping_rows(swc, morphology):
        axon_rows = swc[swc["type"] == 2].copy()
        if axon_rows.empty:
            return pd.DataFrame(columns=NeuliteNetwork.SWC_ID_MAP_COLUMNS)

        axon_rows["morphology"] = morphology
        axon_rows["converted_cid"] = pd.NA
        axon_rows["status"] = "axon_removed"
        return axon_rows[NeuliteNetwork.SWC_ID_MAP_COLUMNS]

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

        base_columns = ["id", "type", "x", "y", "z", "r", "parent"]
        final_data = {
            'id': new_df['new_id'].to_numpy(),
            'type': new_df['type'].to_numpy(),
            'x': new_df['x'].to_numpy(),
            'y': new_df['y'].to_numpy(),
            'z': new_df['z'].to_numpy(),
            'r': new_df['r'].to_numpy(),
            'parent': new_df['new_parent'].to_numpy()
        }
        for column in [c for c in swc_df.columns if c not in base_columns]:
            final_data[column] = new_df[column].to_numpy()
        final_df = pd.DataFrame(final_data)

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
        for column in [c for c in df.columns if c not in new_axon_df.columns]:
            new_axon_df[column] = pd.NA
        if "zero_based_id" in new_axon_df.columns:
            new_axon_df["zero_based_id"] = new_axon_df["id"]
        if "status" in new_axon_df.columns:
            new_axon_df["status"] = "synthetic_axon"

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

        # Update circuit_config.json with node files
        self.update_circuit_config()

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
        node_types_property = self._node_types_properties[model_id]
        if "morphology" not in node_types_property:
            raise ValueError(f"morphology is required for biophysical node_type_id={model_id}")

        model_morph = node_types_property["morphology"]
        if self._is_missing_value(model_morph):
            raise ValueError(f"morphology is required for biophysical node_type_id={model_id}")

        # swc file path
        swc_file = os.path.join(self.morphologies_dir, model_morph)

        if not os.path.exists(swc_file):
            raise FileNotFoundError(f"swc file does not exist: {swc_file}")

        # use converted swc file path
        converted_swc_file = os.path.join(self.neulite_dir, self.swc_dir, model_morph)

        if os.path.exists(converted_swc_file):
            # if converted swc file exists, use it
            data_morph = pd.read_csv(converted_swc_file, delimiter=" ", comment="#", header=None)
            n_comp = data_morph.shape[0]
        else:
            # if converted swc file does not exist, use original swc file (fallback)
            data_morph = pd.read_csv(swc_file, delimiter=" ", comment="#", header=None)
            n_comp = data_morph.shape[0]

        swc_file_for_out = os.path.join(self.swc_dir, model_morph)
        return n_comp, swc_file_for_out

    @staticmethod
    def _is_missing_value(value):
        if value is None:
            return True
        if isinstance(value, str):
            return value.strip() in ("", "NULL", "None", "nan", "NaN")
        try:
            return bool(pd.isna(value))
        except (TypeError, ValueError):
            return False

    @staticmethod
    def _property_or_null(properties, name):
        value = properties.get(name, 'NULL')
        return value if value is not None else 'NULL'

    def _required_edge_type_value(self, edge_type, name):
        value = edge_type.get(name)
        edge_type_id = edge_type.get("edge_type_id", "unknown")
        if self._is_missing_value(value):
            raise ValueError(f"{name} is required for edge_type_id={edge_type_id}")
        return value

    def _edge_type_row(self, edge_type):
        base_params = [
            self._required_edge_type_value(edge_type, column_name)
            for column_name in self.EDGE_TYPE_COLUMNS[:4]
        ]

        tau1, tau2, erev = self._extract_synapse_params(
            edge_type.get('dynamics_params'),
            edge_type_id=edge_type.get('edge_type_id'),
        )
        base_params.extend([tau2, tau1, erev])
        base_params.extend([
            self._property_or_null(edge_type, 'weight_function'),
            self._property_or_null(edge_type, 'weight_sigma'),
        ])
        return base_params

    def _edge_types_dataframe(self, src_network, trg_network, edge_type_ids=None):
        edge_type_id_set = None
        if edge_type_ids is not None:
            edge_type_id_set = {int(edge_type_id) for edge_type_id in edge_type_ids}

        matching_et = [
            c.edge_type_properties for c in self._connection_maps
            if c.source_network_name == src_network and c.target_network_name == trg_network
            and (
                edge_type_id_set is None
                or int(c.edge_type_properties.get("edge_type_id")) in edge_type_id_set
            )
        ]
        edge_types_data = pd.DataFrame(
            [self._edge_type_row(edge_type) for edge_type in matching_et],
            columns=self.EDGE_TYPE_COLUMNS,
        )
        return self._clamp_edge_type_delays(edge_types_data)

    @staticmethod
    def _edge_type_ids_from_collator(merged_edges):
        edge_type_ids = set()
        for chunk_id, idx_beg, idx_end in merged_edges.itr_chunks():
            edge_type_ids.update(int(edge_type_id) for edge_type_id in merged_edges.get_edge_type_ids(chunk_id))
        return edge_type_ids

    @staticmethod
    def _clamp_edge_type_delays(edge_types_data):
        edge_types_data = edge_types_data.copy()
        if edge_types_data.empty:
            return edge_types_data

        non_null_mask = edge_types_data["delay"] != 'NULL'
        if not non_null_mask.any():
            return edge_types_data

        # Neulite 1.1.0 supports floating-point synaptic delays, but the
        # kernel still defines 1.0 ms as the minimum supported delay.
        numeric_delays = pd.to_numeric(edge_types_data.loc[non_null_mask, "delay"])
        clamped_delays = numeric_delays.clip(lower=1.0)
        clamp_mask = numeric_delays < 1.0

        for edge_type_id, delay in edge_types_data.loc[
            numeric_delays[clamp_mask].index,
            ["edge_type_id", "delay"],
        ].itertuples(index=False):
            logger.warning(
                f"Clamping edge delay below Neulite minimum: edge_type_id={edge_type_id}, delay={delay} -> 1.0"
            )

        edge_types_data.loc[non_null_mask, "delay"] = clamped_delays
        return edge_types_data

    @staticmethod
    def _node_property_or_nan(node, property_name):
        try:
            return node[property_name]
        except (KeyError, TypeError):
            return np.nan

    def _node_properties_dataframe(self):
        node_tuning = {}
        node_property_rows = []

        for node in self.nodes():
            node_id = int(node.node_id)
            node_props = {
                "node_id": node_id,
                "node_type_id": int(node.node_type_id),
            }
            for prop_name in ("tuning_angle", "x", "z"):
                node_props[prop_name] = self._node_property_or_nan(node, prop_name)
            try:
                if not pd.isna(node_props["tuning_angle"]):
                    node_tuning[node_id] = node_props["tuning_angle"]
            except (TypeError, ValueError):
                pass
            node_property_rows.append(node_props)

        node_properties = pd.DataFrame(
            node_property_rows,
            columns=self.NODE_PROPERTY_COLUMNS,
        )
        return node_properties, node_tuning

    @staticmethod
    def _write_edges_h5_from_collator(merged_edges, edges_file_name, src_network, trg_network, pop_name,
                                      index_by=('target_node_id', 'source_node_id'), compression='gzip',
                                      sort_by='target_node_id', sort_on_disk=False, edges_file_name_final=None):
        edges_file_name_final = edges_file_name if edges_file_name_final is None else edges_file_name_final

        with h5py.File(edges_file_name, 'w') as hf:
            add_hdf5_attrs(hf)
            pop_grp = hf.create_group('/edges/{}'.format(pop_name))

            n_total_conns = merged_edges.n_total_edges
            pop_grp.create_dataset('source_node_id', (n_total_conns,), dtype='uint64', compression=compression)
            pop_grp['source_node_id'].attrs['node_population'] = src_network
            pop_grp.create_dataset('target_node_id', (n_total_conns,), dtype='uint64', compression=compression)
            pop_grp['target_node_id'].attrs['node_population'] = trg_network
            pop_grp.create_dataset('edge_group_id', (n_total_conns,), dtype='uint16', compression=compression)
            pop_grp.create_dataset('edge_group_index', (n_total_conns,), dtype='uint32', compression=compression)
            pop_grp.create_dataset('edge_type_id', (n_total_conns,), dtype='uint32', compression=compression)

            for group_id in merged_edges.group_ids:
                model_grp = pop_grp.create_group(str(group_id))
                for prop_mdata in merged_edges.get_group_metadata(group_id):
                    model_grp.create_dataset(
                        prop_mdata['name'],
                        shape=prop_mdata['dim'],
                        dtype=prop_mdata['type'],
                        compression=compression,
                    )

            for chunk_id, idx_beg, idx_end in merged_edges.itr_chunks():
                pop_grp['source_node_id'][idx_beg:idx_end] = merged_edges.get_source_node_ids(chunk_id)
                pop_grp['target_node_id'][idx_beg:idx_end] = merged_edges.get_target_node_ids(chunk_id)
                pop_grp['edge_type_id'][idx_beg:idx_end] = merged_edges.get_edge_type_ids(chunk_id)
                pop_grp['edge_group_id'][idx_beg:idx_end] = merged_edges.get_edge_group_ids(chunk_id)
                pop_grp['edge_group_index'][idx_beg:idx_end] = merged_edges.get_edge_group_indices(chunk_id)

                for group_id, prop_name, grp_idx_beg, grp_idx_end in merged_edges.get_group_data(chunk_id):
                    prop_array = merged_edges.get_group_property(prop_name, group_id, chunk_id)
                    pop_grp[str(group_id)][prop_name][grp_idx_beg:grp_idx_end] = prop_array

        if sort_on_disk:
            logger.debug('Sorting {} by {} to {}'.format(edges_file_name, sort_by, edges_file_name_final))
            sort_edges(
                input_edges_path=edges_file_name,
                output_edges_path=edges_file_name_final,
                edges_population='/edges/{}'.format(pop_name),
                sort_by=sort_by,
                compression=compression,
            )
            try:
                logger.debug('Deleting intermediate edges file {}.'.format(edges_file_name))
                os.remove(edges_file_name)
            except OSError as e:
                logger.warning('Unable to remove intermediate edges file {}'.format(edges_file_name))

        if index_by:
            index_by = index_by if isinstance(index_by, (list, tuple)) else [index_by]
            for index_type in index_by:
                logger.debug('Creating index {}'.format(index_type))
                create_index_in_memory(
                    edges_file=edges_file_name_final,
                    edges_population='/edges/{}'.format(pop_name),
                    index_type=index_type,
                    compression=compression,
                )

    @staticmethod
    def _connection_chunk_nsyns(merged_edges, chunk_id, n_edges):
        nsyns = np.ones(n_edges, dtype=np.uint32)
        if n_edges == 0:
            return nsyns

        edge_group_ids = np.asarray(merged_edges.get_edge_group_ids(chunk_id))
        edge_group_indices = np.asarray(merged_edges.get_edge_group_indices(chunk_id))

        for group_id in np.unique(edge_group_ids):
            group_mask = edge_group_ids == group_id
            if not np.any(group_mask):
                continue

            try:
                group_nsyns = np.asarray(merged_edges.get_group_property('nsyns', int(group_id), chunk_id))
            except KeyError:
                # Groups with per-synapse properties are already represented as one row per synapse.
                continue

            group_indices = edge_group_indices[group_mask]
            local_indices = (group_indices - group_indices.min()).astype(np.int64)
            nsyns[group_mask] = group_nsyns[local_indices]

        return nsyns

    @staticmethod
    def _connection_chunk_group_property(merged_edges, chunk_id, prop_name, n_edges):
        values = np.full(n_edges, np.nan, dtype=object)
        if n_edges == 0:
            return values

        edge_group_ids = np.asarray(merged_edges.get_edge_group_ids(chunk_id))
        edge_group_indices = np.asarray(merged_edges.get_edge_group_indices(chunk_id))

        for group_id in np.unique(edge_group_ids):
            group_mask = edge_group_ids == group_id
            if not np.any(group_mask):
                continue

            try:
                group_values = np.asarray(merged_edges.get_group_property(prop_name, int(group_id), chunk_id))
            except KeyError:
                continue

            group_indices = edge_group_indices[group_mask]
            local_indices = (group_indices - group_indices.min()).astype(np.int64)
            values[group_mask] = group_values[local_indices]

        return values

    @staticmethod
    def _connection_chunk_dataframe(merged_edges, chunk_id, edge_types_data):
        source_ids = merged_edges.get_source_node_ids(chunk_id)
        target_ids = merged_edges.get_target_node_ids(chunk_id)
        edge_type_ids = merged_edges.get_edge_type_ids(chunk_id)

        logger.info("Loading edge data for CSV...")
        edge_type_id_df = pd.DataFrame(edge_type_ids, columns=["edge_type_id"])
        source_node_id_df = pd.DataFrame(source_ids, columns=["pre nid"])
        target_node_id_df = pd.DataFrame(target_ids, columns=["post nid"])

        nsyns_array = NeuliteNetwork._connection_chunk_nsyns(merged_edges, chunk_id, len(source_ids))
        nsyns_data = pd.DataFrame(nsyns_array, columns=["nsyns"], dtype="object")

        connections = pd.merge(edge_type_id_df, edge_types_data, on='edge_type_id', how='left')
        connections = pd.concat([connections, nsyns_data], axis=1)
        legacy_post_cids = NeuliteNetwork._connection_chunk_group_property(
            merged_edges,
            chunk_id,
            "post_cid",
            len(source_ids),
        )
        if not pd.Series(legacy_post_cids).map(NeuliteNetwork._is_missing_value).all():
            raise ValueError("Edge property post_cid is not supported; use neulite_post_cid")

        preselected_post_cids = NeuliteNetwork._connection_chunk_group_property(
            merged_edges,
            chunk_id,
            "neulite_post_cid",
            len(source_ids),
        )
        missing_mask = pd.Series(preselected_post_cids).map(NeuliteNetwork._is_missing_value)
        if not missing_mask.all():
            if missing_mask.any():
                raise ValueError(
                    "Preselected post cid group property coverage is partial; "
                    f"neulite_post_cid is missing for {int(missing_mask.sum())} edge row(s)"
                )
            connections["post cid"] = pd.Series(preselected_post_cids, dtype="object")
        connections.insert(0, "pre nid", source_node_id_df)
        connections.insert(1, "post nid", target_node_id_df)
        return connections

    def _prepare_connection_dataframe(self, connections):
        connections = connections.copy()
        connections["target_sections"] = connections["target_sections"].apply(self._str_to_list)
        connections["target_sections"] = connections["target_sections"].apply(self._replace_strings)
        connections["ei"] = ""
        connections["nsyns"] = connections["nsyns"].apply(self._to_list)
        return connections.explode("nsyns").reset_index(drop=True)

    def _connections_dataframe_from_collator(self, merged_edges, edge_types_data):
        connections_list = [
            self._connection_chunk_dataframe(merged_edges, chunk_id, edge_types_data)
            for chunk_id, idx_beg, idx_end in merged_edges.itr_chunks()
        ]

        logger.debug(f"before concat: {len(connections_list)=}")
        connections = pd.concat(connections_list)
        logger.debug(f"after concat: {len(connections)=}")

        connections = self._prepare_connection_dataframe(connections)
        logger.info(f"Finished merging edge types data.")
        return connections

    def _assign_post_cids_before_distribution(self, connections, h5_node_type):
        file_cache = {}
        converted_swc_path = os.path.join(self.neulite_dir, self.swc_dir)
        NeuliteNetwork._validate_preassigned_post_cids(connections)
        if "post cid" in connections.columns:
            return connections

        connections = self._post_cid_sampling_order(connections)
        connections.loc[:, "post cid"] = self._assign_post_cids(
            connections,
            file_cache=file_cache,
            h5_node_type=h5_node_type,
            node_types=self._node_types_properties,
            morphologies_path=self.morphologies_dir,
            converted_morphologies_path=converted_swc_path,
        )
        return connections

    def _post_cid_sampling_order(self, connections):
        connections = connections.copy()
        if mpi_size <= 1 or connections.empty or "edge_type_id" not in connections.columns:
            return connections

        def edge_type_order(edge_type_id):
            edge_type_id = int(edge_type_id)
            return self._edge_type_order_by_id.get(edge_type_id, edge_type_id)

        ordered = connections.copy()
        ordered["_bionetlite_add_order"] = ordered["edge_type_id"].map(edge_type_order)
        ordered["_bionetlite_original_order"] = np.arange(len(ordered))
        ordered = ordered.sort_values(
            ["_bionetlite_add_order", "pre nid", "post nid", "_bionetlite_original_order"],
            kind="mergesort",
        )

        # Match BMTK's single-rank target_node_id sort so post-cid sampling
        # keeps the historical single-rank random sequence under MPI.
        target_order = np.argsort(ordered["post nid"].to_numpy(), kind="quicksort")
        ordered = ordered.iloc[target_order]
        return ordered.drop(columns=["_bionetlite_add_order", "_bionetlite_original_order"])

    def _serialized_connection_chunks(self, connections, node_tuning):
        chunk_size = (len(connections) // mpi_size) + 1
        chunks = [connections.iloc[i*chunk_size:(i+1)*chunk_size] for i in range(mpi_size)]
        logger.debug(f"{mpi_size=}, {chunk_size=}, {len(chunks)=}")
        return [self._serialize_data((chunk, node_tuning)) for chunk in chunks]

    def _write_neulite_connection_csv(self, processed_chunks, src_network, trg_network, sort_output=True):
        valid_chunks = [chunk for chunk in processed_chunks if chunk is not None]

        if not valid_chunks:
            return None

        output_data = pd.concat(valid_chunks)
        if sort_output:
            output_data = output_data.sort_values(["#pre nid", "post nid", "post cid"])

        output_path = os.path.join(self.neulite_dir, f"{src_network}_{trg_network}_connection.csv")
        output_data.to_csv(output_path, mode='w', header=True, index=False)
        return output_path

    def _write_processed_connection_chunks(self, processed_chunks, src_network, trg_network, sort_output=True):
        output_path = self._write_neulite_connection_csv(
            processed_chunks,
            src_network,
            trg_network,
            sort_output=sort_output,
        )

        if output_path is None:
            logger.warning(f"No valid edges to save for {src_network} -> {trg_network}")
            return None

        logger.info(f"Created Neulite connection file: {output_path}")
        return output_path

    @staticmethod
    def _source_ei_series(chunk, node_type_lookup, node_types_data):
        source_node_type_ids = chunk["pre nid"].map(node_type_lookup)
        missing_source_nodes = chunk.loc[source_node_type_ids.isna(), "pre nid"].drop_duplicates().tolist()[:5]
        if missing_source_nodes:
            raise ValueError(f"source node_id not found in local biophysical nodes: {missing_source_nodes}")

        source_node_type_ids = source_node_type_ids.astype(int)
        missing_node_types = sorted(set(source_node_type_ids) - set(node_types_data))[:5]
        if missing_node_types:
            raise ValueError(f"source node_type_id(s) not found: {missing_node_types}")

        node_type_ei = {
            int(node_type_id): props.get("ei")
            for node_type_id, props in node_types_data.items()
        }
        source_ei = source_node_type_ids.map(node_type_ei)
        missing_ei_mask = source_ei.map(NeuliteNetwork._is_missing_value)
        missing_ei = source_node_type_ids[missing_ei_mask].drop_duplicates().tolist()[:5]
        if missing_ei:
            raise ValueError(f"ei is required for source node_type_id(s): {missing_ei}")

        return source_ei

    def _distribute_connection_chunk(self, serialized_chunks):
        if mpi_size == 1:
            (serialized_chunk,) = serialized_chunks
            return self._deserialize_data(serialized_chunk)

        if comm is None:
            raise RuntimeError("MPI communicator is not available")

        if mpi_rank == 0:
            for i, serialized_chunk in enumerate(serialized_chunks):
                if i == 0:
                    chunk, node_tuning = self._deserialize_data(serialized_chunk)
                else:
                    comm.send(serialized_chunk, dest=i)
            return chunk, node_tuning

        serialized_chunk = comm.recv(source=0)
        logger.debug(f"chunk is received:rank{mpi_rank}")
        return self._deserialize_data(serialized_chunk)

    @staticmethod
    def _broadcast_node_type_dataframe(h5_node_type):
        if mpi_size == 1:
            return pd.DataFrame.from_dict(h5_node_type)

        if comm is None:
            raise RuntimeError("MPI communicator is not available")

        try:
            h5_node_type = comm.bcast(h5_node_type, root=0)
            return pd.DataFrame.from_dict(h5_node_type)
        except Exception as e:
            logger.error(f"Error during bcast: {e}")
            comm.Abort()
            raise

    def _collect_processed_chunks(self, processed_chunk):
        if mpi_size == 1:
            return [processed_chunk]

        if comm is None:
            raise RuntimeError("MPI communicator is not available")

        if mpi_rank == 0:
            processed_list = [processed_chunk]
            for i in range(1, mpi_size):
                serialized_chunk = comm.recv(source=i)
                processed_chunk = self._deserialize_data(serialized_chunk)
                logger.debug(f"processed_chunk is received:rank{mpi_rank}")
                processed_list.append(processed_chunk)
            return processed_list

        comm.send(self._serialize_data(processed_chunk), dest=0)
        return None

    def _save_edges(self, edges_file_name, src_network, trg_network, pop_name=None, sort_by='target_node_id',
                    index_by=('target_node_id', 'source_node_id'), compression='gzip',
                    sort_on_disk=False, mode='w', **opts):
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

        if src_network != self.name:
            if mpi_rank == 0:
                logger.info(f"Skipping Neulite edge output for external source network: {src_network} -> {trg_network}")
                logger.info("_save_edges end")
            barrier()
            self.update_circuit_config()
            barrier()
            return

        # Try to sort before writing file, If edges are split across ranks/files for MPI/size issues then we need to
        # write to disk first then sort the hdf5 file
        edges_file_name_final = edges_file_name
        if sort_by:
            if merged_edges.can_sort and not sort_on_disk:
                merged_edges.sort(sort_by=sort_by)
            else:
                sort_on_disk = True
                edges_file_name_final = edges_file_name
                edges_file_basename = os.path.basename(edges_file_name)
                edges_file_dirname = os.path.dirname(edges_file_name)
                edges_file_name = os.path.join(edges_file_dirname, '.unsorted.{}'.format(edges_file_basename))
                if mpi_rank == 0:
                    logger.debug('Unable to sort edges in memory, will temporarly save to {}'.format(edges_file_name) +
                                 ' before sorting hdf5 file.')
        barrier()

        start_time = time.time()
        if mpi_rank == 0:
            pop_name = '{}_to_{}'.format(src_network, trg_network) if pop_name is None else pop_name

            logger.debug('Saving {} edges to disk'.format(n_total_conns))
            self._write_edges_h5_from_collator(
                merged_edges,
                edges_file_name,
                src_network,
                trg_network,
                pop_name,
                index_by=index_by,
                compression=compression,
                sort_by=sort_by,
                sort_on_disk=sort_on_disk,
                edges_file_name_final=edges_file_name_final,
            )

            # neulite start
            start_time = time.time()
            logger.info("Starting the processing...")

            logger.info("Loading edge types data...")
            edge_type_ids = self._edge_type_ids_from_collator(merged_edges)
            edge_types_data = self._edge_types_dataframe(
                src_network,
                trg_network,
                edge_type_ids=edge_type_ids,
            )
            logger.debug(f"{len(edge_types_data)=}")

            h5_node_type, node_tuning = self._node_properties_dataframe()

            merged_connections = self._connections_dataframe_from_collator(merged_edges, edge_types_data)
            merged_connections = self._assign_post_cids_before_distribution(merged_connections, h5_node_type)
            # for MPI
            h5_node_type = h5_node_type.to_dict()
            # for MPI - include node_tuning with each chunk to avoid extra bcast
            serialized_chunks = self._serialized_connection_chunks(merged_connections, node_tuning)
            logger.info("Processing chunks...")
            # neulite end

        else:
            h5_node_type = None
            serialized_chunks = None

        barrier()
        del merged_edges

        # neulite start
        # Use "send" to prevent overflow.
        chunk, node_tuning = self._distribute_connection_chunk(serialized_chunks)
        # neulite end

        h5_node_type = self._broadcast_node_type_dataframe(h5_node_type)

        barrier()

        # neulite start
        file_cache = {}
        # Calculate path to converted (preprocessed) SWC files
        converted_swc_path = os.path.join(self.neulite_dir, self.swc_dir)
        processed_chunk = self._process_chunk(
            chunk,
            file_cache,
            h5_node_type,
            self._node_types_properties,
            self.morphologies_dir,
            converted_morphologies_path=converted_swc_path,
            node_tuning=node_tuning,
            enable_direction_rules=self.enable_direction_rules,
        )
        end_time = time.time()
        logger.info(f"Processing completed in {end_time - start_time} seconds: from rank{mpi_rank}")

        barrier()

        # Collect processed chunks from all ranks
        processed_list = self._collect_processed_chunks(processed_chunk)
        # neulite end

        barrier()

        if mpi_rank == 0:
            self._write_processed_connection_chunks(
                processed_list,
                src_network,
                trg_network,
                sort_output=bool(sort_by),
            )

            logger.info("_save_edges end")
            logger.debug('Saving completed.')

        # Update circuit_config.json with network files after edges are created
        # This resolves the timing issue where circuit_config.json is created
        # before network files exist during module import
        barrier()
        self.update_circuit_config()
        barrier()

    def _extract_synapse_params(self, dynamics_params, edge_type_id=None):
        """Extract tau1, tau2, and erev from dynamics_params
        
        :param dynamics_params: JSON file path or dictionary containing synapse parameters
        :param edge_type_id: edge_type_id used only for error messages
        :return: tuple (tau1, tau2, erev)
        """
        edge_context = f" for edge_type_id={edge_type_id}" if edge_type_id is not None else ""
        if self._is_missing_value(dynamics_params):
            raise ValueError(f"dynamics_params is required{edge_context}")
        
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
            elif isinstance(dynamics_params, dict):
                # if dictionary
                params = dynamics_params
            else:
                raise TypeError(
                    f"dynamics_params must be a JSON path or dict{edge_context}; "
                    f"got {type(dynamics_params).__name__}"
                )

            missing_params = [
                name for name in ("tau1", "tau2", "erev")
                if name not in params or self._is_missing_value(params.get(name))
            ]
            if missing_params:
                raise ValueError(
                    f"dynamics_params is missing required synapse parameter(s){edge_context}: "
                    f"{', '.join(missing_params)}"
                )
            
            tau1 = params['tau1']
            tau2 = params['tau2']
            erev = params['erev']
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
    def _read_morphology_for_cid(morph_path, file_cache, morphologies_path, converted_morphologies_path=None):
        if morph_path in file_cache:
            return file_cache[morph_path]

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
            return morph_data
        except Exception as e:
            logger.error(f"Error reading file {morph_path_full}: {e}")
            raise

    @staticmethod
    def _node_type_lookup(h5_node_type):
        if "node_id" in h5_node_type.columns:
            return pd.Series(h5_node_type["node_type_id"].values, index=h5_node_type["node_id"].values)
        return h5_node_type["node_type_id"]

    @staticmethod
    def _get_post_node_type(post_node_id, h5_node_type=None, node_type_lookup=None):
        if node_type_lookup is None:
            node_type_lookup = NeuliteNetwork._node_type_lookup(h5_node_type)
        return node_type_lookup.loc[post_node_id]

    @staticmethod
    def _section_cache_key(sections):
        if isinstance(sections, np.ndarray):
            return tuple(sections.tolist())
        if isinstance(sections, list):
            return tuple(sections)
        if isinstance(sections, tuple):
            return sections
        return sections

    @staticmethod
    def _post_cid_candidates(post_node_type, sections, file_cache, node_types, morphologies_path, converted_morphologies_path=None):
        if post_node_type not in node_types:
            raise ValueError(f"node_type not found: {post_node_type}")
        if node_types[post_node_type].get("morphology") is None:
            raise ValueError(f"morphology is not set for node_type: {post_node_type}")

        morph_path = node_types[post_node_type].get("morphology")
        morph_data = NeuliteNetwork._read_morphology_for_cid(
            morph_path,
            file_cache,
            morphologies_path,
            converted_morphologies_path=converted_morphologies_path,
        )
        return morph_data.loc[morph_data["type"].isin(sections), "id"].to_numpy()

    @staticmethod
    def _sample_post_cid(candidate_ids):
        return np.random.choice(candidate_ids, size=1, replace=False)[0]

    @staticmethod
    def _assign_post_cids(chunk, file_cache, h5_node_type, node_types, morphologies_path,
                          converted_morphologies_path=None, node_type_lookup=None):
        if node_type_lookup is None:
            node_type_lookup = NeuliteNetwork._node_type_lookup(h5_node_type)

        candidates_cache = {}
        post_cids = []

        for post_node_id, sections in zip(chunk["post nid"].to_numpy(), chunk["target_sections"].to_numpy()):
            try:
                post_node_type = NeuliteNetwork._get_post_node_type(
                    post_node_id,
                    node_type_lookup=node_type_lookup,
                )
                section_key = NeuliteNetwork._section_cache_key(sections)
                cache_key = (post_node_type, section_key)
                if cache_key not in candidates_cache:
                    candidates_cache[cache_key] = NeuliteNetwork._post_cid_candidates(
                        post_node_type,
                        section_key,
                        file_cache,
                        node_types,
                        morphologies_path,
                        converted_morphologies_path=converted_morphologies_path,
                    )
                post_cids.append(NeuliteNetwork._sample_post_cid(candidates_cache[cache_key]))
            except Exception as e:
                logger.error(f"Error _get_cid: {e}")
                logger.error({"post nid": post_node_id, "target_sections": sections})
                raise

        return pd.Series(post_cids, index=chunk.index)

    @staticmethod
    def _has_complete_post_cids(chunk):
        return "post cid" in chunk.columns and not chunk["post cid"].map(NeuliteNetwork._is_missing_value).any()

    @staticmethod
    def _validate_preassigned_post_cids(chunk):
        if "post cid" not in chunk.columns:
            return

        missing_mask = chunk["post cid"].map(NeuliteNetwork._is_missing_value)
        if missing_mask.any():
            raise ValueError(
                "post cid column has partial missing values; "
                f"{int(missing_mask.sum())} row(s) are missing preassigned post cid"
            )

    @staticmethod
    def _validate_post_cid_values(chunk, file_cache, h5_node_type, node_types, morphologies_path,
                                  converted_morphologies_path=None, node_type_lookup=None):
        if "post cid" not in chunk.columns:
            return
        if node_type_lookup is None:
            node_type_lookup = NeuliteNetwork._node_type_lookup(h5_node_type)

        morphology_ids_cache = {}
        invalid_rows = []
        for row_index, (post_node_id, post_cid) in enumerate(
            zip(chunk["post nid"].to_numpy(), chunk["post cid"].to_numpy())
        ):
            post_node_type = NeuliteNetwork._get_post_node_type(
                post_node_id,
                node_type_lookup=node_type_lookup,
            )
            morph_path = node_types[post_node_type].get("morphology")
            if morph_path not in morphology_ids_cache:
                morph_data = NeuliteNetwork._read_morphology_for_cid(
                    morph_path,
                    file_cache,
                    morphologies_path,
                    converted_morphologies_path=converted_morphologies_path,
                )
                morphology_ids_cache[morph_path] = set(morph_data["id"].astype(int).tolist())

            try:
                post_cid_float = float(post_cid)
                post_cid_int = int(post_cid_float)
            except (TypeError, ValueError):
                invalid_rows.append((row_index, post_node_id, post_cid, morph_path))
                continue

            if post_cid_float != post_cid_int or post_cid_int not in morphology_ids_cache[morph_path]:
                invalid_rows.append((row_index, post_node_id, post_cid, morph_path))

        if invalid_rows:
            row_index, post_node_id, post_cid, morph_path = invalid_rows[0]
            raise ValueError(
                "preassigned post cid is outside target morphology compartments; "
                f"row={row_index}, post_nid={post_node_id}, post_cid={post_cid}, morphology={morph_path}"
            )

    @staticmethod
    def _get_cid(row, file_cache, h5_node_type, node_types, morphologies_path, converted_morphologies_path=None):
        try:
            post_node_id = row["post nid"]
            post_node_type = NeuliteNetwork._get_post_node_type(post_node_id, h5_node_type=h5_node_type)
            sections = row["target_sections"]
            candidates = NeuliteNetwork._post_cid_candidates(
                post_node_type,
                sections,
                file_cache,
                node_types,
                morphologies_path,
                converted_morphologies_path=converted_morphologies_path,
            )
            return NeuliteNetwork._sample_post_cid(candidates)
        except Exception as e:
            logger.error(f"Error _get_cid: {e}")
            logger.error(row)
            raise

    def _is_biophysical_source(self, node_id, h5_node_type_df):
        """Check if source node_id belongs to this network (and thus is biophysical).

        NeuliteNetwork only accepts biophysical nodes, so if the node_id exists
        in h5_node_type (this network's nodes), it's biophysical.
        External network nodes (like virtual lgn) won't be in h5_node_type.
        """
        try:
            # h5_node_type_df has 'node_id' column
            return node_id in h5_node_type_df['node_id'].values
        except Exception:
            return False

    @staticmethod
    def _process_chunk(chunk, file_cache, h5_node_type, node_types_data, morphologies_path,
                       converted_morphologies_path=None, node_tuning=None,
                       enable_direction_rules=False):
        logger.debug(f"Start processing chunk...")
        chunk = chunk.copy()
        logger.debug(f"End coppy chunk...")

        # Build node_id -> node_type_id mapping for vectorized lookups
        nid_to_ntid = pd.Series(h5_node_type["node_type_id"].values, index=h5_node_type["node_id"].values)

        has_preassigned_post_cids = "post cid" in chunk.columns
        NeuliteNetwork._validate_preassigned_post_cids(chunk)

        if not has_preassigned_post_cids:
            chunk.loc[:, "post cid"] = NeuliteNetwork._assign_post_cids(
                chunk,
                file_cache=file_cache,
                h5_node_type=h5_node_type,
                node_types=node_types_data,
                morphologies_path=morphologies_path,
                converted_morphologies_path=converted_morphologies_path,
                node_type_lookup=nid_to_ntid,
            )
        else:
            NeuliteNetwork._validate_post_cid_values(
                chunk,
                file_cache=file_cache,
                h5_node_type=h5_node_type,
                node_types=node_types_data,
                morphologies_path=morphologies_path,
                converted_morphologies_path=converted_morphologies_path,
                node_type_lookup=nid_to_ntid,
            )
        logger.debug(f"End _get_cid chunk...")
        # Vectorized ei lookup: pre nid -> node_type_id -> ei value
        chunk.loc[:, "ei"] = NeuliteNetwork._source_ei_series(chunk, nid_to_ntid, node_types_data)
        logger.debug(f"End _get_ei chunk...")

        # DirectionRule_* remains no-direction by default; truly unknown names fail fast.
        chunk = apply_weight_functions(
            chunk,
            node_tuning=node_tuning,
            node_props=h5_node_type,
            strict=True,
            enable_direction_rules=enable_direction_rules,
        )

        # Filter edges to keep only biophysical source connections (vectorized)
        # Note: Target is always biophysical since NeuliteNetwork only accepts biophysical nodes
        bio_ntids = {ntid for ntid, props in node_types_data.items() if props.get("model_type") == "biophysical"}
        pre_biophysical_mask = chunk["pre nid"].map(nid_to_ntid).isin(bio_ntids)

        filtered_count = (~pre_biophysical_mask).sum()
        if filtered_count > 0:
            logger.info(f"Filtering out {filtered_count} edges from non-biophysical source nodes")
            chunk = chunk[pre_biophysical_mask]

        if len(chunk) == 0:
            logger.info("All edges were filtered out")
            return None

        # Drop columns not needed in output, including syn_weight (replaced by weight)
        cols_to_drop = ["nsyns", "target_sections", "syn_weight", "weight_function", "weight_sigma"]
        cols_to_drop = [c for c in cols_to_drop if c in chunk.columns]
        chunk.drop(cols_to_drop, axis=1, inplace=True)

        chunk = chunk.rename(columns={"pre nid": "#pre nid", "ei": "e/i", "tau2": "tau_decay", "tau1": "tau_rise"})
        chunk = chunk[NeuliteNetwork.CONNECTION_CSV_HEADER]
        logger.debug(f"Finished processing chunk.")
        return chunk

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
            "#define ITR_BIT_WIDTH ( 8U )",
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

    def __init__(self, name, adaptor_cls=NeuliteNetwork, convert_morphologies=True, convert_ion_channels=True,
                 simulation_config=None, generate_config_h=True,
                 enable_direction_rules=False,
                 preselect_synapse_locations=False,
                 preselect_synapse_seed=None,
                 write_precise_target_debug_properties=True,
                 **network_props):
        self.adaptor = adaptor_cls(name, convert_morphologies=convert_morphologies, convert_ion_channels=convert_ion_channels,
                                   simulation_config=simulation_config,
                                   generate_config_h=generate_config_h,
                                   enable_direction_rules=enable_direction_rules,
                                   preselect_synapse_locations=preselect_synapse_locations,
                                   preselect_synapse_seed=preselect_synapse_seed,
                                   write_precise_target_debug_properties=write_precise_target_debug_properties,
                                   **network_props)

    def save_edges(self, edges_file_name=None, edge_types_file_name=None, output_dir='.', src_network=None,
                   trg_network=None, name=None, force_build=True, force_overwrite=False,
                   compression='gzip', **opts):
        self.adaptor.save_edges(
            edges_file_name=edges_file_name,
            edge_types_file_name=edge_types_file_name,
            output_dir=output_dir,
            src_network=src_network,
            trg_network=trg_network,
            name=name,
            force_build=force_build,
            force_overwrite=force_overwrite,
            compression=compression,
            **opts,
        )
