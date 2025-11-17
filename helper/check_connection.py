#!/usr/bin/env python3
"""
Connection Validator

Required CSV Format:

Population CSV (e.g., V1_population.csv):
  Header: Required (# n_cell, n_comp, name, swc_file, ion_file)
  Column 0: # n_cell - Number of neurons of this type
  Column 1: n_comp   - Number of compartments per neuron
  Column 2: name     - Cell type name
  Column 3: swc_file - Morphology file path
  Column 4: ion_file - Ion channel file path
  
Connection CSV (e.g., V1_V1_connection.csv):
  Header: Optional (auto-detected)
  Column 0: #pre nid   - Presynaptic neuron ID (0-based)
  Column 1: post nid   - Postsynaptic neuron ID (0-based)
  Column 2: post cid   - Postsynaptic compartment ID (0-based)
  Column 3: weight     - Synaptic weight
  Column 4: tau_decay  - Decay time constant
  Column 5: tau_rise   - Rise time constant
  Column 6: erev       - Reversal potential
  Column 7: delay      - Synaptic delay
  Column 8: e/i        - Connection type (e=excitatory, i=inhibitory)

Usage:
    python check_connection.py population.csv connection.csv
"""
import sys
import pandas as pd

def check_connection(pop_path, conn_path):
    try:
        # Load population and calculate cumulative neuron IDs
        pop_map = pd.read_csv(pop_path)
        if '# n_cell' in pop_map.columns:
            pop_map["nid"] = pop_map['# n_cell'].cumsum()
        else:
            pop_map["nid"] = pop_map['#n_cell'].cumsum()
        
        # Load connection (auto-detect header)
        # Check if first row is numeric to determine if header exists
        first_line = pd.read_csv(conn_path, nrows=1, header=None)
        try:
            int(first_line.iloc[0, 0])
            # First column is numeric -> no header
            conn = pd.read_csv(
                conn_path,
                header=None,
                names=['#pre nid', 'post nid', 'post cid', 'weight', 'tau_decay', 'tau_rise', 'erev', 'delay', 'e/i']
            )
        except (ValueError, TypeError):
            # First column is not numeric -> has header
            conn = pd.read_csv(conn_path)

        # Validate connections
        error_found = False
        error_count = 0
        for idx, row in conn.iterrows():
            post_nid = int(row['post nid'])
            post_cid = int(row['post cid'])
            
            pop_rows = pop_map[pop_map['nid'] > post_nid]
            if len(pop_rows) == 0:
                print(f"Error: nid={post_nid} exceeds total (neuron ID does not exist)")
                error_found = True
                error_count += 1
                continue
                
            n_comp = pop_rows.iloc[0]['n_comp']
            if post_cid >= n_comp:
                print(f"Error: nid={post_nid}, cid={post_cid} >= n_comp={n_comp} (compartment ID out of range)")
                error_found = True
                error_count += 1
        
        if not error_found:
            print(f"✓ All connections valid")
            return True
        else:
            print(f"\n✗ Found {error_count} invalid connections")
            return False
        
    except Exception as e:
        print(f"Error: {e}")
        return False


def main():
    if len(sys.argv) != 3:
        print("Usage: python check_connection.py <population_csv> <connection_csv>")
        sys.exit(1)
    
    if check_connection(sys.argv[1], sys.argv[2]):
        sys.exit(0)
    else:
        sys.exit(1)


if __name__ == "__main__":
    main()