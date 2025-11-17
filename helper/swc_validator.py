#!/usr/bin/env python3
import sys
import csv
from collections import defaultdict
import os
import glob
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA

def validate_swc_file(swc_file_path):
    swc_filename = os.path.basename(swc_file_path)

    try:
        with open(swc_file_path, 'r') as file:
            lines = file.readlines()
    except IOError:
        print(f"Error: Cannot open SWC file '{swc_file_path}'.")
        return False

    compartments = []
    point_ids = set()
    xyz_set = set()
    children = defaultdict(list)
    data_line_count = 0
    axon_count = 0

    for line_num, line in enumerate(lines, 1):
        if line.startswith('#'):
            continue  # Skip comment lines

        parts = line.strip().split()

        if len(parts) != 7:
            print(f"Error: Invalid format at line {line_num}. Seven space-separated elements are required.")
            return False

        try:
            id, type, x, y, z, radius, pid = int(parts[0]), int(parts[1]), float(parts[2]), float(parts[3]), float(parts[4]), float(parts[5]), int(parts[6])
        except ValueError:
            print(f"Error: Invalid numeric value at line {line_num}.")
            return False

        data_line_count += 1

        # Check if id starts from 0
        if data_line_count == 1 and id != 0:
            print(f"Error: The ID of the first data line (line number: {line_num}) is not 0.")
            return False

        # Check id continuity
        if id != data_line_count - 1:
            print(f"Error: ID {id} at line {line_num} is not consecutive. Expected ID: {data_line_count - 1}")
            return False

        # Check id uniqueness
        if id in point_ids:
            print(f"Error: ID {id} at line {line_num} is duplicated.")
            return False
        point_ids.add(id)

        # Check xyz uniqueness
        xyz = (x, y, z)
        if xyz in xyz_set:
            print(f"Error: Coordinates {xyz} at line {line_num} are duplicated.")
            return False
        xyz_set.add(xyz)

        # Check pid validity
        if pid != -1 and pid >= id:
            print(f"Error: Parent ID {pid} at line {line_num} is invalid.")
            return False

        # Check soma
        if id == 0 and (pid != -1 or type != 1):
            print(f"Error: ID 0 (soma, line number: {line_num}) must have parent ID -1 and type 1.")
            return False

        if pid == -1 and id != 0:
            print(f"Error: ID {id} at line {line_num} has parent ID -1, which is only allowed for soma (ID 0).")
            return False

        # Check axon (type 2)
        if type == 2:
            axon_count = axon_count + 1
            
            if axon_count > 2:
                print(f"Error: More than 2 axon compartments exist at line {line_num}. In the perisomatic model, existing axons should be removed and replaced with 2 connected compartments.")

        compartments.append((id, type, x, y, z, radius, pid))
        children[pid].append(id)
    
    print(f"Axon count: {axon_count}")
    if axon_count != 2:
        print(f"Error: Invalid number of axon compartments. In the perisomatic model, existing axons should be removed and replaced with 2 connected compartments.")
        return False

    # Check DFS sorting and existence of subtrees
    visited = set()
    stack = [0]  # Start from soma
    sorted_ids = []

    while stack:
        current = stack.pop()
        if current not in visited:
            visited.add(current)
            sorted_ids.append(current)
            stack.extend(reversed(children[current]))

    if sorted_ids != list(range(len(compartments))):
        print("Error: Compartments are not sorted by DFS or subtrees exist.")
        return False

    print(f"Validation successful: SWC file {swc_file_path} meets all conditions.")
    return True

def check_duplicate_coordinates(df):
    """
    Check if duplicate coordinates exist in the DataFrame.

    Args:
        df: pandas DataFrame with 'x', 'y', 'z' columns

    Returns:
        True if duplicates exist, False otherwise
    """
    xyz_set = set()
    for _, row in df.iterrows():
        xyz = (row['x'], row['y'], row['z'])
        if xyz in xyz_set:
            return True
        xyz_set.add(xyz)
    return False

def convert_to_zero_based(swc_df):
    """
    Convert SWC IDs to zero-based indexing if they start from 1.

    Args:
        swc_df: pandas DataFrame with SWC data

    Returns:
        Modified DataFrame with zero-based IDs
    """
    if swc_df['id'][0] == 1:
        swc_df['id'] = swc_df['id'] - 1
        swc_df['parent'] = swc_df['parent'].mask(swc_df['parent'] > 0, swc_df['parent']-1)
    return swc_df

def sort_ids_by_dfs_recursive(graph, current, visited, sorted_nodes):
    """
    Recursively sort IDs using depth-first search.

    Args:
        graph: Dictionary mapping parent ID to list of child IDs
        current: Current node ID being visited
        visited: Set of already visited node IDs
        sorted_nodes: List to store sorted node IDs

    Returns:
        List of sorted node IDs
    """
    visited.add(current)
    sorted_nodes.append(current)

    for child in graph.get(current, []):
        if child not in visited:
            sort_ids_by_dfs_recursive(graph, child, visited, sorted_nodes)

    return sorted_nodes

def sort_by_dfs(swc_df):
    """
    Sort SWC data by DFS order and reassign IDs.

    Args:
        swc_df: pandas DataFrame with SWC data

    Returns:
        New DataFrame with DFS-sorted IDs
    """
    # Create graph: key=parent, value=list of IDs
    graph = swc_df[swc_df["parent"] != -1].groupby("parent")["id"].apply(list).to_dict()

    # Get root node
    root = swc_df[swc_df['parent'] == -1]['id'].values[0]

    # Sort IDs by DFS
    sorted_ids = sort_ids_by_dfs_recursive(graph, root, set(), [])

    # Create ID mapping
    id_mapping = {old_id: new_id for new_id, old_id in enumerate(sorted_ids, 0)}

    # Create new DataFrame
    new_df = swc_df.copy()
    new_df['new_id'] = new_df['id'].map(id_mapping)
    new_df['new_parent'] = new_df['parent'].map(lambda x: id_mapping.get(x, -1))

    # Sort by new ID
    new_df = new_df.sort_values('new_id')

    # Create final DataFrame
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

def get_axon_direction(df):
    """
    Calculate axon direction using PCA and return new axon segment coordinates.

    Args:
        df: pandas DataFrame with SWC data

    Returns:
        numpy array with shape (2, 3) containing coordinates for two new axon segments
    """
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

def fix_axon_perisomatic(df):
    """
    Replace existing axon with two new perisomatic axon compartments.

    Args:
        df: pandas DataFrame with SWC data

    Returns:
        New DataFrame with axon replaced
    """
    # Check if axon exists
    if 2 not in df['type'].values:
        raise Exception('There is no axonal reconstruction in swc file.')

    # Get new axon coordinates
    axon_coords = get_axon_direction(df)

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

def rebuild_swc_file(input_path, output_path):
    """
    Rebuild a single SWC file by applying preprocessing steps.

    Args:
        input_path: Path to input SWC file
        output_path: Path to output SWC file

    Returns:
        True on success, False on failure
    """
    try:
        print(f"Rebuilding: {input_path} -> {output_path}")

        # Read file
        swc = pd.read_csv(
            input_path,
            names=['id', 'type', 'x', 'y', 'z', 'r', 'parent'],
            comment="#",
            delimiter=" ",
            dtype={'id': int, 'parent': int}
        )

        # Check for duplicate coordinates before preprocessing
        if check_duplicate_coordinates(swc):
            print(f"Warning: Duplicate coordinates detected in {input_path}")
            print("         These will NOT be fixed during rebuild.")

        # Apply preprocessing steps
        swc = convert_to_zero_based(swc)
        swc = fix_axon_perisomatic(swc)
        swc = sort_by_dfs(swc)

        # Check for duplicate coordinates after preprocessing
        if check_duplicate_coordinates(swc):
            print("Warning: Duplicate coordinates still present after rebuild.")

        # Save to output file
        swc = swc.rename(columns={'id': '#id'})
        swc.to_csv(output_path, index=False, sep=" ", float_format='%.6f')

        print(f"Successfully rebuilt: {output_path}")
        return True

    except Exception as e:
        print(f"Error rebuilding {input_path}: {e}")
        return False

def rebuild_directory(src_dir, dst_dir):
    """
    Rebuild all SWC files in a directory.

    Args:
        src_dir: Source directory containing SWC files
        dst_dir: Destination directory for rebuilt files

    Returns:
        True if all files succeed, False if any fail
    """
    swc_files = glob.glob(os.path.join(src_dir, "*.swc"))

    if not swc_files:
        print(f"No SWC files found in {src_dir}")
        return False

    print(f"Found {len(swc_files)} SWC files to rebuild")

    # Create destination directory if it doesn't exist
    os.makedirs(dst_dir, exist_ok=True)

    for file in swc_files:
        output_file = os.path.join(dst_dir, os.path.basename(file))
        if not rebuild_swc_file(file, output_file):
            print(f"Failed to rebuild {file}")
            return False

    print(f"\nAll SWC files rebuilt successfully to {dst_dir}!")
    return True

def validate_directory(directory_path):
    """
    Validate all SWC files in a directory.
    
    Args:
        directory_path: Path to directory containing SWC files
    """
    swc_files = glob.glob(os.path.join(directory_path, "*.swc"))
    
    if not swc_files:
        print(f"No SWC files found in {directory_path}")
        return False
    
    print(f"Found {len(swc_files)} SWC files to validate")
    
    all_valid = True
    for file in swc_files:
        print(f"\nValidating: {os.path.basename(file)}")
        res = validate_swc_file(file)
        if not res:
            all_valid = False
            print(f"Validation failed for {file}")
            break
    
    if all_valid:
        print("\nAll SWC files validated successfully!")
    
    return all_valid

def main():
    if len(sys.argv) < 2:
        print("Usage:")
        print("  Validate single file:    python swc_validator.py <swc_file>")
        print("  Rebuild  single file:    python swc_validator.py <swc_file_src> <swc_file_dst>")
        print("  Validate directory:      python swc_validator.py -d <directory>")
        print("  Rebuild  directory:      python swc_validator.py -d <directory_src> <directory_dst>")
        sys.exit(1)

    if sys.argv[1] == "-d":
        # Directory mode
        if len(sys.argv) == 3:
            # Validate only
            if validate_directory(sys.argv[2]):
                sys.exit(0)
            else:
                sys.exit(1)
        elif len(sys.argv) == 4:
            # Rebuild
            if rebuild_directory(sys.argv[2], sys.argv[3]):
                sys.exit(0)
            else:
                sys.exit(1)
        else:
            print("Error: Invalid number of arguments for directory mode")
            sys.exit(1)
    else:
        # File mode
        if len(sys.argv) == 2:
            # Validate only
            if validate_swc_file(sys.argv[1]):
                sys.exit(0)
            else:
                sys.exit(1)
        elif len(sys.argv) == 3:
            # Rebuild
            if rebuild_swc_file(sys.argv[1], sys.argv[2]):
                sys.exit(0)
            else:
                sys.exit(1)
        else:
            print("Error: Invalid number of arguments")
            sys.exit(1)

if __name__ == "__main__":
    main()