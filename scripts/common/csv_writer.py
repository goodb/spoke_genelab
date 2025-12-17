"""
CSV output utilities for Neo4j bulk import.

Provides functions for saving node and relationship DataFrames
to CSV files compatible with the Neo4j bulk import tool.
"""

import glob
import os
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Union

import pandas as pd

from .graph_builder import convert_lists_in_dataframe


def save_nodes_to_csv(
    df: pd.DataFrame,
    node_name: str,
    output_dir: Union[str, Path],
    id_column: str = "identifier",
    date_suffix: bool = True,
) -> Path:
    """
    Save a DataFrame of nodes to a CSV file for Neo4j import.

    Args:
        df: DataFrame containing node data
        node_name: Name of the node type (e.g., "Study", "Gene")
        output_dir: Directory to write the CSV file
        id_column: Column to use as node identifier
        date_suffix: Whether to add date suffix to filename

    Returns:
        Path to the saved CSV file
    """
    if df.empty:
        print(f"Warning: No {node_name} nodes to save")
        return None

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Remove previous versions
    for old_file in glob.glob(str(output_dir / f"{node_name}_*.csv")):
        os.remove(old_file)

    # Convert lists to pipe-delimited strings
    df = convert_lists_in_dataframe(df.copy())

    # Deduplicate by identifier
    if id_column in df.columns:
        df = df.drop_duplicates(subset=id_column)

    # Generate filename
    if date_suffix:
        date_str = datetime.today().strftime("%Y-%m-%d")
        filename = f"{node_name}_{date_str}.csv"
    else:
        filename = f"{node_name}.csv"

    file_path = output_dir / filename
    df.to_csv(file_path, index=False)

    print(f"Saved {len(df)} {node_name} nodes to {file_path}")
    return file_path


def save_relationships_to_csv(
    df: pd.DataFrame,
    rel_name: str,
    output_dir: Union[str, Path],
    from_column: str = "from",
    to_column: str = "to",
    date_suffix: bool = True,
) -> Path:
    """
    Save a DataFrame of relationships to a CSV file for Neo4j import.

    Args:
        df: DataFrame containing relationship data
        rel_name: Name of the relationship type (e.g., "Study-PERFORMED-Assay")
        output_dir: Directory to write the CSV file
        from_column: Column for source node ID
        to_column: Column for target node ID
        date_suffix: Whether to add date suffix to filename

    Returns:
        Path to the saved CSV file
    """
    if df.empty:
        print(f"Warning: No {rel_name} relationships to save")
        return None

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Remove previous versions
    for old_file in glob.glob(str(output_dir / f"{rel_name}_*.csv")):
        os.remove(old_file)

    # Convert lists to pipe-delimited strings
    df = convert_lists_in_dataframe(df.copy())

    # Deduplicate by from/to pair
    if from_column in df.columns and to_column in df.columns:
        df = df.drop_duplicates(subset=[from_column, to_column])

    # Generate filename
    if date_suffix:
        date_str = datetime.today().strftime("%Y-%m-%d")
        filename = f"{rel_name}_{date_str}.csv"
    else:
        filename = f"{rel_name}.csv"

    file_path = output_dir / filename
    df.to_csv(file_path, index=False)

    print(f"Saved {len(df)} {rel_name} relationships to {file_path}")
    return file_path


def save_graph_to_csv(
    nodes: Dict[str, pd.DataFrame],
    relationships: Dict[str, pd.DataFrame],
    node_dir: Union[str, Path],
    rel_dir: Union[str, Path],
    date_suffix: bool = True,
) -> Dict[str, List[Path]]:
    """
    Save all nodes and relationships to CSV files.

    Args:
        nodes: Dictionary of node type -> DataFrame
        relationships: Dictionary of relationship type -> DataFrame
        node_dir: Directory for node CSV files
        rel_dir: Directory for relationship CSV files
        date_suffix: Whether to add date suffix to filenames

    Returns:
        Dictionary with 'nodes' and 'relationships' keys containing lists of saved paths
    """
    saved_paths = {"nodes": [], "relationships": []}

    # Save nodes
    for node_type, df in nodes.items():
        path = save_nodes_to_csv(df, node_type, node_dir, date_suffix=date_suffix)
        if path:
            saved_paths["nodes"].append(path)

    # Save relationships
    for rel_type, df in relationships.items():
        path = save_relationships_to_csv(df, rel_type, rel_dir, date_suffix=date_suffix)
        if path:
            saved_paths["relationships"].append(path)

    return saved_paths


def validate_csv_for_neo4j(
    file_path: Union[str, Path],
    is_node: bool = True,
    id_column: str = "identifier",
    from_column: str = "from",
    to_column: str = "to",
) -> bool:
    """
    Validate a CSV file for Neo4j bulk import.

    Args:
        file_path: Path to CSV file
        is_node: Whether this is a node file (vs relationship)
        id_column: Expected identifier column for nodes
        from_column: Expected from column for relationships
        to_column: Expected to column for relationships

    Returns:
        True if valid, raises ValueError otherwise
    """
    df = pd.read_csv(file_path)

    if is_node:
        if id_column not in df.columns:
            raise ValueError(f"Node file missing '{id_column}' column: {file_path}")
        if df[id_column].isna().any():
            raise ValueError(f"Node file has null identifiers: {file_path}")
        if df[id_column].duplicated().any():
            raise ValueError(f"Node file has duplicate identifiers: {file_path}")
    else:
        if from_column not in df.columns:
            raise ValueError(f"Relationship file missing '{from_column}' column: {file_path}")
        if to_column not in df.columns:
            raise ValueError(f"Relationship file missing '{to_column}' column: {file_path}")
        if df[from_column].isna().any() or df[to_column].isna().any():
            raise ValueError(f"Relationship file has null endpoints: {file_path}")

    return True


def generate_neo4j_header(
    metadata_file: Union[str, Path],
    node_or_rel_name: str,
    is_node: bool = True,
) -> str:
    """
    Generate Neo4j bulk import header from metadata CSV.

    Args:
        metadata_file: Path to metadata CSV file
        node_or_rel_name: Name of node or relationship type
        is_node: Whether this is a node (vs relationship)

    Returns:
        Neo4j bulk import header string
    """
    meta = pd.read_csv(metadata_file)

    headers = []
    for _, row in meta.iterrows():
        prop = row["property"]
        prop_type = row["type"]

        # Map types to Neo4j format
        type_map = {
            "string": "string",
            "int": "int",
            "float": "float",
            "boolean": "boolean",
            "string[]": "string[]",
        }

        neo4j_type = type_map.get(prop_type, "string")

        if is_node:
            if prop == "identifier":
                headers.append(f"identifier:ID({node_or_rel_name}-ID)")
            else:
                headers.append(f"{prop}:{neo4j_type}")
        else:
            if prop == "from":
                # Extract source node type from relationship name
                from_node = node_or_rel_name.split("-")[0]
                headers.append(f"from:START_ID({from_node}-ID)")
            elif prop == "to":
                # Extract target node type from relationship name
                to_node = node_or_rel_name.split("-")[-1]
                headers.append(f"to:END_ID({to_node}-ID)")
            else:
                headers.append(f"{prop}:{neo4j_type}")

    return ",".join(headers)
