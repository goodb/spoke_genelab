"""
Graph building utilities for creating nodes and relationships.

Provides common functions for creating DataFrame-based representations
of graph nodes and relationships that can be output to Neo4j CSV or RDF.
"""

import hashlib
import json
from dataclasses import dataclass, field
from datetime import datetime
from typing import Any, Dict, List, Optional, Set

import pandas as pd


@dataclass
class NodeSchema:
    """Schema definition for a node type."""

    name: str
    id_property: str = "identifier"
    required_properties: List[str] = field(default_factory=list)
    optional_properties: List[str] = field(default_factory=list)

    @property
    def all_properties(self) -> List[str]:
        return [self.id_property] + self.required_properties + self.optional_properties


@dataclass
class RelationshipSchema:
    """Schema definition for a relationship type."""

    name: str
    from_node: str
    to_node: str
    from_property: str = "from"
    to_property: str = "to"
    properties: List[str] = field(default_factory=list)


# Standard node schemas for SPOKE-GeneLab
NODE_SCHEMAS = {
    "Mission": NodeSchema(
        name="Mission",
        required_properties=["name"],
        optional_properties=["flight_program", "space_program", "start_date", "end_date"],
    ),
    "Study": NodeSchema(
        name="Study",
        required_properties=["name"],
        optional_properties=[
            "project_title", "project_type", "description",
            "organism", "taxonomy",
        ],
    ),
    "Assay": NodeSchema(
        name="Assay",
        required_properties=["name", "technology", "measurement"],
        optional_properties=[
            "factors_1", "factors_2", "material_1", "material_2",
            "material_id_1", "material_id_2", "factor_space_1", "factor_space_2",
        ],
    ),
    "MGene": NodeSchema(
        name="MGene",
        required_properties=["symbol"],
        optional_properties=["name", "organism", "taxonomy"],
    ),
    "Gene": NodeSchema(
        name="Gene",
        required_properties=[],
        optional_properties=["symbol", "name"],
    ),
    "Anatomy": NodeSchema(
        name="Anatomy",
        required_properties=[],
        optional_properties=["name"],
    ),
    "CellType": NodeSchema(
        name="CellType",
        required_properties=[],
        optional_properties=["name"],
    ),
    "PathwayEnrichment": NodeSchema(
        name="PathwayEnrichment",
        required_properties=["enrichment_type", "contrast_id"],
        optional_properties=[
            "name", "p_value", "adj_p_value", "effect_size",
            "genes_total", "genes_significant",
        ],
    ),
    "GOTerm": NodeSchema(
        name="GOTerm",
        required_properties=[],
        optional_properties=["name", "category"],
    ),
    "ReactomePathway": NodeSchema(
        name="ReactomePathway",
        required_properties=[],
        optional_properties=["name"],
    ),
    "InterProDomain": NodeSchema(
        name="InterProDomain",
        required_properties=[],
        optional_properties=["name"],
    ),
}


class GraphBuilder:
    """
    Builder class for constructing graph nodes and relationships.

    Tracks created nodes and relationships to ensure uniqueness
    and enable batch output to CSV or RDF.
    """

    def __init__(self):
        self._nodes: Dict[str, pd.DataFrame] = {}
        self._relationships: Dict[str, pd.DataFrame] = {}
        self._node_ids: Dict[str, Set[str]] = {}

    def add_nodes(
        self,
        node_type: str,
        data: pd.DataFrame,
        id_column: str = "identifier",
    ) -> pd.DataFrame:
        """
        Add nodes of a given type from a DataFrame.

        Args:
            node_type: Type of node (e.g., "Study", "Gene")
            data: DataFrame with node data
            id_column: Column to use as node identifier

        Returns:
            DataFrame of added nodes (deduplicated)
        """
        if data.empty:
            return data

        # Ensure identifier column exists
        if id_column != "identifier" and "identifier" not in data.columns:
            data = data.rename(columns={id_column: "identifier"})

        # Deduplicate by identifier
        data = data.drop_duplicates(subset="identifier")

        # Track node IDs
        if node_type not in self._node_ids:
            self._node_ids[node_type] = set()

        new_ids = set(data["identifier"].astype(str))
        existing_ids = self._node_ids[node_type]

        # Filter to only new nodes
        new_nodes = data[~data["identifier"].astype(str).isin(existing_ids)]

        if not new_nodes.empty:
            self._node_ids[node_type].update(new_ids)

            if node_type in self._nodes:
                self._nodes[node_type] = pd.concat(
                    [self._nodes[node_type], new_nodes],
                    ignore_index=True,
                )
            else:
                self._nodes[node_type] = new_nodes.copy()

        return new_nodes

    def add_relationships(
        self,
        rel_type: str,
        data: pd.DataFrame,
        from_column: str = "from",
        to_column: str = "to",
    ) -> pd.DataFrame:
        """
        Add relationships of a given type from a DataFrame.

        Args:
            rel_type: Type of relationship (e.g., "CONDUCTED", "PERFORMED")
            data: DataFrame with relationship data
            from_column: Column for source node ID
            to_column: Column for target node ID

        Returns:
            DataFrame of added relationships (deduplicated)
        """
        if data.empty:
            return data

        # Standardize column names
        if from_column != "from":
            data = data.rename(columns={from_column: "from"})
        if to_column != "to":
            data = data.rename(columns={to_column: "to"})

        # Deduplicate by from/to pair
        data = data.drop_duplicates(subset=["from", "to"])

        if rel_type in self._relationships:
            # Merge and deduplicate
            combined = pd.concat(
                [self._relationships[rel_type], data],
                ignore_index=True,
            )
            self._relationships[rel_type] = combined.drop_duplicates(
                subset=["from", "to"]
            )
        else:
            self._relationships[rel_type] = data.copy()

        return data

    def get_nodes(self, node_type: str) -> pd.DataFrame:
        """Get all nodes of a given type."""
        return self._nodes.get(node_type, pd.DataFrame())

    def get_relationships(self, rel_type: str) -> pd.DataFrame:
        """Get all relationships of a given type."""
        return self._relationships.get(rel_type, pd.DataFrame())

    def get_all_nodes(self) -> Dict[str, pd.DataFrame]:
        """Get all nodes by type."""
        return self._nodes.copy()

    def get_all_relationships(self) -> Dict[str, pd.DataFrame]:
        """Get all relationships by type."""
        return self._relationships.copy()

    def node_exists(self, node_type: str, identifier: str) -> bool:
        """Check if a node with given identifier exists."""
        if node_type not in self._node_ids:
            return False
        return str(identifier) in self._node_ids[node_type]

    def get_stats(self) -> Dict[str, int]:
        """Get statistics about the current graph."""
        stats = {}
        for node_type, df in self._nodes.items():
            stats[f"nodes_{node_type}"] = len(df)
        for rel_type, df in self._relationships.items():
            stats[f"relationships_{rel_type}"] = len(df)
        return stats

    def clear(self):
        """Clear all nodes and relationships."""
        self._nodes.clear()
        self._relationships.clear()
        self._node_ids.clear()


def generate_hash_id(data: Dict[str, Any]) -> str:
    """
    Generate a deterministic hash ID from a dictionary.

    Args:
        data: Dictionary of values to hash

    Returns:
        MD5 hash string
    """
    json_str = json.dumps(data, sort_keys=True)
    return hashlib.md5(json_str.encode()).hexdigest()


def list_to_string(values: List[Any], separator: str = "|") -> str:
    """
    Convert a list to a separator-delimited string.

    Args:
        values: List of values
        separator: Delimiter (default "|" for Neo4j)

    Returns:
        Delimited string
    """
    if not values:
        return ""
    return separator.join(str(v) for v in values)


def convert_lists_in_dataframe(df: pd.DataFrame, separator: str = "|") -> pd.DataFrame:
    """
    Convert all list-type columns in a DataFrame to delimited strings.

    Args:
        df: DataFrame to convert
        separator: Delimiter for list values

    Returns:
        DataFrame with lists converted to strings
    """
    df = df.copy()
    for col in df.columns:
        if df[col].apply(lambda x: isinstance(x, list)).any():
            df[col] = df[col].apply(
                lambda x: separator.join(str(v) for v in x) if isinstance(x, list) else x
            )
    return df


def create_study_node(
    identifier: str,
    name: str,
    **kwargs,
) -> Dict[str, Any]:
    """Create a Study node dictionary."""
    node = {
        "identifier": identifier,
        "name": name,
    }
    node.update(kwargs)
    return node


def create_assay_node(
    identifier: str,
    name: str,
    technology: str,
    measurement: str,
    **kwargs,
) -> Dict[str, Any]:
    """Create an Assay node dictionary."""
    node = {
        "identifier": identifier,
        "name": name,
        "technology": technology,
        "measurement": measurement,
    }
    node.update(kwargs)
    return node


def create_gene_node(
    identifier: str,
    symbol: str = "",
    name: str = "",
    organism: str = "",
    taxonomy: str = "",
) -> Dict[str, Any]:
    """Create a Gene or MGene node dictionary."""
    return {
        "identifier": identifier,
        "symbol": symbol,
        "name": name,
        "organism": organism,
        "taxonomy": taxonomy,
    }


def create_pathway_enrichment_node(
    identifier: str,
    enrichment_type: str,
    contrast_id: str,
    name: str = "",
    p_value: float = None,
    adj_p_value: float = None,
    effect_size: float = None,
    genes_total: int = None,
    genes_significant: int = None,
) -> Dict[str, Any]:
    """Create a PathwayEnrichment node dictionary."""
    node = {
        "identifier": identifier,
        "enrichment_type": enrichment_type,
        "contrast_id": contrast_id,
        "name": name,
    }
    if p_value is not None:
        node["p_value"] = p_value
    if adj_p_value is not None:
        node["adj_p_value"] = adj_p_value
    if effect_size is not None:
        node["effect_size"] = effect_size
    if genes_total is not None:
        node["genes_total"] = genes_total
    if genes_significant is not None:
        node["genes_significant"] = genes_significant
    return node
