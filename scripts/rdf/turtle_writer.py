"""
RDF Turtle writer.

Generates RDF Turtle (.ttl) files from graph node and relationship DataFrames
following the Biolink Model ontology.
"""

import hashlib
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

import pandas as pd
from rdflib import Graph, Literal, URIRef, BNode
from rdflib.namespace import RDF, RDFS, XSD

from .rdf_config import (
    NAMESPACES,
    BIOLINK,
    SPOKEGENELAB,
    get_xsd_type,
    create_node_uri,
    create_uri_from_ontology_uri,
)
from .biolink_mapper import (
    get_biolink_class,
    get_biolink_predicate,
    get_property_predicate,
    is_reified_relationship,
    get_association_class,
    parse_relationship_type,
)


class TurtleWriter:
    """
    RDF Turtle file writer following Biolink Model.

    Creates RDF triples from node and relationship DataFrames
    and serializes them to Turtle format.
    """

    def __init__(self, output_path: Optional[Union[str, Path]] = None):
        """
        Initialize the Turtle writer.

        Args:
            output_path: Optional path for output file
        """
        self.graph = Graph()
        self.output_path = Path(output_path) if output_path else None
        self._bind_namespaces()

    def _bind_namespaces(self):
        """Bind all standard namespaces to the graph."""
        for prefix, namespace in NAMESPACES.items():
            self.graph.bind(prefix, namespace)

    def add_node(
        self,
        node_type: str,
        identifier: str,
        properties: Dict[str, Any],
    ) -> URIRef:
        """
        Add a node to the graph as RDF triples.

        Args:
            node_type: Type of node (e.g., "Study", "Gene")
            identifier: Node identifier (can be a full ontology URI)
            properties: Dictionary of property name -> value

        Returns:
            URIRef of the created node
        """
        # Check if identifier is already a full URI (e.g., ontology URI)
        if identifier.startswith("http://") or identifier.startswith("https://"):
            node_uri = create_uri_from_ontology_uri(identifier)
        else:
            # Create node URI using namespace
            node_uri = create_node_uri(node_type, identifier)

        # Add type triple
        biolink_class = get_biolink_class(
            node_type,
            category=properties.get("category"),
        )
        self.graph.add((node_uri, RDF.type, biolink_class))

        # Add identifier
        self.graph.add((
            node_uri,
            BIOLINK.id,
            Literal(str(identifier), datatype=XSD.string),
        ))

        # Add other properties
        for prop_name, value in properties.items():
            if prop_name == "identifier":
                continue  # Already handled
            if value is None or (isinstance(value, float) and pd.isna(value)):
                continue

            predicate = get_property_predicate(prop_name, node_type)

            # Handle lists
            if isinstance(value, list):
                for item in value:
                    self.graph.add((
                        node_uri,
                        predicate,
                        Literal(str(item), datatype=XSD.string),
                    ))
            else:
                # Determine datatype
                if isinstance(value, bool):
                    literal = Literal(value, datatype=XSD.boolean)
                elif isinstance(value, int):
                    literal = Literal(value, datatype=XSD.integer)
                elif isinstance(value, float):
                    literal = Literal(value, datatype=XSD.float)
                else:
                    literal = Literal(str(value), datatype=XSD.string)

                self.graph.add((node_uri, predicate, literal))

        return node_uri

    def add_relationship(
        self,
        rel_type: str,
        from_type: str,
        from_id: str,
        to_type: str,
        to_id: str,
        properties: Optional[Dict[str, Any]] = None,
    ) -> Optional[URIRef]:
        """
        Add a relationship to the graph.

        For relationships with properties, uses Biolink Association pattern
        (reified edges).

        Args:
            rel_type: Relationship type
            from_type: Source node type
            from_id: Source node identifier (can be a full URI)
            to_type: Target node type
            to_id: Target node identifier (can be a full URI)
            properties: Optional relationship properties

        Returns:
            URIRef of association node if reified, None otherwise
        """
        # Handle full URIs for from_id
        if from_id.startswith("http://") or from_id.startswith("https://"):
            from_uri = create_uri_from_ontology_uri(from_id)
        else:
            from_uri = create_node_uri(from_type, from_id)

        # Handle full URIs for to_id
        if to_id.startswith("http://") or to_id.startswith("https://"):
            to_uri = create_uri_from_ontology_uri(to_id)
        else:
            to_uri = create_node_uri(to_type, to_id)

        predicate = get_biolink_predicate(rel_type)

        # Check if we need to reify this relationship
        if properties and is_reified_relationship(rel_type):
            return self._add_reified_relationship(
                rel_type, from_uri, to_uri, predicate, properties
            )
        else:
            # Simple triple
            self.graph.add((from_uri, predicate, to_uri))
            return None

    def _add_reified_relationship(
        self,
        rel_type: str,
        from_uri: URIRef,
        to_uri: URIRef,
        predicate: URIRef,
        properties: Dict[str, Any],
    ) -> URIRef:
        """
        Add a reified relationship as a Biolink Association.

        Args:
            rel_type: Relationship type
            from_uri: Source node URI
            to_uri: Target node URI
            predicate: Relationship predicate
            properties: Relationship properties

        Returns:
            URIRef of the Association node
        """
        # Generate unique association ID
        assoc_id = hashlib.md5(
            f"{from_uri}_{to_uri}_{rel_type}".encode()
        ).hexdigest()[:12]

        assoc_uri = SPOKEGENELAB[f"Association/{assoc_id}"]

        # Add Association type
        assoc_class = get_association_class(rel_type)
        self.graph.add((assoc_uri, RDF.type, assoc_class))

        # Add subject, predicate, object
        self.graph.add((assoc_uri, BIOLINK.subject, from_uri))
        self.graph.add((assoc_uri, BIOLINK.predicate, predicate))
        self.graph.add((assoc_uri, BIOLINK.object, to_uri))

        # Add properties
        for prop_name, value in properties.items():
            if prop_name in ("from", "to"):
                continue
            if value is None or (isinstance(value, float) and pd.isna(value)):
                continue

            prop_uri = SPOKEGENELAB[prop_name]

            if isinstance(value, (int, float)):
                literal = Literal(value, datatype=XSD.float)
            else:
                literal = Literal(str(value), datatype=XSD.string)

            self.graph.add((assoc_uri, prop_uri, literal))

        return assoc_uri

    def add_nodes_from_dataframe(
        self,
        df: pd.DataFrame,
        node_type: str,
        id_column: str = "identifier",
    ):
        """
        Add multiple nodes from a DataFrame.

        Args:
            df: DataFrame with node data
            node_type: Type of nodes
            id_column: Column to use as identifier
        """
        for _, row in df.iterrows():
            identifier = row[id_column]
            properties = row.to_dict()
            self.add_node(node_type, identifier, properties)

    def add_relationships_from_dataframe(
        self,
        df: pd.DataFrame,
        rel_type: str,
        from_column: str = "from",
        to_column: str = "to",
    ):
        """
        Add multiple relationships from a DataFrame.

        Args:
            df: DataFrame with relationship data
            rel_type: Relationship type
            from_column: Column for source IDs
            to_column: Column for target IDs
        """
        # Parse relationship type to get node types
        from_type, _, to_type = parse_relationship_type(rel_type)

        # Get property columns (excluding from/to)
        prop_columns = [c for c in df.columns if c not in (from_column, to_column)]

        # Check if this relationship should be reified (has properties or is DE/methylation)
        should_reify = len(prop_columns) > 0 or is_reified_relationship(rel_type)

        for _, row in df.iterrows():
            from_id = row[from_column]
            to_id = row[to_column]

            # Collect properties
            properties = {c: row[c] for c in prop_columns} if prop_columns else None

            if should_reify and properties:
                # Create reified association
                self._add_reified_relationship_from_row(
                    rel_type, from_type, from_id, to_type, to_id, properties
                )
            else:
                self.add_relationship(
                    rel_type,
                    from_type, from_id,
                    to_type, to_id,
                    properties,
                )

    def _add_reified_relationship_from_row(
        self,
        rel_type: str,
        from_type: str,
        from_id: str,
        to_type: str,
        to_id: str,
        properties: Dict[str, Any],
    ) -> URIRef:
        """
        Create a reified relationship as a Biolink Association.

        This method ensures proper handling of DE relationships with log2fc/p-value.
        """
        # Handle full URIs for from_id
        if from_id.startswith("http://") or from_id.startswith("https://"):
            from_uri = create_uri_from_ontology_uri(from_id)
        else:
            from_uri = create_node_uri(from_type, from_id)

        # Handle full URIs for to_id
        if to_id.startswith("http://") or to_id.startswith("https://"):
            to_uri = create_uri_from_ontology_uri(to_id)
        else:
            to_uri = create_node_uri(to_type, to_id)

        predicate = get_biolink_predicate(rel_type)

        # Generate unique association ID
        assoc_id = hashlib.md5(
            f"{from_id}_{to_id}_{rel_type}".encode()
        ).hexdigest()[:12]

        assoc_uri = SPOKEGENELAB[f"Association/{assoc_id}"]

        # Add Association type
        assoc_class = get_association_class(rel_type)
        self.graph.add((assoc_uri, RDF.type, assoc_class))

        # Add subject, predicate, object
        self.graph.add((assoc_uri, BIOLINK.subject, from_uri))
        self.graph.add((assoc_uri, BIOLINK.predicate, predicate))
        self.graph.add((assoc_uri, BIOLINK.object, to_uri))

        # Add properties
        for prop_name, value in properties.items():
            if value is None or (isinstance(value, float) and pd.isna(value)):
                continue

            prop_uri = SPOKEGENELAB[prop_name]

            if isinstance(value, bool):
                literal = Literal(value, datatype=XSD.boolean)
            elif isinstance(value, int):
                literal = Literal(value, datatype=XSD.integer)
            elif isinstance(value, float):
                literal = Literal(value, datatype=XSD.float)
            else:
                literal = Literal(str(value), datatype=XSD.string)

            self.graph.add((assoc_uri, prop_uri, literal))

        return assoc_uri

    def serialize(self, format: str = "turtle") -> str:
        """
        Serialize the graph to a string.

        Args:
            format: Output format (turtle, xml, n3, etc.)

        Returns:
            Serialized RDF string
        """
        return self.graph.serialize(format=format)

    def write(self, output_path: Optional[Union[str, Path]] = None):
        """
        Write the graph to a file.

        Args:
            output_path: Path to output file (uses self.output_path if not provided)
        """
        path = Path(output_path) if output_path else self.output_path
        if not path:
            raise ValueError("No output path specified")

        path.parent.mkdir(parents=True, exist_ok=True)
        self.graph.serialize(destination=str(path), format="turtle")
        print(f"Wrote {len(self.graph)} triples to {path}")

    def get_triple_count(self) -> int:
        """Get the number of triples in the graph."""
        return len(self.graph)

    def query(self, sparql_query: str) -> List[Dict]:
        """
        Execute a SPARQL query on the graph.

        Args:
            sparql_query: SPARQL query string

        Returns:
            List of result dictionaries
        """
        results = []
        for row in self.graph.query(sparql_query):
            results.append({str(var): str(val) for var, val in zip(row.labels, row)})
        return results


def write_graph_to_turtle(
    nodes: Dict[str, pd.DataFrame],
    relationships: Dict[str, pd.DataFrame],
    output_path: Union[str, Path],
) -> TurtleWriter:
    """
    Write all graph data to a Turtle file.

    Args:
        nodes: Dictionary of node type -> DataFrame
        relationships: Dictionary of relationship type -> DataFrame
        output_path: Path to output file

    Returns:
        TurtleWriter instance with the graph
    """
    writer = TurtleWriter(output_path)

    # Add nodes
    for node_type, df in nodes.items():
        if not df.empty:
            print(f"Adding {len(df)} {node_type} nodes...")
            writer.add_nodes_from_dataframe(df, node_type)

    # Add relationships
    for rel_type, df in relationships.items():
        if not df.empty:
            print(f"Adding {len(df)} {rel_type} relationships...")
            writer.add_relationships_from_dataframe(df, rel_type)

    # Write to file
    writer.write()

    return writer


def validate_turtle_syntax(turtle_file: Union[str, Path]) -> bool:
    """
    Validate that a Turtle file is syntactically correct.

    Args:
        turtle_file: Path to Turtle file

    Returns:
        True if valid, raises exception otherwise
    """
    g = Graph()
    g.parse(str(turtle_file), format="turtle")
    print(f"Validated {len(g)} triples in {turtle_file}")
    return True
