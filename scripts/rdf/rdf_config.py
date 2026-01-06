"""
RDF namespace and configuration definitions.

Defines standard namespaces for RDF output following Biolink Model conventions.
"""

import re
from urllib.parse import quote

from rdflib import Namespace, URIRef
from rdflib.namespace import RDF, RDFS, XSD, OWL

# Standard ontology namespaces
BIOLINK = Namespace("https://w3id.org/biolink/vocab/")
OBO = Namespace("http://purl.obolibrary.org/obo/")

# Project-specific namespace
SPOKEGENELAB = Namespace("https://spoke.ucsf.edu/genelab/")

# External identifier namespaces
NCBIGENE = Namespace("https://www.ncbi.nlm.nih.gov/gene/")
NCBITAXON = Namespace("http://purl.obolibrary.org/obo/NCBITaxon_")
ENSEMBL = Namespace("http://identifiers.org/ensembl/")
UBERON = Namespace("http://purl.obolibrary.org/obo/UBERON_")
CL = Namespace("http://purl.obolibrary.org/obo/CL_")
GO = Namespace("http://purl.obolibrary.org/obo/GO_")
REACTOME = Namespace("https://reactome.org/content/detail/")
INTERPRO = Namespace("https://www.ebi.ac.uk/interpro/entry/InterPro/")

# Ontology namespaces for characteristics/factors
MONDO = Namespace("http://purl.obolibrary.org/obo/MONDO_")
PATO = Namespace("http://purl.obolibrary.org/obo/PATO_")
HANCESTRO = Namespace("http://purl.obolibrary.org/obo/HANCESTRO_")
EFO = Namespace("http://www.ebi.ac.uk/efo/EFO_")

# Data source namespaces
OSDR = Namespace("https://osdr.nasa.gov/bio/repo/data/studies/")
GEA = Namespace("https://www.ebi.ac.uk/gxa/experiments/")

# All namespaces for binding
NAMESPACES = {
    "biolink": BIOLINK,
    "obo": OBO,
    "spokegenelab": SPOKEGENELAB,
    "ncbigene": NCBIGENE,
    "ncbitaxon": NCBITAXON,
    "ensembl": ENSEMBL,
    "uberon": UBERON,
    "cl": CL,
    "go": GO,
    "reactome": REACTOME,
    "interpro": INTERPRO,
    "mondo": MONDO,
    "pato": PATO,
    "hancestro": HANCESTRO,
    "efo": EFO,
    "osdr": OSDR,
    "gea": GEA,
    "rdf": RDF,
    "rdfs": RDFS,
    "xsd": XSD,
    "owl": OWL,
}

# Property type to XSD datatype mapping
PROPERTY_XSD_TYPES = {
    "string": XSD.string,
    "int": XSD.integer,
    "integer": XSD.integer,
    "float": XSD.float,
    "double": XSD.double,
    "boolean": XSD.boolean,
    "date": XSD.date,
    "datetime": XSD.dateTime,
    "string[]": XSD.string,  # Lists handled separately
}


def get_xsd_type(type_str: str):
    """
    Get XSD datatype for a property type string.

    Args:
        type_str: Type string (e.g., "string", "int", "float")

    Returns:
        XSD datatype URIRef
    """
    return PROPERTY_XSD_TYPES.get(type_str.lower(), XSD.string)


def sanitize_uri_identifier(identifier: str) -> str:
    """
    Sanitize an identifier for use in a URI.

    Removes or encodes characters that are invalid in URIs.

    Args:
        identifier: Raw identifier string

    Returns:
        URI-safe identifier string
    """
    id_str = str(identifier)

    # Replace common problematic characters with underscores
    # These are characters that have special meaning or are invalid in URIs
    replacements = {
        " ": "_",
        ":": "_",
        ">": "_",
        "<": "_",
        '"': "_",
        "'": "_",
        "|": "_",
        "\\": "_",
        "^": "_",
        "`": "_",
        "{": "_",
        "}": "_",
        "[": "_",
        "]": "_",
        "#": "_",
        "%": "_",
        "?": "_",
        "&": "_",
        "=": "_",
        "+": "_",
        "@": "_",
        "$": "_",
        ",": "_",
        ";": "_",
        "!": "_",
        "*": "_",
        "(": "_",
        ")": "_",
    }

    for char, replacement in replacements.items():
        id_str = id_str.replace(char, replacement)

    # Collapse multiple underscores into one
    id_str = re.sub(r"_+", "_", id_str)

    # Remove leading/trailing underscores
    id_str = id_str.strip("_")

    # If empty after cleaning, use a placeholder
    if not id_str:
        id_str = "unknown"

    return id_str


def create_uri(namespace: Namespace, identifier: str) -> URIRef:
    """
    Create a URI from a namespace and identifier.

    Args:
        namespace: RDFLib Namespace
        identifier: Local identifier

    Returns:
        URIRef
    """
    clean_id = sanitize_uri_identifier(identifier)
    return namespace[clean_id]


def get_namespace_for_node_type(node_type: str) -> Namespace:
    """
    Get the appropriate namespace for a node type.

    Args:
        node_type: Node type name

    Returns:
        Namespace for that node type
    """
    namespace_map = {
        "Study": SPOKEGENELAB,
        "Mission": SPOKEGENELAB,
        "Assay": SPOKEGENELAB,
        "MGene": NCBIGENE,
        "Gene": NCBIGENE,
        "Anatomy": UBERON,
        "CellType": CL,
        "PathwayEnrichment": SPOKEGENELAB,
        "GOTerm": GO,
        "ReactomePathway": REACTOME,
        "InterProDomain": INTERPRO,
        "MethylationRegion": SPOKEGENELAB,
        # Characteristic/factor node types - use OBO namespace as base
        "Disease": MONDO,
        "Sex": PATO,
        "DevelopmentalStage": EFO,
        "EthnicGroup": HANCESTRO,
        "OrganismStatus": PATO,
    }
    return namespace_map.get(node_type, SPOKEGENELAB)


# Mapping of URI prefixes to namespace objects
URI_PREFIX_TO_NAMESPACE = {
    "http://purl.obolibrary.org/obo/MONDO_": MONDO,
    "http://purl.obolibrary.org/obo/PATO_": PATO,
    "http://purl.obolibrary.org/obo/HANCESTRO_": HANCESTRO,
    "http://www.ebi.ac.uk/efo/EFO_": EFO,
    "http://purl.obolibrary.org/obo/UBERON_": UBERON,
    "http://purl.obolibrary.org/obo/CL_": CL,
    "http://purl.obolibrary.org/obo/NCBITaxon_": NCBITAXON,
    "http://purl.obolibrary.org/obo/GO_": GO,
}


def parse_ontology_uri(uri: str) -> tuple:
    """
    Parse an ontology URI to extract the namespace prefix and local ID.

    Args:
        uri: Full ontology URI (e.g., 'http://purl.obolibrary.org/obo/MONDO_0004975')

    Returns:
        Tuple of (prefix, local_id) e.g., ('MONDO', '0004975')
        Returns (None, None) if URI format not recognized
    """
    if not uri:
        return None, None

    # Try to match known prefixes
    for prefix_uri, namespace in URI_PREFIX_TO_NAMESPACE.items():
        if uri.startswith(prefix_uri):
            local_id = uri[len(prefix_uri):]
            # Extract prefix name from URI (e.g., 'MONDO' from 'http://.../MONDO_')
            prefix = prefix_uri.rstrip("_").split("/")[-1].split("_")[0]
            return prefix, local_id

    return None, None


def create_uri_from_ontology_uri(ontology_uri: str) -> URIRef:
    """
    Create an RDFLib URIRef directly from an ontology URI.

    Args:
        ontology_uri: Full ontology URI string

    Returns:
        URIRef of the ontology term
    """
    return URIRef(ontology_uri)


def format_identifier_for_namespace(identifier: str, node_type: str) -> str:
    """
    Format an identifier appropriately for its namespace.

    For example, GO:0000001 should become just 0000001 for the GO namespace.

    Args:
        identifier: Raw identifier
        node_type: Node type

    Returns:
        Formatted identifier
    """
    id_str = str(identifier)

    # Handle prefixed identifiers
    if node_type == "GOTerm" and id_str.startswith("GO:"):
        return id_str.replace("GO:", "")
    if node_type == "Anatomy" and id_str.startswith("UBERON:"):
        return id_str.replace("UBERON:", "")
    if node_type == "CellType" and id_str.startswith("CL:"):
        return id_str.replace("CL:", "")

    # Handle Reactome IDs (R-MMU-12345 format)
    if node_type == "ReactomePathway":
        return id_str

    # Handle InterPro IDs (IPR000001 format)
    if node_type == "InterProDomain":
        return id_str

    return id_str


def create_node_uri(node_type: str, identifier: str) -> URIRef:
    """
    Create a URI for a node.

    Args:
        node_type: Type of node
        identifier: Node identifier

    Returns:
        URIRef for the node
    """
    namespace = get_namespace_for_node_type(node_type)
    formatted_id = format_identifier_for_namespace(identifier, node_type)
    return create_uri(namespace, formatted_id)
