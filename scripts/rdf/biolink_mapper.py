"""
Biolink Model mapping utilities.

Maps SPOKE-GeneLab node types and relationship types to Biolink Model classes
and predicates for semantic interoperability.
"""

from typing import Dict, Optional, Tuple

from rdflib import URIRef

from .rdf_config import BIOLINK, SPOKEGENELAB


# Node type to Biolink class mapping
BIOLINK_NODE_CLASSES = {
    "Mission": BIOLINK.Activity,
    "Study": BIOLINK.Study,
    "Assay": BIOLINK.Assay,
    "MGene": BIOLINK.Gene,
    "Gene": BIOLINK.Gene,
    "Anatomy": BIOLINK.AnatomicalEntity,
    "CellType": BIOLINK.Cell,
    "MethylationRegion": BIOLINK.GenomicEntity,
    "PathwayEnrichment": BIOLINK.Association,
    "GOTerm": BIOLINK.BiologicalProcess,  # Default, may be overridden
    "ReactomePathway": BIOLINK.Pathway,
    "InterProDomain": BIOLINK.ProteinDomain,
    # Characteristic/factor node types
    "Disease": BIOLINK.Disease,
    "Sex": BIOLINK.BiologicalSex,
    "DevelopmentalStage": BIOLINK.LifeStage,
    "EthnicGroup": BIOLINK.PopulationOfIndividualOrganisms,
    "OrganismStatus": BIOLINK.Attribute,
}

# GO category to Biolink class mapping
GO_CATEGORY_CLASSES = {
    "biological_process": BIOLINK.BiologicalProcess,
    "molecular_function": BIOLINK.MolecularActivity,
    "cellular_component": BIOLINK.CellularComponent,
    "BP": BIOLINK.BiologicalProcess,
    "MF": BIOLINK.MolecularActivity,
    "CC": BIOLINK.CellularComponent,
}

# Relationship type to Biolink predicate mapping
BIOLINK_PREDICATES = {
    # Mission/Study relationships
    "CONDUCTED": BIOLINK.has_part,
    "Mission-CONDUCTED_MIcS-Study": BIOLINK.has_part,

    # Study/Assay relationships
    "PERFORMED": BIOLINK.has_output,
    "Study-PERFORMED_SpAS-Assay": BIOLINK.has_output,

    # Assay/Gene relationships (differential expression)
    "MEASURED_DIFFERENTIAL_EXPRESSION": BIOLINK.affects_expression_of,
    "Assay-MEASURED_DIFFERENTIAL_EXPRESSION_ASmMG-MGene": BIOLINK.affects_expression_of,

    # Assay/Methylation relationships
    "MEASURED_DIFFERENTIAL_METHYLATION": BIOLINK.affects,
    "Assay-MEASURED_DIFFERENTIAL_METHYLATION_ASmMR-MethylationRegion": BIOLINK.affects,

    # Ortholog relationships
    "IS_ORTHOLOG": BIOLINK.orthologous_to,
    "MGene-IS_ORTHOLOG_MGiG-Gene": BIOLINK.orthologous_to,

    # Investigation relationships
    "INVESTIGATED": BIOLINK.has_input,
    "Assay-INVESTIGATED_ASiA-Anatomy": BIOLINK.has_input,
    "Assay-INVESTIGATED_ASiCT-CellType": BIOLINK.has_input,

    # Methylation relationships
    "METHYLATED_IN": BIOLINK.overlaps,
    "MGene-METHYLATED_IN_MGmMR-MethylationRegion": BIOLINK.overlaps,

    # Enrichment relationships
    "HAS_ENRICHMENT": BIOLINK.has_output,
    "Assay-HAS_ENRICHMENT_AShPE-PathwayEnrichment": BIOLINK.has_output,

    # Pathway/Term relationships
    "ENRICHES": BIOLINK.participates_in,
    "PathwayEnrichment-ENRICHES_PEeGO-GOTerm": BIOLINK.participates_in,
    "PathwayEnrichment-ENRICHES_PEeRP-ReactomePathway": BIOLINK.participates_in,
    "PathwayEnrichment-ENRICHES_PEeIP-InterProDomain": BIOLINK.participates_in,

    # Study-Characteristic relationships
    "HAS_DISEASE": BIOLINK.studies,
    "Study-HAS_DISEASE-Disease": BIOLINK.studies,
    "HAS_SEX": BIOLINK.has_attribute,
    "Study-HAS_SEX-Sex": BIOLINK.has_attribute,
    "HAS_DEVELOPMENTAL_STAGE": BIOLINK.has_attribute,
    "Study-HAS_DEVELOPMENTAL_STAGE-DevelopmentalStage": BIOLINK.has_attribute,
    "HAS_ETHNIC_GROUP": BIOLINK.has_attribute,
    "Study-HAS_ETHNIC_GROUP-EthnicGroup": BIOLINK.has_attribute,
    "HAS_ORGANISM_STATUS": BIOLINK.has_attribute,
    "Study-HAS_ORGANISM_STATUS-OrganismStatus": BIOLINK.has_attribute,
    # Assay-Characteristic relationships (characteristics vary per assay group)
    "Assay-HAS_ATTRIBUTE-CellType": BIOLINK.has_attribute,
    "Assay-HAS_ATTRIBUTE-Anatomy": BIOLINK.has_attribute,
    "Assay-HAS_INPUT-Disease": BIOLINK.has_input,
}


def get_biolink_class(node_type: str, category: Optional[str] = None) -> URIRef:
    """
    Get the Biolink class URI for a node type.

    Args:
        node_type: Node type name (e.g., "Study", "Gene")
        category: Optional category for GO terms

    Returns:
        Biolink class URIRef
    """
    # Handle GO terms with category
    if node_type == "GOTerm" and category:
        return GO_CATEGORY_CLASSES.get(category, BIOLINK.BiologicalProcess)

    return BIOLINK_NODE_CLASSES.get(node_type, BIOLINK.NamedThing)


def get_biolink_predicate(relationship_type: str) -> URIRef:
    """
    Get the Biolink predicate URI for a relationship type.

    Args:
        relationship_type: Relationship type name

    Returns:
        Biolink predicate URIRef
    """
    # Try exact match first
    if relationship_type in BIOLINK_PREDICATES:
        return BIOLINK_PREDICATES[relationship_type]

    # Try to extract base relationship type
    for key, value in BIOLINK_PREDICATES.items():
        if key in relationship_type:
            return value

    # Default to related_to
    return BIOLINK.related_to


def get_property_predicate(property_name: str, node_type: str) -> URIRef:
    """
    Get the predicate URI for a node property.

    Args:
        property_name: Property name
        node_type: Node type (for context)

    Returns:
        Predicate URIRef
    """
    # Standard Biolink properties
    biolink_properties = {
        "identifier": BIOLINK.id,
        "name": BIOLINK.name,
        "description": BIOLINK.description,
        "symbol": BIOLINK.symbol,
        "category": BIOLINK.category,
        "organism": BIOLINK.in_taxon,
        "taxonomy": BIOLINK.in_taxon,
    }

    if property_name in biolink_properties:
        return biolink_properties[property_name]

    # Use custom namespace for non-standard properties
    return SPOKEGENELAB[property_name]


def parse_relationship_type(rel_type: str) -> Tuple[str, str, str]:
    """
    Parse a relationship type string to extract source, predicate, and target.

    Example: "Study-PERFORMED_SpAS-Assay" -> ("Study", "PERFORMED_SpAS", "Assay")

    Args:
        rel_type: Full relationship type string

    Returns:
        Tuple of (source_type, predicate_name, target_type)
    """
    parts = rel_type.split("-")

    if len(parts) >= 3:
        source = parts[0]
        target = parts[-1]
        predicate = "-".join(parts[1:-1])
        return source, predicate, target

    # Fallback
    return "", rel_type, ""


def is_reified_relationship(rel_type: str) -> bool:
    """
    Check if a relationship type should be reified (has properties).

    Relationships with properties like log2fc and p-value should be
    represented as Biolink Associations.

    Args:
        rel_type: Relationship type name

    Returns:
        True if the relationship should be reified
    """
    reified_types = {
        "MEASURED_DIFFERENTIAL_EXPRESSION",
        "MEASURED_DIFFERENTIAL_METHYLATION",
        "Assay-MEASURED_DIFFERENTIAL_EXPRESSION_ASmMG-MGene",
        "Assay-MEASURED_DIFFERENTIAL_METHYLATION_ASmMR-MethylationRegion",
    }
    return any(rt in rel_type for rt in reified_types)


def get_association_class(rel_type: str) -> URIRef:
    """
    Get the Biolink Association class for a reified relationship.

    Args:
        rel_type: Relationship type name

    Returns:
        Biolink Association class URIRef
    """
    if "EXPRESSION" in rel_type.upper():
        return BIOLINK.GeneExpressionMixin
    if "METHYLATION" in rel_type.upper():
        return BIOLINK.Association

    return BIOLINK.Association
