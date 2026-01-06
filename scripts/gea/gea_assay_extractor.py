"""
GEA Assay extractor.

Extracts Assay nodes from GEA experiments based on contrasts defined
in the configuration.xml file. Also extracts Anatomy and CellType
proxy nodes from sample metadata.
"""

import hashlib
import json
import sys
from pathlib import Path
from typing import Dict, List, Optional, Set

import pandas as pd

from .gea_parser import GEAExperiment, AssayGroup, Contrast

# Add notebooks directory to path for ontology_mapper
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "notebooks"))

try:
    from ontology_mapper import map_ontology
except ImportError:
    map_ontology = None


def generate_assay_id(experiment_accession: str, contrast_id: str) -> str:
    """
    Generate a unique assay identifier.

    Args:
        experiment_accession: Experiment accession (e.g., "E-GEOD-5305")
        contrast_id: Contrast identifier (e.g., "g1_g3")

    Returns:
        Assay identifier
    """
    return f"{experiment_accession}-{contrast_id}"


def extract_assay_nodes(experiment: GEAExperiment) -> pd.DataFrame:
    """
    Extract Assay nodes from a GEA experiment.

    Each contrast in the experiment becomes an Assay node.

    Args:
        experiment: Parsed GEAExperiment object

    Returns:
        DataFrame with Assay node data
    """
    assays = []

    for contrast in experiment.contrasts:
        # Get assay group details
        ref_group = experiment.assay_groups.get(contrast.reference_group_id)
        test_group = experiment.assay_groups.get(contrast.test_group_id)

        ref_label = ref_group.label if ref_group else ""
        test_label = test_group.label if test_group else ""

        # Parse factors from labels (format: "condition1; condition2")
        ref_factors = [f.strip() for f in ref_label.split(";")] if ref_label else []
        test_factors = [f.strip() for f in test_label.split(";")] if test_label else []

        assay = {
            "identifier": generate_assay_id(experiment.accession, contrast.id),
            "name": contrast.name,
            "study_id": experiment.accession,
            "contrast_id": contrast.id,
            "technology": "DNA microarray",  # GEA experiments are typically microarray
            "measurement": "transcription profiling",
            "array_design": contrast.array_design,
            "reference_group_id": contrast.reference_group_id,
            "reference_group_label": ref_label,
            "test_group_id": contrast.test_group_id,
            "test_group_label": test_label,
            "factors_1": ref_factors,
            "factors_2": test_factors,
        }

        assays.append(assay)

    return pd.DataFrame(assays)


def extract_materials_from_samples(
    experiment: GEAExperiment,
) -> pd.DataFrame:
    """
    Extract unique material types (tissues, cell types) from sample metadata.

    Args:
        experiment: Parsed GEAExperiment object

    Returns:
        DataFrame with unique materials
    """
    materials = set()

    for sample in experiment.samples:
        # Look for cell type or tissue characteristics
        for key, value in sample.items():
            if any(term in key.lower() for term in ["cell type", "tissue", "material", "organ"]):
                if value and pd.notna(value):
                    materials.add(str(value))

    return pd.DataFrame({"material": list(materials)})


def map_materials_to_ontology(
    materials: pd.DataFrame,
    apikey: Optional[str] = None,
) -> pd.DataFrame:
    """
    Map material terms to ontology identifiers (UBERON, CL).

    Args:
        materials: DataFrame with 'material' column
        apikey: BioPortal API key

    Returns:
        DataFrame with material, material_name, material_id columns
    """
    if map_ontology is None or apikey is None:
        # Return unmapped materials
        materials["material_name"] = materials["material"]
        materials["material_id"] = ""
        return materials

    try:
        # Use UBERON for anatomical entities and CL for cell types
        mapped = map_ontology(
            materials,
            column="material",
            ontologies=["UBERON", "CL"],
            apikey=apikey,
        )
        return mapped
    except Exception as e:
        print(f"Warning: Ontology mapping failed: {e}")
        materials["material_name"] = materials["material"]
        materials["material_id"] = ""
        return materials


def assign_materials_to_assays(
    assays: pd.DataFrame,
    experiment: GEAExperiment,
    mapped_materials: Optional[pd.DataFrame] = None,
) -> pd.DataFrame:
    """
    Assign material/tissue information to assays based on sample metadata.

    Args:
        assays: DataFrame with Assay data
        experiment: Parsed GEAExperiment object
        mapped_materials: Optional DataFrame with mapped materials

    Returns:
        DataFrame with material columns added to assays
    """
    # Create sample to material mapping
    sample_materials = {}
    for sample in experiment.samples:
        sample_id = sample.get("Sample", "")
        for key, value in sample.items():
            if any(term in key.lower() for term in ["cell type", "tissue", "material"]):
                if value and pd.notna(value):
                    sample_materials[sample_id] = str(value)
                    break

    # Create material lookup from mapped materials
    material_lookup = {}
    if mapped_materials is not None and not mapped_materials.empty:
        for _, row in mapped_materials.iterrows():
            material_lookup[row["material"]] = {
                "name": row.get("material_name", row["material"]),
                "id": row.get("material_id", ""),
            }

    # Assign materials to assays
    def get_material_for_group(group_id: str) -> tuple:
        group = experiment.assay_groups.get(group_id)
        if not group:
            return "", ""

        # Get material from first sample in group
        for sample_id in group.samples:
            if sample_id in sample_materials:
                material = sample_materials[sample_id]
                if material in material_lookup:
                    return material_lookup[material]["name"], material_lookup[material]["id"]
                return material, ""
        return "", ""

    assays["material_1"], assays["material_id_1"] = zip(
        *assays["reference_group_id"].apply(get_material_for_group)
    )
    assays["material_2"], assays["material_id_2"] = zip(
        *assays["test_group_id"].apply(get_material_for_group)
    )

    return assays


def create_anatomy_nodes(assays: pd.DataFrame) -> pd.DataFrame:
    """
    Create Anatomy proxy nodes from assay material IDs.

    Args:
        assays: DataFrame with Assay data including material_id columns

    Returns:
        DataFrame with Anatomy node data
    """
    anatomy_ids = set()

    for col in ["material_id_1", "material_id_2"]:
        if col in assays.columns:
            ids = assays[col].dropna().unique()
            for id_val in ids:
                if id_val and "UBERON" in str(id_val):
                    anatomy_ids.add(str(id_val))

    if anatomy_ids:
        return pd.DataFrame({
            "identifier": list(anatomy_ids),
        })
    return pd.DataFrame(columns=["identifier"])


def create_celltype_nodes(assays: pd.DataFrame) -> pd.DataFrame:
    """
    Create CellType proxy nodes from assay material IDs.

    Args:
        assays: DataFrame with Assay data including material_id columns

    Returns:
        DataFrame with CellType node data
    """
    celltype_ids = set()

    for col in ["material_id_1", "material_id_2"]:
        if col in assays.columns:
            ids = assays[col].dropna().unique()
            for id_val in ids:
                if id_val and "CL" in str(id_val):
                    celltype_ids.add(str(id_val))

    if celltype_ids:
        return pd.DataFrame({
            "identifier": list(celltype_ids),
        })
    return pd.DataFrame(columns=["identifier"])


def create_study_assay_relationships(assays: pd.DataFrame) -> pd.DataFrame:
    """
    Create Study-PERFORMED-Assay relationships.

    Args:
        assays: DataFrame with Assay data

    Returns:
        DataFrame with relationship data
    """
    if "study_id" not in assays.columns:
        return pd.DataFrame(columns=["from", "to"])

    relationships = assays[["study_id", "identifier"]].copy()
    relationships = relationships.rename(columns={
        "study_id": "from",
        "identifier": "to",
    })
    return relationships.drop_duplicates()


def create_assay_anatomy_relationships(assays: pd.DataFrame) -> pd.DataFrame:
    """
    Create Assay-INVESTIGATED-Anatomy relationships.

    Args:
        assays: DataFrame with Assay data

    Returns:
        DataFrame with relationship data
    """
    relationships = []

    for col in ["material_id_1", "material_id_2"]:
        if col in assays.columns:
            for _, row in assays.iterrows():
                material_id = row.get(col, "")
                if material_id and "UBERON" in str(material_id):
                    relationships.append({
                        "from": row["identifier"],
                        "to": str(material_id),
                    })

    if relationships:
        return pd.DataFrame(relationships).drop_duplicates()
    return pd.DataFrame(columns=["from", "to"])


def create_assay_celltype_relationships(assays: pd.DataFrame) -> pd.DataFrame:
    """
    Create Assay-INVESTIGATED-CellType relationships.

    Args:
        assays: DataFrame with Assay data

    Returns:
        DataFrame with relationship data
    """
    relationships = []

    for col in ["material_id_1", "material_id_2"]:
        if col in assays.columns:
            for _, row in assays.iterrows():
                material_id = row.get(col, "")
                if material_id and "CL" in str(material_id):
                    relationships.append({
                        "from": row["identifier"],
                        "to": str(material_id),
                    })

    if relationships:
        return pd.DataFrame(relationships).drop_duplicates()
    return pd.DataFrame(columns=["from", "to"])


# =============================================================================
# SDRF Characteristic/Factor Extraction Functions
# =============================================================================

# Mapping of SDRF annotation names to node types and their expected ontology prefixes
CHARACTERISTIC_TYPE_MAPPING = {
    "disease": {"node_type": "Disease", "uri_prefixes": ["MONDO", "PATO"]},
    "sex": {"node_type": "Sex", "uri_prefixes": ["PATO"]},
    "developmental stage": {"node_type": "DevelopmentalStage", "uri_prefixes": ["EFO"]},
    "ethnic group": {"node_type": "EthnicGroup", "uri_prefixes": ["HANCESTRO"]},
    "organism status": {"node_type": "OrganismStatus", "uri_prefixes": ["PATO"]},
    "organism part": {"node_type": "Anatomy", "uri_prefixes": ["UBERON"]},
    "cell type": {"node_type": "CellType", "uri_prefixes": ["CL"]},
}

# Characteristics to skip (sample-specific identifiers, not ontology terms)
SKIP_CHARACTERISTICS = {"individual", "age", "organism"}


def extract_sdrf_characteristics(experiment: GEAExperiment) -> Dict[str, List[Dict]]:
    """
    Extract all characteristics and factors from SDRF data.

    Parses the wide-format SDRF data (stored in experiment.samples) to extract
    unique characteristic and factor values with their ontology URIs.

    Args:
        experiment: Parsed GEAExperiment object

    Returns:
        Dictionary mapping characteristic names to lists of unique values with URIs.
        Example: {
            "disease": [
                {"value": "Alzheimers disease", "uri": "http://...MONDO_0004975", "type": "characteristic"},
                {"value": "normal", "uri": "http://...PATO_0000461", "type": "characteristic"},
            ],
            ...
        }
    """
    characteristics: Dict[str, Set[tuple]] = {}

    for sample in experiment.samples:
        for key, value in sample.items():
            if not value or pd.isna(value):
                continue

            # Parse column name to extract annotation type and name
            # Format: "characteristic:disease" or "factor:organism part"
            if ":" not in key or key.endswith("_URI"):
                continue

            parts = key.split(":", 1)
            if len(parts) != 2:
                continue

            annot_type, annot_name = parts
            annot_name = annot_name.strip()

            # Skip certain characteristics
            if annot_name.lower() in SKIP_CHARACTERISTICS:
                continue

            # Only process characteristics and factors
            if annot_type not in ("characteristic", "factor"):
                continue

            # Look for corresponding URI - the URI key is just "{annot_name}_URI"
            # (not prefixed with characteristic: or factor:)
            uri_key = f"{annot_name}_URI"
            uri = sample.get(uri_key, "")
            if pd.isna(uri):
                uri = ""

            # Store unique (value, uri, type) tuples
            if annot_name not in characteristics:
                characteristics[annot_name] = set()
            characteristics[annot_name].add((str(value), str(uri), annot_type))

    # Convert to list format
    result = {}
    for annot_name, values in characteristics.items():
        result[annot_name] = [
            {"value": v, "uri": u, "type": t}
            for v, u, t in values
        ]

    return result


def extract_factor_values_per_assay_group(
    experiment: GEAExperiment,
) -> Dict[str, Dict[str, List[Dict]]]:
    """
    Extract factor values for each assay group.

    Factors vary between assay groups, so we need to map them per group.

    Args:
        experiment: Parsed GEAExperiment object

    Returns:
        Dictionary mapping group_id to factor values.
        Example: {
            "g1": {"disease": [{"value": "normal", "uri": "..."}], ...},
            "g2": {"disease": [{"value": "Alzheimers disease", "uri": "..."}], ...},
        }
    """
    group_factors: Dict[str, Dict[str, Set[tuple]]] = {}

    # Build sample to group mapping
    sample_to_group = {}
    for group_id, group in experiment.assay_groups.items():
        for sample_id in group.samples:
            sample_to_group[sample_id] = group_id

    # Extract factors per group
    for sample in experiment.samples:
        sample_id = sample.get("Sample", "")
        group_id = sample_to_group.get(sample_id)
        if not group_id:
            continue

        if group_id not in group_factors:
            group_factors[group_id] = {}

        for key, value in sample.items():
            if not value or pd.isna(value):
                continue

            # Only process factors
            if not key.startswith("factor:") or key.endswith("_URI"):
                continue

            annot_name = key.split(":", 1)[1].strip()

            # Look for corresponding URI - the URI key is just "{annot_name}_URI"
            uri_key = f"{annot_name}_URI"
            uri = sample.get(uri_key, "")
            if pd.isna(uri):
                uri = ""

            if annot_name not in group_factors[group_id]:
                group_factors[group_id][annot_name] = set()
            group_factors[group_id][annot_name].add((str(value), str(uri)))

    # Convert to list format
    result = {}
    for group_id, factors in group_factors.items():
        result[group_id] = {}
        for factor_name, values in factors.items():
            result[group_id][factor_name] = [
                {"value": v, "uri": u} for v, u in values
            ]

    return result


def extract_characteristics_per_assay_group(
    experiment: GEAExperiment,
) -> Dict[str, Dict[str, List[Dict]]]:
    """
    Extract characteristic values for each assay group.

    Characteristics may vary between assay groups (e.g., cell type, organism part),
    so we need to map them per group to link to the correct Assays.

    Args:
        experiment: Parsed GEAExperiment object

    Returns:
        Dictionary mapping group_id to characteristic values.
        Example: {
            "g1": {"cell type": [{"value": "pyramidal neuron", "uri": "..."}], ...},
            "g2": {"cell type": [{"value": "layer III neuron", "uri": ""}], ...},
        }
    """
    group_characteristics: Dict[str, Dict[str, Set[tuple]]] = {}

    # Build sample to group mapping
    sample_to_group = {}
    for group_id, group in experiment.assay_groups.items():
        for sample_id in group.samples:
            sample_to_group[sample_id] = group_id

    # Extract characteristics per group
    for sample in experiment.samples:
        sample_id = sample.get("Sample", "")
        group_id = sample_to_group.get(sample_id)
        if not group_id:
            continue

        if group_id not in group_characteristics:
            group_characteristics[group_id] = {}

        for key, value in sample.items():
            if not value or pd.isna(value):
                continue

            # Only process characteristics
            if not key.startswith("characteristic:") or key.endswith("_URI"):
                continue

            annot_name = key.split(":", 1)[1].strip()

            # Skip certain characteristics (sample-specific identifiers)
            if annot_name.lower() in SKIP_CHARACTERISTICS:
                continue

            # Look for corresponding URI
            uri_key = f"{annot_name}_URI"
            uri = sample.get(uri_key, "")
            if pd.isna(uri):
                uri = ""

            if annot_name not in group_characteristics[group_id]:
                group_characteristics[group_id][annot_name] = set()
            group_characteristics[group_id][annot_name].add((str(value), str(uri)))

    # Convert to list format
    result = {}
    for group_id, chars in group_characteristics.items():
        result[group_id] = {}
        for char_name, values in chars.items():
            result[group_id][char_name] = [
                {"value": v, "uri": u} for v, u in values
            ]

    return result


def create_characteristic_nodes(
    characteristics: Dict[str, List[Dict]],
    node_type: str,
    characteristic_name: str,
) -> pd.DataFrame:
    """
    Create nodes for a specific characteristic type.

    Args:
        characteristics: Output from extract_sdrf_characteristics
        node_type: The node type to create (e.g., "Disease", "Sex")
        characteristic_name: The SDRF annotation name (e.g., "disease", "sex")

    Returns:
        DataFrame with node data (identifier, name, uri columns)
    """
    if characteristic_name not in characteristics:
        return pd.DataFrame(columns=["identifier", "name", "uri"])

    values = characteristics[characteristic_name]
    nodes = []

    for item in values:
        uri = item.get("uri", "")
        value = item.get("value", "")

        # Handle cases where URI field contains multiple URIs separated by spaces
        if uri and " " in uri and uri.startswith("http"):
            # Split multiple URIs and create a node for each
            uris = uri.split()
            for single_uri in uris:
                if single_uri.startswith("http"):
                    nodes.append({
                        "identifier": single_uri,
                        "name": value,
                        "uri": single_uri,
                    })
        elif uri:
            # Single URI - use as identifier
            nodes.append({
                "identifier": uri,
                "name": value,
                "uri": uri,
            })
        else:
            # Generate a spokegenelab URI for values without ontology URIs
            # Sanitize the value for use in URI
            safe_value = value.replace(" ", "_").replace("'", "")
            identifier = f"https://spoke.ucsf.edu/genelab/{node_type}/{safe_value}"
            nodes.append({
                "identifier": identifier,
                "name": value,
                "uri": uri,
            })

    if nodes:
        df = pd.DataFrame(nodes).drop_duplicates(subset=["identifier"])
        return df
    return pd.DataFrame(columns=["identifier", "name", "uri"])


def create_disease_nodes(characteristics: Dict[str, List[Dict]]) -> pd.DataFrame:
    """Create Disease nodes from characteristics."""
    return create_characteristic_nodes(characteristics, "Disease", "disease")


def create_sex_nodes(characteristics: Dict[str, List[Dict]]) -> pd.DataFrame:
    """Create Sex nodes from characteristics."""
    return create_characteristic_nodes(characteristics, "Sex", "sex")


def create_developmental_stage_nodes(characteristics: Dict[str, List[Dict]]) -> pd.DataFrame:
    """Create DevelopmentalStage nodes from characteristics."""
    return create_characteristic_nodes(characteristics, "DevelopmentalStage", "developmental stage")


def create_ethnic_group_nodes(characteristics: Dict[str, List[Dict]]) -> pd.DataFrame:
    """Create EthnicGroup nodes from characteristics."""
    return create_characteristic_nodes(characteristics, "EthnicGroup", "ethnic group")


def create_organism_status_nodes(characteristics: Dict[str, List[Dict]]) -> pd.DataFrame:
    """Create OrganismStatus nodes from characteristics."""
    return create_characteristic_nodes(characteristics, "OrganismStatus", "organism status")


def create_anatomy_nodes_from_sdrf(characteristics: Dict[str, List[Dict]]) -> pd.DataFrame:
    """Create Anatomy nodes from SDRF organism part characteristics."""
    return create_characteristic_nodes(characteristics, "Anatomy", "organism part")


def create_celltype_nodes_from_sdrf(characteristics: Dict[str, List[Dict]]) -> pd.DataFrame:
    """Create CellType nodes from SDRF cell type characteristics."""
    return create_characteristic_nodes(characteristics, "CellType", "cell type")


def create_study_disease_relationships(
    study_id: str,
    characteristics: Dict[str, List[Dict]],
) -> pd.DataFrame:
    """
    Create Study-HAS_DISEASE-Disease relationships.

    Args:
        study_id: Study identifier (experiment accession)
        characteristics: Output from extract_sdrf_characteristics

    Returns:
        DataFrame with relationship data
    """
    if "disease" not in characteristics:
        return pd.DataFrame(columns=["from", "to"])

    relationships = []
    for item in characteristics["disease"]:
        uri = item.get("uri", "")
        value = item.get("value", "")

        if uri:
            to_id = uri
        else:
            safe_value = value.replace(" ", "_").replace("'", "")
            to_id = f"https://spoke.ucsf.edu/genelab/Disease/{safe_value}"

        relationships.append({"from": study_id, "to": to_id})

    if relationships:
        return pd.DataFrame(relationships).drop_duplicates()
    return pd.DataFrame(columns=["from", "to"])


def create_study_characteristic_relationships(
    study_id: str,
    characteristics: Dict[str, List[Dict]],
    characteristic_name: str,
    node_type: str,
) -> pd.DataFrame:
    """
    Create Study-HAS_*-* relationships for a characteristic type.

    Args:
        study_id: Study identifier
        characteristics: Output from extract_sdrf_characteristics
        characteristic_name: SDRF annotation name
        node_type: Node type for generating fallback URIs

    Returns:
        DataFrame with relationship data
    """
    if characteristic_name not in characteristics:
        return pd.DataFrame(columns=["from", "to"])

    relationships = []
    for item in characteristics[characteristic_name]:
        uri = item.get("uri", "")
        value = item.get("value", "")

        # Handle cases where URI field contains multiple URIs separated by spaces
        if uri and " " in uri and uri.startswith("http"):
            uris = uri.split()
            for single_uri in uris:
                if single_uri.startswith("http"):
                    relationships.append({"from": study_id, "to": single_uri})
        elif uri:
            relationships.append({"from": study_id, "to": uri})
        else:
            safe_value = value.replace(" ", "_").replace("'", "")
            to_id = f"https://spoke.ucsf.edu/genelab/{node_type}/{safe_value}"
            relationships.append({"from": study_id, "to": to_id})

    if relationships:
        return pd.DataFrame(relationships).drop_duplicates()
    return pd.DataFrame(columns=["from", "to"])


def create_study_sex_relationships(
    study_id: str,
    characteristics: Dict[str, List[Dict]],
) -> pd.DataFrame:
    """Create Study-HAS_SEX-Sex relationships."""
    return create_study_characteristic_relationships(
        study_id, characteristics, "sex", "Sex"
    )


def create_study_developmental_stage_relationships(
    study_id: str,
    characteristics: Dict[str, List[Dict]],
) -> pd.DataFrame:
    """Create Study-HAS_DEVELOPMENTAL_STAGE-DevelopmentalStage relationships."""
    return create_study_characteristic_relationships(
        study_id, characteristics, "developmental stage", "DevelopmentalStage"
    )


def create_study_ethnic_group_relationships(
    study_id: str,
    characteristics: Dict[str, List[Dict]],
) -> pd.DataFrame:
    """Create Study-HAS_ETHNIC_GROUP-EthnicGroup relationships."""
    return create_study_characteristic_relationships(
        study_id, characteristics, "ethnic group", "EthnicGroup"
    )


def create_study_organism_status_relationships(
    study_id: str,
    characteristics: Dict[str, List[Dict]],
) -> pd.DataFrame:
    """Create Study-HAS_ORGANISM_STATUS-OrganismStatus relationships."""
    return create_study_characteristic_relationships(
        study_id, characteristics, "organism status", "OrganismStatus"
    )


def create_assay_factor_relationships(
    assays: pd.DataFrame,
    group_factors: Dict[str, Dict[str, List[Dict]]],
    factor_name: str,
    node_type: str,
) -> pd.DataFrame:
    """
    Create Assay-HAS_INPUT-* relationships for factors.

    Factors are linked to assays since they vary between experimental groups.

    Args:
        assays: DataFrame with Assay data
        group_factors: Output from extract_factor_values_per_assay_group
        factor_name: Factor annotation name (e.g., "disease", "organism part")
        node_type: Node type for the factor

    Returns:
        DataFrame with relationship data
    """
    relationships = []

    for _, assay in assays.iterrows():
        assay_id = assay["identifier"]

        # Get factors for both reference and test groups
        for group_col in ["reference_group_id", "test_group_id"]:
            group_id = assay.get(group_col, "")
            if not group_id or group_id not in group_factors:
                continue

            factors = group_factors[group_id].get(factor_name, [])
            for item in factors:
                uri = item.get("uri", "")
                value = item.get("value", "")

                # Handle cases where URI field contains multiple URIs separated by spaces
                if uri and " " in uri and uri.startswith("http"):
                    uris = uri.split()
                    for single_uri in uris:
                        if single_uri.startswith("http"):
                            relationships.append({"from": assay_id, "to": single_uri})
                elif uri:
                    relationships.append({"from": assay_id, "to": uri})
                else:
                    safe_value = value.replace(" ", "_").replace("'", "")
                    to_id = f"https://spoke.ucsf.edu/genelab/{node_type}/{safe_value}"
                    relationships.append({"from": assay_id, "to": to_id})

    if relationships:
        return pd.DataFrame(relationships).drop_duplicates()
    return pd.DataFrame(columns=["from", "to"])


def create_assay_anatomy_relationships_from_sdrf(
    assays: pd.DataFrame,
    group_factors: Dict[str, Dict[str, List[Dict]]],
) -> pd.DataFrame:
    """Create Assay-INVESTIGATED-Anatomy relationships from SDRF factors."""
    return create_assay_factor_relationships(
        assays, group_factors, "organism part", "Anatomy"
    )


def create_assay_celltype_relationships_from_sdrf(
    assays: pd.DataFrame,
    group_factors: Dict[str, Dict[str, List[Dict]]],
) -> pd.DataFrame:
    """Create Assay-INVESTIGATED-CellType relationships from SDRF factors."""
    return create_assay_factor_relationships(
        assays, group_factors, "cell type", "CellType"
    )


def create_assay_disease_relationships(
    assays: pd.DataFrame,
    group_factors: Dict[str, Dict[str, List[Dict]]],
) -> pd.DataFrame:
    """Create Assay-HAS_INPUT-Disease relationships from SDRF factors."""
    return create_assay_factor_relationships(
        assays, group_factors, "disease", "Disease"
    )


def create_study_anatomy_relationships(
    study_id: str,
    characteristics: Dict[str, List[Dict]],
) -> pd.DataFrame:
    """Create Study-HAS_ANATOMY-Anatomy relationships for anatomy characteristics."""
    return create_study_characteristic_relationships(
        study_id, characteristics, "organism part", "Anatomy"
    )


def create_study_celltype_relationships(
    study_id: str,
    characteristics: Dict[str, List[Dict]],
) -> pd.DataFrame:
    """Create Study-HAS_CELLTYPE-CellType relationships for cell type characteristics."""
    return create_study_characteristic_relationships(
        study_id, characteristics, "cell type", "CellType"
    )


def create_assay_characteristic_relationships(
    assays: pd.DataFrame,
    group_characteristics: Dict[str, Dict[str, List[Dict]]],
    characteristic_name: str,
    node_type: str,
) -> pd.DataFrame:
    """
    Create Assay-HAS_INPUT-* relationships for characteristics.

    Characteristics are linked to assays based on the samples in each assay group.

    Args:
        assays: DataFrame with Assay data
        group_characteristics: Output from extract_characteristics_per_assay_group
        characteristic_name: Characteristic annotation name (e.g., "cell type", "organism part")
        node_type: Node type for the characteristic

    Returns:
        DataFrame with relationship data
    """
    relationships = []

    for _, assay in assays.iterrows():
        assay_id = assay["identifier"]

        # Get characteristics for both reference and test groups
        for group_col in ["reference_group_id", "test_group_id"]:
            group_id = assay.get(group_col, "")
            if not group_id or group_id not in group_characteristics:
                continue

            chars = group_characteristics[group_id].get(characteristic_name, [])
            for item in chars:
                uri = item.get("uri", "")
                value = item.get("value", "")

                # Handle cases where URI field contains multiple URIs separated by spaces
                if uri and " " in uri and uri.startswith("http"):
                    uris = uri.split()
                    for single_uri in uris:
                        if single_uri.startswith("http"):
                            relationships.append({"from": assay_id, "to": single_uri})
                elif uri:
                    relationships.append({"from": assay_id, "to": uri})
                else:
                    safe_value = value.replace(" ", "_").replace("'", "")
                    to_id = f"https://spoke.ucsf.edu/genelab/{node_type}/{safe_value}"
                    relationships.append({"from": assay_id, "to": to_id})

    if relationships:
        return pd.DataFrame(relationships).drop_duplicates()
    return pd.DataFrame(columns=["from", "to"])


def create_assay_celltype_relationships_from_characteristics(
    assays: pd.DataFrame,
    group_characteristics: Dict[str, Dict[str, List[Dict]]],
) -> pd.DataFrame:
    """Create Assay-HAS_INPUT-CellType relationships from SDRF characteristics."""
    return create_assay_characteristic_relationships(
        assays, group_characteristics, "cell type", "CellType"
    )


def create_assay_anatomy_relationships_from_characteristics(
    assays: pd.DataFrame,
    group_characteristics: Dict[str, Dict[str, List[Dict]]],
) -> pd.DataFrame:
    """Create Assay-HAS_INPUT-Anatomy relationships from SDRF characteristics."""
    return create_assay_characteristic_relationships(
        assays, group_characteristics, "organism part", "Anatomy"
    )
