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
