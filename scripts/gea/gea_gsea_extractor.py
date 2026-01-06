"""
GEA GSEA extractor.

Extracts PathwayEnrichment nodes from GEA Gene Set Enrichment Analysis results.
Creates nodes for enriched GO terms, Reactome pathways, and InterPro domains.
"""

import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd

from .gea_parser import (
    GEAExperiment,
    parse_gsea_file,
    get_gsea_enrichment_type,
    get_contrast_from_gsea_filename,
)


def extract_gsea_results(
    experiment: GEAExperiment,
    p_value_threshold: float = 0.01,
    max_terms_per_type: int = 20,
) -> pd.DataFrame:
    """
    Extract all GSEA enrichment results from a GEA experiment.

    Args:
        experiment: Parsed GEAExperiment object
        p_value_threshold: Adjusted p-value threshold for significance (default 0.01)
        max_terms_per_type: Maximum number of enriched terms per enrichment type per contrast (default 20)

    Returns:
        DataFrame with GSEA results including enrichment type and contrast
    """
    all_results = []

    for gsea_file in experiment.gsea_files:
        filename = Path(gsea_file).name
        enrichment_type = get_gsea_enrichment_type(filename)
        contrast_id = get_contrast_from_gsea_filename(filename)

        if enrichment_type == "unknown" or not contrast_id:
            continue

        try:
            df = parse_gsea_file(gsea_file)

            # Add metadata columns
            df["enrichment_type"] = enrichment_type
            df["contrast_id"] = contrast_id
            df["experiment_accession"] = experiment.accession
            df["source_file"] = filename

            # Filter by p-value threshold
            if "adj_p_value" in df.columns:
                df = df[df["adj_p_value"] <= p_value_threshold]
                # Sort by p-value and limit to top N terms
                df = df.sort_values("adj_p_value").head(max_terms_per_type)

            if not df.empty:
                all_results.append(df)

        except Exception as e:
            print(f"Warning: Failed to parse {gsea_file}: {e}")

    if all_results:
        return pd.concat(all_results, ignore_index=True)

    return pd.DataFrame()


def create_pathway_enrichment_nodes(
    gsea_results: pd.DataFrame,
) -> pd.DataFrame:
    """
    Create PathwayEnrichment nodes from GSEA results.

    Args:
        gsea_results: DataFrame with parsed GSEA results

    Returns:
        DataFrame with PathwayEnrichment node data
    """
    if gsea_results.empty:
        return pd.DataFrame(columns=[
            "identifier", "name", "enrichment_type", "contrast_id",
            "term_id", "p_value", "adj_p_value", "effect_size",
            "genes_total", "genes_significant",
        ])

    enrichments = []

    for _, row in gsea_results.iterrows():
        term_id = row.get("term_id", row.get("Accession", ""))
        term_name = row.get("term_name", row.get("Term", ""))

        # Generate unique identifier
        identifier = f"{row['experiment_accession']}_{row['contrast_id']}_{term_id}"

        enrichment = {
            "identifier": identifier,
            "name": term_name,
            "enrichment_type": row["enrichment_type"],
            "contrast_id": row["contrast_id"],
            "term_id": term_id,
            "experiment_accession": row["experiment_accession"],
        }

        # Add numeric fields if present
        for field in ["p_value", "adj_p_value", "effect_size", "genes_total", "significant_in_set"]:
            if field in row and pd.notna(row[field]):
                if field == "significant_in_set":
                    enrichment["genes_significant"] = int(row[field])
                else:
                    enrichment[field] = row[field]

        enrichments.append(enrichment)

    return pd.DataFrame(enrichments)


def create_go_term_nodes(gsea_results: pd.DataFrame) -> pd.DataFrame:
    """
    Create GOTerm proxy nodes from GO enrichment results.

    Args:
        gsea_results: DataFrame with GSEA results

    Returns:
        DataFrame with GOTerm node data
    """
    go_results = gsea_results[gsea_results["enrichment_type"] == "go"]

    if go_results.empty:
        return pd.DataFrame(columns=["identifier", "name", "category"])

    # Extract unique GO terms
    go_terms = []
    seen = set()

    for _, row in go_results.iterrows():
        term_id = row.get("term_id", row.get("Accession", ""))
        if term_id and term_id not in seen:
            seen.add(term_id)
            term_name = row.get("term_name", row.get("Term", ""))

            # Try to determine GO category from term ID prefix or name
            category = determine_go_category(term_id, term_name)

            go_terms.append({
                "identifier": term_id,
                "name": term_name,
                "category": category,
            })

    return pd.DataFrame(go_terms)


def determine_go_category(term_id: str, term_name: str) -> str:
    """
    Attempt to determine GO category (BP, MF, CC) from term.

    Args:
        term_id: GO term identifier
        term_name: GO term name

    Returns:
        Category string or empty string if unknown
    """
    # Common keywords that hint at category
    bp_keywords = ["process", "regulation", "response", "development", "metabolism"]
    mf_keywords = ["activity", "binding", "transporter", "receptor", "enzyme"]
    cc_keywords = ["membrane", "complex", "organelle", "cytoplasm", "nucleus"]

    name_lower = term_name.lower()

    for kw in bp_keywords:
        if kw in name_lower:
            return "biological_process"
    for kw in mf_keywords:
        if kw in name_lower:
            return "molecular_function"
    for kw in cc_keywords:
        if kw in name_lower:
            return "cellular_component"

    return ""


def create_reactome_pathway_nodes(gsea_results: pd.DataFrame) -> pd.DataFrame:
    """
    Create ReactomePathway proxy nodes from Reactome enrichment results.

    Args:
        gsea_results: DataFrame with GSEA results

    Returns:
        DataFrame with ReactomePathway node data
    """
    reactome_results = gsea_results[gsea_results["enrichment_type"] == "reactome"]

    if reactome_results.empty:
        return pd.DataFrame(columns=["identifier", "name"])

    # Extract unique Reactome pathways
    pathways = []
    seen = set()

    for _, row in reactome_results.iterrows():
        term_id = row.get("term_id", row.get("Accession", ""))
        if term_id and term_id not in seen:
            seen.add(term_id)
            pathways.append({
                "identifier": term_id,
                "name": row.get("term_name", row.get("Term", "")),
            })

    return pd.DataFrame(pathways)


def create_interpro_domain_nodes(gsea_results: pd.DataFrame) -> pd.DataFrame:
    """
    Create InterProDomain proxy nodes from InterPro enrichment results.

    Args:
        gsea_results: DataFrame with GSEA results

    Returns:
        DataFrame with InterProDomain node data
    """
    interpro_results = gsea_results[gsea_results["enrichment_type"] == "interpro"]

    if interpro_results.empty:
        return pd.DataFrame(columns=["identifier", "name"])

    # Extract unique InterPro domains
    domains = []
    seen = set()

    for _, row in interpro_results.iterrows():
        term_id = row.get("term_id", row.get("Accession", ""))
        if term_id and term_id not in seen:
            seen.add(term_id)
            domains.append({
                "identifier": term_id,
                "name": row.get("term_name", row.get("Term", "")),
            })

    return pd.DataFrame(domains)


def create_assay_enrichment_relationships(
    enrichment_nodes: pd.DataFrame,
) -> pd.DataFrame:
    """
    Create Assay-HAS_ENRICHMENT-PathwayEnrichment relationships.

    Args:
        enrichment_nodes: DataFrame with PathwayEnrichment node data

    Returns:
        DataFrame with relationship data
    """
    if enrichment_nodes.empty:
        return pd.DataFrame(columns=["from", "to"])

    relationships = []

    for _, row in enrichment_nodes.iterrows():
        # Assay ID is experiment_accession-contrast_id
        assay_id = f"{row['experiment_accession']}-{row['contrast_id']}"
        relationships.append({
            "from": assay_id,
            "to": row["identifier"],
        })

    return pd.DataFrame(relationships).drop_duplicates()


def create_enrichment_term_relationships(
    enrichment_nodes: pd.DataFrame,
) -> Dict[str, pd.DataFrame]:
    """
    Create relationships from PathwayEnrichment to term proxy nodes.

    Args:
        enrichment_nodes: DataFrame with PathwayEnrichment node data

    Returns:
        Dictionary with relationship DataFrames for each enrichment type
    """
    relationships = {
        "PathwayEnrichment-ENRICHES_PEeGO-GOTerm": [],
        "PathwayEnrichment-ENRICHES_PEeRP-ReactomePathway": [],
        "PathwayEnrichment-ENRICHES_PEeIP-InterProDomain": [],
    }

    rel_type_map = {
        "go": "PathwayEnrichment-ENRICHES_PEeGO-GOTerm",
        "reactome": "PathwayEnrichment-ENRICHES_PEeRP-ReactomePathway",
        "interpro": "PathwayEnrichment-ENRICHES_PEeIP-InterProDomain",
    }

    for _, row in enrichment_nodes.iterrows():
        rel_type = rel_type_map.get(row["enrichment_type"])
        if rel_type and row.get("term_id"):
            relationships[rel_type].append({
                "from": row["identifier"],
                "to": row["term_id"],
            })

    # Convert to DataFrames
    result = {}
    for rel_type, rels in relationships.items():
        if rels:
            result[rel_type] = pd.DataFrame(rels).drop_duplicates()
        else:
            result[rel_type] = pd.DataFrame(columns=["from", "to"])

    return result


def get_upregulated_genes_by_pathway(
    experiment: GEAExperiment,
    pathway_id: str,
    log2fc_threshold: float = 1.0,
    p_value_threshold: float = 0.05,
) -> pd.DataFrame:
    """
    Get genes that are upregulated and enriched in a specific pathway.

    This is useful for queries like "what genes are upregulated in this dataset
    that are involved in [pathway X]?"

    Args:
        experiment: Parsed GEAExperiment object
        pathway_id: Pathway/GO term identifier to filter by
        log2fc_threshold: Minimum log2 fold change for upregulation
        p_value_threshold: Maximum p-value for significance

    Returns:
        DataFrame with upregulated genes in the pathway
    """
    from .gea_gene_extractor import extract_differential_expression

    # Get DE data
    de_data = extract_differential_expression(experiment, p_value_threshold)

    # Filter for upregulated genes
    upregulated = de_data[de_data["log2fc"] >= log2fc_threshold]

    # TODO: Cross-reference with pathway gene lists from GSEA
    # This would require parsing the .gsea_list.tsv files

    return upregulated


def extract_all_gsea_data(
    experiment: GEAExperiment,
    p_value_threshold: float = 0.01,
    max_terms_per_type: int = 20,
) -> Dict[str, pd.DataFrame]:
    """
    Extract all GSEA-related nodes and relationships from an experiment.

    Args:
        experiment: Parsed GEAExperiment object
        p_value_threshold: Adjusted p-value threshold (default 0.01)
        max_terms_per_type: Maximum enriched terms per type per contrast (default 20)

    Returns:
        Dictionary with all nodes and relationships DataFrames
    """
    # Extract GSEA results
    gsea_results = extract_gsea_results(experiment, p_value_threshold, max_terms_per_type)

    if gsea_results.empty:
        return {
            "PathwayEnrichment": pd.DataFrame(),
            "GOTerm": pd.DataFrame(),
            "ReactomePathway": pd.DataFrame(),
            "InterProDomain": pd.DataFrame(),
            "Assay-HAS_ENRICHMENT_AShPE-PathwayEnrichment": pd.DataFrame(),
            "PathwayEnrichment-ENRICHES_PEeGO-GOTerm": pd.DataFrame(),
            "PathwayEnrichment-ENRICHES_PEeRP-ReactomePathway": pd.DataFrame(),
            "PathwayEnrichment-ENRICHES_PEeIP-InterProDomain": pd.DataFrame(),
        }

    # Create nodes
    enrichment_nodes = create_pathway_enrichment_nodes(gsea_results)
    go_nodes = create_go_term_nodes(gsea_results)
    reactome_nodes = create_reactome_pathway_nodes(gsea_results)
    interpro_nodes = create_interpro_domain_nodes(gsea_results)

    # Create relationships
    assay_enrichment_rels = create_assay_enrichment_relationships(enrichment_nodes)
    term_rels = create_enrichment_term_relationships(enrichment_nodes)

    return {
        "PathwayEnrichment": enrichment_nodes,
        "GOTerm": go_nodes,
        "ReactomePathway": reactome_nodes,
        "InterProDomain": interpro_nodes,
        "Assay-HAS_ENRICHMENT_AShPE-PathwayEnrichment": assay_enrichment_rels,
        **term_rels,
    }
