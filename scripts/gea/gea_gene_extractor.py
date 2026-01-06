"""
GEA Gene extractor.

Extracts MGene (model organism gene) and Gene (human ortholog) nodes
from GEA differential expression analytics files.
"""

import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd

from .gea_parser import GEAExperiment, parse_analytics_file
from ..common.gene_id_mapper import add_ncbi_gene_ids

# Add notebooks directory to path for ortholog_mapper
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "notebooks"))

try:
    from ortholog_mapper import map_orthologs
except ImportError:
    map_orthologs = None


def extract_genes_from_analytics(
    analytics_file: str,
    organism: str = "",
    taxonomy: str = "",
) -> pd.DataFrame:
    """
    Extract gene information from a GEA analytics file.

    Args:
        analytics_file: Path to analytics TSV file
        organism: Organism name
        taxonomy: NCBI taxonomy ID

    Returns:
        DataFrame with gene information including NCBI and Ensembl IDs
    """
    df = parse_analytics_file(analytics_file)

    # Extract gene columns
    genes = df[["gene_id", "gene_name"]].drop_duplicates()

    # Filter out NA/empty gene IDs
    genes = genes[genes["gene_id"].notna() & (genes["gene_id"] != "")]

    # Add organism info
    genes["organism"] = organism
    genes["taxonomy"] = taxonomy

    # Rename to match schema
    genes = genes.rename(columns={
        "gene_id": "ensembl_id",
        "gene_name": "symbol",
    })

    # Add empty name column (GEA doesn't typically have full gene names)
    genes["name"] = ""

    # For human genes, map Ensembl IDs to NCBI gene IDs
    if taxonomy == "9606":
        genes = add_ncbi_gene_ids(genes, ensembl_col="ensembl_id", ncbi_col="ncbi_gene_id")
        # Use NCBI gene ID as identifier when available, otherwise use Ensembl ID
        genes["identifier"] = genes["ncbi_gene_id"].fillna(genes["ensembl_id"])
        genes["id_source"] = genes["ncbi_gene_id"].apply(
            lambda x: "NCBIGene" if pd.notna(x) else "Ensembl"
        )
    else:
        # For non-human, use Ensembl ID as identifier
        genes["identifier"] = genes["ensembl_id"]
        genes["ncbi_gene_id"] = None
        genes["id_source"] = "Ensembl"

    return genes[["identifier", "ensembl_id", "ncbi_gene_id", "symbol", "name", "organism", "taxonomy", "id_source"]]


def extract_genes_from_experiment(experiment: GEAExperiment) -> pd.DataFrame:
    """
    Extract all unique genes from a GEA experiment.

    Args:
        experiment: Parsed GEAExperiment object

    Returns:
        DataFrame with MGene node data
    """
    all_genes = []

    taxonomy = experiment.taxonomy_id
    if not taxonomy:
        from .gea_parser import get_organism_taxonomy
        taxonomy = get_organism_taxonomy(experiment.organism)

    for analytics_file in experiment.analytics_files:
        genes = extract_genes_from_analytics(
            analytics_file,
            organism=experiment.organism,
            taxonomy=taxonomy,
        )
        all_genes.append(genes)

    if all_genes:
        combined = pd.concat(all_genes, ignore_index=True)
        # Deduplicate by gene ID
        combined = combined.drop_duplicates(subset="identifier")
        return combined

    return pd.DataFrame(columns=["identifier", "symbol", "name", "organism", "taxonomy"])


def map_to_human_orthologs(
    mgenes: pd.DataFrame,
    ortholog_dbs: Optional[List[str]] = None,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Map model organism genes to human orthologs.

    Args:
        mgenes: DataFrame with MGene data (must have 'identifier' and 'taxonomy' columns)
        ortholog_dbs: List of ortholog databases to use (default: ["JAX", "Ensembl"])

    Returns:
        Tuple of (Gene nodes DataFrame, MGene-IS_ORTHOLOG-Gene relationships DataFrame)
    """
    if ortholog_dbs is None:
        ortholog_dbs = ["JAX", "Ensembl"]

    if map_orthologs is None:
        print("Warning: ortholog_mapper not available. Skipping ortholog mapping.")
        return pd.DataFrame(columns=["identifier"]), pd.DataFrame(columns=["from", "to"])

    # Filter to non-human genes
    human_taxid = "9606"
    non_human = mgenes[mgenes["taxonomy"] != human_taxid].copy()

    if non_human.empty:
        # All genes are human, no mapping needed
        genes = mgenes[["identifier"]].drop_duplicates()
        return genes, pd.DataFrame(columns=["from", "to"])

    # Map orthologs using the existing mapper
    try:
        mapped = map_orthologs(non_human, ortholog_dbs)

        # Extract unique human gene IDs
        human_genes = mapped["human_entrezid"].dropna().unique()
        gene_nodes = pd.DataFrame({"identifier": human_genes})

        # Create relationships
        relationships = mapped[["identifier", "human_entrezid"]].dropna()
        relationships = relationships.rename(columns={
            "identifier": "from",
            "human_entrezid": "to",
        })
        relationships = relationships.drop_duplicates()

        return gene_nodes, relationships

    except Exception as e:
        print(f"Warning: Ortholog mapping failed: {e}")
        return pd.DataFrame(columns=["identifier"]), pd.DataFrame(columns=["from", "to"])


def create_mgene_nodes(experiment: GEAExperiment) -> pd.DataFrame:
    """
    Create MGene (model organism gene) nodes from a GEA experiment.

    Args:
        experiment: Parsed GEAExperiment object

    Returns:
        DataFrame with MGene node data
    """
    return extract_genes_from_experiment(experiment)


def create_gene_and_ortholog_data(
    mgenes: pd.DataFrame,
    ortholog_dbs: Optional[List[str]] = None,
) -> Dict[str, pd.DataFrame]:
    """
    Create Gene nodes and ortholog relationships from MGene data.

    Args:
        mgenes: DataFrame with MGene data
        ortholog_dbs: Ortholog databases to use

    Returns:
        Dictionary with 'Gene' nodes and 'MGene-IS_ORTHOLOG-Gene' relationships
    """
    genes, orthologs = map_to_human_orthologs(mgenes, ortholog_dbs)

    return {
        "Gene": genes,
        "MGene-IS_ORTHOLOG_MGiG-Gene": orthologs,
    }


def extract_differential_expression(
    experiment: GEAExperiment,
    p_value_threshold: float = 0.01,
    max_genes_per_assay: int = 200,
) -> pd.DataFrame:
    """
    Extract differential expression data for creating Assay-MGene relationships.

    Args:
        experiment: Parsed GEAExperiment object
        p_value_threshold: Adjusted p-value threshold for significance (default 0.01)
        max_genes_per_assay: Maximum number of DE genes to include per assay (default 200)

    Returns:
        DataFrame with columns: from, to, log2fc, adj_p_value
    """
    all_de = []

    # Get taxonomy for gene ID mapping
    taxonomy = experiment.taxonomy_id
    if not taxonomy:
        from .gea_parser import get_organism_taxonomy
        taxonomy = get_organism_taxonomy(experiment.organism)

    # Build Ensembl -> NCBI gene ID mapping for human genes
    ensembl_to_identifier = {}
    if taxonomy == "9606":
        from ..common.gene_id_mapper import get_ensembl_to_ncbi_map
        ncbi_map = get_ensembl_to_ncbi_map()
        # Map to NCBI ID if available, otherwise use Ensembl ID
        ensembl_to_identifier = {k: v if v else k for k, v in ncbi_map.items()}

    for analytics_file in experiment.analytics_files:
        df = parse_analytics_file(analytics_file)

        # Get contrast columns
        for contrast in experiment.contrasts:
            contrast_id = contrast.id

            # Column names in analytics file
            pval_col = f"{contrast_id}.p-value"
            log2fc_col = f"{contrast_id}.log2foldchange"

            if pval_col not in df.columns or log2fc_col not in df.columns:
                continue

            # Select and filter
            de_data = df[["gene_id", pval_col, log2fc_col]].copy()
            de_data = de_data.rename(columns={
                pval_col: "p_value",
                log2fc_col: "log2fc",
            })

            # Filter by p-value threshold
            de_data = de_data[de_data["p_value"].notna()]
            de_data = de_data[de_data["p_value"] <= p_value_threshold]

            if de_data.empty:
                continue

            # Sort by p-value and limit to top N genes
            de_data = de_data.sort_values("p_value").head(max_genes_per_assay)

            # Create assay ID from experiment accession and contrast
            assay_id = f"{experiment.accession}-{contrast_id}"
            de_data["assay_id"] = assay_id

            # Map gene IDs to NCBI IDs for human genes
            if taxonomy == "9606" and ensembl_to_identifier:
                # Use NCBI ID if available, otherwise keep Ensembl ID
                de_data["to"] = de_data["gene_id"].map(
                    lambda x: ensembl_to_identifier.get(x, x)
                )
            else:
                de_data["to"] = de_data["gene_id"]

            de_data["from"] = assay_id

            all_de.append(de_data[["from", "to", "log2fc", "p_value"]])

    if all_de:
        combined = pd.concat(all_de, ignore_index=True)
        # Rename p_value to adj_p_value to match schema
        combined = combined.rename(columns={"p_value": "adj_p_value"})
        return combined

    return pd.DataFrame(columns=["from", "to", "log2fc", "adj_p_value"])
