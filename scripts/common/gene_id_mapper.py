"""
Gene ID mapping utilities.

Maps between Ensembl and NCBI gene identifiers using HGNC data.
"""

import pandas as pd
from typing import Dict, Optional
from functools import lru_cache


HGNC_URL = "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"


@lru_cache(maxsize=1)
def load_hgnc_mapping() -> pd.DataFrame:
    """
    Load HGNC gene mapping data.

    Returns:
        DataFrame with ensembl_gene_id, entrez_id, and symbol columns
    """
    print("Loading HGNC gene ID mappings...")
    try:
        df = pd.read_csv(
            HGNC_URL,
            sep="\t",
            usecols=["symbol", "ensembl_gene_id", "entrez_id"],
            dtype=str,
        )
        # Drop rows without both IDs
        df = df.dropna(subset=["ensembl_gene_id", "entrez_id"])
        print(f"  Loaded {len(df)} gene mappings")
        return df
    except Exception as e:
        print(f"Warning: Failed to load HGNC mappings: {e}")
        return pd.DataFrame(columns=["symbol", "ensembl_gene_id", "entrez_id"])


def get_ensembl_to_ncbi_map() -> Dict[str, str]:
    """
    Get a dictionary mapping Ensembl gene IDs to NCBI gene IDs.

    Returns:
        Dict mapping Ensembl ID -> NCBI gene ID
    """
    df = load_hgnc_mapping()
    return dict(zip(df["ensembl_gene_id"], df["entrez_id"]))


def get_ncbi_to_ensembl_map() -> Dict[str, str]:
    """
    Get a dictionary mapping NCBI gene IDs to Ensembl gene IDs.

    Returns:
        Dict mapping NCBI gene ID -> Ensembl ID
    """
    df = load_hgnc_mapping()
    return dict(zip(df["entrez_id"], df["ensembl_gene_id"]))


def map_ensembl_to_ncbi(ensembl_ids: pd.Series) -> pd.DataFrame:
    """
    Map Ensembl gene IDs to NCBI gene IDs.

    Args:
        ensembl_ids: Series of Ensembl gene IDs

    Returns:
        DataFrame with columns: ensembl_id, ncbi_gene_id
    """
    mapping = get_ensembl_to_ncbi_map()

    result = pd.DataFrame({"ensembl_id": ensembl_ids})
    result["ncbi_gene_id"] = result["ensembl_id"].map(mapping)

    mapped_count = result["ncbi_gene_id"].notna().sum()
    total_count = len(result)
    print(f"  Mapped {mapped_count}/{total_count} Ensembl IDs to NCBI gene IDs")

    return result


def add_ncbi_gene_ids(
    df: pd.DataFrame,
    ensembl_col: str = "identifier",
    ncbi_col: str = "ncbi_gene_id",
) -> pd.DataFrame:
    """
    Add NCBI gene IDs to a DataFrame with Ensembl IDs.

    Args:
        df: DataFrame with Ensembl gene IDs
        ensembl_col: Column name containing Ensembl IDs
        ncbi_col: Column name for NCBI gene IDs (will be added)

    Returns:
        DataFrame with ncbi_gene_id column added
    """
    mapping = get_ensembl_to_ncbi_map()

    df = df.copy()
    df[ncbi_col] = df[ensembl_col].map(mapping)

    return df
