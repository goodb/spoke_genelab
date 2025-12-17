"""
Configuration management for SPOKE-GeneLab pipelines.

Handles loading configuration from .env files and YAML config files,
and provides a unified Config object for all pipeline modules.
"""

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional

from dotenv import load_dotenv

# Default paths relative to project root
DEFAULT_DATA_DIR = "data"
DEFAULT_OUTPUT_DIR = "output"


@dataclass
class Config:
    """Configuration settings for SPOKE-GeneLab pipelines."""

    # Knowledge graph version
    kg_version: str = "v0.2.0"

    # Directory paths
    project_root: Path = field(default_factory=lambda: Path(__file__).parent.parent.parent)
    data_dir: Path = field(default_factory=lambda: Path(DEFAULT_DATA_DIR))
    output_dir: Path = field(default_factory=lambda: Path(DEFAULT_OUTPUT_DIR))

    # Neo4j settings
    neo4j_home: Optional[str] = None
    neo4j_bin: Optional[str] = None
    neo4j_database: str = "spoke-genelab"
    neo4j_metadata_dir: Optional[Path] = None
    neo4j_data_dir: Optional[Path] = None

    # API settings
    osdr_api_root: str = "https://visualization.osdr.nasa.gov/biodata/api/v2/"
    bioportal_api_key: Optional[str] = None

    # Processing settings
    adj_p_value_threshold: float = 0.1
    q_value_threshold: float = 0.1

    # Supported organisms (taxid -> name)
    organisms: Dict[str, str] = field(default_factory=lambda: {
        "10090": "Mus musculus",
        "10116": "Rattus norvegicus",
        "7955": "Danio rerio",
        "9606": "Homo sapiens",
        "6239": "Caenorhabditis elegans",
        "7227": "Drosophila melanogaster",
    })

    # Technology types to process
    technology_types: List[str] = field(default_factory=lambda: [
        "RNA Sequencing (RNA-Seq)",
        "DNA microarray",
    ])

    # RDF output settings
    rdf_output_format: str = "turtle"
    rdf_base_uri: str = "https://spoke.ucsf.edu/genelab/"

    def __post_init__(self):
        """Initialize derived paths after dataclass initialization."""
        # Ensure paths are Path objects
        if isinstance(self.project_root, str):
            self.project_root = Path(self.project_root)
        if isinstance(self.data_dir, str):
            self.data_dir = Path(self.data_dir)
        if isinstance(self.output_dir, str):
            self.output_dir = Path(self.output_dir)

        # Set up Neo4j paths if not specified
        if self.neo4j_metadata_dir is None:
            self.neo4j_metadata_dir = self.project_root / "kg" / self.kg_version / "metadata"
        if self.neo4j_data_dir is None:
            self.neo4j_data_dir = self.project_root / "kg" / self.kg_version / "data"

    @property
    def node_output_dir(self) -> Path:
        """Directory for node CSV files."""
        return self.neo4j_data_dir / "nodes"

    @property
    def relationship_output_dir(self) -> Path:
        """Directory for relationship CSV files."""
        return self.neo4j_data_dir / "relationships"

    @property
    def rdf_output_dir(self) -> Path:
        """Directory for RDF output files."""
        return self.output_dir / "rdf"

    def ensure_directories(self):
        """Create all required output directories."""
        dirs = [
            self.data_dir,
            self.output_dir,
            self.node_output_dir,
            self.relationship_output_dir,
            self.rdf_output_dir,
        ]
        for d in dirs:
            d.mkdir(parents=True, exist_ok=True)

    def validate(self):
        """Validate configuration settings."""
        errors = []

        if self.bioportal_api_key is None:
            print("Warning: BIOPORTAL_API_KEY not set. Ontology mapping will be disabled.")

        if not self.neo4j_metadata_dir.exists():
            errors.append(f"NEO4J_METADATA directory does not exist: {self.neo4j_metadata_dir}")

        if errors:
            raise ValueError("Configuration validation failed:\n" + "\n".join(errors))

        return True


def load_config(
    env_file: Optional[str] = None,
    config_file: Optional[str] = None
) -> Config:
    """
    Load configuration from environment variables and optional config file.

    Args:
        env_file: Path to .env file. If None, looks for .env in project root.
        config_file: Path to YAML config file (optional).

    Returns:
        Config object with all settings.
    """
    # Find project root
    project_root = Path(__file__).parent.parent.parent

    # Load .env file
    if env_file is None:
        env_file = project_root / ".env"

    if Path(env_file).exists():
        load_dotenv(env_file, override=True)

    # Build config from environment
    config = Config(
        project_root=project_root,
        kg_version=os.getenv("KG_VERSION", "v0.2.0"),
        data_dir=Path(os.getenv("DATA_DIR", project_root / "data")),
        output_dir=Path(os.getenv("OUTPUT_DIR", project_root / "output")),
        neo4j_home=os.getenv("NEO4J_HOME"),
        neo4j_bin=os.getenv("NEO4J_BIN"),
        neo4j_database=os.getenv("NEO4J_DATABASE", "spoke-genelab"),
        bioportal_api_key=os.getenv("BIOPORTAL_API_KEY"),
        adj_p_value_threshold=float(os.getenv("ADJ_P_VALUE_THRESHOLD", "0.1")),
        q_value_threshold=float(os.getenv("Q_VALUE_THRESHOLD", "0.1")),
    )

    # Override with NEO4J_METADATA and NEO4J_DATA if set
    if os.getenv("NEO4J_METADATA"):
        config.neo4j_metadata_dir = Path(os.getenv("NEO4J_METADATA"))
    if os.getenv("NEO4J_DATA"):
        config.neo4j_data_dir = Path(os.getenv("NEO4J_DATA"))

    # Load YAML config file if provided
    if config_file and Path(config_file).exists():
        import yaml
        with open(config_file) as f:
            yaml_config = yaml.safe_load(f)

        # Override config with YAML values
        for key, value in yaml_config.items():
            if hasattr(config, key):
                setattr(config, key, value)

    return config


def get_output_dirs(config: Config) -> Dict[str, Path]:
    """
    Get dictionary of output directories for pipeline use.

    Args:
        config: Config object

    Returns:
        Dict with 'nodes', 'relationships', and 'rdf' keys
    """
    config.ensure_directories()
    return {
        "nodes": config.node_output_dir,
        "relationships": config.relationship_output_dir,
        "rdf": config.rdf_output_dir,
    }
