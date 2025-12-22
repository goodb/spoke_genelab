"""
GEA FTP Downloader - Incremental download from EBI Gene Expression Atlas.

Downloads experiment data from the EBI FTP server, filtering to only include
files needed by the GEA processing pipeline. Optionally generates RDF output
after each successful download.

Usage:
    python -m scripts.gea.gea_downloader [options]

Options:
    --data-dir PATH       Override GEA_DATA_DIR from .env
    --prefix PREFIX       Experiment prefix filter (default: E-GEOD)
    --experiment E-XXX    Download single experiment only
    --max-size MB         Max file size in MB (default: 10)
    --dry-run             List files without downloading
    --rdf                 Generate RDF after each successful download
    --output-dir PATH     Output directory for RDF files (default: ./output)

Examples:
    python -m scripts.gea.gea_downloader --prefix E-GEOD
    python -m scripts.gea.gea_downloader --experiment E-GEOD-5305 --rdf
    python -m scripts.gea.gea_downloader --prefix E-MTAB --max-size 5 --rdf
"""

import argparse
import ftplib
import json
import logging
import os
import re
import sys
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Optional, Set, List, Dict, Any

from dotenv import load_dotenv

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger(__name__)

# FTP server configuration
FTP_HOST = "ftp.ebi.ac.uk"
FTP_PATH = "/pub/databases/microarray/data/atlas/experiments"
DEFAULT_PREFIX = "E-GEOD"
DEFAULT_MAX_SIZE_MB = 10.0

# File extensions to skip (images, plots, etc.)
SKIP_EXTENSIONS = {
    ".png",
    ".eps",
    ".jpg",
    ".jpeg",
    ".gif",
    ".pdf",
    ".bedgraph",
    ".rdata",
}

# Patterns in filenames to skip (unprocessed versions)
SKIP_PATTERNS = [".undecorated", ".unrounded"]

# File patterns we want to download (from gea_parser.py)
WANTED_PATTERNS = [
    r"\.idf\.txt$",  # Investigation metadata
    r"\.condensed-sdrf\.tsv$",  # Sample metadata
    r"\.sdrf\.tsv$",  # Alternative sample format
    r"-configuration\.xml$",  # Assay groups/contrasts
    r"-analytics\.tsv$",  # Differential expression
    r"\.go\.gsea\.tsv$",  # GO enrichment
    r"\.reactome\.gsea\.tsv$",  # Reactome enrichment
    r"\.interpro\.gsea\.tsv$",  # InterPro enrichment
    r"normalized-expressions\.tsv$",  # Expression data
]


@dataclass
class DownloadState:
    """Track download progress for resume capability."""

    completed_experiments: Set[str] = field(default_factory=set)
    rdf_completed: Set[str] = field(default_factory=set)
    in_progress: Optional[str] = None
    completed_files: Set[str] = field(default_factory=set)
    failed_files: Set[str] = field(default_factory=set)
    rdf_failed: Set[str] = field(default_factory=set)

    def save(self, path: Path) -> None:
        """Save state to JSON file."""
        data = {
            "completed_experiments": list(self.completed_experiments),
            "rdf_completed": list(self.rdf_completed),
            "in_progress": self.in_progress,
            "completed_files": list(self.completed_files),
            "failed_files": list(self.failed_files),
            "rdf_failed": list(self.rdf_failed),
        }
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "w") as f:
            json.dump(data, f, indent=2)

    @classmethod
    def load(cls, path: Path) -> "DownloadState":
        """Load state from JSON file."""
        if not path.exists():
            return cls()

        try:
            with open(path) as f:
                data = json.load(f)
            return cls(
                completed_experiments=set(data.get("completed_experiments", [])),
                rdf_completed=set(data.get("rdf_completed", [])),
                in_progress=data.get("in_progress"),
                completed_files=set(data.get("completed_files", [])),
                failed_files=set(data.get("failed_files", [])),
                rdf_failed=set(data.get("rdf_failed", [])),
            )
        except (json.JSONDecodeError, KeyError) as e:
            logger.warning(f"Could not load state file: {e}. Starting fresh.")
            return cls()


@dataclass
class ExperimentResult:
    """Result of processing a single experiment."""
    accession: str
    download_success: bool = False
    download_files: int = 0
    rdf_success: Optional[bool] = None  # None if RDF not attempted
    rdf_triples: int = 0
    error: Optional[str] = None


class GEADownloader:
    """Downloads GEA experiment data from EBI FTP server."""

    def __init__(
        self,
        data_dir: Path,
        prefix: str = DEFAULT_PREFIX,
        max_size_mb: float = DEFAULT_MAX_SIZE_MB,
        dry_run: bool = False,
        generate_rdf: bool = False,
        output_dir: Optional[Path] = None,
    ):
        self.data_dir = Path(data_dir)
        self.prefix = prefix
        self.max_size = int(max_size_mb * 1024 * 1024)
        self.dry_run = dry_run
        self.generate_rdf = generate_rdf
        self.output_dir = Path(output_dir) if output_dir else Path("./output")
        self.state_file = self.data_dir / ".download_state.json"
        self.ftp: Optional[ftplib.FTP] = None
        self.state = DownloadState.load(self.state_file)

        # Compile wanted patterns for efficiency
        self._wanted_regex = [re.compile(p, re.IGNORECASE) for p in WANTED_PATTERNS]

    def connect(self) -> None:
        """Establish FTP connection."""
        logger.info(f"Connecting to {FTP_HOST}...")
        self.ftp = ftplib.FTP(FTP_HOST)
        self.ftp.login()  # Anonymous login
        self.ftp.cwd(FTP_PATH)
        logger.info("Connected successfully.")

    def disconnect(self) -> None:
        """Close FTP connection."""
        if self.ftp:
            try:
                self.ftp.quit()
            except Exception:
                pass
            self.ftp = None

    def list_experiments(self) -> List[str]:
        """List all experiment directories matching the prefix."""
        if not self.ftp:
            raise RuntimeError("Not connected to FTP server")

        logger.info(f"Listing experiments with prefix '{self.prefix}'...")
        experiments = []

        # Get directory listing
        lines = []
        self.ftp.dir(lines.append)

        for line in lines:
            # Parse directory listing (format varies but name is usually last)
            parts = line.split()
            if not parts:
                continue

            name = parts[-1]

            # Check if it's a directory and matches prefix
            if line.startswith("d") and name.startswith(self.prefix):
                experiments.append(name)

        experiments.sort()
        logger.info(f"Found {len(experiments)} experiments matching '{self.prefix}'")
        return experiments

    def list_files(self, experiment: str) -> List[Dict[str, Any]]:
        """List files in an experiment directory with their sizes."""
        if not self.ftp:
            raise RuntimeError("Not connected to FTP server")

        files = []

        try:
            # Change to experiment directory
            self.ftp.cwd(f"{FTP_PATH}/{experiment}")

            # Get file listing with sizes using MLSD if available
            try:
                for name, facts in self.ftp.mlsd():
                    if facts.get("type") == "file":
                        size = int(facts.get("size", 0))
                        files.append({"name": name, "size": size})
            except ftplib.error_perm:
                # Fallback to LIST if MLSD not supported
                lines = []
                self.ftp.dir(lines.append)

                for line in lines:
                    parts = line.split()
                    if len(parts) >= 5 and not line.startswith("d"):
                        # Standard Unix ls -l format: size is 5th column, name is last
                        try:
                            size = int(parts[4])
                            name = parts[-1]
                            files.append({"name": name, "size": size})
                        except (ValueError, IndexError):
                            continue

        except ftplib.error_perm as e:
            logger.warning(f"Could not list files in {experiment}: {e}")

        return files

    def should_download(self, filename: str, size: int) -> bool:
        """Determine if a file should be downloaded based on filters."""
        filename_lower = filename.lower()

        # Check file extension
        for ext in SKIP_EXTENSIONS:
            if filename_lower.endswith(ext):
                return False

        # Check skip patterns
        for pattern in SKIP_PATTERNS:
            if pattern in filename_lower:
                return False

        # Check file size
        if size > self.max_size:
            return False

        # Check if it matches any wanted pattern
        for regex in self._wanted_regex:
            if regex.search(filename):
                return True

        # If no wanted pattern matches, skip
        return False

    def download_file(
        self, experiment: str, filename: str, size: int
    ) -> bool:
        """Download a single file with resume support."""
        if not self.ftp:
            raise RuntimeError("Not connected to FTP server")

        # Create local directory with -gea suffix (expected by gea_pipeline.py)
        local_dir_name = f"{experiment}-gea" if not experiment.endswith("-gea") else experiment
        local_dir = self.data_dir / local_dir_name
        local_dir.mkdir(parents=True, exist_ok=True)
        local_path = local_dir / filename

        remote_path = f"{FTP_PATH}/{experiment}/{filename}"
        file_key = f"{experiment}/{filename}"

        # Check if already completed
        if file_key in self.state.completed_files:
            logger.debug(f"Skipping {file_key} (already downloaded)")
            return True

        # Check if local file exists with correct size
        if local_path.exists():
            local_size = local_path.stat().st_size
            if local_size == size:
                logger.debug(f"Skipping {file_key} (already exists with correct size)")
                self.state.completed_files.add(file_key)
                return True
            elif local_size < size:
                # Resume partial download
                logger.info(f"Resuming {file_key} from byte {local_size}")
                mode = "ab"
                rest = local_size
            else:
                # Local file is larger (corrupted?), re-download
                logger.warning(f"Re-downloading {file_key} (local size > remote)")
                mode = "wb"
                rest = None
        else:
            mode = "wb"
            rest = None

        if self.dry_run:
            logger.info(f"[DRY-RUN] Would download: {file_key} ({size:,} bytes)")
            return True

        # Download the file
        try:
            self.ftp.cwd(f"{FTP_PATH}/{experiment}")

            with open(local_path, mode) as f:
                if rest:
                    self.ftp.retrbinary(f"RETR {filename}", f.write, rest=rest)
                else:
                    self.ftp.retrbinary(f"RETR {filename}", f.write)

            # Verify size
            if local_path.stat().st_size == size:
                self.state.completed_files.add(file_key)
                logger.info(f"Downloaded: {file_key}")
                return True
            else:
                logger.warning(f"Size mismatch after download: {file_key}")
                return False

        except ftplib.error_perm as e:
            logger.error(f"FTP error downloading {file_key}: {e}")
            self.state.failed_files.add(file_key)
            return False
        except Exception as e:
            logger.error(f"Error downloading {file_key}: {e}")
            self.state.failed_files.add(file_key)
            return False

    def download_experiment(self, accession: str) -> tuple[bool, int]:
        """Download all needed files for a single experiment.

        Returns:
            Tuple of (success, files_downloaded)
        """
        # Check if already completed
        if accession in self.state.completed_experiments:
            logger.debug(f"Skipping download for {accession} (already completed)")
            # Count existing files
            local_dir = self.data_dir / f"{accession}-gea"
            existing_files = len(list(local_dir.glob("*"))) if local_dir.exists() else 0
            return True, existing_files

        # Mark as in progress
        self.state.in_progress = accession
        self.state.save(self.state_file)

        # List files in experiment
        files = self.list_files(accession)

        # Filter and download
        downloaded = 0
        skipped = 0
        failed = 0

        for file_info in files:
            filename = file_info["name"]
            size = file_info["size"]

            if self.should_download(filename, size):
                if self.download_file(accession, filename, size):
                    downloaded += 1
                else:
                    failed += 1
            else:
                skipped += 1

        # Mark as completed if no failures
        if failed == 0:
            self.state.completed_experiments.add(accession)
            self.state.in_progress = None

        self.state.save(self.state_file)
        return failed == 0, downloaded

    def generate_rdf_for_experiment(self, accession: str) -> tuple[bool, int]:
        """Generate RDF for a single experiment.

        Returns:
            Tuple of (success, num_triples)
        """
        # Check if already completed
        if accession in self.state.rdf_completed:
            logger.debug(f"Skipping RDF for {accession} (already completed)")
            return True, 0

        local_dir = self.data_dir / f"{accession}-gea"
        if not local_dir.exists():
            logger.warning(f"Cannot generate RDF: directory not found: {local_dir}")
            return False, 0

        try:
            # Import here to avoid circular imports
            from .gea_pipeline import process_gea_experiment
            from ..rdf.turtle_writer import write_graph_to_turtle

            # Process the experiment
            result = process_gea_experiment(
                local_dir,
                p_value_threshold=0.1,
                include_gsea=True,
                include_orthologs=True,
            )

            if result.errors:
                for error in result.errors:
                    logger.warning(f"  Pipeline warning: {error}")

            # Write RDF - use accession-specific file
            self.output_dir.mkdir(parents=True, exist_ok=True)
            rdf_file = self.output_dir / f"{accession}.ttl"

            writer = write_graph_to_turtle(
                result.nodes,
                result.relationships,
                rdf_file,
            )
            num_triples = writer.get_triple_count()

            # Mark as completed
            self.state.rdf_completed.add(accession)
            self.state.save(self.state_file)

            return True, num_triples

        except Exception as e:
            logger.error(f"RDF generation failed for {accession}: {e}")
            self.state.rdf_failed.add(accession)
            self.state.save(self.state_file)
            return False, 0

    def process_experiment(self, accession: str) -> ExperimentResult:
        """Process a single experiment: download and optionally generate RDF.

        Returns:
            ExperimentResult with download and RDF status
        """
        result = ExperimentResult(accession=accession)

        # Download
        try:
            download_success, files_downloaded = self.download_experiment(accession)
            result.download_success = download_success
            result.download_files = files_downloaded
        except Exception as e:
            result.download_success = False
            result.error = f"Download error: {e}"
            return result

        # Generate RDF if requested and download succeeded
        if self.generate_rdf and result.download_success and not self.dry_run:
            try:
                rdf_success, num_triples = self.generate_rdf_for_experiment(accession)
                result.rdf_success = rdf_success
                result.rdf_triples = num_triples
            except Exception as e:
                result.rdf_success = False
                result.error = f"RDF error: {e}"

        return result

    def run(self, single_experiment: Optional[str] = None) -> List[ExperimentResult]:
        """Main download orchestration.

        Returns:
            List of ExperimentResult for each processed experiment
        """
        results: List[ExperimentResult] = []

        # Ensure data directory exists
        self.data_dir.mkdir(parents=True, exist_ok=True)

        try:
            self.connect()

            if single_experiment:
                # Download single experiment
                experiments = [single_experiment]
            else:
                # Get all experiments matching prefix
                experiments = self.list_experiments()

            # Filter out already completed (both download and RDF if applicable)
            if self.generate_rdf:
                remaining = [
                    e for e in experiments
                    if e not in self.state.completed_experiments
                    or e not in self.state.rdf_completed
                ]
            else:
                remaining = [e for e in experiments if e not in self.state.completed_experiments]

            logger.info(
                f"Processing {len(remaining)} experiments "
                f"({len(self.state.completed_experiments)} downloads completed"
                + (f", {len(self.state.rdf_completed)} RDFs completed)" if self.generate_rdf else ")")
            )

            # Process each experiment
            for i, experiment in enumerate(remaining, 1):
                logger.info(f"=" * 60)
                logger.info(f"[{i}/{len(remaining)}] Processing: {experiment}")

                try:
                    result = self.process_experiment(experiment)
                    results.append(result)

                    # Log result summary
                    if result.download_success:
                        dl_msg = f"Download: SUCCESS ({result.download_files} files)"
                    else:
                        dl_msg = f"Download: FAILED"

                    if result.rdf_success is None:
                        rdf_msg = ""
                    elif result.rdf_success:
                        rdf_msg = f" | RDF: SUCCESS ({result.rdf_triples:,} triples)"
                    else:
                        rdf_msg = f" | RDF: FAILED"

                    logger.info(f"  {dl_msg}{rdf_msg}")

                    if result.error:
                        logger.error(f"  Error: {result.error}")

                except Exception as e:
                    logger.error(f"  Unexpected error processing {experiment}: {e}")
                    results.append(ExperimentResult(
                        accession=experiment,
                        download_success=False,
                        error=str(e),
                    ))
                    # Reconnect in case of connection issues
                    try:
                        self.disconnect()
                        self.connect()
                    except Exception:
                        pass

            # Final summary
            logger.info(f"=" * 60)
            logger.info("SUMMARY")
            dl_success = sum(1 for r in results if r.download_success)
            dl_failed = sum(1 for r in results if not r.download_success)
            total_files = sum(r.download_files for r in results)
            logger.info(f"Downloads: {dl_success} succeeded, {dl_failed} failed, {total_files} total files")

            if self.generate_rdf:
                rdf_success = sum(1 for r in results if r.rdf_success is True)
                rdf_failed = sum(1 for r in results if r.rdf_success is False)
                total_triples = sum(r.rdf_triples for r in results)
                logger.info(f"RDF: {rdf_success} succeeded, {rdf_failed} failed, {total_triples:,} total triples")

        finally:
            self.disconnect()

        return results


def main():
    """CLI entry point."""
    # Load .env file
    project_root = Path(__file__).parent.parent.parent
    env_file = project_root / ".env"
    if env_file.exists():
        load_dotenv(env_file)

    # Default directories from environment
    default_data_dir = os.getenv(
        "GEA_DATA_DIR", "/Users/bgood/Documents/Scripps/gea"
    )
    default_rdf_output_dir = os.getenv(
        "RDF_OUTPUT_DIR", "/Users/bgood/Documents/Scripps/rdf"
    )

    parser = argparse.ArgumentParser(
        description="Download GEA experiment data from EBI FTP server.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    %(prog)s --prefix E-GEOD
    %(prog)s --experiment E-GEOD-5305 --rdf
    %(prog)s --prefix E-MTAB --max-size 5 --rdf
    %(prog)s --dry-run --prefix E-GEOD
        """,
    )

    parser.add_argument(
        "--data-dir",
        type=Path,
        default=default_data_dir,
        help=f"Directory to store downloaded data (default: {default_data_dir})",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=default_rdf_output_dir,
        help=f"Output directory for RDF files (default: {default_rdf_output_dir})",
    )
    parser.add_argument(
        "--prefix",
        type=str,
        default=DEFAULT_PREFIX,
        help=f"Experiment prefix filter (default: {DEFAULT_PREFIX})",
    )
    parser.add_argument(
        "--experiment",
        type=str,
        help="Download a single experiment only (e.g., E-GEOD-5305)",
    )
    parser.add_argument(
        "--max-size",
        type=float,
        default=DEFAULT_MAX_SIZE_MB,
        help=f"Max file size in MB (default: {DEFAULT_MAX_SIZE_MB})",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="List files without downloading",
    )
    parser.add_argument(
        "--rdf",
        action="store_true",
        help="Generate RDF after each successful download",
    )
    parser.add_argument(
        "--verbose",
        "-v",
        action="store_true",
        help="Enable verbose output",
    )

    args = parser.parse_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    logger.info(f"GEA Downloader")
    logger.info(f"  Data directory: {args.data_dir}")
    logger.info(f"  Prefix filter: {args.prefix}")
    logger.info(f"  Max file size: {args.max_size} MB")
    if args.rdf:
        logger.info(f"  RDF output: {args.output_dir}/")
    if args.dry_run:
        logger.info("  Mode: DRY RUN (no files will be downloaded)")

    downloader = GEADownloader(
        data_dir=args.data_dir,
        prefix=args.prefix,
        max_size_mb=args.max_size,
        dry_run=args.dry_run,
        generate_rdf=args.rdf,
        output_dir=args.output_dir,
    )

    try:
        downloader.run(single_experiment=args.experiment)
    except KeyboardInterrupt:
        logger.info("Download interrupted by user. Progress has been saved.")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Download failed: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
