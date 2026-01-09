"""Command-line interface for asap_to_kite."""

import logging
from importlib.metadata import version

import click

from asap_to_kite.core.config import ProcessingConfig
from asap_to_kite.core.exceptions import AsapToKiteError
from asap_to_kite.core.processor import process_fastq_folder

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)


@click.command()
@click.version_option(version=version("asap-to-kite"))
@click.option(
    "--fastq-folder",
    "-ff",
    required=True,
    help="Path(s) to folder(s) created by mkfastq or bcl2fastq. Multiple paths can be comma-separated.",
)
@click.option(
    "--sample-prefix",
    "-sp",
    required=True,
    help="Prefix(es) of the filenames of fastq_folder to select. Multiple prefixes can be comma-separated.",
)
@click.option(
    "--output-folder",
    "-of",
    default="output",
    help="Path to the output folder where results will be saved. Default is 'output'.",
)
@click.option(
    "--output-name",
    "-on",
    default="asap2kite",
    help="Unique run ID used to name the output files. Default is 'asap2kite'.",
)
@click.option(
    "--totalseq-conjugation",
    "-tc",
    default="TotalSeq-A",
    type=click.Choice(["TotalSeq-A", "TotalSeq-B"], case_sensitive=False),
    help="Antibody conjugation type. Options are 'TotalSeq-A' (default) or 'TotalSeq-B'.",
)
@click.option(
    "--cores",
    "-c",
    default=4,
    type=int,
    help="Number of CPU cores for parallel processing. Default is 4.",
)
@click.option(
    "--nreads",
    "-nr",
    default=10000000,
    type=int,
    help="Maximum number of reads to process in one iteration. Default is 10,000,000.",
)
@click.option(
    "--no-rc-R2",
    "-nrr2",
    is_flag=True,
    default=False,
    help="Disable reverse complement of R2 (barcode). Default is to perform reverse complement.",
)
@click.option("--example", is_flag=True, help="Show example usage.")
def cli(
    fastq_folder: str,
    sample_prefix: str,
    output_folder: str,
    output_name: str,
    totalseq_conjugation: str,
    cores: int,
    nreads: int,
    no_rc_r2: bool,
    example: bool,
) -> None:
    """Convert ASAP-seq FASTQ files for processing with kite (kallisto | bustools)."""
    if example:
        print("\nExample usage:")
        print(
            "asap-to-kite -ff /path/to/fastq_folder -sp sample_prefix -of output_folder -on output_id\n"
        )
        return

    logger.info(f"asap-to-kite version {version('asap-to-kite')}")

    # Log the configuration details
    config_details = {
        "FASTQ folder": fastq_folder,
        "Sample name": sample_prefix,
        "Output folder": output_folder,
        "Output name": output_name,
        "Number of cores": cores,
        "Number of reads to process at a time": nreads,
        "Reverse complement": not no_rc_r2,
        "TotalSeq format": totalseq_conjugation,
    }

    logger.info("Configuration for analysis:")
    for key, value in config_details.items():
        logger.info(f"{key:40}: {value}")

    try:
        from asap_to_kite.utils.file_utils import parse_directories

        # Parse input directories and find FASTQ files
        r1_files = parse_directories(fastq_folder, sample_prefix)
        logger.info(f"Processing {len(r1_files)} FASTQ samples:")
        for r in r1_files:
            logger.info(f"  {r.replace('_R1_001.fastq.gz', '')}")

        # Create processing configuration
        config = ProcessingConfig(
            n_cores=cores,
            batch_size=nreads,
            rc_r2=not no_rc_r2,
            conjugation=totalseq_conjugation,
        )

        # Process FASTQ files
        stats = process_fastq_folder(r1_files, output_name, output_folder, config)
        logger.info(f"Complete: {stats.reads_processed:,} reads, {stats.files_processed} files")

    except AsapToKiteError as e:
        logger.error(f"Error: {e}")
        raise SystemExit(1) from None
    except Exception as e:
        logger.error(f"Unexpected: {e}")
        raise SystemExit(1) from e


def main() -> None:
    """Entry point."""
    cli()


if __name__ == "__main__":
    main()
