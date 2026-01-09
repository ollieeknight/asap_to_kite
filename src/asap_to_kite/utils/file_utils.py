"""Utility functions for file parsing and validation."""

import glob
import logging
import os
import re

from asap_to_kite.core.exceptions import InvalidInputError, ValidationError

logger = logging.getLogger(__name__)


def validate_fastq_files(r1_file: str) -> tuple[str, str]:
    """Validate R2 and R3 files exist for R1."""
    if not os.path.exists(r1_file):
        raise InvalidInputError(f"File not found: {r1_file}")

    if not r1_file.endswith(".fastq.gz"):
        raise InvalidInputError(f"Not a gzipped FASTQ: {r1_file}")

    r2_file = r1_file.replace("_R1_001.fastq.gz", "_R2_001.fastq.gz")
    r3_file = r1_file.replace("_R1_001.fastq.gz", "_R3_001.fastq.gz")

    if not os.path.exists(r2_file):
        raise ValidationError(f"Missing R2: {r2_file}")

    if not os.path.exists(r3_file):
        raise ValidationError(f"Missing R3: {r3_file}")

    if os.path.getsize(r1_file) == 0:
        raise ValidationError(f"Empty file: {r1_file}")

    return r2_file, r3_file


def parse_directories(folder_of_fastq_folder: str, sample_name: str) -> list[str]:
    """Find matching R1 FASTQ files in directories."""
    list_folders = folder_of_fastq_folder.split(",")
    list_samples = sample_name.split(",")

    all_r1s = []

    for path_to_folder in list_folders:
        path_to_folder = path_to_folder.strip()

        if not os.path.exists(path_to_folder):
            raise InvalidInputError(f"Folder not found: {path_to_folder}")

        if not os.path.isdir(path_to_folder):
            raise InvalidInputError(f"Not a directory: {path_to_folder}")

        for sample_name_one in list_samples:
            sample_name_one = sample_name_one.strip()
            matching_r1s = glob.glob(f"{path_to_folder}/*{sample_name_one}*_R1_001.fastq.gz")

            for r1_file in matching_r1s:
                if not re.search(r"L00\d", r1_file):
                    logger.warning(f"Skipping (no lane number): {r1_file}")
                    continue

                try:
                    validate_fastq_files(r1_file)
                    all_r1s.append(r1_file)
                except ValidationError as e:
                    logger.error(str(e))
                    raise

    if not all_r1s:
        raise InvalidInputError(
            f"No FASTQ files found in {folder_of_fastq_folder} " f"matching {sample_name}"
        )

    return sorted(all_r1s)
