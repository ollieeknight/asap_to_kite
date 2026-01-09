"""Configuration for asap_to_kite."""

from dataclasses import dataclass
from typing import Literal

ConjugationType = Literal["TotalSeq-A", "TotalSeq-B"]


@dataclass
class ProcessingConfig:
    """FASTQ processing configuration."""

    n_cores: int = 4
    batch_size: int = 10_000_000
    rc_r2: bool = True
    conjugation: ConjugationType = "TotalSeq-A"

    def __post_init__(self):
        if self.n_cores < 1:
            raise ValueError(f"Need at least 1 core, got {self.n_cores}")
        if self.batch_size < 1000:
            raise ValueError(f"Batch size too small: {self.batch_size} (minimum 1000)")
        if self.conjugation not in ["TotalSeq-A", "TotalSeq-B"]:
            raise ValueError(f"Unknown conjugation: {self.conjugation}")


@dataclass
class ConversionStats:
    """Conversion statistics."""

    total_reads: int = 0
    reads_processed: int = 0
    files_processed: int = 0
    batches_processed: int = 0

    def __str__(self):
        return (
            f"Files: {self.files_processed}, "
            f"Reads: {self.total_reads:,}, "
            f"Batches: {self.batches_processed}"
        )
