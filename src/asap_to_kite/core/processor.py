"""ASAP-seq to kite conversion."""

import gc
import gzip
import logging
import os
import sys
from concurrent.futures import ProcessPoolExecutor

from Bio.Seq import Seq
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from tqdm import tqdm

from asap_to_kite.core.config import ConversionStats, ProcessingConfig
from asap_to_kite.core.exceptions import ProcessingError

logger = logging.getLogger(__name__)

# Use fork on Linux for better performance, spawn on macOS/Windows
MP_CONTEXT = "fork" if sys.platform.startswith("linux") else "spawn"


def batch_iterator(iterator, batch_size: int):
    """Yield batches from FASTQ iterator."""
    batch = []
    for entry in iterator:
        batch.append(entry)
        if len(batch) >= batch_size:
            yield batch
            batch = []
    if batch:
        yield batch


def asap_to_kite_conversion_batch(
    batch_reads: list[tuple[tuple[str, str, str], tuple[str, str, str], tuple[str, str, str]]],
    config: ProcessingConfig,
) -> list[tuple[str, str]]:
    """Convert batch of ASAP-seq reads to kite format."""
    results = []

    for list_read1, list_read2, list_read3 in batch_reads:
        title1, sequence1, quality1 = list_read1
        title2, sequence2, quality2 = list_read2
        title3, sequence3, quality3 = list_read3

        # Reverse complement R2 if requested
        if config.rc_r2:
            sequence2 = str(Seq(sequence2).reverse_complement())
            quality2 = quality2[::-1]

        # Convert based on conjugation type
        if config.conjugation == "TotalSeq-A":
            # TotalSeq-A: cell barcode (R2) + UMI (first 10bp of R1)
            new_sequence1 = sequence2 + sequence1[:10]
            new_sequence2 = sequence3
            new_quality1 = quality2 + quality1[:10]
            new_quality2 = quality3
        elif config.conjugation == "TotalSeq-B":
            # TotalSeq-B: cell barcode (R2) + feature barcode parts
            new_sequence1 = sequence2 + sequence3[:10] + sequence3[25:34]
            new_sequence2 = sequence3[10:25]
            new_quality1 = quality2 + quality3[:10] + quality3[25:34]
            new_quality2 = quality3[10:25]
        else:
            raise ValueError(f"Unknown conjugation type: {config.conjugation}")

        out_fq1 = format_read(title1, new_sequence1, new_quality1)
        out_fq2 = format_read(title2, new_sequence2, new_quality2)

        results.append((out_fq1, out_fq2))

    return results


def format_read(title: str, sequence: str, quality: str) -> str:
    """Format read in FASTQ format."""
    return f"@{title}\n{sequence}\n+\n{quality}\n"


def process_fastq_folder(
    r1s_for_analysis: list[str],
    output_name: str,
    output_folder: str,
    config: ProcessingConfig,
) -> ConversionStats:
    """Process FASTQ files and convert to kite format."""
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    outfq1_file = os.path.join(output_folder, f"{output_name}_R1.fastq.gz")
    outfq2_file = os.path.join(output_folder, f"{output_name}_R2.fastq.gz")

    logger.info(f"Output: {outfq1_file} and {outfq2_file}")
    logger.info(f"Using {config.n_cores} cores, batch size {config.batch_size:,}")

    stats = ConversionStats()

    try:
        with gzip.open(outfq1_file, "wt") as out_f1, gzip.open(outfq2_file, "wt") as out_f2:
            for r1_file in r1s_for_analysis:
                r2_file = r1_file.replace("_R1_001.fastq.gz", "_R2_001.fastq.gz")
                r3_file = r1_file.replace("_R1_001.fastq.gz", "_R3_001.fastq.gz")

                logger.info(f"Processing: {os.path.basename(r1_file)}")

                file_stats = _process_single_fastq(
                    r1_file, r2_file, r3_file, out_f1, out_f2, config
                )

                stats.total_reads += file_stats.total_reads
                stats.reads_processed += file_stats.reads_processed
                stats.batches_processed += file_stats.batches_processed
                stats.files_processed += 1
                gc.collect()

        logger.info(f"Done: {stats}")
        return stats

    except Exception as e:
        raise ProcessingError(f"Failed to process {os.path.basename(r1_file)}: {e}") from e


def _process_single_fastq(
    r1_file: str,
    r2_file: str,
    r3_file: str,
    out_f1,
    out_f2,
    config: ProcessingConfig,
) -> ConversionStats:
    """Process a single FASTQ file triplet."""
    stats = ConversionStats()

    # Open all three FASTQ files
    with (
        gzip.open(r1_file, "rt") as f1,
        gzip.open(r2_file, "rt") as f2,
        gzip.open(r3_file, "rt") as f3,
    ):
        it1 = FastqGeneralIterator(f1)
        it2 = FastqGeneralIterator(f2)
        it3 = FastqGeneralIterator(f3)

        # Batch the iterators
        batches_r1 = batch_iterator(it1, config.batch_size)
        batches_r2 = batch_iterator(it2, config.batch_size)
        batches_r3 = batch_iterator(it3, config.batch_size)

        # Process batches with progress bar
        with tqdm(desc="Processing reads", unit="reads", unit_scale=True) as pbar:
            if config.n_cores == 1:
                # Sequential processing
                for batch_r1, batch_r2, batch_r3 in zip(
                    batches_r1, batches_r2, batches_r3, strict=True
                ):
                    batch_data = list(zip(batch_r1, batch_r2, batch_r3, strict=True))
                    results = asap_to_kite_conversion_batch(batch_data, config)

                    # Write results
                    for r1_out, r2_out in results:
                        out_f1.write(r1_out)
                        out_f2.write(r2_out)

                    stats.reads_processed += len(results)
                    stats.total_reads += len(batch_data)
                    stats.batches_processed += 1
                    pbar.update(len(batch_data))
            else:
                # Parallel processing with ProcessPoolExecutor
                _process_parallel(
                    batches_r1, batches_r2, batches_r3, out_f1, out_f2, config, stats, pbar
                )

    return stats


def _process_parallel(
    batches_r1, batches_r2, batches_r3, out_f1, out_f2, config, stats, pbar
) -> None:
    """Process batches in parallel."""
    import multiprocessing as mp

    with ProcessPoolExecutor(
        max_workers=config.n_cores, mp_context=mp.get_context(MP_CONTEXT)
    ) as executor:
        for batch_r1, batch_r2, batch_r3 in zip(batches_r1, batches_r2, batches_r3, strict=True):
            batch_data = list(zip(batch_r1, batch_r2, batch_r3, strict=True))

            # Submit batch for processing
            future = executor.submit(asap_to_kite_conversion_batch, batch_data, config)

            try:
                results = future.result()

                # Write results
                for r1_out, r2_out in results:
                    out_f1.write(r1_out)
                    out_f2.write(r2_out)

                stats.reads_processed += len(results)
                stats.total_reads += len(batch_data)
                stats.batches_processed += 1
                pbar.update(len(batch_data))

            except Exception as e:
                raise ProcessingError(f"Batch failed: {e}") from e
