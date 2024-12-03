import sys
import re
import shutil
import os
import glob
import gzip
import numpy as np
import pandas as pd
import logging
import click
from multiprocessing import Pool
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqIO.QualityIO import FastqGeneralIterator

# Setup logger
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

logger.info("ASAP-to-kite, version 3")

@click.command()
@click.version_option()
@click.option('--fastq-folder', '-ff', required=True, help="Path(s) to folder(s) created by mkfastq or bcl2fastq. Multiple paths can be comma-separated.")
@click.option('--sample-prefix', '-sp', required=True, help="Prefix(es) of the filenames of fastq_folder to select. Multiple prefixes can be comma-separated.")
@click.option('--output-folder', '-of', default="output", help="Path to the output folder where results will be saved. Default is 'output'.")
@click.option('--output-name', '-on', default="asap2kite", help="Unique run ID used to name the output files. Default is 'asap2kite'.")
@click.option('--totalseq-conjugation', '-tc', default="TotalSeq-A", help="Antibody conjugation type. Options are 'TotalSeq-A' (default) or 'TotalSeq-B'.")
@click.option('--cores', '-c', default=4, help="Number of CPU cores for parallel processing. Default is 4.")
@click.option('--nreads', '-nr', default=10000000, help="Maximum number of reads to process in one iteration. Default is 10,000,000.")
@click.option('--no-rc-R2', '-nrr2', is_flag=True, default=False, help="Disable reverse complement of R2 (barcode). Default is to perform reverse complement.")
@click.option('--example', is_flag=True, help="Show example usage.")

def main(fastq_folder, sample_prefix, output_folder, output_name, totalseq_conjugation, cores, nreads, no_rc_r2, example):
    if example:
        print("Example usage:")
        print("python asap_to_kite.py -ff /path/to/fastq_folder -sp sample_prefix -of output_folder -on output_id -tc TotalSeq-A -c 4 -nr 10000000 -nrr2")
        return

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

    # Ensure the output folder exists
    if output_folder and not os.path.exists(output_folder):
        os.makedirs(output_folder)

    R1s_for_analysis = parse_directories(fastq_folder, sample_prefix)

    logger.info("Processing these fastq samples: ")
    for r in R1s_for_analysis:
        logger.info(r.replace("_R1_001.fastq.gz", ""))

    try:
        process_fastq_folder(R1s_for_analysis, output_name, output_folder, nreads, cores, not no_rc_r2, totalseq_conjugation)
        logger.info("Conversion completed successfully")
    except Exception as e:
        logger.error(f"An error occurred during processing: {e}")

def parse_directories(folder_of_fastq_folder, sample_name):
    list_folders = folder_of_fastq_folder.split(",")
    list_samples = sample_name.split(",")

    all_R1s = []
    for path_to_folder in list_folders:
        for sample_name_one in list_samples:
            matching_R1s = glob.glob(f"{path_to_folder}/*{sample_name_one}*_R1_001.fastq.gz")
            for R1file in matching_R1s:
                L00x = re.search(r'L00\d', R1file)
                if L00x:
                    R2file = R1file.replace("_R1_001.fastq.gz", f"_R2_001.fastq.gz")
                    R3file = R1file.replace("_R1_001.fastq.gz", f"_R3_001.fastq.gz")
                    if os.path.exists(R2file) and os.path.exists(R3file):
                        all_R1s.append(R1file)
    return all_R1s

def batch_iterator(iterator, batch_size):
    entry = True
    while entry:
        batch = []
        while len(batch) < batch_size:
            try:
                entry = next(iterator)
            except StopIteration:
                entry = None
            if entry is None:
                break
            batch.append(entry)
        if batch:
            yield batch

def process_fastq_folder(R1s_for_analysis, out, output_folder, n_reads, n_cpu, rc_R2, conjugation):
    outfq1file = os.path.join(output_folder, f"{out}_R1.fastq.gz")
    outfq2file = os.path.join(output_folder, f"{out}_R2.fastq.gz")
    with gzip.open(outfq1file, "wt") as out_f1, gzip.open(outfq2file, "wt") as out_f2:
        for R1file in R1s_for_analysis:
            R2file = R1file.replace("_R1_001.fastq.gz", "_R2_001.fastq.gz")
            R3file = R1file.replace("_R1_001.fastq.gz", "_R3_001.fastq.gz")

            it1 = batch_iterator(FastqGeneralIterator(gzip.open(R1file, "rt")), n_reads)
            it2 = batch_iterator(FastqGeneralIterator(gzip.open(R2file, "rt")), n_reads)
            it3 = batch_iterator(FastqGeneralIterator(gzip.open(R3file, "rt")), n_reads)

            for batch_R1 in it1:
                batch_R2 = next(it2)
                batch_R3 = next(it3)

                with Pool(processes=n_cpu) as pool:
                    pm = pool.map(asap_to_kite_conversion_wrapper, zip(batch_R1, batch_R2, batch_R3, [rc_R2]*len(batch_R1), [conjugation]*len(batch_R1)))

                fq_data = list(map(''.join, zip(*[item.pop(0) for item in pm])))
                out_f1.writelines(fq_data[0])
                out_f2.writelines(fq_data[1])

def asap_to_kite_conversion(listRead1, listRead2, listRead3, rc_R2, conjugation):
    title1, sequence1, quality1 = listRead1
    title2, sequence2, quality2 = listRead2
    title3, sequence3, quality3 = listRead3

    if rc_R2:
        sequence2 = str(Seq(sequence2).reverse_complement())
        quality2 = quality2[::-1]

    if conjugation == "TotalSeq-A":
        new_sequence1 = sequence2 + sequence1[:10]
        new_sequence2 = sequence3
        new_quality1 = quality2 + quality1[:10]
        new_quality2 = quality3
    elif conjugation == "TotalSeq-B":
        new_sequence1 = sequence2 + sequence3[:10] + sequence3[25:34]
        new_sequence2 = sequence3[10:25]
        new_quality1 = quality2 + quality3[:10] + quality3[25:34]
        new_quality2 = quality3[10:25]

    out_fq1 = formatRead(title1, new_sequence1, new_quality1)
    out_fq2 = formatRead(title2, new_sequence2, new_quality2)

    return [[out_fq1, out_fq2]]

def asap_to_kite_conversion_wrapper(args):
    return asap_to_kite_conversion(*args)

def formatRead(title, sequence, quality):
    return f"@{title}\n{sequence}\n+\n{quality}\n"

if __name__ == "__main__":
    main()