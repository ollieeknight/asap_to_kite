# ASAP to kite

This script was originally developed by [Caleb Lareau](https://www.mskcc.org/research/ski/labs/caleb-lareau), but this repo is maintained by [Oliver Knight](https://immunologie.charite.de/metas/person/person/address_detail/oliver_knight/).  


## Installation

### From source
```bash
git clone https://github.com/ollieeknight/asap_to_kite
cd asap_to_kite
pip install -e .
```

## About ASAP-seq

This tool converts ASAP-seq FASTQ files (R1, R2, R3) into the format required by kite/kallisto for feature counting:

- **TotalSeq-A**: Combines cell barcode (R2) with UMI (first 10bp of R1)
- **TotalSeq-B**: Rearranges feature barcode components from R3

## Sample use cases

### One sample, one directory

The most basic use case is when we have one library sequenced once. From the demultiplexing,
we should see files that look like this:

```sh
git clone https://github.com/ollieeknight/asap_to_kite

cd asap_to_kite
bash
asap-to-kite

tests/data1/
├── test1_S1_L003_R1_001.fastq.gz
├── test1_S1_L003_R2_001.fastq.gz
├── test1_S1_L003_R3_001.fastq.gz

├── test1_S1_L004_R1_001.fastq.gz
├── test1_S1_L004_R2_001.fastq.gz
├── test1_S1_L004_R3_001.fastq.gz

├── test2_S1_L001_R1_001.fastq.gz
├── test2_S1_L001_R2_001.fastq.gz
├── test2_S1_L001_R3_001.fastq.gz

├── test2_S1_L002_R1_001.fastq.gz
├── test2_S1_L002_R2_001.fastq.gz
└── test2_S1_L002_R3_001.fastq.gz
```

We can process these fastqs in one line:
```bash
python asap_to_kite.py --fastq-folder tests/data1 --sample-prefix test1 --output-folder single_sample --output-name test1 --cores 2
```

### One sample, multiple directories 

If multiple sequencing rounds are performed, we can supply all sequencing libraries as a comma-separated list:

```bash
asap-to-kite --fastq-folder tests/data1,tests/data2 --sample-prefix test2,test2 --output-folder single_sample --output-name test1 --cores 2
```

## Important

This code works for one biological sample at a time. If multiple samples are supplied in the command line execution, then they will be merged (under the assumption that they were called different things). Execute the code sequentially for each sample in the event of multiple biological samples. 

## Citation

If you use this tool in your research, please cite:

Lareau, C.A., Ludwig, L.S., Muus, C. et al. Massively parallel single-cell mitochondrial DNA genotyping and chromatin profiling. Nat Biotechnol 39, 451–461 (2021). https://doi.org/10.1038/s41587-020-0645-6