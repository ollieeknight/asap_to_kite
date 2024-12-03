# ASAP to kite
This script was originally developed by [Caleb Lareau](https://www.mskcc.org/research/ski/labs/caleb-lareau), but this repo is maintained by [Oliver Knight](https://immunologie.charite.de/metas/person/person/address_detail/oliver_knight/).  

It describes a script (`asap_to_kite.py`) which converts FASTQ files generated from ASAP-seq for downstream processing with kite (kallisto | bustools). 

## About

## Sample use cases

### One sample, one directory

The most basic use case is when we have one library sequenced once. From the demultiplexing,
we should see files that look like this:

```sh
git clone https://github.com/ollieeknight/asap_to_kite

cd asap_to_kite

tree tests/data1/

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

```
python asap_to_kite.py --fastq-folder tests/data1 --sample-prefix test1 --output-folder single_sample --output-name test1 --cores 2

2024-12-03 10:34:06,537 - INFO - ASAP-to-kite, version 3
2024-12-03 10:34:06,537 - INFO - Configuration for analysis:
2024-12-03 10:34:06,538 - INFO - FASTQ folder                            : tests/data1
2024-12-03 10:34:06,538 - INFO - Sample name                             : test1
2024-12-03 10:34:06,538 - INFO - Output folder                           : single_sample
2024-12-03 10:34:06,538 - INFO - Output name                             : test1
2024-12-03 10:34:06,538 - INFO - Number of cores                         : 2
2024-12-03 10:34:06,538 - INFO - Number of reads to process at a time    : 10000000
2024-12-03 10:34:06,538 - INFO - Reverse complement                      : True
2024-12-03 10:34:06,538 - INFO - TotalSeq format                         : TotalSeq-A
2024-12-03 10:34:06,551 - INFO - Processing these fastq samples:
2024-12-03 10:34:06,551 - INFO - tests/data1/test1_S1_L003
2024-12-03 10:34:06,551 - INFO - tests/data1/test1_S1_L004
2024-12-03 10:34:11,940 - INFO - Conversion completed successfully
```

### One sample, multiple directories 

If multiple sequencing rounds are performed, we can supply all sequencing libraries as a comma-separated list:

```
python asap_to_kite.py --fastq-folder tests/data1,tests/data2 --sample-prefix test2,test2 --output-folder single_sample --output-name test1 --cores 2
2024-12-03 10:38:29,545 - INFO - ASAP-to-kite, version 3
2024-12-03 10:38:29,546 - INFO - Configuration for analysis:
2024-12-03 10:38:29,546 - INFO - FASTQ folder                            : tests/data1,tests/data2
2024-12-03 10:38:29,546 - INFO - Sample name                             : test2,test2
2024-12-03 10:38:29,546 - INFO - Output folder                           : single_sample
2024-12-03 10:38:29,546 - INFO - Output name                             : test1
2024-12-03 10:38:29,546 - INFO - Number of cores                         : 2
2024-12-03 10:38:29,546 - INFO - Number of reads to process at a time    : 10000000
2024-12-03 10:38:29,546 - INFO - Reverse complement                      : True
2024-12-03 10:38:29,546 - INFO - TotalSeq format                         : TotalSeq-A
2024-12-03 10:38:29,547 - INFO - Processing these fastq samples:
2024-12-03 10:38:29,547 - INFO - tests/data1/test2_S1_L002
2024-12-03 10:38:29,547 - INFO - tests/data1/test2_S1_L001
2024-12-03 10:38:29,547 - INFO - tests/data1/test2_S1_L002
2024-12-03 10:38:29,547 - INFO - tests/data1/test2_S1_L001
2024-12-03 10:38:29,547 - INFO - tests/data2/test2_S1_L003
2024-12-03 10:38:29,547 - INFO - tests/data2/test2_S1_L003
2024-12-03 10:38:45,450 - INFO - Conversion completed successfully
```

## Important

This code works for one biological sample at a time. If multiple samples are supplied in the command line execution, then they will be merged (under the assumption that they were called different things). Execute the code sequentially for each sample in the event of multiple biological samples. 

## Options

```
python asap_to_kite.py --help
```

yields

```
Usage: asap_to_kite.py [OPTIONS]

Options:
  --fastq-folder, -ff TEXT        Path(s) to folder(s) created by mkfastq or
                                  bcl2fastq. Multiple paths can be comma-
                                  separated.  [required]
  --sample-prefix, -sp TEXT       Prefix(es) of the filenames of fastq_folder
                                  to select. Multiple prefixes can be comma-
                                  separated.  [required]
  --output-folder, -of TEXT       Path to the output folder where results will
                                  be saved. Default is 'output'.
  --output-name, -on TEXT         Unique run ID used to name the output files.
                                  Default is 'asap2kite'.
  --totalseq-conjugation, -tc TEXT
                                  Antibody conjugation type. Options are
                                  'TotalSeq-A' (default) or 'TotalSeq-B'.
  --cores, -c INTEGER             Number of CPU cores for parallel processing.
                                  Default is 4.
  --nreads, -nr INTEGER           Maximum number of reads to process in one
                                  iteration. Default is 10,000,000.
  --no-rc-R2, -nrr2               Disable reverse complement of R2 (barcode).
                                  Default is to perform reverse complement.
  --example                       Show example usage.
  --help                          Show this message and exit.
```
