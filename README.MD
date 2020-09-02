# Core Genome Align Pipeline in Nextflow

### This project is still undercontruction.
This pipeline takes a set of fastq files, assembles them (with Spades), and align the core genome against again a set of reference sequences using Prokka/Roary/Fasttree

### Dependencies:
1. [Nextflow](https://www.nextflow.io/)
2. [Singularity](https://sylabs.io/docs/) - if you are running this on a cluster, please install version 3.6 or later. 

### Usage: 
Test the full pipeline with prepackaged data:
```
nextflow run coregenome_align_nf/main.nf 
```

Use the configuration file to run on the SGE system in SCBS
```
nextflow run coregenome_align_nf/main.nf -c coregenome_align_nf/cdc_sge.config
```

Put in your own parameters and get log tracing example:
```
nextflow run coregenome_align_nf/main.nf --in="coregenome_align_nf/fastq/*{1,2}.fq.gz" --refseq="/coregenome_align_nf/refseq/*.fna" --out="test" --sizefilter=3000 -with-report report.html
```


