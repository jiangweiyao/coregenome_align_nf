#!/usr/bin/env nextflow


Channel.fromFilePairs( params.in ).into { fastq_files }

println """\
         This workflow finds all the paired files matching your glob pattern and prints them. 
         Use it to check which files match to troubleshoot or prevent unexpected surprises. 
         Usage: nextflow run <path>/filepair_finder.nf --in="<glob>"

         Common glob patterns:
         "<path>/*{1,2}.fastq.gz" for paired SRA files
         "<path>/*{R1_001,R2_001}.fastq.gz" for paired Illumina files

         Finding files in subfolders of a directory:
         "<path>/**/*{R1_001,R2_001}.fastq.gz"

         File Pair Finder for    NEXTFLOW      PIPELINE   
         =====================================================
         input reads (--in)                  : ${params.in}

         """
         .stripIndent()

fastq_files.subscribe { println it }
