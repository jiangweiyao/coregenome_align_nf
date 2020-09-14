#!/usr/bin/env nextflow


Channel.fromFilePairs( params.in ).into { fastq_files; fastq_files2; fastq_files3 }
Channel.fromPath( params.refseq ).set {refseq_files}

adapters = file(params.adapters)
phix = file(params.phix)
mash_genome_db = file(params.genome_db)

mash_parser = file("$baseDir/bin/mash_screen_parser.py")

println """\
         CORE GENOME ALIGN     NEXTFLOW      PIPELINE
         =====================================================
         input reads (--in)                  : ${params.in}
         reference sequences (--refseq)      : ${params.refseq}
         outdir (--out)                      : ${params.out}
         Mash Genome Reference (--genome_db) : ${params.genome_db}
         Max threads per process (--thread)  : ${params.thread}
         size cutoff (--sizefilter)          : ${params.sizefilter}
         """
         .stripIndent()


mash_genome_file = file("$baseDir/db/refseq.genomes.k21s1000.msh")
if(!mash_genome_file.exists()){
    println("Mash genome reference missing. Downloading...")
    mash_genome_file = file('https://gembox.cbcb.umd.edu/mash/refseq.genomes.k21s1000.msh')
    mash_genome_file.copyTo("$baseDir/db/refseq.genomes.k21s1000.msh")
}

//fastq_files.mix(fasta_files).subscribe onNext: { println it }, onComplete: { println 'Done' }

//fasta_files.subscribe {println it.simpleName}

process fastqc {
    
    //errorStrategy 'ignore'
    publishDir params.out, pattern: "*.html", mode: 'copy', overwrite: true

    input:
    set val(name), file(fastq) from fastq_files2
 
    output:
    file "*_fastqc.{zip,html}" into qc_files

    """
    fastqc -q ${fastq}
    """
}

process multiqc {

    errorStrategy 'ignore'
    publishDir params.out, mode: 'copy', overwrite: true

    input:
    file reports from qc_files.collect().ifEmpty([])

    output:
    path "multiqc_report.html" into multiqc_output

    """
    multiqc $reports
    """
}


process mash_screen_genome {

    //errorStrategy 'ignore'
    //publishDir params.out, mode: 'copy', overwrite: true
    memory '8 GB'

    input:
    tuple val(name), file(fastq) from fastq_files3

    output:
    tuple val(name), path("*_pathogen_id_raw.out") into mash_screen_genome_out

    """
    cat ${fastq} | mash screen ${mash_genome_db} - | sort -gr > ${name}_pathogen_id_raw.out
    """
}

process tabulate_mash_genome {

    //errorStrategy 'ignore'
    publishDir params.out, mode: 'copy', overwrite: true

    input:
    tuple val(name), file(table) from mash_screen_genome_out

    output:
    path "*_pathogen_id.out" into tabulate_mash_genome_out

    """
    python3 ${mash_parser} -i ${table} -o ${name}_pathogen_id.out -n ${params.mashtopn}
    """
}

process prokka_fasta {

    //errorStrategy 'ignore'
    //publishDir params.out, mode: 'copy', overwrite: true

    cpus params.thread

    input:
    file(fasta) from refseq_files

    output:
    path("*.gff") into prokka_fasta_output

    """
    prokka --cpus ${params.thread} --outdir ${fasta.simpleName}_prokka --prefix ${fasta.simpleName} ${fasta}
    cp ${fasta.simpleName}_prokka/${fasta.simpleName}.gff ./
    """
}

process bbmap_adapter_trimming {

    //errorStrategy 'ignore'
    //publishDir params.out, mode: 'copy', overwrite: true

    input:
    tuple val(name), file(fastq) from fastq_files

    output:
    tuple val(name), file("*_clean{1,2}.fq.gz") into trimmed_fastq

    """
    bbduk.sh -Xmx1g in1=${fastq[0]} in2=${fastq[1]} out1=int1.fq.gz out2=int2.fq.gz ref=${phix} k=31 hdist=1 t=1
    bbduk.sh -Xmx1g in1=int1.fq.gz in2=int2.fq.gz out1=${name}_clean1.fq.gz out2=${name}_clean2.fq.gz ref=${adapters} ktrim=r k=23 mink=11 hdist=1 t=1 tpe tbo
    """
}

process assembly {

    //errorStrategy 'ignore'
    //publishDir params.out, mode: 'copy', overwrite: true
    errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }
    maxRetries 3

    memory { 4.GB * params.thread * task.attempt * task.attempt }
    cpus params.thread

    input:
    tuple val(name), file(fastq) from trimmed_fastq

    output:
    tuple val(name), path("*_scaffolds.fasta") into assembly_output

    """
    spades.py -1 ${fastq[0]} -2 ${fastq[1]} -o ${name} -t ${params.thread} -m \$((4 * $params.thread * $task.attempt * $task.attempt))
    cp ${name}/scaffolds.fasta ${name}_scaffolds.fasta
    """

}

process bbmap_size_filter {

    //errorStrategy 'ignore'
    publishDir params.out, mode: 'copy', overwrite: true

    input:
    tuple val(name), file(assembly) from assembly_output

    output:
    tuple val(name), path("*_scaffolds_filtered.fasta") into assembly_filter_output, assembly_filter_output2

    """
    reformat.sh in=${assembly} out=${assembly.simpleName}_filtered.fasta minlength=${params.sizefilter}
    """
}

process busco {

    errorStrategy 'ignore'
    publishDir params.out, mode: 'copy', overwrite: true
    memory '8 GB'

    input:
    tuple val(name), file(assembly) from assembly_filter_output2

    output:
    path("*/short_summary*.txt") into busco_output

    """
    busco --auto-lineage-prok -f -m geno -o ${name}_busco -i ${assembly} -c 1
    """
}


process prokka_assembly {

    //errorStrategy 'ignore'
    //publishDir params.out, mode: 'copy', overwrite: true

    cpus params.thread

    input:
    tuple val(name), file(assembly) from assembly_filter_output

    output:
    tuple val(name), path("*.gff") into prokka_assembly_output

    """
    prokka --cpus ${params.thread} --outdir ${name}_prokka --prefix ${name} ${assembly}
    cp ${name}_prokka/${name}.gff ./
    """
}

process roary {

    //errorStrategy 'ignore'
    publishDir params.out, mode: 'copy', overwrite: true

    memory { 2.GB * params.roarythread }
    cpus params.roarythread

    input:
    file gff from prokka_fasta_output.mix(prokka_assembly_output).collect()

    output:
    path("*") into roary_output

    """
    roary -f . -e -n -v -r *.gff -p ${params.roarythread}
    """
}

process raxmlng {

    //errorStrategy 'ignore'
    publishDir params.out, mode: 'copy', overwrite: true

    memory '16 GB'
    cpus params.thread

    input:
    file roary_output from roary_output

    output:
    path("core_gene.raxml*") into raxml_output

    """
    raxml-ng --msa core_gene_alignment.aln --model GTR+G --prefix core_gene --threads ${params.thread} --seed 1234 
    """
}
