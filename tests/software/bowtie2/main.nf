#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BOWTIE2_BUILD } from '../../../software/bowtie2/build/main.nf' addParams( options: ['args': '--seed 1'] )
include { BOWTIE2_ALIGN } from '../../../software/bowtie2/align/main.nf'  addParams( options: ['args': ''] )

process grep_bam_header {

    conda (params.enable_conda ? 'bioconda::samtools=1.10' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container 'https://depot.galaxyproject.org/singularity/mulled-v2-ffbf83a6b0ab6ec567a336cf349b80637135bca3:9e14e16c284d6860574cf5b624bbc44c793cb024-0'
    } else {
        container 'quay.io/biocontainers/mulled-v2-ffbf83a6b0ab6ec567a336cf349b80637135bca3:9e14e16c284d6860574cf5b624bbc44c793cb024-0'
    }

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path('*.sam'), emit: bam

    script:
    //def prefix    = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    samtools view $bam > ${meta.id}.sam
    """
}

workflow test_bowtie2_build {
    fasta = file("${launchDir}/tests/data/fasta/E_coli/NC_010473.fa", checkIfExists: true)
    BOWTIE2_BUILD ( fasta )
}

workflow test_bowtie2_alignment_single_end {

    fasta = file("${launchDir}/tests/data/fasta/E_coli/NC_010473.fa", checkIfExists: true)
    BOWTIE2_BUILD ( fasta )

    def input = []
    input = [ [ id:'test', single_end:true ], // meta map
              [ file("${launchDir}/tests/data/fastq/rna/test_R1.fastq.gz", checkIfExists: true) ] ]
    BOWTIE2_ALIGN ( input, BOWTIE2_BUILD.out.index )
    grep_bam_header (BOWTIE2_ALIGN.out.bam)
}

workflow test_bowtie2_alignment_paired_end {
    fasta = file("${launchDir}/tests/data/fasta/E_coli/NC_010473.fa", checkIfExists: true)
    BOWTIE2_BUILD ( fasta )

    def input = []
    input = [ [ id:'test', single_end:false ], // meta map
              [ file("${launchDir}/tests/data/fastq/rna/test_R1.fastq.gz", checkIfExists: true),
                file("${launchDir}/tests/data/fastq/rna/test_R2.fastq.gz", checkIfExists: true) ] ]
    BOWTIE2_ALIGN ( input, BOWTIE2_BUILD.out.index )
    grep_bam_header (BOWTIE2_ALIGN.out.bam)
}
