#!/usr/bin/env nextflow

// Step 1: FastQC
process FastQC {
    input:
    set set sample_id, file(reads) from read_pairs // read_pairs is channel to pass paired end reads

    output:
    path "fastqc/*" into fastqc_results

    script:
    """
    mkdir -p fastqc
    fastqc -o fastqc ${reads.join(' ')}

    """
}

// Step 2: Trim Galore
process TrimGalore {
    input:
    set sample_id, file(reads) from read_pairs

    output:
    set sample_id, file("*.fq.gz") into trimmed_reads // pushing the trimmed reads into a trimmed reads channels

    script:
    """
    trim_galore --paired --gzip --fastqc ${reads[0]} ${reads[1]}
    """
}

// Step 3: Building taxus genome index

process BuildHisat2Index {
    input:
    file taxus_genome from taxus_genome_ch

    output:
    path "taxus_genome.*.ht2" into hisat2_index

    script:
    """
    hisat2-build ${taxus_genome} taxus_genome_index_prefix
    """
}

// Step 4: Read alignment against Taxus genome

process HISAT2Align {

    input:
    set sample_id, file(reads) from trimmed_reads
    path hisat2_index

    output:
    set sample_id, file("${sample_id}.bam") into bam_files

    script:
    """
    hisat2 -x taxus_genome_index_prefix -1 ${reads[0]} -2 ${reads[1]} | \
    samtools sort -o ${sample_id}.bam
    """
}

// Step 5: featureCounts

process FeatureCounts {
    input:
    set sample_id, file(bam) from bam_files

    output:
    set sample_id, file("${sample_id}_counts.txt") into count_files

    script:
    """
    featureCounts -a ${params.annotation} -o ${sample_id}_counts.txt ${bam}
    """
}

workflow {

    // Create initial channels from raw data
    read_pairs = Channel.fromFilePairs(params.fastq, flat: true)         // Paired-end reads
    taxus_genome_ch = Channel.fromPath("ref/taxus_genome.fa")            // Reference genome fasta

    // Run processes in order

    // Step 1: FastQC on raw reads
    FastQC(read_pairs.map { sample_id, files -> tuple(sample_id, files.flatten()) })

    // Step 2: TrimGalore
    TrimGalore(read_pairs)

    // Step 3: Build index (runs once)
    BuildHisat2Index(taxus_genome_ch)

    // Step 4: Align reads
    HISAT2Align(trimmed_reads, hisat2_index)

    // Step 5: Count features
    FeatureCounts(bam_files)
}
