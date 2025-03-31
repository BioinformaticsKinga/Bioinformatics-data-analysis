// Define input parameters and set the appropriate paths for your files
params.reads = 'zygmukin_reads/*.fastq.gz'
params.outdir = 'zygmukin_results'
params.adapters = 'TruSeq3-PE.fa'  // Adapter sequences for trimming

// Process 1: FastQC (Quality Control)
process fastqc {
    container 'biocontainers/fastqc:v0.11.9'
    input:
    path reads

    output:
    path 'fastqc_reports'

    script:
    """
    fastqc -o fastqc_reports $reads
    """
}

// Process 2: Trimmomatic (Read Trimming)
process trimmomatic {
    container 'biocontainers/trimmomatic:v0.39'
    input:
    path reads
    path params.adapters

    output:
    path 'trimmed_reads'

    script:
    """
    trimmomatic PE -threads 4 ${reads[0]} ${reads[1]} \
        trimmed_reads/1P.fastq.gz trimmed_reads/1U.fastq.gz \
        trimmed_reads/2P.fastq.gz trimmed_reads/2U.fastq.gz \
        ILLUMINACLIP:${params.adapters}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    """
}

// Process 3: Trinity Assembly (De novo assembly)
process trinity_assembly {
    container 'trinityrnaseq/trinityrnaseq:2.13.2'
    input:
    path reads

    output:
    path 'trinity_output'

    script:
    """
    Trinity --seqType fq --max_memory 50G --left ${reads[0]} --right ${reads[1]} --CPU 4 --output trinity_output
    """
}

// Process 4: RSEM (RNA-Seq Quantification)
process rsem_analysis {
    container 'biocontainers/rsem:v1.3.1'
    input:
    path 'trinity_output/Trinity.fasta'
    path 'trimmed_reads/1P.fastq.gz'
    path 'trimmed_reads/2P.fastq.gz'

    output:
    path 'rsem_results'

    script:
    """
    rsem-calculate-expression --paired-end \
        --bowtie2 --alignments \
        --output-genome-bam --estimate-rspd \
        --strand-specific rsem_results \
        trimmed_reads/1P.fastq.gz trimmed_reads/2P.fastq.gz trinity_output/Trinity.fasta
    """
}

// Process 5: TransDecoder (ORF prediction)
process transdecoder_prediction {
    container 'biocontainers/transdecoder:v5.5.0'
    input:
    path 'trinity_output/Trinity.fasta'

    output:
    path 'transdecoder_output'

    script:
    """
    TransDecoder.LongOrfs -t trinity_output/Trinity.fasta
    TransDecoder.Predict -t trinity_output/Trinity.fasta
    mv trinity_output/Trinity.fasta.transdecoder.* transdecoder_output
    """
}

// Process 6: Trinotate Annotation
process trinnotate_annotation {
    container 'trinotate/trinotate:latest'
    input:
    path 'trinity_output/Trinity.fasta'
    path 'rsem_results/counts.txt'

    output:
    path 'trinotate_results'

    script:
    """
    ./trinotate_annotation.sh trinity_output/Trinity.fasta rsem_results/counts.txt trinotate_results
    """
}

// Process 7: KOBAS Annotation (Functional Annotation)
process kobas_annotation {
    container 'agbase/kobas:3.0.3_3'
    input:
    path 'transdecoder_output/Trinity.fasta.transdecoder.pep'

    output:
    path 'kobas_results'

    script:
    """
    kobas-annotate -i transdecoder_output/Trinity.fasta.transdecoder.pep -s ath -t fasta:nuc -o kobas_results
    """
}

// Process 8: Differential Expression Analysis (DESeq2, Limma-voom, edgeR)
process deseq2_analysis {
    container 'biocontainers/r-ver:v3.6.1'
    input:
    path 'rsem_results/counts.txt'

    output:
    path 'deseq2_results.csv'

    script:
    """
    Rscript deseq2_analysis.R rsem_results/counts.txt deseq2_results.csv
    """
}

process limma_voom_analysis {
    container 'biocontainers/r-ver:v3.6.1'
    input:
    path 'rsem_results/counts.txt'

    output:
    path 'limma_voom_results.csv'

    script:
    """
    Rscript limma_voom_analysis.R rsem_results/counts.txt limma_voom_results.csv
    """
}

process edger_analysis {
    container 'biocontainers/r-ver:v3.6.1'
    input:
    path 'rsem_results/counts.txt'

    output:
    path 'edger_results.csv'

    script:
    """
    Rscript edger_analysis.R rsem_results/counts.txt edger_results.csv
    """
}

// Define workflow execution
workflow {
    reads = params.reads

    qc_results = fastqc(reads)  
    trimmed_reads = trimmomatic(qc_results.out)  
    assembled_transcriptome = trinity_assembly(trimmed_reads.out)

    quantified_reads = rsem_analysis(assembled_transcriptome.out, trimmed_reads.out)
    predicted_orfs = transdecoder_prediction(assembled_transcriptome.out)
    annotation_results = trinnotate_annotation(assembled_transcriptome.out, quantified_reads.out)
    functional_annotation = kobas_annotation(predicted_orfs.out)

    de_analysis1 = deseq2_analysis(quantified_reads.out)
    de_analysis2 = limma_voom_analysis(quantified_reads.out)
    de_analysis3 = edger_analysis(quantified_reads.out)
}
