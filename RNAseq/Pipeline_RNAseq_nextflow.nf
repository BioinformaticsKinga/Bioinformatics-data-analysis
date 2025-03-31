// Define input parameters and set the appropriate paths for your files
params.reads = 'zygmukin_reads/*.fastq.gz' // Path to raw sequencing data (FASTQ files)
params.outdir = 'zygmukin_results' // Directory to store results
params.adapters = 'TruSeq3-PE.fa' // Adapter sequences for trimming 

// Process 1: FastQC (Quality Control)
// Perform quality control on the raw reads to assess their quality 
process fastqc {
    container 'biocontainers/fastqc:v0.11.9' // Docker image for FastQC tool
    input:
    path reads // Input FASTQ files

    output:
    path 'fastqc_reports' // Directory to store FastQC output

    script:
    """
    fastqc -o fastqc_reports $reads
    """
}

// Process 2: Trimmomatic (Read Trimming)
// Trim the raw reads to remove low-quality bases and adapter sequences using Trimmomatic
process trimmomatic {
    container 'biocontainers/trimmomatic:v0.39' // Docker image for Trimmomatic
    input:
    path reads // Input paired-end reads
    path params.adapters // Adapter sequences for trimming

    output:
    path 'trimmed_reads' // Output directory for trimmed reads

    script:
    """
    trimmomatic PE -threads 4 ${reads[0]} ${reads[1]} \
        trimmed_reads/1P.fastq.gz trimmed_reads/1U.fastq.gz \
        trimmed_reads/2P.fastq.gz trimmed_reads/2U.fastq.gz \
        ILLUMINACLIP:${params.adapters}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    """
}

// Process 3: Trinity Assembly (De novo assembly)
// Perform a de novo assembly of the transcriptome using Trinity, an RNA-Seq assembly tool
process trinity_assembly {
    container 'trinityrnaseq/trinityrnaseq:2.13.2' // Docker image for Trinity RNA-Seq assembly tool
    input:
    path reads // Input trimmed reads

    output:
    path 'trinity_output' // Output directory for Trinity assembly

    script:
    """
    Trinity --seqType fq --max_memory 50G --left ${reads[0]} --right ${reads[1]} --CPU 4 --output trinity_output
    """
}

// Process 4: RSEM (RNA-Seq Quantification)
// Use RSEM to quantify gene expression levels from the assembled transcriptome
process rsem_analysis {
    container 'biocontainers/rsem:v1.3.1' // Docker image for RSEM (RNA-Seq expression tool)
    input:
    path 'trinity_output/Trinity.fasta' // Assembled transcriptome (FASTA file)
    path 'trimmed_reads/1P.fastq.gz' // Trimmed read pair 1
    path 'trimmed_reads/2P.fastq.gz' // Trimmed read pair 2

    output:
    path 'rsem_results' // Directory to store RSEM results

    script:
    """
    rsem-calculate-expression --paired-end \
        --bowtie2 --alignments \
        --output-genome-bam --estimate-rspd \
        --strand-specific rsem_results \
        trimmed_reads/1P.fastq.gz trimmed_reads/2P.fastq.gz trinity_output/Trinity.fasta
    """
}

// Process 5: TransDecoder (ORF Prediction)
// Predict open reading frames (ORFs) from the assembled transcriptome using TransDecoder
process transdecoder_prediction {
    container 'biocontainers/transdecoder:v5.5.0' // Docker image for TransDecoder
    input:
    path 'trinity_output/Trinity.fasta' // Assembled transcriptome (FASTA file)

    output:
    path 'transdecoder_output' // Directory to store TransDecoder results

    script:
    """
    TransDecoder.LongOrfs -t trinity_output/Trinity.fasta
    TransDecoder.Predict -t trinity_output/Trinity.fasta
    mv trinity_output/Trinity.fasta.transdecoder.* transdecoder_output
    """
}

// Process 6: Trinotate Annotation
// Annotate the predicted proteins using Trinotate, integrating information from various sources (e.g., Gene Ontology, protein domains)
process trinnotate_annotation {
    container 'trinotate/trinotate:latest' // Docker image for Trinotate annotation tool
    input:
    path 'trinity_output/Trinity.fasta' // Assembled transcriptome (FASTA file)
    path 'rsem_results/counts.txt' // RSEM quantification file

    output:
    path 'trinotate_results' // Directory to store Trinotate annotation results

    script:
    """
    ./trinotate_annotation.sh trinity_output/Trinity.fasta rsem_results/counts.txt trinotate_results
    """
}

// Process 7: KOBAS Annotation (Functional Annotation)
// Perform functional annotation using KOBAS, a tool for annotating proteins with functional information
process kobas_annotation {
    container 'agbase/kobas:3.0.3_3' // Docker image for KOBAS annotation tool
    input:
    path 'transdecoder_output/Trinity.fasta.transdecoder.pep' // Predicted proteins (peptide sequences)

    output:
    path 'kobas_results' // Directory to store KOBAS annotation results

    script:
    """
    kobas-annotate -i transdecoder_output/Trinity.fasta.transdecoder.pep -s ath -t fasta:nuc -o kobas_results
    """
}

// Process 8: Differential Expression Analysis (DESeq2, Limma-voom, edgeR)
// Perform differential gene expression analysis using three different methods: DESeq2, Limma-voom, and edgeR
process deseq2_analysis {
    container 'biocontainers/r-ver:v3.6.1' // Docker image for R (DESeq2 analysis)
    input:
    path 'rsem_results/counts.txt' // RSEM quantification file

    output:
    path 'deseq2_results.csv' // Output DESeq2 results

    script:
    """
    Rscript deseq2_analysis.R rsem_results/counts.txt deseq2_results.csv
    """
}

process limma_voom_analysis {
    container 'biocontainers/r-ver:v3.6.1' // Docker image for R (Limma-voom analysis)
    input:
    path 'rsem_results/counts.txt' // RSEM quantification file

    output:
    path 'limma_voom_results.csv' // Output Limma-voom results

    script:
    """
    Rscript limma_voom_analysis.R rsem_results/counts.txt limma_voom_results.csv
    """
}

process edger_analysis {
    container 'biocontainers/r-ver:v3.6.1' // Docker image for R (edgeR analysis)
    input:
    path 'rsem_results/counts.txt' // RSEM quantification file

    output:
    path 'edger_results.csv' // Output edgeR results

    script:
    """
    Rscript edger_analysis.R rsem_results/counts.txt edger_results.csv
    """
}

// Define workflow execution
workflow {
    reads = params.reads // Load raw sequencing data

    // Quality control step
    qc_results = fastqc(reads)  

    // Trimming step to clean up reads
    trimmed_reads = trimmomatic(qc_results.out)  

    // Perform de novo assembly of the transcriptome
    assembled_transcriptome = trinity_assembly(trimmed_reads.out)

    // Quantify gene expression using RSEM
    quantified_reads = rsem_analysis(assembled_transcriptome.out, trimmed_reads.out)

    // Predict ORFs using TransDecoder
    predicted_orfs = transdecoder_prediction(assembled_transcriptome.out)

    // Annotate the assembled transcriptome using Trinotate
    annotation_results = trinnotate_annotation(assembled_transcriptome.out, quantified_reads.out)

    // Perform functional annotation using KOBAS
    functional_annotation = kobas_annotation(predicted_orfs.out)

    // Perform differential expression analysis using multiple methods
    de_analysis1 = deseq2_analysis(quantified_reads.out)
    de_analysis2 = limma_voom_analysis(quantified_reads.out)
    de_analysis3 = edger_analysis(quantified_reads.out)
}
