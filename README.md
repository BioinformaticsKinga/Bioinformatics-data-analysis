# RNA-Seq de novo Assembly Pipeline

# Project Overview
This project outlines a comprehensive RNA-seq analysis pipeline designed for de novo transcriptome assembly and downstream functional annotation, as well as differential expression analysis. By integrating key bioinformatics tools, this pipeline provides a streamlined approach to analyze raw RNA-seq data, from quality control to gene expression quantification and functional insights.

# Key Objectives
Ensure high-quality RNA-seq data through preprocessing and quality control.

Assemble RNA-seq reads into full transcript sequences using de novo methods.

Quantify gene expression levels to assess transcript abundance.

Identify coding sequences within the transcriptome to predict open reading frames (ORFs).

Annotate transcripts to provide functional insights related to gene function and biological pathways.

Identify differentially expressed genes across different experimental conditions.

# Pipeline Stages and Tools
1. Quality Control (FastQC)
Raw RNA-seq reads are assessed for quality using FastQC, which helps identify potential issues such as adapter contamination, base quality, or sequence duplication. This step ensures that the data is suitable for downstream processing.

2. Read Trimming (Trimmomatic)
Trimmomatic is used to remove adapter sequences and trim low-quality bases from the RNA-seq reads, enhancing the quality of the data for transcript assembly.

3. Transcriptome Assembly (Trinity)
The cleaned RNA-seq reads are assembled into transcript sequences using Trinity, a widely used de novo transcriptome assembler. This step reconstructs full-length transcripts from short sequencing reads without the need for a reference genome.

4. Expression Quantification (RSEM)
Once the transcriptome is assembled, RSEM is employed to quantify the abundance of each transcript. This tool calculates the expression levels of the assembled transcripts, enabling the identification of highly expressed genes.

5. ORF Prediction (TransDecoder)
To identify protein-coding regions within the transcriptome, TransDecoder is used to predict open reading frames (ORFs). This step helps identify which transcripts have the potential to encode functional proteins.

6. Functional Annotation (Trinotate & KOBAS)
The assembled transcripts are annotated using Trinotate and KOBAS, which provide functional insights into gene functions, associated pathways, and possible biological roles. This annotation step links the transcript data to known gene functions, pathways, and potential diseases.

7. Differential Expression Analysis (DESeq2, Limma-voom, edgeR)
Differential expression analysis is conducted using multiple tools: DESeq2, Limma-voom, and edgeR. These tools are used to detect genes that are differentially expressed across different experimental conditions, identifying those that are upregulated or downregulated under specific treatments or conditions.

# Results and Findings
Data quality was ensured through rigorous quality control and trimming, which improved the reliability of subsequent analyses.

Transcript assembly generated a comprehensive set of full-length transcript sequences, offering a detailed representation of the transcriptome.

Gene expression was quantified, providing insights into the relative abundance of different transcripts across conditions.

Coding sequences were predicted, identifying potential protein-coding genes within the transcriptome.

Functional annotation provided valuable insights into the biological roles of the identified transcripts, linking them to known pathways and functions.

Differential expression analysis revealed genes with significant changes in expression between experimental conditions, highlighting potential candidates for further study.

# Technologies and Tools Used
Programming Languages: Python, R

Quality Control & Preprocessing: FastQC, Trimmomatic

Transcriptome Assembly: Trinity

Expression Quantification: RSEM

ORF Prediction: TransDecoder

Functional Annotation: Trinotate, KOBAS

Differential Expression Analysis: DESeq2, Limma-voom, edgeR


