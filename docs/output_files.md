## Small RNA mapping

| File name | Descrpition |
| --------- | ----------- |
| `snakemake/sequential_mapping.snakemake` | Snakefile for sequential mapping. Required by snakemake/mapping_small.snakemake |
| `${output_dir}/cutadapt/${sample_id}.fastq` | Reads with adaptor trimmed |
| `${output_dir}/tbam/${sample_id}/${rna_type}.bam` | BAM files in transcript coordinates |
| `${output_dir}/gbam/${sample_id}/${rna_type}.bam` | BAM files in genome coordinates |
| `${output_dir}/unmapped/${sample_id}/${rna_type}.fa.gz` | Unmapped reads in each step |
| `${output_dir}/fastqc/${sample_id}_fastqc.html` | FastQC report file |
| `${output_dir}/summary/fastqc.html` | Summary report for FastQC (HTML) |
| `${output_dir}/summary/fastqc.txt`  | Summary table for FastQC |
| `${output_dir}/summary/fastqc.ipynb` | Summary report for FastQC (Jupyter notebook) |
| `${output_dir}/summary/read_counts.txt` | Summary table for read counts |
| `${output_dir}/stats/mapped_read_length_by_sample/${sample_id}` | Length distribution of mapped reads |

## Long RNA mapping

| File name | Descrpition |
| --------- | ----------- |
| `${output_dir}/cutadapt/${sample_id}.fastq` | Reads with adaptor trimmed |
| `${output_dir}/bam/${sample_id}/rRNA.bam` | BAM files for reads mapped to rRNA |
| `${output_dir}/bam/${sample_id}/genome.bam` | BAM files for reads mapped to genome |
| `${output_dir}/bam/${sample_id}/remove_duplicates.bam` | BAM files for reads after removing duplicates |
| `${output_dir}/bam/${sample_id}/circRNA.bam` | BAM files for reads after removing duplicates |
| `${output_dir}/unmapped/${sample_id}/${map_step}_1.fa.gz` | Unmapped reads in each step |
| `${output_dir}/fastqc/${sample_id}_fastqc.html` | FastQC report file |
| `${output_dir}/summary/read_counts.txt` | Summary table for read counts |
| `${output_dir}/stats/mapped_read_length_by_sample/${sample_id}` | Length distribution of mapped reads |
| `${output_dir}/stats/mapped_insert_size_by_sample/${sample_id}` | Length distribution of mapped reads |

## Count matrix (small RNA-seq)

| File name | Descrpition |
| --------- | ----------- |
| `${output_dir}/count_matrix/transcript.txt` | Count matrix of transcripts |
| `${output_dir}/count_matrix/htseq.txt` | Count matrix of genes generated using HTSeq-count |
| `${output_dir}/count_matrix/featurecounts.txt` | Count matrix of genes generated using featureCounts |
| `${output_dir}/counts_by_biotype/${count_method}/${sample_id}/${rna_type}` | Gene/transcript counts generated using a feature counting tool |

## Long RNA domains

| File name | Descrpition |
| --------- | ----------- |
| `${output_dir}/domain_counts/${bin_size}/${pvalue}/${sample_id}.bed` | Read counts in long RNA domains (BED format with read counts in Column 5 |
| `${output_dir}/count_matrix/domain_${pvalue}.txt` | Read count matrix of long RNA domains |
| `${output_dir}/domains/${bin_size}/${pvalue}.bed` | Long RNA domain locations |
| `${output_dir}/domains_recurrence/${bin_size}/${pvalue}.bed` | Recurrence of long RNA domains among samples (Column 5) |

## Matrix processing

| File name | Description |
| --------- | ----------- |
| `${output_dir}/normalized_matrix/${normalization_method}.${imputation_method}.${batch_removal_method}.txt` | 
| `${output_dir}/matrix_processing/normalization/${normalization_method}.txt` |
| `${output_dir}/matrix_processing/imputation/${normalization_method}.${imputation_method}.txt` |
| `${output_dir}/matrix_processing/batch_removal/${batch_removal_method}.${batch_index}.txt` |

## Differential expression

| File name | Description |
| --------- | ----------- |
| `${output_dir}/differential_expression/${count_method}/${compare_group}/${diffexp_method}.txt` | 

