# IGV Genome Browser

## Reference genome configuration file

A template for reference configuration can be found in `templates/igv/config/genome.yaml`:

```yaml
genome: hg38
reference:
    id: hg38
    name: hg38
    fastaURL: genome/hg38/fasta/genome.fa
    indexURL: genome/hg38/fasta/genome.fa.fai
    cytobandURL: genome/hg38/igv/cytoBandIdeo.txt
tracks:
  GENCODE_V27:
      name: GENCODE_V27
      type: annotation
      format: bed
      url: genome/hg38/bed/gencode.bed
      indexURL: genome/hg38/bed/gencode.bed.idx
      displayMode: "EXPANDED"
      searchable: true
      visibilityWindow:  300000000
      height: 100
      show: true
  long_RNA_gene:
      name: long_RNA_gene
      type: annotation
      format: bed
      url: genome/hg38/bed/long_RNA.gene.bed
      indexURL: genome/hg38/bed/long_RNA.gene.bed.idx
      displayMode: "EXPANDED"
      searchable: true
      visibilityWindow:  300000000
      height: 100
      show: true
```

Two keys are required: `genome` and `reference`. The annotation tracks can be provided in BED, genePred, genePredExt, GTF or GFF format.

## Custom reference genome from FASTA file

### Human rRNA
```bash
# map NR_* ids to gene names
tr '|' $'\t' < "genome/hg38/source/refSeq_rRNA.ids.txt" > genome/hg38/source/refSeq_rRNA.gene_names.txt
# create reference
bin/create_igv.py create_reference --genome rRNA --name 'Human (rRNA)' \
    --gene-names genome/hg38/source/refSeq_rRNA.gene_names.txt \
    --fasta genome/hg38/fasta/rRNA.fa -o genome/hg38/igv/rRNA
```

The `create_reference` command generates a directory names `genome/hg38/igv/rRNA` that contains the following files:

| Filename | Description |
| -------- | ----------- |
| `reference.fa` | Reference genome sequences |
| `reference.fa.fai` | FASTA index for `reference.fa` |
| `config.yaml` | Track configuration file for creating IGV web page |
| `annotation.bed` | Annotation for each sequence in FASTA file |
| `annotation.genePred` | Annotation in genePred format |
| `cytoband.txt` | Cytoband file |


## Generate IGV HTML

### Configure web server

Setup a web server using Apache or other HTTP engines and set the base URL:
```bash
base_url="http://example.com/igv"
```

The directory structure of IGV should be like:
```
/genome
    hg38/
        igv/
        fasta/
            genome.fa
            genome.fai
        bed/
            gencode.bed
            gencode.bed.idx
            long_RNA.bed
            long_RNA.bed.idx
/igv
    config/
        ${dataset}_${map_step}.yaml
    html/
        ${dataset}_${map_step}.html
    
```

### Transcriptomic BigWig files (small RNA)
```bash
bin/create_igv.py generate_config \
    --sample-classes data/${dataset}/sample_classes.txt \
    --bigwig-pattern "${dataset}/{sample_id}.transcriptome.{strand}.bigWig" \
    --base-url "$base_url" \
    --max-samples-per-class 10 \
    --reference "templates/igv/config/genome.yaml" \
    -o genome_browser/config/${dataset}_transcriptome.yaml
bin/create_igv.py render -i templates/igv/main.html \
    -c genome_browser/config/${dataset}_transcriptome.yaml \
    -o genome_browser/igv/${dataset}_transcriptome.html
```

### BigWig files on custom reference genome (Long RNA)
```bash
map_step="rRNA"
bin/create_igv.py generate_config --locus "$locus" \
    --sample-classes data/${dataset}/sample_classes.txt \
    --bigwig-pattern "${dataset}/{sample_id}.${map_step}_rmdup.{strand}.bigWig" \
    --base-url "$base_url" \
    --strand "+" \
    --max-samples-per-class 10 \
    --reference "genome/hg38/igv/${map_step}/config.yaml" \
    -o genome_browser/config/${dataset}_${map_step}.yaml
```

## Create feature database

Create a feature database for searching genomic locus with feature names from a list of annotation files:

```bash
[ -d "igv/database" ] || mkdir -p "igv/database"
bin/web_server.py --build-database --genome hg38 \
    -i genome/hg38/gtf/gencode.gtf \
    -i genome/hg38/gtf/gencode_tRNA.gtf \
    -i genome/hg38/gff3/refseq.gff3 \
    -i genome/hg38/gtf/mitranscriptome_lncRNA.gtf \
    -i genome/hg38/gtf/mitranscriptome_tucp.gtf \
    -i genome/hg38/gff3/miRBase.gff3 \
    -o igv/database/hg38.pkl
```

## Start web server

Start a web server that listens on port 5000
```bash
bin/web_server.py --host 0.0.0.0 --port 5000 -i igv/database/hg38.pkl
```

Then navigate to `http://<server>:5000/igv/${dataset}_genome.html` to visit the genome browser.