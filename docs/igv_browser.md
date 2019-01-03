# IGV Genome Browser

## Templates

```bash
for d in img js css;do
    [ -d "templates/igv/$d" ] || mkdir -p "templates/igv/$d"
done
wget -O templates/igv/img/favicon.ico 'https://igv.org/web/img/favicon.ico'
wget -O templates/igv/js/igv.min.js 'https://igv.org/web/release/2.1.0/dist/igv.min.js'
wget -O templates/igv/css/bootstrap.min.css 'https://maxcdn.bootstrapcdn.com/bootstrap/4.0.0/css/bootstrap.min.css'
wget -O templates/igv/js/bootstrap.min.js 'https://maxcdn.bootstrapcdn.com/bootstrap/4.0.0/js/bootstrap.min.js'
wget -O templates/igv/js/jquery-3.2.1.slim.min.js 'https://code.jquery.com/jquery-3.2.1.slim.min.js'
wget -O templates/igv/js/popper.min.js 'https://cdnjs.cloudflare.com/ajax/libs/popper.js/1.12.9/umd/popper.min.js'
```

## Genome reference from IGV

```bash
wget -O genome/hg38/source/igv.hg38.genome http://s3.amazonaws.com/igv.broadinstitute.org/genomes/hg38.genome
unzip genome/hg38/source/igv.hg38.genome -d genome/hg38/igv
```

## Generate IGV html

```bash
bin/create_igv.py -i templates/igv/main.html -g hg38 --genome-dir genome/hg38/ \
    --track genome/hg38/bed/long_RNA.intron.bed \
    --track genome/hg38/bed/rmsk.bed \
    --track genome/hg38/bed/promoter.merged.bed \
    --track genome/hg38/bed/enhancer.merged.bed \
    --track output/scirep/bigwig/Sample_1S1.genome.+.bigWig \
    --track output/scirep/bigwig/Sample_1S1.genome.-.bigWig \
    -o igv.html

--track output/scirep/bigwig/Sample_1S1.genome.+.bigWig -o igv.html
```