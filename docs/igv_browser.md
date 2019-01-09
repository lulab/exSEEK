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

## IGV Browser reference

### igv.browser

addMultiLocusPanelWithGenomicStateAtIndex: ƒ (genomicState,index,viewportWidth)
addTrack: ƒ (track)
buildViewportsWithGenomicStateList: ƒ (genomicStateList)
cancelTrackPan: ƒ ()
compressedSession: ƒ ()
dispose: ƒ ()
emptyViewportContainers: ƒ ()
endTrackDrag: ƒ ()
findTracks: ƒ (property,value)
fireEvent: ƒ (eventName,args,thisObj)
getChromosomeLengthBP: ƒ (genome,referenceFrame)
goto: ƒ (chrName,start,end)
hideCenterGuide: ƒ ()
hideCursorGuide: ƒ ()
hideTrackLabels: ƒ ()
isMultiLocusMode: ƒ ()
isMultiLocusWholeGenomeView: ƒ ()
loadGenome: ƒ (idOrConfig,initialLocus)
loadInProgress: ƒ ()
loadSampleInformation: ƒ (url)
loadSession: ƒ (sessionURL,config)
loadTrack: ƒ (config)
loadTrackList: ƒ (configList)
minimumBases: ƒ ()
mouseDownOnViewport: ƒ (e,viewport)
on: ƒ (eventName,fn)
presentAlert: ƒ (alert,$parent,callback)
presentMessageWithCallback: ƒ (message,callback)
presentSplitScreenMultiLocusPanel: ƒ (alignment,genomicState)
removeAllTracks: ƒ (removeSequence)
removeMultiLocusPanelWithGenomicState: ƒ (genomicState,doResize)
removeTrack: ƒ (track)
removeTrackByName: ƒ (name)
renderSVG: ƒ (config)
reorderTracks: ƒ ()
resize: ƒ ()
search: ƒ (string,init)
selectMultiLocusPanelWithGenomicState: ƒ (selectedGenomicState)
sessionURL: ƒ ()
setTrackHeight: ƒ (newHeight)
setTrackLabelName: ƒ (trackView,name)
showCenterGuide: ƒ ()
showCursorGuide: ƒ ()
showTrackLabels: ƒ ()
startTrackDrag: ƒ (trackView)
toJSON: ƒ ()
un: ƒ (eventName,fn)
updateLocusSearchWidget: ƒ (genomicState)
updateTrackDrag: ƒ (dragDestination)
updateUIWithGenomicStateListChange: ƒ (genomicStateList)
updateViews: ƒ (genomicState,views)
updateZoomSlider: ƒ ($slider)
viewportContainerWidth: ƒ ()
viewportWidth: ƒ ()
visibilityChange: ƒ ()
zoom: ƒ (scaleFactor)
zoomIn: ƒ ()
zoomOut: ƒ ()
zoomWithRangePercentage: ƒ (percentage)
zoomWithScaleFactor: ƒ (centerBPOrUndefined,viewportOrUndefined,scaleFactor)

### File extensions

* narrowpeak
* broadpeak
* peaks
* bedgraph
* wig
* gff3
* gff
* gtf
* fusionjuncspan
* refflat
* seg
* aed
* bed
* vcf
* bb
* bigbed
* bw
* bigwig
* bam
* tdf
* refgene
* genepred
* genepredext
* bedpe
* bp
* snp
* rmsk


## Track configuration reference
```yaml
track1:
    type: wig
    format: bigwig
    height: 200
    min: 0
    max: 100
    color: 'rgb(0,0,0)'
    autoScale: true,
    logScale: false
```

```bash
long_RNA:
    url: genome/hg38/genePred/long_RNA.genePred
    type: annotation
    format: genepred
    height: 200

{
    echo 'reference:
    id: hg38
    fastaURL: genome/hg38/fasta/genome.fa
    indexURL: genome/hg38/fasta/genome.fa.fai
    cytobandURL: genome/hg38/igv/cytoBandIdeo.txt
Refseq Genes:
    type: annotation
    format: refGene
    url: genome/hg38/igv/refGene.txt
    visibilityWindow: 300000000
    height: 80
GENCODE_V27:
    type: annotation
    format: bed
    url: genome/hg38/igv/gencode.bed
    indexURL: genome/hg38/igv/gencode.bed.idx
    displayMode: "EXPANDED"
    searchable: true
    visibilityWindow: 300000000
    height: 200'
    for sample_id in $(head -n 50 data/scirep/sample_ids.txt);do
        printf "${sample_id}(+):
    url: scirep/${sample_id}.transcriptome.+.bigWig
    type: wig
    format: bigwig
    height: 25
${sample_id}(-):
    url: scirep/${sample_id}.transcriptome.-.bigWig
    type: wig
    format: bigwig
    height: 25
"
    done
} > genome_browser/config/scirep.yaml
bin/create_igv.py -i templates/igv/main.html -g hg38 \
    --tracklist genome_browser/config/scirep.yaml \
    --base-url http://166.111.156.12:10080/genomebrowse \
    -o genome_browser/igv/scirep.html
sync-bigwig
```

```bash
bin/create_igv.py -i templates/igv/main.html -g hg38 --genome-dir genome/hg38/ \
    --track genome/hg38/bed/long_RNA.intron.bed \
    --track genome/hg38/bed/rmsk.bed \
    --track genome/hg38/bed/promoter.merged.bed \
    --track genome/hg38/bed/enhancer.merged.bed \
    --track output/scirep/bigwig/Sample_1S1.genome.+.bigWig \
    --track output/scirep/bigwig/Sample_1S1.genome.-.bigWig \
    -o igv.html

bin/create_igv.py -i templates/igv/main.html -g hg38 --genome-dir genome/hg38/ \
    --tracklist tmp/tracklist.yaml \
    -o igv.html
--track output/scirep/bigwig/Sample_1S1.genome.+.bigWig -o igv.html

```

## Reference tracks
```bash
gtfToGenePred -geneNameAsName2 -genePredExt genome/hg38/gtf/gencode.gtf genome/hg38/genePred/gencode.genePredExt
gffread -E --bed genome/hg38/gtf/gencode.gtf -o /dev/stdout | bedtools sort > genome/hg38/igv/gencode.bed
./tools/IGVTools/igvtools index genome/hg38/igv/gencode.bed
```

## Domains
```bash
gffread --bed -o genome/hg38/bed/long_RNA.bed genome/hg38/gtf/long_RNA.gtf
gffread --bed -o genome/hg38/bed/gencode_tRNA.bed genome/hg38/gtf/gencode_tRNA.gtf

[ -d "output/scirep/gdomains/20/05.bed" ] || mkdir -p output/scirep/gdomains/20
{
    bin/tbed2gbed <(cat genome/hg38/bed/long_RNA.bed genome/hg38/bed/gencode_tRNA.bed) <(grep -v '^chr' output/scirep/domains/20/05.bed) /dev/stdout
    awk 'BEGIN{OFS="\t";FS="\t"}/^chr/{print $0,0,0,0,1,$3-$2,0}' output/scirep/domains/20/05.bed
} | bedtools sort > output/scirep/gdomains/20/05.bed
```

## Custom reference genome

* IGV genome repository: (https://s3.amazonaws.com/igv.org.genomes/genomes.json)
* Cytoband schema: [UCSC](http://rohsdb.cmb.usc.edu/GBshape/cgi-bin/hgTables?db=hg19&hgta_group=map&hgta_track=cytoBand&hgta_table=cytoBand&hgta_doSchema=describe+table+schema)

```bash
genome="rRNA"
fastaURL="genome/hg38/fasta/${genome}.fa"
fastaIndexURL="genome/hg38/fasta/${genome}.fa.fai"
cytobandURL="genome/hg38/cytoband/${genome}.txt"
refURL="genome/hg38/genePred/${genome}.genePred"
# convert NR_ ids to rRNA names
# gene_id|gene_name
rRNAIdFile="genome/hg38/source/refSeq_rRNA.ids.txt"

awk 'BEGIN{OFS="\t";FS="\t"}{print $1,0,$2,"p"NR".1","gneg"}' "$fastaIndexURL" > "$cytobandURL"
bedToGenePred <(awk 'BEGIN{OFS="\t";FS="\t"}NR==FNR{split($0,a,"|");name[a[1]]=a[2];next}{print $1,0,$2,name[$1],0,"+"}' "$rRNAIdFile" "$fastaIndexURL") "$refURL"

```

```bash
genome="circRNA"
fastaURL="genome/hg38/fasta/${genome}.fa"
fastaIndexURL="genome/hg38/fasta/${genome}.fa.fai"
cytobandURL="genome/hg38/cytoband/${genome}.txt"
refURL="genome/hg38/genePred/${genome}.genePred"

awk 'BEGIN{OFS="\t";FS="\t"}{print $1,0,$2,"p"NR".1","gneg"}' "$fastaIndexURL" > "$cytobandURL"
bedToGenePred <(awk 'BEGIN{OFS="\t";FS="\t"}{print $1,0,$2,$1,0,"+"}' "$fastaIndexURL") "$refURL"
```
