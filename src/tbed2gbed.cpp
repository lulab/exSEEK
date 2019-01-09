#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <sstream>
#include <fstream>
#include <stdexcept>
#include <cstdlib>
#include <algorithm>
#include <stdint.h>

#include "utils.h"
#include "formats.h"

bool tbed2gbed_record(const Bed12Record& anno, const Bed6Record& tbed, Bed12Record& gbed)
{
    if(tbed.chrom != anno.name)
        throw std::invalid_argument("different chromosome in tbed and anno: " + tbed.chrom + ", " + anno.name);
    int64_t tChromStart = tbed.chromStart;
    int64_t tChromEnd = tbed.chromEnd;
    if(anno.strand == '-')
    {
        int64_t tx_length = 0;
        for(int i = 0; i < anno.blockCount; i ++)
            tx_length += anno.blockSizes[i];
        tChromStart = tx_length - tbed.chromEnd;
        tChromEnd = tx_length - tbed.chromStart;
    }
    gbed.chrom = anno.chrom;
    //gbed.name = tbed.name + "|" + anno.name;
    gbed.name = tbed.name;
    gbed.score = tbed.score;
    gbed.strand = anno.strand;
    gbed.thickStart = 0;
    gbed.thickEnd = 0;
    gbed.itemRgb = 0;

    int64_t tBlockStart = 0;
    // relative to gene
    int64_t gChromStart = 0;
    int start_exon = 0;
    for(start_exon = 0; start_exon < anno.blockCount; start_exon ++)
    {
        if((tBlockStart <= tChromStart) && (tChromStart < (tBlockStart + anno.blockSizes[start_exon])))
        {
            gChromStart = anno.blockStarts[start_exon] + tChromStart - tBlockStart;
            gbed.chromStart = anno.blockStarts[start_exon] + anno.chromStart + tChromStart - tBlockStart;
            break;
        }
        tBlockStart += anno.blockSizes[start_exon];
    }
    for(int end_exon = start_exon; end_exon < anno.blockCount; end_exon ++)
    {
        gbed.chromEnd = anno.blockStarts[end_exon] + anno.chromStart + tChromEnd - tBlockStart;
        int64_t blockStart = std::max(tChromStart, tBlockStart) - tBlockStart + anno.blockStarts[end_exon] - gChromStart;
        int64_t blockEnd = std::min(tChromEnd, tBlockStart + int64_t(anno.blockSizes[end_exon])) - tBlockStart + anno.blockStarts[end_exon] - gChromStart;
        gbed.blockSizes.push_back(blockEnd - blockStart);
        gbed.blockStarts.push_back(blockStart);

        if((tBlockStart < tChromEnd) && (tChromEnd <= (tBlockStart + anno.blockSizes[end_exon])))
            break;
        tBlockStart += anno.blockSizes[end_exon];
    }
    gbed.blockCount = gbed.blockSizes.size();
    return true;
}

void tbed2gbed(const std::string& tbed_file, 
    const std::string& gbed_file, 
    const std::string& anno_file)
{
    std::map<std::string, Bed12Record> anno = read_bed12(anno_file);
    std::ifstream tbed(tbed_file.c_str());
    if(!tbed)
        throw std::runtime_error("cannot open tbed file: " + tbed_file);
    std::ofstream gbed(gbed_file.c_str());
    if(!gbed)
        throw std::runtime_error("cannot open gbed file: " + gbed_file);
    for(std::string line; std::getline(tbed, line); )
    {
        if(!line.empty())
        {
            Bed6Record tbed_record(line);
            Bed12Record gbed_record;
            std::map<std::string, Bed12Record>::iterator it = anno.find(tbed_record.chrom);
            if(it != anno.end())
            {
                if(tbed2gbed_record(it->second, tbed_record, gbed_record))
                    gbed << gbed_record.to_string() << std::endl;
            }
            else
                std::cerr << "Warning: cannot find transcript_id: " << tbed_record.chrom << std::endl;
        }
    }
    tbed.close();
    gbed.close();
}

int main(int argc, char** argv)
{
    if(argc != 4)
    {
        std::cerr << "Usage: " << argv[0] << " anno_file tbed_file gbed_file" << std::endl;
        std::cerr << "    anno_file: input transcript annotation in BED12 format\n"
                  << "        (may use gffread --bed to convert GTF to BED12 format)" << std::endl;
        std::cerr << "    tbed_file: input intervals in transcript coordinates (BED format)" << std::endl;
        std::cerr << "    gbed_file: output intervals in genomic coordinates (BED12 format)" << std::endl;
        std::exit(1);
    }
    std::string anno_file(argv[1]);
    std::string tbed_file(argv[2]);
    std::string gbed_file(argv[3]);

    tbed2gbed(tbed_file, gbed_file, anno_file);
    
    return 0;
}