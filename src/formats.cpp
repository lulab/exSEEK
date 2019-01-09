#include <sstream>
#include <map>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <stdexcept>

#include "formats.h"
#include "utils.h"

Bed12Record::Bed12Record(const std::string& s)
{
    std::vector<std::string> c = str_split(s, '\t');
    if(c.size() < 12)
        throw std::invalid_argument("expect at least 12 columns in BED12 format: " + s);
    chrom = c[0];
    chromStart = std::atoi(c[1].c_str());
    chromEnd = std::atoi(c[2].c_str());
    name = c[3];
    score = c[4];
    strand = c[5][0];
    blockCount = std::atoi(c[9].c_str());

    std::vector<std::string> blockSizesStr = str_split(c[10], ',');
    blockSizes.resize(blockCount);
    for(size_t i = 0; i < blockCount; i ++)
        blockSizes[i] = std::atoi(blockSizesStr[i].c_str());
    
    std::vector<std::string> blockStartsStr = str_split(c[11], ',');
    blockStarts.resize(blockCount);
    for(size_t i = 0; i < blockCount; i ++)
        blockStarts[i] = std::atoi(blockStartsStr[i].c_str());
    
    if(c.size() > 12)
    {
        extras.resize(c.size() - 12);
        for(size_t i = 12; i < c.size(); i ++)
            extras[i - 12] = c[i];
    }
}

std::string Bed12Record::to_string()
{
    std::ostringstream os;
    os << chrom << '\t'
       << chromStart << '\t'
       << chromEnd << '\t'
       << name << '\t'
       << score << '\t'
       << strand << '\t'
       << thickStart << '\t'
       << thickEnd << '\t'
       << itemRgb << '\t'
       << blockCount << '\t'
       << vector_to_string(blockSizes, ',') << '\t'
       << vector_to_string(blockStarts, ',')
       << vector_to_string(extras, '\t');
    return os.str();
}

std::map<std::string, Bed12Record> read_bed12(const std::string& filename)
{
    std::ifstream fin(filename.c_str());
    if(!fin)
        throw std::runtime_error("cannot open the bed12 file: " + filename);
    std::map<std::string, Bed12Record> bed12;
    for(std::string line; std::getline(fin, line); )
    {
        if(!line.empty())
        {
            Bed12Record record(line);
            bed12[record.name] = record;
        }
    }
    fin.close();
    return bed12;
}

Bed6Record::Bed6Record(const std::string& s)
{
    std::vector<std::string> c = str_split(s, '\t');
    if(c.size() < 6)
        throw std::invalid_argument("expect at least 6 columns in BED6 format: " + s);
    chrom = c[0];
    chromStart = std::atoi(c[1].c_str());
    chromEnd = std::atoi(c[2].c_str());
    name = c[3];
    score = c[4];
    strand = c[5][0];
    if(c.size() > 6)
    {
        extras.resize(c.size() - 6);
        for(size_t i = 6; i < c.size(); i ++)
            extras[i - 6] = c[i];
    }
}

std::string Bed6Record::to_string()
{
    std::ostringstream os;
    os << chrom << '\t'
       << chromStart << '\t'
       << chromEnd << '\t'
       << name << '\t'
       << score << '\t'
       << strand << '\t'
       << vector_to_string(extras, '\t');
    return os.str();
}