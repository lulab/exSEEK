#ifndef __FORMATS_H__
#define __FORMATS_H__

#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <map>

struct Bed12Record
{
    // fields
    std::string chrom;
    int64_t chromStart;
    int64_t chromEnd;
    std::string name;
    std::string score;
    char strand;
    int thickStart;
    int thickEnd;
    int itemRgb;
    int64_t blockCount;
    std::vector<int64_t> blockSizes;
    std::vector<int64_t> blockStarts;
    std::vector<std::string> extras;
    // methods
    Bed12Record() {}
    Bed12Record(const std::string& s);
    std::string to_string();
};

std::map<std::string, Bed12Record> read_bed12(const std::string& filename);

struct Bed6Record
{
    // fields
    std::string chrom;
    int64_t chromStart;
    int64_t chromEnd;
    std::string name;
    std::string score;
    char strand;
    std::vector<std::string> extras;
    // methods
    Bed6Record() {}
    Bed6Record(const std::string& s);
    std::string to_string();
};



#endif