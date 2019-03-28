#include <iostream>
#include <list>
#include <string>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <vector>
#include <stdexcept>
#include <limits>
#include <fstream>
#include <map>
#include <cstring>
#include <errno.h>
using namespace std;

typedef long int GenomicPosition;

class IOError: public std::exception
{
    std::string message;
    std::string filename;
public:
    explicit IOError(const std::string message, const std::string filename)
        : message(message), filename(filename) {}
    virtual ~IOError() throw () {}
    virtual const char* what()
    {
        return (message + " " + filename + " :" + std::string(std::strerror(errno))).c_str();
    }
};

vector<string> split_string(const std::string& str, char delimiter)
{
    std::istringstream iss(str);
    vector<string> tokens;
    string token;
    while(std::getline(iss, token, delimiter))
        tokens.push_back(token);
    return tokens;
}


std::map<std::string, GenomicPosition> read_chrom_sizes(const std::string& filename)
{
    std::map<std::string, GenomicPosition> chrom_sizes;
    ifstream chrom_sizes_f(filename.c_str());
    if(!chrom_sizes_f) {
        throw IOError("cannot open chrom sizes file", filename);
    }
    string line;
    while(!chrom_sizes_f.eof()){
    std::getline(chrom_sizes_f, line);
    if(line.empty())
        continue;
    istringstream is(line);
    string chrom, size;
    is >> chrom >> size;
    if(!is)
        throw std::invalid_argument("invalid format in chrom sizes file");
    chrom_sizes[chrom] = atoi(size.c_str());
    }
    chrom_sizes_f.close();
    return chrom_sizes;
}

struct GenomicInterval
{
    string chrom;
    GenomicPosition start;
    GenomicPosition end;
    string name;
    string score;
    char strand;

    GenomicInterval(const string& chrom = "", GenomicPosition start = -1, GenomicPosition end = -1,
        const std::string& name = "X", const std::string& score = "0", char strand = '+') 
        : chrom(chrom), start(start), end(end), name(name), score(score), strand(strand) {}

    static GenomicInterval from_bed(const std::string& s, int min_columns = 3)
    {
        GenomicInterval interval;
        vector<string> c = split_string(s, '\t');
        if(c.size() < min_columns)
            throw std::invalid_argument("not enough columns found in BED file");
        interval.chrom = c[0];
        interval.start = std::strtol(c[1].c_str(), NULL, 10);
        interval.end = std::strtol(c[2].c_str(), NULL, 10);
        if(c.size() >= 6)
        {
            interval.name = c[3];
            interval.score = c[4];
            interval.strand = c[5][0];
        }
        return interval;
    }
};

struct Interval
{
    GenomicPosition start;
    GenomicPosition end;
    Interval(GenomicPosition start = -1, GenomicPosition end = -1) : start(start), end(end) {}
};

void output_reads(std::ostream& out, const std::list<Interval>& reads)
{
    bool first = true;
    for(std::list<Interval>::const_iterator it = reads.begin(); it != reads.end(); ++it)
    {
        if(first)
            first = false;
        else
            out << " ";
        out << '(' << it->start << ", " << it->end << ')';
    }
    out << '\n';
}

inline bool is_empty_interval(const Interval& x)
{
    return x.start >= x.end;
}

void output_coverage(std::ostream& out, std::list<Interval>& reads, std::string& chrom,
     GenomicPosition start, GenomicPosition end, GenomicPosition bin_size, char strand)
{
    GenomicPosition position = start;
    while((position < end) && !reads.empty())
    {
        //output_reads(out, reads);
        GenomicPosition next_position = std::min(position + bin_size, end);
        out << chrom << '\t' << position << '\t' << next_position
                    << "\tX\t" << reads.size() << '\t' << strand << '\n';
        position = next_position;
        for(std::list<Interval>::iterator it = reads.begin(); it != reads.end(); ++it)
            it->start = position;
        reads.remove_if(is_empty_interval);
    }
}

void bin_coverage(std::istream& bed, std::ostream& out, 
    const std::map<string, GenomicPosition>& chrom_sizes, GenomicPosition bin_size)
{
    string line;
    string chrom;
    GenomicPosition chrom_size = 0;
    std::list<Interval> reads_pos, reads_neg;
    GenomicPosition gbin = 0;
    while(std::getline(bed, line))
    {
        //std::cout << line << std::endl;
        GenomicInterval interval = GenomicInterval::from_bed(line, 3);
        //std::cout << interval.chrom << '\t' << interval.start << '\t' << interval.end << std::endl;
        if(interval.chrom != chrom)
        {
            output_coverage(out, reads_pos, chrom, gbin*bin_size, chrom_size, bin_size, '+');
            output_coverage(out, reads_neg, chrom, gbin*bin_size, chrom_size, bin_size, '-');
            gbin = 0;
            chrom = interval.chrom;
            std::map<string, GenomicPosition>::const_iterator it = chrom_sizes.find(chrom);
            if(it != chrom_sizes.end())
                chrom_size = it->second;
            else
                throw std::invalid_argument(std::string("unknown chromosome found in reads: ") + chrom);
        }

        GenomicPosition rbin = interval.start/bin_size;
        // count a new bin
        if(rbin > gbin)
        {
            // output read count for previous bin
            output_coverage(out, reads_pos, chrom, gbin*bin_size, std::min(rbin*bin_size, chrom_size), bin_size, '+');
            output_coverage(out, reads_neg, chrom, gbin*bin_size, std::min(rbin*bin_size, chrom_size), bin_size, '-');
            gbin = rbin;
        }
        // not sorted
        else if(rbin < gbin)
        {
            throw std::invalid_argument("BED is not sorted");
        }
        // add read to list
        if(interval.strand == '+')
            reads_pos.push_back(Interval(interval.start, interval.end));
        else if(interval.strand == '-')
            reads_neg.push_back(Interval(interval.start, interval.end));
    }
}

int main(int argc, char** argv)
{
    if(argc != 5)
    {
        std::cerr << "Usage: " << argv[0] << " bed chrom_sizes output bin_size" << std::endl;
        exit(1);
    }
    std::string bed_file(argv[1]);
    std::string chrom_sizes_file(argv[2]);
    std::string output_file(argv[3]);
    GenomicPosition bin_size = std::strtol(argv[4], NULL, 10);

    std::ios::sync_with_stdio(false);
    std::cin.tie(NULL);
    
    std::ifstream bedf;
    if(bed_file != "-")
        bedf.open(bed_file.c_str());
    std::istream bed(bed_file == "-" ? std::cin.rdbuf() : bedf.rdbuf());

    std::map<std::string, GenomicPosition> chrom_sizes = read_chrom_sizes(chrom_sizes_file);

    std::ofstream outf;
    if(output_file != "-")
        outf.open(output_file.c_str());
    std::ostream output(output_file == "-" ? std::cout.rdbuf() : outf.rdbuf());
    
    bin_coverage(bed, output, chrom_sizes, bin_size);
    return 0;
}