#include <string>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <cmath>
#include <stdint.h>
using namespace std;

#define LINE_BUFFER_SIZE 1024*64
#define OUTPUT_BUFFER_SIZE 1024*1024

int nucleotide_table[256];

static void init_nucleotide_table()
{
    memset(nucleotide_table, 0, sizeof(int)*256);
    nucleotide_table[int('A')] = 0;
    nucleotide_table[int('T')] = 1;
    nucleotide_table[int('U')] = 1;
    nucleotide_table[int('C')] = 2;
    nucleotide_table[int('G')] = 3;
    nucleotide_table[int('N')] = 0;
}


struct SequenceComplexity
{
    int length;
    double entropy;
    double dust;
    SequenceComplexity(char* seq);
};

SequenceComplexity::SequenceComplexity(char* seq)
{
    length = strlen(seq);
    int counts[64] = {};
    for(int i = 2; i < length; i ++)
    {
        if((seq[i] == 0) || (seq[i] == '\r') || (seq[i] == '\n'))
        {
            length = i;
            seq[i] = 0;
            break;
        }
        int hash = 0;
        hash += nucleotide_table[seq[i - 2]];
        hash += nucleotide_table[seq[i - 1]]*4;
        hash += nucleotide_table[seq[i - 0]]*16;
        //printf("%d\n", length);
        //printf("hash[%d](%c%c%c) = %d\n", i, seq[i - 2], seq[i - 1], seq[i - 0], hash);
        counts[hash] ++;
    }
    entropy = 0;
    if(length > 3)
    {
        double r = 1.0/(length - 2);
        for(int i = 0; i < 64; i ++)
        {
            if(counts[i] > 0)
                entropy -= counts[i]*r*log(counts[i]*r);
        }
        entropy = entropy/log(length - 2)*100;
    }
    else
        entropy = 0.0;

    dust = 0;
    if(length > 3)
    {
        for(int i = 0; i < 64; i ++)
        {
            if(counts[i] > 0)
                dust += 0.5*counts[i]*(counts[i] - 1);
        }
        dust  = dust/(length - 2)/(length - 3)*200.0;
    }
    else
        dust = 100.0;
}

int main(int argc, char** argv)
{
    init_nucleotide_table();

    char* seq = new char[LINE_BUFFER_SIZE];

    size_t record_length = 0;
    int64_t lineno = 0;
    while(!feof(stdin))
    {
        lineno ++;
        char* line = fgets(seq, LINE_BUFFER_SIZE, stdin);
        if(!line)
            break;
        if(lineno%4 == 2)
        {
            SequenceComplexity c(seq);
            printf("%s\t%d\t%d\t%d\n", seq, int(c.length), int(c.entropy), int(c.dust));
        }
    }
    delete[] seq;

    return 0;
}