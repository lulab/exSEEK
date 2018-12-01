#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <cctype>
#include "generic_stream.h"

using namespace std;

#define LINE_BUFFER_SIZE 1024*64
#define OUTPUT_BUFFER_SIZE 1024*1024

int main(int argc, char** argv)
{
    char* buffer = new char[LINE_BUFFER_SIZE];
    map<string, int> counts;
    string cur_seq;

    size_t record_length = 0;
    int64_t lineno = 0;
    while(!feof(stdin))
    {
        lineno ++;
        char* line = fgets(buffer, LINE_BUFFER_SIZE, stdin);
        if(!line)
            break;
        if(lineno%4 == 2)
        {
            cur_seq = buffer;
            int i = cur_seq.size() - 1;
            while((i >= 0) && isspace(cur_seq[i]))
                i --;
            cur_seq.resize(i + 1);
            counts[cur_seq] ++;
        }
    }

    FILE* fout = fdopen(fileno(stdout), "wb");
    char* output_buffer = new char[OUTPUT_BUFFER_SIZE];
    setbuffer(fout, output_buffer, OUTPUT_BUFFER_SIZE);
    size_t readnum = 1;
    for(map<string, int>::iterator it = counts.begin(); it != counts.end(); ++it)
    {
        //fprintf(fout, ">%u|%d\n", readnum, it->second);
        //fprintf(fout, "%s\n", it->first.c_str());
        fprintf(fout, "%d\t%s\n", it->second, it->first.c_str());
        readnum ++;
    }
    delete[] buffer;

    return 0;
}