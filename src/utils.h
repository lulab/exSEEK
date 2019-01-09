#ifndef __UTILS_H__
#define __UTILS_H__

#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <cstring>

inline std::vector<std::string> str_split(const std::string& s, char delim='\t')
{
    std::istringstream ss(s);
    std::string token;
    std::vector<std::string> a;
    while (std::getline(ss, token, delim)) {
        a.push_back(token);
    }
    return a;
}

template <class T>
inline std::string vector_to_string(const std::vector<T>& a, char delim)
{
    std::ostringstream ss;
    for(size_t i = 0; i < a.size(); i ++)
    {
        if(i > 0)
            ss << delim;
        ss << a[i];
    }
    return ss.str();
}

inline bool next_fasta(std::istream& in, std::string& id, std::string& comment, std::string& seq)
{
    std::string line;
    if(!std::getline(in, line))
        return false;
    if(line[0] == '>')
    {
        size_t pos = 0;
        for(pos = 1; pos < line.size(); pos ++)
        {
            if(std::isspace(line[pos]))
                break;
        }
        if(pos == 0)
            throw std::invalid_argument("empty FASTA header found");
        id = line.substr(1, pos - 1);
        if(pos < (line.size() - 1))
            comment = line.substr(pos + 1);
    }
    else
        throw std::invalid_argument("FASTA header not found");
    seq = "";
    while(std::getline(in, line))
    {
        if(line.empty())
            throw std::invalid_argument("empty line found");
        if(line[0] == '>')
            break;
    }
    return true;
}

#endif // utils.h