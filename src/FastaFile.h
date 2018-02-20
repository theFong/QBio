//
//  Fasta.hpp
//  pa3
//
//  Created by Alec Fong on 2/18/18.
//
//
#pragma once

#include <string>
#include <stdio.h>
#include <unordered_map>

class FastaFile
{
    void Read(const std::string& source);
    std::string mHeader;
    std::string mSequence;
    
    public:
    FastaFile(const std::string& source);
    std::unordered_map<char, int> AminoAcidCount();
    std::string GetHeader();
    std::string GetSequence();
    void Write(const std::unordered_map<char, int> &count);
};

