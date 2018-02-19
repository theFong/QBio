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

class FastaFile
{
    void Read(const std::string& source);
    std::string header;
    std::string sequence;
    
    public:
    FastaFile(const std::string& source);
    FastaFile(const std::string& source1, const std::string& source2);
    int AminoAcidCount();
    std::string GetHeader();
    std::string GetSequence();
};

