//
//  Nwa.hpp
//  pa3
//
//  Created by Alec Fong on 2/19/18.
//
//
#pragma once

#include <stdio.h>
#include "FastaFile.h"

class Nwa
{
    FastaFile ff1;
    FastaFile ff2;
    short score;
    std::string seq1;
    std::string seq2;
    
public:
    Nwa(const FastaFile &ff1,const FastaFile &ff2);
    void SequenceAlign();
    void Write();
};
