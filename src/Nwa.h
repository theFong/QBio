//
//  Nwa.hpp
//  pa3
//
//  Created by Alec Fong on 2/19/18.
//
//
#pragma once

#include <cstdio>
#include "FastaFile.h"

class Nwa
{
    FastaFile mFf1;
    FastaFile mFf2;
    short mScore;
    std::string mSeq1;
    std::string mSeq2;
    
public:
    Nwa(const FastaFile &ff1,const FastaFile &ff2);
    void SequenceAlign();
    void Write();
};
