//
//  Fasta.cpp
//  pa3
//
//  Created by Alec Fong on 2/18/18.
//
//

#include "FastaFile.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>

FastaFile::FastaFile(const std::string& source)
{
    Read(source);
}

std::string FastaFile::GetSequence()
{
    return sequence;
}

std::string FastaFile::GetHeader()
{
    return header;
}

void FastaFile::Read(const std::string& source)
{
    // Requires <fstream>
    std::ifstream::pos_type size;
    // Open the file for input, in binary mode, and start ATE (at the end)
    std::ifstream file (source, std::ios::in|std::ios::ate);
    if (file.is_open())
    {
        size = file.tellg(); // Save the size of the file
        file.seekg(0, std::ios::beg); // Seek back to start of file
        std::getline(file, this->header);
        this->sequence.reserve(static_cast<unsigned int>(size) - this->header.length());
        std::string buf;
        while(std::getline(file,buf))
        {
            this->sequence.append(buf);
        }
        file.close();
        
    }
}


int FastaFile::AminoAcidCount()
{
    int stateMachine[24][4] =
    {
    //St| T | C | A | G |
    //------------------------
    //0 | 0 | 0 | 1 | 0 | // Looking for start codon
        { 0,  0,  1,  0 },
    //1 | 2 | 0 | 0 | 0 |
        { 2,  0,  1,  0 },
    //2 | 0 | 0 | 0 | M | // M and start
        { 0,  0,  1, 'M'},
    //3 | 4 | 9 |14 |19 | // first letter of amino acid
        { 4,  9, 14, 19 },
    //4                   // T??
        { 5,  6,  7,  8 },
    //5                   // TT?
        {'F','F','L','L'},
    //6                   // TC?
        {'S','S','S','S'},
    //7                   // TA?
        {'Y','Y', 0, 0 },
    //8                   // TG?
        {'C','C', 0, 'W'},
    //9                   // C??
        {10, 11, 12, 13 },
    //10                  // CT?
        {'L','L','L','L'},
    //11                  // CC?
        {'P','P','P','P'},
    //12                  // CA?
        {'H','H','Q','Q'},
    //13                  // CG?
        {'R','R','R','R'},
    //14                  // A??
        { 15, 16, 17, 18},
    //15                  // AT?
        {'I','I','I','M'},
    //16                  // AC?
        {'T','T','T','T'},
    //17                  // AA?
        {'N','N','K','K'},
    //18                  // AG?
        {'S','S','R','R'},
    //19                  // G??
        { 20, 21, 22, 23},
    //20                  // GT?
        {'V','V','V','V'},
    //21                  // GC?
        {'A','A','A','A'},
    //22                  // GA?
        {'D','D','E','E'},
    //23                  // GG?
        {'G','G','G','G'}
    };
}

