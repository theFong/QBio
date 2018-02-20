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
#include <unordered_map>
#include <numeric>
#include <vector>
#include <algorithm>
#include <map>

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
        this->header.erase(0,1);
        this->sequence.reserve(static_cast<unsigned int>(size) - this->header.length());
        std::string buf;
        while(std::getline(file,buf))
        {
            this->sequence.append(buf);
        }
        file.close();
        
    }
}


std::unordered_map<char, int> FastaFile::AminoAcidCount()
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
    
    std::unordered_map<char, int> aminoCount =
    {
        {'A',0},
        {'R',0},
        {'N',0},
        {'D',0},
        {'C',0},
        {'Q',0},
        {'E',0},
        {'G',0},
        {'H',0},
        {'I',0},
        {'L',0},
        {'K',0},
        {'M',0},
        {'F',0},
        {'P',0},
        {'S',0},
        {'T',0},
        {'W',0},
        {'Y',0},
        {'V',0}
    };
    
    std::unordered_map<char, int> baseToInd =
    {
        {'T',0},
        {'C',1},
        {'A',2},
        {'G',3},
    };
    
    char state = 0;
    for (char c : sequence) {
        char nextInstr = stateMachine[int(state)][baseToInd[c]];
        
        // change state
        if(nextInstr < 24)
        {
            state = nextInstr;
        } else
        {
            state = 3;
            ++aminoCount[nextInstr];
        }
    }
    
    return aminoCount;
}

void FastaFile::WriteAminoCount(const std::unordered_map<char, int> &count)
{
    std::unordered_map<char, std::string> stringify =
    {
        {'A',"(A) Alanine: "},
        {'C',"(C) Cysteine: "},
        {'D',"(D) Aspartic acid: "},
        {'E',"(E) Glutamic acid: "},
        {'F',"(F) Phenylalanine: "},
        {'G',"(G) Glycine: "},
        {'H',"(H) Histidine: "},
        {'I',"(I) Isoleucine: "},
        {'K',"(K) Lysine: "},
        {'L',"(L) Leucine: "},
        {'M',"(M) Methionine: "},
        {'N',"(N) Asparagine: "},
        {'P',"(P) Proline: "},
        {'Q',"(Q) Glutamine: "},
        {'R',"(R) Arginine: "},
        {'S',"(S) Serine: "},
        {'T',"(T) Threonine: "},
        {'V',"(V) Valine: "},
        {'W',"(W) Tryptophan: "},
        {'Y',"(Y) Tyrosine: "}
    };
    
    auto countCopy = count;
    int sum = std::accumulate(countCopy.begin(), countCopy.end(),0,
                              [](const auto &a,const auto &b)
                              {
                                  return a + b.second;
                              });

    std::map<char,int> ordered(countCopy.begin(), countCopy.end());
    
    std::ofstream out("amino.txt", std::ios::out|std::ios::trunc);
    if (out.is_open())
    {
        out << header << std::endl;
        out << "Total amino acids produced: " << sum << std::endl;
        for (auto i = ordered.begin(); i!= ordered.end(); ++i)
        {
            out << stringify[(*i).first] << (*i).second << ((std::next(i) != ordered.end()) ? "\n" : "");
        }
        
    }
    
    
}

