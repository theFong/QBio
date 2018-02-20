//
//  Nwa.cpp
//  pa3
//
//  Created by Alec Fong on 2/19/18.
//
//

#include <numeric>
#include <vector>
#include <algorithm>
#include <map>
#include "Nwa.h"
#include<iostream>
#include <sstream>
#include <fstream>

Nwa::Nwa(const FastaFile &ff1,const FastaFile &ff2) : ff1(ff1),  ff2(ff2)
{
    // check some stuff?
}


void Nwa::SequenceAlign()
{
    
    enum Direction : char { ABOVELEFT,LEFT,ABOVE };
    
    std::vector<std::vector<short>> scoreTable(ff1.GetSequence().length()+1,std::vector<short>(ff2.GetSequence().length()+1,0));
    std::vector<std::vector<char>> directionTable(ff1.GetSequence().length()+1,std::vector<char>(ff2.GetSequence().length()+1,-1));
    
    // intialize row 0 and col 0 to not have to worry about corner case
    // rows
    for (int r = 0; r < scoreTable.size(); ++r)
    {
        scoreTable[r][0] = -r;
        directionTable[r][0] = ABOVE;
    }
    // cols
    for (int c = 0; c < scoreTable[0].size(); ++c)
    {
        scoreTable[0][c] = -c;
        directionTable[0][c] = LEFT;
    }
    
    // main loop
    for (int r = 1; r < scoreTable.size(); ++r)
    {
        for (int c = 1; c < scoreTable[0].size(); ++c)
        {
            // scores indices correlate with Direction enum
            short scores[3];
            bool isMatch = ff1.GetSequence()[r-1] == ff2.GetSequence()[c-1];
            // Score of neighbor cell to the above left + (match [1] or mismatch score [-1])
            scores[ABOVELEFT] = scoreTable[r-1][c-1] + (isMatch ? 1 : -1);
            // Score of neighbor cell to the left + gap score [-1]
            scores[LEFT] = scoreTable[r][c-1] + -1;
            // Score of neighbor cell above + gap score [-1]
            scores[ABOVE] = scoreTable[r-1][c] + -1;
            // get max of values
            Direction maxDir = scores[ABOVELEFT] >= scores[LEFT] ? ABOVELEFT : LEFT;
            maxDir = scores[maxDir] >= scores[ABOVE] ? maxDir : ABOVE;
            // set traceback direction
            directionTable[r][c] = maxDir;
            // set score table
            scoreTable[r][c] = scores[maxDir];
        }
    }
    this->score = scoreTable[ff1.GetSequence().size()][ff2.GetSequence().size()];
    
    int maxSize = std::max(ff1.GetSequence().size(), ff2.GetSequence().size());
    seq1.clear();
    seq1.reserve(maxSize);
    seq2.clear();
    seq2.reserve(maxSize);
    
    // traceback
    int r = ff1.GetSequence().size();
    int c = ff2.GetSequence().size();
    while(((r > 0) || (c > 0)))
    {
        switch (directionTable[r][c]) {
            case ABOVELEFT:
                seq1.push_back(ff1.GetSequence()[r-1]);
                seq2.push_back(ff2.GetSequence()[c-1]);
                --r;
                --c;
                break;
                
            case LEFT:
                seq1.push_back('_');
                seq2.push_back(ff2.GetSequence()[c-1]);
                --c;
                break;
            
            case ABOVE:
                seq1.push_back(ff1.GetSequence()[r-1]);
                seq2.push_back('_');
                --r;
                break;
                
            default:
                // do nothing
                break;
        }
    }
    
    std::reverse(seq1.begin(), seq1.end());
    std::reverse(seq2.begin(), seq2.end());
    
}

void Nwa::Write()
{
    std::ofstream out("match.result", std::ios::out|std::ios::trunc);
    if (out.is_open())
    {
        out << "A: " << ff1.GetHeader() << std::endl;
        out << "B: " << ff2.GetHeader() << std::endl;
        out << "Score: " << this->score << std::endl;
        
        enum Line : char { SEQ1, SEQ2, SIML };
        Line toggle = SEQ2
        ;
        
        int seq1Count = 0;
        int seq2Count = 0;
        int simlCount = 0;
        
        for (int i = 0; i < seq1.size()*3; ++i) {
            bool over = seq1Count >= seq1.size() || seq2Count >= seq2.size() || simlCount >= seq1.size();
            if((i % 70) ==  0 || over)
            {
                if (over)
                {
                    if (seq1Count >= seq1.size())
                    {
                        seq1Count = 0;
                    }
                    else if (seq2Count >= seq2.size())
                    {
                        seq2Count =  0;
                    }
                    else if(simlCount >= seq1.size())
                    {
                        simlCount = 0;
                    }
                }
                // set toggle
                switch (toggle) {
                    case SEQ1:
                        toggle = SIML;
                        break;
                        
                    case SEQ2:
                        toggle = SEQ1;
                        break;
                        
                    case SIML:
                        toggle = SEQ2;
                        break;
                        
                    default:
                        break;
                }

                out << std::endl;
            }
            switch (toggle) {
                case SEQ1:
                    out << seq1[seq1Count];
                    ++seq1Count;
                    break;
                    
                case SEQ2:
                    out << seq2[seq2Count];
                    ++seq2Count;
                    break;
                    
                case SIML:
                    out << (seq1[simlCount] == seq2[simlCount] ? '|' : ' ');
                    ++simlCount;
                    break;
                    
                default:
                    // do nothing
                    break;
            }
            
        }
        out << std::endl << std::endl;
    }
}


//
//for(int i=0; i<ff1.GetSequence().size()+1; i++)    //This loops on the rows.
//{
//    for(int j=0; j<ff2.GetSequence().size()+1; j++) //This loops on the columns
//    {
//        std::cout << scoreTable[i][j]  << "  ";
//        }
//        std::cout << std::endl;
//        }
//        
//        for(int i=0; i<ff1.GetSequence().size()+1; i++)    //This loops on the rows.
//        {
//            for(int j=0; j<ff2.GetSequence().size()+1; j++) //This loops on the columns
//            {
//                printf("%x ",directionTable[i][j]);
//            }
//            std::cout << std::endl;
//        }
