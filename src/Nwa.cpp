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
            bool isMatch = ff1.GetSequence()[r] == ff2.GetSequence()[c];
            // Score of neighbor cell to the above left + (match [1] or mismatch score [-1])
            scores[ABOVELEFT] = scoreTable[r-1][c-1] + (isMatch ? 1 : -1);
            // Score of neighbor cell to the left + gap score [-1]
            scores[LEFT] = scoreTable[r][c-1] + -1;
            // Score of neighbor cell above + gap score [-1]
            scores[ABOVE] = scoreTable[r-1][c] + -1;
            // get max of values
            char maxDir = std::distance(scores, std::max_element(scores, scores + 3));
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
    
}


//for(int i=0; i<ff1.GetSequence().size()+1; i++)    //This loops on the rows.
//{
//    for(int j=0; j<ff2.GetSequence().size()+1; j++) //This loops on the columns
//    {
//        std::cout << scoreTable[i][j]  << "  ";
//        }
//        std::cout << std::endl;
//        }

