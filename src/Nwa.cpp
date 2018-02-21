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

Nwa::Nwa(const FastaFile &ff1,const FastaFile &ff2) : mFf1(ff1),  mFf2(ff2), mScore(0)
{
    // check some stuff?
}


void Nwa::SequenceAlign()
{
    
    enum Direction : char { ABOVELEFT,LEFT,ABOVE };
    
    std::vector<std::vector<short>> scoreTable(mFf2.GetSequence().length()+1,std::vector<short>(mFf1.GetSequence().length()+1,0));
    std::vector<std::vector<char>> directionTable(mFf2.GetSequence().length()+1,std::vector<char>(mFf1.GetSequence().length()+1,-1));
    
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
    
//     main loop
    for (int r = 1; r < scoreTable.size(); ++r)
    {
        for (int c = 1; c < scoreTable[0].size(); ++c)
        {
            // scores indices correlate with Direction enum
            short scores[3];
            bool isMatch = mFf2.GetSequence()[r-1] == mFf1.GetSequence()[c-1];
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
    
    this->mScore= scoreTable[mFf2.GetSequence().size()][mFf1.GetSequence().size()];
    
    int maxSize = std::max(mFf2.GetSequence().size(), mFf1.GetSequence().size());
    mSeq1.clear();
    mSeq1.reserve(maxSize);
    mSeq2.clear();
    mSeq2.reserve(maxSize);

    // traceback
    int r = mFf2.GetSequence().size();
    int c = mFf1.GetSequence().size();
    while(((r > 0) || (c > 0)))
    {
        switch (directionTable[r][c]) {
            case ABOVELEFT:
                mSeq2.push_back(mFf2.GetSequence()[r-1]);
                mSeq1.push_back(mFf1.GetSequence()[c-1]);
                --r;
                --c;
                break;
                
            case LEFT:
                mSeq2.push_back('_');
                mSeq1.push_back(mFf1.GetSequence()[c-1]);
                --c;
                break;
            
            case ABOVE:
                mSeq2.push_back(mFf2.GetSequence()[r-1]);
                mSeq1.push_back('_');
                --r;
                break;
                
            default:
                // do nothing
                break;
        }
    }
    
    std::reverse(mSeq1.begin(), mSeq1.end());
    std::reverse(mSeq2.begin(), mSeq2.end());
    
}

//void Nwa::Write()
//{
//    std::ofstream out("match.result", std::ios::out|std::ios::trunc);
//    if (out.is_open())
//    {
//        out << "A: " << mFf1.GetHeader() << std::endl;
//        out << "B: " << mFf2.GetHeader() << std::endl;
//        out << "Score: " << this->mScore;
//        
//        enum Line : char { SEQ1, SEQ2, SIML };
//        // will change to SEQ1
//        Line toggle = SEQ2;
//        
//        int seq1Count = 0;
//        int seq2Count = 0;
//        int simlCount = 0;
//        
//        for (int i = 0; i < mSeq1.size()*3; ++i) {
//            bool over = seq1Count >= mSeq1.size() || seq2Count >= mSeq2.size() || simlCount >= mSeq1.size();
//            if((i % 70) ==  0 || over)
//            {
//                if (over)
//                {
//                    if (seq1Count >= mSeq1.size())
//                    {
//                        seq1Count = 0;
//                    }
//                    else if (seq2Count >= mSeq2.size())
//                    {
//                        seq2Count =  0;
//                    }
//                    else if(simlCount >= mSeq1.size())
//                    {
//                        simlCount = 0;
//                    }
//                }
//                // set toggle
//                switch (toggle) {
//                    case SEQ1:
//                        toggle = SIML;
//                        break;
//                        
//                    case SEQ2:
//                        toggle = SEQ1;
//                        out << std::endl;
//                        break;
//                        
//                    case SIML:
//                        toggle = SEQ2;
//                        break;
//                        
//                    default:
//                        break;
//                }
//
//                out << std::endl;
//            }
//            switch (toggle) {
//                case SEQ1:
//                    out << mSeq1[seq1Count];
//                    ++seq1Count;
//                    break;
//                    
//                case SEQ2:
//                    out << mSeq2[seq2Count];
//                    ++seq2Count;
//                    break;
//                    
//                case SIML:
//                    out << (mSeq1[simlCount] == mSeq2[simlCount] ? '|' : ' ');
//                    ++simlCount;
//                    break;
//                    
//                default:
//                    // do nothing
//                    break;
//            }
//            
//        }
//        out << std::endl << std::endl;
//    }
//}

std::string Nwa::BuildMatch(const std::string &a,const std::string &b)
{
    std::string matchStr = "";
    for(int i = 0; i < a.length(); i++)
    {
        if(a[i] == b[i])
        {
            matchStr.push_back('|');
        }
        else
        {
            matchStr.push_back(' ');
        }
    }
    return matchStr;
}

void Nwa::Write()
{
    std::ofstream out ("match.result", std::ios::out|std::ios::trunc);
    if (out.is_open())
    {
        // write out header
        out << "A: " <<  mFf1.GetHeader() << std::endl;
        out << "B: " <<  mFf2.GetHeader() << std::endl;
        out << "Score: " << mScore;
        out << std::endl << std::endl;
        std::string line1 = "";
        std::string matchLine = "";
        std::string line2 = "";
        
        int i = 0;
        while(i < mSeq1.length())
        {
            int lineLen = 70;
            if(i + 70 > mSeq1.length())
            {
                lineLen = mSeq1.length()  - i;
            }
            line1 = mSeq1.substr(i, lineLen);
            line2 = mSeq2.substr(i, lineLen);
            matchLine = BuildMatch(line1,line2);
            
            out << line1 << std::endl;
            out << matchLine << std::endl;
            out << line2 << std::endl;
            out << std::endl;
            
            i = i + lineLen;
        }
        out.close();
    }
}


//
//for(int i=0; i<mFf1.GetSequence().size()+1; i++)    //This loops on the rows.
//{
//    for(int j=0; j<mFf2.GetSequence().size()+1; j++) //This loops on the columns
//    {
//        std::cout << scoreTable[i][j]  << "  ";
//        }
//        std::cout << std::endl;
//        }
//        
//        for(int i=0; i<mFf1.GetSequence().size()+1; i++)    //This loops on the rows.
//        {
//            for(int j=0; j<mFf2.GetSequence().size()+1; j++) //This loops on the columns
//            {
//                printf("%x ",directionTable[i][j]);
//            }
//            std::cout << std::endl;
//        }

