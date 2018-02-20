#include "SrcMain.h"
#include <iostream>
#include "FastaFile.h"
#include "Nwa.h"

void ProcessCommandArgs(int argc, const char* argv[])
{
    FastaFile ff1(argv[1]);
	if(argc == 2)
    {
        // amino count
        auto res = ff1.AminoAcidCount();
        ff1.Write(res);
    }
    else if(argc == 3)
    {
        // sequence alignment
        FastaFile ff2(argv[2]);
        Nwa geneticSimilarty(ff1,ff2);
        geneticSimilarty.SequenceAlign();
        geneticSimilarty.Write();
    }
}
