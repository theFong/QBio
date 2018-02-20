#include "SrcMain.h"
#include <iostream>
#include "FastaFile.h"

void ProcessCommandArgs(int argc, const char* argv[])
{
    FastaFile ff(argv[1]);
	if(argc == 2)
    {
        auto res = ff.AminoAcidCount();
        ff.WriteAminoCount(res);
    }
    else if(argc == 3)
    {
        // sequence alignment
    }
}
