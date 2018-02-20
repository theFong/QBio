# If you create new headers/cpp files, add them to these list!
set(HEADER_FILES
	SrcMain.h
    FastaFile.h
    Nwa.h
)

set(SOURCE_FILES
	SrcMain.cpp
    FastaFile.cpp
    Nwa.cpp
)

# Don't change this
add_library(src ${SOURCE_FILES} ${HEADER_FILES})
