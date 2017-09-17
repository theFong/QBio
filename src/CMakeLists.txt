# If you create new headers/cpp files, add them to these list!
set(HEADER_FILES
	SrcMain.h
)

set(SOURCE_FILES
	SrcMain.cpp
)

# Don't change this
add_library(src ${SOURCE_FILES} ${HEADER_FILES})
