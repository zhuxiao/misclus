## Process this file with automake to produce Makefile.in
AUTOMAKE_OPTIONS=foreign 
ACLOCAL_AMFLAGS = -I m4
##CXX=$(g++ std=c++11)
CXXFLAGS=-O0 -g -Wall -fmessage-length=0 -std=c++11
nobase_include_HEADERS=*.h

bin_PROGRAMS=misclus
misclus_SOURCES=alnDataLoader.cpp Base.cpp Block.cpp covLoader.cpp InsertSizeEst.cpp util_header.cpp
misclus_SOURCES+=events.cpp Paras.cpp misEval_main.cpp util.cpp RefSeqLoader.cpp Evaluate.cpp miswork.cpp
misclus_CXXFLAGS=$(AM_CXXFLAGS)

##lib_LIBRARIES=libparser.a
##libparser_a_SOURCES=alnDataLoader.cpp Base.cpp Block.cpp covLoader.cpp 
##libparser_a_SOURCES+=events.cpp Paras.cpp parser_main.cpp util.cpp RefSeqLoader.cpp


##lib_LTLIBRARIES=libmiseval.la
##libmiseval_la_SOURCES=$(miseval_SOURCES)
##libmiseval_la_LDFLAGS=-release 0.1.0
##alnDataLoader.cpp Base.cpp Block.cpp covLoader.cpp 
##libparser_la_SOURCES+=events.cpp Paras.cpp parser_main.cpp util.cpp RefSeqLoader.cpp
