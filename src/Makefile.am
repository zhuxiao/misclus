## Process this file with automake to produce Makefile.in
AUTOMAKE_OPTIONS=foreign 
ACLOCAL_AMFLAGS = -I m4
##CXX=$(g++ std=c++11)
CXXFLAGS=-O0 -g -Wall -fmessage-length=0 -std=c++11
nobase_include_HEADERS=*.h

bin_PROGRAMS=misclus
misclus_SOURCES=alnDataLoader.cpp Base.cpp Block.cpp covLoader.cpp InsertSizeEst.cpp util_header.cpp
misclus_SOURCES+=events.cpp Paras.cpp misclus_main.cpp util.cpp RefSeqLoader.cpp Evaluate.cpp miswork.cpp
misclus_CXXFLAGS=$(AM_CXXFLAGS)
