#!/bin/sh
	mkdir bin
	cd src/
	libtoolize
	aclocal
	autoheader
	automake --add-missing
	automake
	autoconf
	./configure
	make
	mv -f misclas_bin ../bin
	cd ../
	
