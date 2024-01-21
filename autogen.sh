#!/bin/sh
	if [ ! -d "bin" ]; then
		mkdir bin
	fi
	cd src/
	libtoolize
	aclocal
	autoheader
	automake --add-missing
	automake
	autoconf
	./configure
	make
	make mostlyclean-compile
	mv -f misclus ../bin
	cd ../
