#!/bin/sh
	cd src/
	libtoolize
	aclocal
	autoheader
	automake --add-missing
	automake
	autoconf
	./configure
	make
	mv -f misclas ../bin
	cd ../
	ln -f -s $PWD/bin/misclas.py $PWD/misclas
