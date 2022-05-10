#include <iostream>
#include <string>
#include <htslib/faidx.h>

#include "Base.h"
#include "Evaluate.h"
#include "Paras.h"

using namespace std;

int main(int arg,char **argv){

	Paras *para = new Paras();
	para->refFile=argv[1];
	para->inBamFile=argv[2];
	para->checkBamFile();
	faidx_t *fai = fai_load(para->refFile.c_str());

	Evaluate eva(argv,fai,para);
	eva.extracfeature();
	delete para;
	fai_destroy(fai);

	return 0;
}
