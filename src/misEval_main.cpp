#include <iostream>
#include <string>
#include <htslib/faidx.h>

#include "Base.h"
#include "Evaluate.h"
#include "Paras.h"

using namespace std;

int main(int argc,char **argv){

	Paras *para = new Paras(argc, argv);

	Evaluate eva(para);
	eva.extracfeature();

	eva.clusterfile();
	eva.analysfile();
	delete para;

	return 0;
}
