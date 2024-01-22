#include <iostream>
#include <string>
#include <htslib/faidx.h>

#include "Base.h"
#include "misCluster.h"
#include "Paras.h"

using namespace std;

int main(int argc,char **argv){

	Paras *para = new Paras(argc, argv);

	misCluster mis_cluster(para);

	cout << "Feature extraction ..." << endl;
	mis_cluster.extracfeature();

	cout << "Misassembly clustering ..." << endl;
	mis_cluster.clusterfile();

	cout << "Cluster information analyzing ..." << endl;
	mis_cluster.analysfile();

	cout << "All works are completed." << endl;

	delete para;

	return 0;
}
