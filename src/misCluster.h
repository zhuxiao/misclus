#ifndef MISCLUSTER_H_
#define MISCLUSTER_H_

#include <iostream>
#include <cstring>
#include <random>
#include <fstream>
#include <algorithm>
#include <iomanip>
#include <sys/stat.h>

#include "Block.h"
#include "covLoader.h"
#include "InsertSizeEst.h"
#include "util.h"
#include "miswork.h"

using namespace std;

#define MAX_PAIR_NUM		200000
#define MAX_ISIZE			2000

class misCluster{
	public:
		bool covflag, abstrandflag, abisizeflag, abmateflag, IDCflag;
//		bool indelflag;
		int32_t minRegsize, num_threads;
//		double exRegFold, minCovFold, maxCovFold, minLocRatio, minStrandRatio, minisizeRatio, IsizeSdevFold, minMateRatio;
		double exRegFold, IsizeSdevFold;
		string outdir;
		string cov_filename = "cov.csv";
		string indel_filename= "indel.csv";
		string strand_filename = "strand.csv";
		string insert_filename = "insert.csv";
		string mate_filename = "mate.csv";
		string IDC_filename = "IDC.csv";

		string cov_cluster_filename, indel_cluster_filename, strand_cluster_filename, insert_cluster_filename, mate_cluster_filename, IDC_cluster_filename;
		string cluster_stat_filename, cluster_detail_filename;
		string final_result_filename;

		int64_t regions_len;//restore the number of vector<>region
		double meaninsertsize;
		double randomCoef;
		double mininsertsize, maxinsertsize;

		Paras *paras;
		vector<indicators> scores;//indicator[0...5]coverage，indel，clip，insertsize，reversed，mate
		vector<kind> result;

		vector<string> inregions;
		vector<region> evaregions, regions, chrregions, ranregion, fairegion;//evaluated region，auxiliary region，merged region，ranregion，chromosal region
		vector<region> abmate, abstrand, abisize;//misworkExtra       mark clump
		faidx_t * fai;
		string fainame;
		Base* basearr;
		vector<bam1_t*> alnDataVec;
		vector<vector<int32_t>> abnormalpos;
		string regionfile, bamFile;

		vector<miswork*> miswork_vec;
		pthread_mutex_t mtx_misworkDone_num;
		int32_t miswork_num, misworkDone_num;

		misCluster(Paras *paras);
		virtual ~misCluster();
		void init();
		void getinregion(string &filename);
		void splitstring(string &origin, char flag, string &front, string &after);
		int64_t getChrLen(string &chrname);
		void getEvaRegions(string &fai);
		static bool compareContigLen(region &a,region &b);
		bool isRegExist(region &reg, vector <region> &vec);
		void getRandRegion();
		Base* initBaseinfo(region &r);
		int32_t loadAlnData();
		void extracfeature();
		void getMisworkVec();
		void clearMisworkVec();
		void outputfile();
		void clusterfile();
		vector<vector<double>> formatfile(string& filename);
		void analysfile();
};

#endif /* MISCLUSTER_H_ */
