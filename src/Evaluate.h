#ifndef EVALUATE_H_
#define EVALUATE_H_

#include <iostream>
#include <cstring>
#include <random>
#include <fstream>
#include <numeric>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <sys/stat.h>

#include "Block.h"
#include "covLoader.h"
#include "InsertSizeEst.h"
#include "util.h"
#include "miswork.h"

#include <thread>

using namespace std;

#define MAX_PAIR_NUM		200000
#define MAX_ISIZE			2000

class Evaluate{
	public:
		bool covflag = true, indelflag = true, abstrandflag = true, abisizeflag = true, abmateflag = true;
		int minRegsize, num_threads;
		double exRegFold = 2, minCovFold = 0.5, maxCovFold = 2, minLocRatio = 0.2, minStrandRatio = 0.3, minisizeRatio = 0.3, IsizeSdevFold = 3, minMateRatio = 0.3;

		string outdir;


		string cov_filename = "cov.csv";
		string indel_filename= "indel.csv";
		string strand_filename = "strand.csv";
		string insert_filename = "insert.csv";
		string mate_filename = "mate.csv";

//		string cov_filename, indel_filename, strand_filename, insert_filename, mate_filename;

//		string cov_filename_path =outdir +"/"+ cov_filename;
//		string indel_filename_path =outdir +"/" + indel_filename ;
//		string strand_filename_path =outdir +"/" + strand_filename;
//		string insert_filename_path =outdir +"/" + insert_filename;
//		string mate_filename_path =outdir +"/" + mate_filename;

		string cov_cluster_filename, indel_cluster_filename, strand_cluster_filename, insert_cluster_filename, mate_cluster_filename;
		string cluster_stat_filename, cluster_detail_filename;
		string final_result_filename;

		long int regions_len;//restore the number of vector<>region
		double meaninsertsize, mininsertsize, maxinsertsize, meancov, mincov, maxcov,chimeriCoef=1, randomCoef;

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
		vector<vector<int>> abnormalpos;
		string regionfile, bamFile;

		vector<miswork*> miswork_vec;
		pthread_mutex_t mtx_misworkDone_num;
		int32_t miswork_num, misworkDone_num;
		

		Evaluate(Paras *paras);
		virtual ~Evaluate();
		void init();
		void getinregion(string &filename);
		void splitstring(string &origin, char flag, string &front, string &after);
		long int getChrLen(string &chrname);
		void getEvaRegions(string &fai);
		static bool compareContigLen(region &a,region &b);
		bool isRegExist(region &reg, vector <region> &vec);
		void getranregion();
		Base* initBaseinfo(region &r);
		int loadAlnData();
		bool isabstrandRead(bam1_t *b);
		void getreadsMarkregion(region &r);
		void getMultiReads(region &r);
		void preData(region &r);
		void clearData();
		double *abnormalreads(region &r);
		double* getabreadsreatio(region &r);//caculate indicators in region
//		void getindicator(uint64_t s);
		void extracfeature();
		void getMisworkVec();
		void clearMisworkVec();
		void outputfile();
		void clusterfile();
		vector<vector<double>> formatfile(string& filename);
		void analysfile();
};

#endif /* EVALUATE_H_ */
