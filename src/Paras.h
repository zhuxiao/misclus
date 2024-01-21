#ifndef SRC_PARAS_H_
#define SRC_PARAS_H_

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <getopt.h>
#include <unistd.h>

#include "structures.h"

using namespace std;

// program variables
#define PROG_NAME					"misClus"
#define PROG_DESC					"Misassembly Clustering"
#define PROG_VERSION				"0.5.0"

#define MIN_REG_SIZE				200
#define EXREG_FOLD					2
#define MIN_COV_FOLD				(0.5f)
#define MAX_COV_FOLD				2
#define MIN_LOCAL_RATIO				(0.2f)
#define MIN_STRAND_RATIO			(0.3f)
#define MIN_SIZE_RATIO				(0.3f)
#define ISIZE_DEV_FOLD				3
#define MIN_MATE_RATIO				(0.5f)
#define RAND_NORM_REG_PERCENT		(0.2f)

#define OUT_DIR						"output"



#define SIZE_EST_OP					0
#define NUM_EST_OP					1
#define SNV_EST_OP					2

#define LARGE_INDEL_SIZE_FACTOR		3


// default parameter values
#define BLOCKSIZE  					1000000
#define SLIDESIZE  					500
#define ASSEM_SLIDE_SIZE			5000

#define MIN_INDEL_EVENT_SIZE		2
#define MIN_INDEL_EVENT_NUM			5
#define MIN_SV_SIZE_USR				2

#define MIN_CLIP_READS_NUM_THRES	7	// to be parameterized

#define MAX_CLIP_REG_SIZE			10000  // to be parameterized

#define SIZE_PERCENTILE_EST			(0.95f)
#define NUM_PERCENTILE_EST			(0.99995f)
#define AUX_ARR_SIZE				1001

#define MAX_REG_SUM_SIZE_EST		500000


#define ASSEMBLY_SUCCESS			"ASS_SUCCESS"
#define ASSEMBLY_FAILURE			"ASS_FAILURE"
#define ALIGN_SUCCESS				"ALN_SUCCESS"
#define ALIGN_FAILURE				"ALN_FAILURE"
#define CALL_SUCCESS				"CALL_SUCCESS"
#define CALL_FAILURE				"CALL_FAILURE"

#define DONE_STR					"DONE"

#define SAMPLED_STR					"SAMPLED"
#define UNSAMPLED_STR				"UNSAMPLED"

#define LIMIT_REG_ALL_STR			"ALL"

#define EXPECTED_COV_ASSEMBLE		(30.0f)
#define MASK_VAL_DEFAULT			0

#define NUM_PARTS_PROGRESS			50
#define NUM_THREADS_PER_WORK		0  // 0: unspecify the limited number of threads for each work


#define MAX_ULTRA_HIGH_COV_THRES	300		// maximal coverage threshold for ultra-high coverage


// program parameters
class Paras
{
	public:
		// user/system defined parameters
		string command, refFile, inBamFile, outFilePrefix; //, canu_version;
		string outDir;
		string regionFile;
		int minRegsize;
		double exRegFold , minmateradio;
		double minCovFold, maxCovFold, minLocRatio, minStrandRatio, minisizeRatio, isizeSdevFold;
		bool cov_flag, indel_flag, abstrand_flag, abisize_flag, abmate_flag;
		double random_norm_reg_percent;
		string sample;
		int32_t num_threads;

		int num_parts_progress;
		vector<simpleReg_t*> limit_reg_vec;
//		bool maskMisAlnRegFlag;
		bool load_from_file_flag;
		bool limit_reg_process_flag = false;	// true for limit regions; default is false for disable limit regions (i.e. process all regions)

	public:
		Paras();
		Paras(int argc, char **argv);
		virtual ~Paras();
//		void initEst();
//		void estimate(size_t op_est);
		void outputParas();
//		void outputEstParas(string info);
		int checkBamFile();

	private:
		void init();
		void show_version();
		//string getCanuVersion();
		int parseParas(int argc, char **argv);
		int parseMisclasParas(int argc, char **argv);
		int parse_long_opt(int32_t option_index, const char *optarg, const struct option *lopts);
		void showUsage();
		size_t estimateSinglePara(size_t *arr, size_t n, double threshold, size_t min_val);
};

#endif /* SRC_PARAS_H_ */
