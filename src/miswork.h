#ifndef MISWORK_H_
#define MISWORK_H_

#include "covLoader.h"
#include "Block.h"

typedef struct regionNode{
	string chrname;
	long int startPos = 0, endPos = 0;
	long int chrlen=0;
}region;
struct warnregion{
	region interval;
	long int startaln, endaln;
};
struct indicators{//
	double covscore[3], indelscore[3], clipscore, insertscore, strandscore, matescore;
};
struct kind{
	int flag[7] = {0, 0, 0, 0, 0, 0, 0};
};
struct sortflag{
	float score;
	int order;
};

class miswork {
public:

	indicators scores = {0};

	vector<region> abmate, abstrand, abisize;//mark clump
	vector<bam1_t*> alnDataVec;//data_loader.lodataAlndata-358
	Base *basearr;
//	int64_t reg_id;
	region reg, evareg;
//	int blankregion, sum, count, nindel, basecov;
//	double *abReadsRatio;
//	double covscore=0, indelscore=0, minLocov;
	int64_t startPos_tmp =0;
//	double meaninsertsize;
	double mininsertsize, maxinsertsize, meancov, mincov, maxcov,chimeriCoef=1, randomCoef = 0.2;
//	uint64_t sequence_num; // nothing used
	float exRegFold = 2, minCovFold = 0.5, maxCovFold = 2, minLocRatio = 0.2, minStrandRatio = 0.3, minisizeRatio = 0.3, IsizeSdevFold = 3, minMateRatio = 0.3;
	faidx_t *fai;
	string bamFile;//regionfile

	miswork(region &r1, region &r2, faidx_t *fai, string bamFile);
	virtual ~miswork();

	void getindicator();
//	indicators getindicator();
	void preData(region &r);
	void clearData();
	double *getabreadsreatio(region &r);//change the *
	void getreadsMarkregion(region &r);
	void getMultiReads(region &r);
	bool isabstrandRead(bam1_t *b);
};

#endif /* MISWORK_H_ */
