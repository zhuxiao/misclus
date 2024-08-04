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
	double covscore[3], indelscore[3], IDCscore, insertscore, strandscore, matescore;
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

	Paras *paras;

	indicators scores = {0};

	vector<region> abmate, abstrand, abisize;//mark clump
	vector<bam1_t*> alnDataVec;//data_loader.lodataAlndata-358
	Base *basearr;
	region reg, evareg;
	faidx_t *fai;

	double chimeriCoef;

	double mininsertsize, maxinsertsize;

	double meancov, mincov, maxcov;

	miswork(region &r1, region &r2, faidx_t *fai, Paras *paras, double mininsertsize, double maxinsertsize);
	virtual ~miswork();

	void getindicator();
	void preData(region &r);
	void clearData();
	double *getabreadsreatio(region &r);//change the *
	void getreadsMarkregion(region &r);
	void getMultiReads(region &r);
	bool isabstrandRead(bam1_t *b);
};

#endif /* MISWORK_H_ */
