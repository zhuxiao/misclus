#include "miswork.h"
#include "covLoader.h"
#include "Evaluate.h"
#include "string.h"
#include "util.h"
miswork::miswork(region &r1, region &r2, faidx_t *fai, string bamFile) {
	this->reg = r1;
	this->evareg = r2;
//	this->reg_id = reg_id;
	this->fai = fai;
	this->bamFile = bamFile;
//	blankregion = 0;
//	sum = 0;
//	count = 0;
//	basecov = 0;
	mininsertsize = 0;
	maxinsertsize = 0;
	meancov = 0;
	mincov = 0;
	maxcov = 0;
//	minLocov = 0;
//	nindel = 0;

	basearr = NULL;
//	abReadsRatio = NULL;

//	covscore=0;
//	indelscore=0;

}

miswork::~miswork() {//	minLocov = 0;
	// TODO Auto-generated destructor stub
}

void miswork::getindicator(){
//indicators miswork::getindicator(){
//	Evaluate preData(evareg);//
	preData(evareg);//
	indicators tmp;//miswork.tmp
	int blankregion=0, sum=0, count=0, nindel, basecov;
	double *abReadsRatio;//
	double covscore=0, indelscore=0, minLocov;
	count = 0;
	int64_t startPos_tmp = evareg.startPos;//
	for(int64_t i=reg.startPos; i<=reg.endPos; i++){ //int64_64
		if(basearr[i-startPos_tmp].coverage.idx_RefBase<4){
			sum = basearr[i-startPos_tmp].coverage.num_bases[5] + sum;
			count++;
		}
	}//if(count == 0
	if(count == 0){//	minLocov = 0;
//		cout << "count is :" << count << endl;
//		meancov = (reg.endPos - reg.startPos)*0.8;
		meancov = 0;
	}else{
		meancov = sum/count;
	}
	mincov = meancov*minCovFold;
	maxcov = meancov*maxCovFold;
	minLocov = meancov*minLocRatio;
//	cout << "regions.at(s).start :" << regions.at(s).startPos << " regions.at(s).end :" << regions.at(s).endPos << endl;

	for(int64_t i = reg.startPos; i <= reg.endPos; i++){		//basearrd:Base *baseArray = new Base[endPos-startPos+1]();
		if(basearr[i-startPos_tmp].coverage.idx_RefBase < 4){				//42217  42418
			basecov = basearr[i-startPos_tmp].coverage.num_bases[5] + basearr[i - startPos_tmp].del_num_from_del_vec;
			nindel = basearr[i-startPos_tmp].insVector.size() + basearr[i-startPos_tmp].del_num_from_del_vec;
			if(basecov > minLocov){
				if(indelscore < double(nindel)/basecov){
					indelscore = double(nindel)/basecov;
				}
			}
		}
	}
//	cout << "evaregions.at(s)start :" << evaregions.at(s).startPos << " evaregions.at(s)end :" << evaregions.at(s).endPos << endl;

	for(int64_t i= evareg.startPos; i <= evareg.endPos; i++){  //41815 42820
//		cout << basearr[i-regions.at(s).startPos] << endl；
		if(basearr[i-startPos_tmp].coverage.idx_RefBase < 4){
			basecov = basearr[i-startPos_tmp].coverage.num_bases[5];
			nindel = basearr[i-startPos_tmp].insVector.size() + basearr[i-startPos_tmp].del_num_from_del_vec;
			if(basecov < mincov or (basecov - 0.5*basearr[i-startPos_tmp].coverage.num_multiReads)>maxcov){
				covscore += 1;
			}
		}else{
			blankregion++;
		}
	}
	abReadsRatio = getabreadsreatio(evareg);
	if(blankregion < (evareg.endPos -evareg.startPos+1)){
		tmp.covscore[0] = chimeriCoef*(covscore)/(evareg.endPos - evareg.startPos-blankregion+1);
	}else{
		tmp.covscore[0] = 0;
	}
	tmp.indelscore[0] = indelscore;
	tmp.matescore = abReadsRatio[0];
	tmp.strandscore = abReadsRatio[1];
	tmp.insertscore = abReadsRatio[2];
	blankregion = 0;
	scores = tmp ;
//	delete abReadsRatio;
	delete[] abReadsRatio;
	clearData();
//	return NULL;
}
void miswork::preData(region &r){
	covLoader cov_loader(r.chrname, r.startPos, r.endPos, fai);
	basearr = cov_loader.initBaseArray();
	alnDataLoader data_loader(r.chrname, r.startPos, r.endPos, bamFile);
	data_loader.loadAlnData(alnDataVec);
	cov_loader.generateBaseCoverage(basearr, alnDataVec);
	getreadsMarkregion(r);
	getMultiReads(r);

}

void miswork::clearData(){
	delete[] basearr;
	basearr=NULL;
	if(!alnDataVec.empty()){
		vector<bam1_t*>::iterator aln;
		for(aln=alnDataVec.begin(); aln!=alnDataVec.end(); aln++)
			bam_destroy1(*aln);
		vector<bam1_t*>().swap(alnDataVec);
	}
	vector<region>().swap(abisize);
	vector<region>().swap(abmate);
	vector<region>().swap(abstrand);
	chimeriCoef = 1;
}

double *miswork::getabreadsreatio(region &r){
	double *score = new double[3];
	double mttmp = 0, sttmp = 0, istmp = 0;
	int markMeancov;
	//caculate abmate
	for (size_t i = 0; i < abmate.size(); i++)
	{
		int count = 0, sum = 0;
		//judge whether coverage is too low
		for (int64_t j = abmate.at(i).startPos; j < abmate.at(i).endPos; j++)
		{
			sum = basearr[j - r.startPos].coverage.num_bases[5] + sum;//basearr problem  evaregions.startPos
			count++;
		}
		markMeancov = sum/count;
		if (markMeancov < mincov) break;
		//caculate score
		int rToal = 0, rAb = 0;
		for (size_t k = 0; k < alnDataVec.size(); k++)
		{
			if((abmate.at(i).startPos < (alnDataVec.at(k)->core.pos + alnDataVec.at(k)->core.l_qseq) and ((alnDataVec.at(k)->core.pos + alnDataVec.at(k)->core.l_qseq) < abmate.at(i).endPos)) or
			(alnDataVec.at(k)->core.pos < abmate.at(i).endPos and (abmate.at(i).endPos < (alnDataVec.at(k)->core.pos+alnDataVec.at(k)->core.l_qseq)))){
				rToal++;
				if(alnDataVec[k]->core.tid != alnDataVec[k]->core.mtid)
				rAb++;//judge whether read is abmate
			}
		}
		if(mttmp < double(rAb)/rToal) mttmp = double(rAb)/rToal;
	}
	score[0] = mttmp;

	//caculate abstrand
	for (size_t i = 0; i < abstrand.size(); i++)
	{
		int count = 0,sum = 0;
		//judge whether coverage is too low
		for (int64_t j = abstrand.at(i).startPos; j < abstrand.at(i).endPos; j++)
		{
			sum = basearr[j - r.startPos].coverage.num_bases[5] + sum;
			count++;
		}
		markMeancov = sum/count;
		if (markMeancov < mincov) break;
		//caculate score
		int rToal = 0, rAb = 0;
		for (size_t k = 0; k < alnDataVec.size(); k++)
		{
			if((abstrand.at(i).startPos < (alnDataVec.at(k)->core.pos + alnDataVec.at(k)->core.l_qseq) and ((alnDataVec.at(k)->core.pos + alnDataVec.at(k)->core.l_qseq) < abstrand.at(i).endPos)) or
			(alnDataVec.at(k)->core.pos < abstrand.at(i).endPos and (abstrand.at(i).endPos < (alnDataVec.at(k)->core.pos + alnDataVec.at(k)->core.l_qseq)))){
				rToal++;
				if(isabstrandRead(alnDataVec.at(k)))
				rAb++;//judge whether read is abstrandread
			}
		}
		if(sttmp < double(rAb)/rToal) sttmp = double(rAb)/rToal;
	}
	score[1] = sttmp;

	//caculate abisize
	for (size_t i = 0; i < abisize.size(); i++)
	{
		int count = 0,sum = 0;
		//judge whether coverage is too low
		for (int64_t j = abisize.at(i).startPos; j < abisize.at(i).endPos; j++)
		{
			sum = basearr[j - r.startPos].coverage.num_bases[5] + sum;
			count++;
		}
		markMeancov = sum/count;
		if (markMeancov < mincov) break;
		//caculate score
		int rToal = 0, rAb = 0;
		for (size_t k = 0; k < alnDataVec.size(); k++)
		{
			if((abisize.at(i).startPos < (alnDataVec.at(k)->core.pos + alnDataVec.at(k)->core.l_qseq) and ((alnDataVec.at(k)->core.pos + alnDataVec.at(k)->core.l_qseq) < abisize.at(i).endPos)) or
			(alnDataVec.at(k)->core.pos < abisize.at(i).endPos and (abisize.at(i).endPos < (alnDataVec.at(k)->core.pos + alnDataVec.at(k)->core.l_qseq)))){
				rToal++;
				if(!((abs(alnDataVec[k]->core.isize) > mininsertsize) and (abs(alnDataVec[k]->core.isize) < maxinsertsize)))
				rAb++;//judge whether read is abisize
			}
		}
		if(istmp < double(rAb)/rToal) istmp = double(rAb)/rToal;
	}
	score[2] = istmp;
	return score;
}


void miswork::getreadsMarkregion(region &r){
	for(uint64_t i=0; i<alnDataVec.size(); i++){

		if(isabstrandRead(alnDataVec.at(i))){
			for(int64_t j = max((int64_t)(alnDataVec[i]->core.pos + 1), (int64_t)r.startPos); j <= min((int64_t)(alnDataVec[i]->core.pos + alnDataVec[i]->core.l_qseq + 1), (int64_t)r.endPos); j++){  //0-based to 1-based needed +1
				basearr[j-r.startPos].coverage.num_reads[1]++;
			}
		}
		if(alnDataVec[i]->core.tid==alnDataVec[i]->core.mtid){

			if(!((abs(alnDataVec[i]->core.isize) > mininsertsize) and (abs(alnDataVec[i]->core.isize) < maxinsertsize))){
				for(int64_t j = max((int64_t)(alnDataVec[i]->core.pos + 1), (int64_t)r.startPos); j <= min((int64_t)(alnDataVec[i]->core.pos+alnDataVec[i]->core.l_qseq + 1), (int64_t)r.endPos); j++){
					basearr[j-r.startPos].coverage.num_reads[2]++;
				}
			}
		}else{

			for(int64_t j = max((int64_t)(alnDataVec[i]->core.pos + 1), (int64_t)r.startPos);j <= min((int64_t)(alnDataVec[i]->core.pos + alnDataVec[i]->core.l_qseq + 1), (int64_t)(r.endPos)); j++){
				basearr[j-r.startPos].coverage.num_reads[0]++;
			}
		}
	}

	bool stflag = false, isflag = false, mtflag = false;
	region sttmp, istmp, mttmp;
	for (int64_t i = r.startPos; i <= r.endPos; i++){
		if(((basearr[i - r.startPos].coverage.num_reads[1]/double(basearr[i - r.startPos].coverage.num_bases[5])) > minStrandRatio) and i!=r.endPos){
			if (!stflag)
			{
				sttmp.chrname = r.chrname;
				sttmp.chrlen = r.chrlen;
				sttmp.startPos = i;
				stflag = true;
			}
		}else{
			if (stflag)
			{
				sttmp.endPos = i;
				abstrand.push_back(sttmp);
				stflag = false;
			}
		}

		if((basearr[i - r.startPos].coverage.num_reads[2]/double(basearr[i - r.startPos].coverage.num_bases[5])) > minisizeRatio and i!=r.endPos){
			if (!isflag)
			{
				istmp.chrname = r.chrname;
				istmp.chrlen = r.chrlen;
				istmp.startPos = i;
				isflag = true;
			}
		}else{
			if (isflag)
			{
				istmp.endPos = i;
				abisize.push_back(istmp);
				isflag = false;
			}
		}

		if((basearr[i - r.startPos].coverage.num_reads[0]/double(basearr[i - r.startPos].coverage.num_bases[5])) > minMateRatio and i!=r.endPos){
			if (!mtflag)
			{
				mttmp.chrname = r.chrname;
				mttmp.chrlen = r.chrlen;
				mttmp.startPos = i;
				mtflag = true;
			}
		}else{
			if (mtflag)
			{
				mttmp.endPos = i;
				abmate.push_back(mttmp);
				mtflag = false;
			}
		}
	}
}

void miswork::getMultiReads(region &r){
	int64_t chimericNum = 0;
	for(uint64_t i=0; i<alnDataVec.size(); i++){
		if((alnDataVec.at(i)->core.flag & BAM_FSUPPLEMENTARY)) chimericNum++;
		if((alnDataVec.at(i)->core.flag & BAM_FSECONDARY)){
			for(int64_t j = max((int64_t)(alnDataVec[i]->core.pos + 1), (int64_t)r.startPos); j <= min((int64_t)(alnDataVec[i]->core.pos + alnDataVec[i]->core.l_qseq + 1), (int64_t)r.endPos); j++){
				basearr[j-r.startPos].coverage.num_multiReads++;
			}
		}

	}
	chimeriCoef = double(chimericNum)/alnDataVec.size();
}

bool miswork::isabstrandRead(bam1_t *b){
	if (b->core.tid == b->core.mtid)
	{
		if(b->core.pos < b->core.mpos){
			if (b->core.flag & BAM_FREVERSE)
			{
				return true;
			}
		}else{
			if(!(b->core.flag & BAM_FREVERSE) and (b->core.flag & BAM_FMREVERSE) )
			{
				return true;
			}
		}
	}
	return false;
}
