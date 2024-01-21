#include "Evaluate.h"

Evaluate::Evaluate(Paras *paras){
	this->paras = paras;
	init();
	cout << "############# Parameters #############" << endl;

	getinregion(regionfile);
	getEvaRegions(fainame);
	getranregion();
	loadAlnData();
	getMisworkVec();
}

Evaluate::~Evaluate(){
	clearMisworkVec();
	fai_destroy(fai);
}

void Evaluate::init(){
	string chrname_tmp, cluster_path, feature_path;

	paras->outputParas();

	fai = fai_load(paras->refFile.c_str());
	regionfile = paras->regionFile;
	bamFile = paras->inBamFile;
	fainame = paras->refFile + ".fai";

	num_threads = paras->num_threads;
	outdir = paras->outDir;

	if (access (outdir.c_str (), 0) == -1){
		mkdir(outdir.c_str(), S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
	}

	cluster_path = outdir + "/"+ "cluster";
	feature_path = outdir + "/"+ "feature";
	mkdir(cluster_path.c_str(), S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
	mkdir(feature_path.c_str(), S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);

	final_result_filename = outdir + "/" + "final_result";
	cov_filename = feature_path +"/"+ cov_filename;
	indel_filename = feature_path +"/" + indel_filename ;
	strand_filename = feature_path +"/" + strand_filename;
	insert_filename = feature_path +"/" + insert_filename;
	mate_filename = feature_path +"/" + mate_filename;

	cov_cluster_filename = cluster_path + "/"+ "cluster_cov.csv";
	indel_cluster_filename = cluster_path + "/" + "cluster_indel.csv";
	strand_cluster_filename = cluster_path + "/" + "cluster_strand.csv";
	insert_cluster_filename = cluster_path + "/" + "cluster_insert.csv";
	mate_cluster_filename = cluster_path + "/" + "cluster_mate.csv";
	cluster_stat_filename = cluster_path + "/" + "cluster_stat";
	cluster_detail_filename = cluster_path + "/" + "cluster_detail";

	randomCoef = paras->random_norm_reg_percent;
	minRegsize = paras->minRegsize;
	exRegFold = paras->exRegFold;
	covflag = paras->cov_flag;
	minCovFold = paras->minCovFold;
	maxCovFold = paras->maxCovFold;
	indelflag = paras->indel_flag;
	minLocRatio = paras->minLocRatio;
	abstrandflag = paras->abstrand_flag;
	minStrandRatio = paras->minStrandRatio;
	abisizeflag = paras->abisize_flag;
	minisizeRatio = paras->minisizeRatio;
	IsizeSdevFold = paras ->isizeSdevFold;
	abmateflag = paras->abmate_flag;
	minMateRatio = paras->minmateradio;

}
void Evaluate::getinregion(string &filename){
	ifstream file;
	file.open(filename, ios::in);
	if(!file.is_open()){
		cout << "file open failed" << endl;
	}
	string buf;
	getline(file, buf);
	getline(file, buf);		//20220510 get rid of titleline
	while(buf.size()){
		inregions.push_back(buf);
		getline(file, buf);
	}
	file.close();
}

void Evaluate::splitstring(string &origin, char flag, string &front, string &after){
	int pos = origin.find(flag);
	front = origin.substr(0, pos);
	after = origin.substr(pos+1, origin.size());
}


long int Evaluate::getChrLen(string &chrname){
	for(uint64_t i=0; i < fairegion.size(); i++){
		if(chrname == fairegion.at(i).chrname){
			return fairegion.at(i).endPos;
		}
	}
	return 0;
}


void Evaluate::getEvaRegions(string &fai){
	fstream in;
	in.open(fai, ios::in);
	if(!in.is_open()) cout << fai << " is not opened" << endl;
	string buf;
	region tmp;
	while(getline(in,buf)){
		vector<string> r;
		r = split(buf,"\t");
		tmp.chrname = r[0];
		tmp.startPos = 1;
		tmp.endPos=stoi(r[1]);
		fairegion.push_back(tmp);
	}
	fairegion.shrink_to_fit();

	string nextchr,startpos,endpos,postmp;
	for(uint64_t i=0; i<inregions.size(); i++){
		splitstring(inregions[i], '\t', tmp.chrname, postmp);
		splitstring(postmp, '\t', startpos, endpos);
		tmp.startPos = atol(startpos.c_str());
		tmp.endPos = atol(endpos.c_str());
		tmp.chrlen = getChrLen(tmp.chrname);
		if(tmp.endPos-tmp.startPos<minRegsize){
			tmp.startPos = max(1l,(atol(startpos.c_str()) + atol(endpos.c_str()))/2 - minRegsize/2);
			tmp.endPos = min(tmp.chrlen,(atol(startpos.c_str()) + atol(endpos.c_str()))/2 + (long int)minRegsize/2);
		}
		regions.push_back(tmp);
	}
	
	// 	tmp.startPos = max(one, (chrregions[i].startPos+chrregions[i].endPos)/2 - 2.5*len);
	// 	tmp.endPos = (chrregions[i].startPos+chrregions[i].endPos)/2 + 2.5*len;

	int len;
	for(uint64_t i=0; i<regions.size(); i++){
		len = regions.at(i).endPos - regions.at(i).startPos + 1;
		tmp.chrname = regions.at(i).chrname;
		if((regions.at(i).startPos - exRegFold*len)>1){
			tmp.startPos = (regions.at(i).startPos - exRegFold*len);
		}else{
			tmp.startPos = 1;
		}
		if((regions.at(i).endPos + exRegFold*len) < regions.at(i).chrlen){
			tmp.endPos = (regions.at(i).endPos + exRegFold*len);
		}else{
			tmp.endPos = regions.at(i).chrlen;
		}
		tmp.chrlen = regions.at(i).chrlen;
		evaregions.push_back(tmp);
	}
}

bool Evaluate::compareContigLen(region &a,region &b){
	return a.endPos > b.endPos;
}

bool Evaluate::isRegExist(region &reg, vector <region> &vec){
	for(uint64_t i=0; i<vec.size();i++){
		if (reg.chrname==vec.at(i).chrname and ((vec.at(i).startPos<reg.endPos and reg.endPos<vec.at(i).endPos) or (reg.startPos<vec.at(i).endPos and vec.at(i).endPos<reg.endPos)))
		{
			return true;
		}
	}
	return false;
}

void Evaluate::getranregion(){
	region tmp;
	int len = evaregions.at(0).endPos - evaregions.at(0).startPos + 1, randomFai;
	int randomNum = int(randomCoef*regions.size());
	long one = 1;
	vector <region> v1;
	v1 = fairegion;
	sort(v1.begin(),v1.end(),compareContigLen);
	srand(time(NULL));

	for(int i=0;i<randomNum;){
        
        
        randomFai = (rand()% (int(0.2*v1.size())));
        tmp.chrname=v1.at(randomFai).chrname;
        tmp.startPos = (rand() % (v1.at(randomFai).endPos-20400));
        tmp.endPos=tmp.startPos+len;
        if (isRegExist(tmp,regions)) continue;
        regions.push_back(tmp);
		long evaStart = tmp.startPos-2*len;
		long evaEnd = tmp.endPos+2*len;

		tmp.startPos = max(one, evaStart);
		tmp.endPos = min(evaEnd,v1.at(randomFai).endPos);
		evaregions.push_back(tmp);
        i++;
    }
	
	// for(uint64_t i=0; i<chrregions.size(); i++){
	// 	tmp.chrname = chrregions[i].chrname;
	// 	tmp.startPos = max(one, (chrregions[i].startPos+chrregions[i].endPos)/2 - 2.5*len);
	// 	tmp.endPos = (chrregions[i].startPos+chrregions[i].endPos)/2 + 2.5*len;
	// 	tmp.num = i;
	// 	regions.push_back(tmp);
	// 	tmp.startPos = max(one, (chrregions[i].startPos+chrregions[i].endPos)/2 - 0.5*len);
	// 	tmp.endPos = (chrregions[i].startPos+chrregions[i].endPos)/2 + len/2;
	// 	evaregions.push_back(tmp);

	// }
	regions.shrink_to_fit();
	evaregions.shrink_to_fit();
}

int Evaluate::loadAlnData(){
	double sdev = 0;
	InsertSizeEst isize(bamFile);
	isize.estInsertSize(meaninsertsize, sdev);
	mininsertsize = meaninsertsize-IsizeSdevFold*sdev;
	maxinsertsize = meaninsertsize+IsizeSdevFold*sdev;

	cout << "minisize: " << mininsertsize << "\tmaxisize: " << maxinsertsize << endl;
	return 0;
}

bool Evaluate::isabstrandRead(bam1_t *b){
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

void Evaluate::getreadsMarkregion(region &r){

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


void Evaluate::getMultiReads(region &r){
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

void Evaluate::preData(region &r){
	covLoader cov_loader(r.chrname, r.startPos, r.endPos, fai);
	basearr = cov_loader.initBaseArray();
	alnDataLoader data_loader(r.chrname, r.startPos, r.endPos, paras->inBamFile);
	data_loader.loadAlnData(alnDataVec);
	cov_loader.generateBaseCoverage(basearr, alnDataVec);
	getreadsMarkregion(r);
	getMultiReads(r);

}

void Evaluate::clearData(){
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


double* Evaluate::getabreadsreatio(region &r){
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

//void Evaluate::getindicator(uint64_t s){
//	preData(evaregions.at(s));//
//	indicators tmp;//miswork.tmp
//	meancov=0;//
//	int blankregion=0, sum=0, count=0, nindel, basecov;//
//	double *abReadsRatio;//
//	double covscore=0, indelscore=0, minLocov;//
//	int64_t startPos_tmp = evaregions.at(s).startPos;//
//	for(int64_t i=regions.at(s).startPos; i<=regions.at(s).endPos; i++){ //int64_64
//		if(basearr[i-startPos_tmp].coverage.idx_RefBase<4){
//			sum = basearr[i-startPos_tmp].coverage.num_bases[5] + sum;
//			count++;
//		}
//	}
//	meancov = sum/count;
//	mincov = meancov*minCovFold;
//	maxcov = meancov*maxCovFold;
//	minLocov = meancov*minLocRatio;
////	cout << "regions.at(s).start :" << regions.at(s).startPos << " regions.at(s).end :" << regions.at(s).endPos << endl;
//
//	for(int64_t i = regions.at(s).startPos; i <= regions.at(s).endPos; i++){		//basearrd:Base *baseArray = new Base[endPos-startPos+1]();
//		if(basearr[i-startPos_tmp].coverage.idx_RefBase < 4){				//42217  42418
//			basecov = basearr[i-startPos_tmp].coverage.num_bases[5] + basearr[i - startPos_tmp].del_num_from_del_vec;
//			nindel = basearr[i-startPos_tmp].insVector.size() + basearr[i-startPos_tmp].del_num_from_del_vec;
//			if(basecov > minLocov){
//				if(indelscore < double(nindel)/basecov){
//					indelscore = double(nindel)/basecov;
//				}
//			}
//		}
//	}
////	cout << "evaregions.at(s)start :" << evaregions.at(s).startPos << " evaregions.at(s)end :" << evaregions.at(s).endPos << endl;//调试
//
//	for(int64_t i=evaregions.at(s).startPos; i <= evaregions.at(s).endPos; i++){  //41815 42820
////		cout << basearr[i-regions.at(s).startPos] << endl；
//		if(basearr[i-startPos_tmp].coverage.idx_RefBase < 4){
//			basecov = basearr[i-startPos_tmp].coverage.num_bases[5];
//			nindel = basearr[i-startPos_tmp].insVector.size() + basearr[i-startPos_tmp].del_num_from_del_vec;
//			if(basecov < mincov or (basecov - 0.5*basearr[i-startPos_tmp].coverage.num_multiReads)>maxcov){
//				covscore += 1;
//			}
//		}else{
//			blankregion++;
//		}
//	}
//	abReadsRatio = getabreadsreatio(evaregions.at(s));
//	if(blankregion < (evaregions.at(s).endPos - evaregions.at(s).startPos+1)){
//		tmp.covscore[0] = chimeriCoef*(covscore)/(evaregions.at(s).endPos - evaregions.at(s).startPos-blankregion+1);
//	}else{
//		tmp.covscore[0] = 0;
//	}
//	tmp.indelscore[0] = indelscore;
//	tmp.matescore = abReadsRatio[0];
//	tmp.strandscore = abReadsRatio[1];
//	tmp.insertscore = abReadsRatio[2];
//	blankregion = 0;
//	scores.push_back(tmp);
//	delete abReadsRatio;
//	clearData();
//
//
//}

//void Evaluate::extracfeature(){
//	int tmp = 0;
//	Time time;
//
//	cout << "############# Start region analysis #############" << endl;
//	cout << "Total regions: " << evaregions.size() << endl;
//	cout << "evaregions.size():" << evaregions.size() << " region.size:" << regions.size() << endl;
//
//	for(uint64_t i=0; i<evaregions.size(); i++){
//		if(i<49000){
//			continue;
//		}
////		cout << i << endl;//test
//		miswork* thread = new miswork (regions.at(i),evaregions.at(i),fai,bamFile);
//		scores.push_back(thread->getindicator());
//		delete thread;
////		getindicator(i);
//		if(100*i/evaregions.size()-tmp >= 1 ){
//			cout << "[" << time.getTime() << "]: " <<"Processed regions: " << i << '/' << evaregions.size();
//			cout << " (" << setprecision(4) << 100*float(i)/evaregions.size() << '%' << ')' << endl;
//			tmp = 100*i/evaregions.size();
//		}
//	}
//	cout << "Processed regions: " << evaregions.size() << '/' << evaregions.size();
//	cout << " (" << setprecision(4) << 100*float(evaregions.size())/evaregions.size() << '%' << ')' << endl;
//	outputfile();
//}

void Evaluate::extracfeature(){
	Time time;
	
	misWork_opt *miswork_opt;
	miswork *mis_work;
	size_t num_work, num_work_percent;	//

	hts_tpool *p = hts_tpool_init(num_threads);//num_threads cin
	hts_tpool_process *q = hts_tpool_process_init(p, num_threads*2, 1);

	num_work = evaregions.size();//to store the numbers of total works;
	num_work_percent = num_work / paras->num_parts_progress;//to store the percentage of work of each process
	if(num_work_percent==0) num_work_percent = 1;

	misworkDone_num = 0;

	pthread_mutex_init(&mtx_misworkDone_num, NULL);//

	cout << "############# Start region analysis #############" << endl;
	cout << "Total regions: " << evaregions.size() << endl;
	cout << "evaregions.size():" << evaregions.size() << " region.size:" << regions.size() << endl;
	
	for(uint64_t i=0; i<evaregions.size(); i++){
//		if(i<49000){
//			continue;
//		}
		mis_work = miswork_vec.at(i);

		miswork_opt = new misWork_opt();
		miswork_opt->mis_work = mis_work;
		miswork_opt->work_id = i;
		miswork_opt->num_work = num_work;
		miswork_opt->num_work_percent = num_work_percent;
		miswork_opt->p_misworkDone_num = &(misworkDone_num);
		miswork_opt->p_mtx_misworkDone_num = &(mtx_misworkDone_num);

		hts_tpool_dispatch(p, q, processSingleMisWork, miswork_opt);
//		hts_tpool_dispatch(p, q, getindicator, call_work_opt);
//		if(100*i/evaregions.size()-tmp >= 1 ){
//			cout << "[" << time.getTime() << "]: " <<"Processed regions: " << i << '/' << evaregions.size();
//			cout << " (" << setprecision(4) << 100*float(i)/evaregions.size() << '%' << ')' << endl;
//			tmp = 100*i/evaregions.size();
//		}
	}
//	cout << "Processed regions: " << evaregions.size() << '/' << evaregions.size();
//	cout << " (" << setprecision(4) << 100*float(evaregions.size())/evaregions.size() << '%' << ')' << endl;

    hts_tpool_process_flush(q);
    hts_tpool_process_destroy(q);
    hts_tpool_destroy(p);

	outputfile();  // test and then will restore
//    clearMisworkVec();
}

void Evaluate::getMisworkVec(){
	miswork* tmp;
	for(size_t i = 0; i < evaregions.size(); i++){
		tmp = new miswork(regions.at(i), evaregions.at(i), fai, bamFile);
		tmp->maxinsertsize = maxinsertsize;		//maxsize chushihua
		tmp->mininsertsize = mininsertsize;
		miswork_vec.push_back(tmp);
	}
}

void Evaluate::clearMisworkVec(){
	for(size_t i = 0;i < miswork_vec.size();i++){
		delete miswork_vec.at(i);
	}
	vector<miswork*>().swap(miswork_vec);
}

void Evaluate::outputfile(){
	string foldername;
	foldername = outdir;

//	if(covflag){
//		string cov = "./" + foldername + "/cov.csv";
//		ofstream covi(cov.c_str());
//		if(!covi.is_open()) cout << "outputfile failed create" << endl;
//		for(uint64_t i=0; i<miswork_vec.size(); i++){  //score[]
//			covi << setprecision(2) << miswork_vec[i]->scores.covscore[0] << endl;
//		}
//		covi.close();
//	}
	if(covflag){
		ofstream covi(cov_filename.c_str());
		if(!covi.is_open()) cout << "outputfile failed create" << endl;
		for(uint64_t i=0; i<miswork_vec.size(); i++){  //score[]
			covi << setprecision(2) << miswork_vec[i]->scores.covscore[0] << endl;
		}
		covi.close();
	}
	if(indelflag){
		ofstream indeli(indel_filename.c_str());
		if(!indeli.is_open()) cout << "outputfile failed create" << endl;
		for(uint64_t i=0; i<miswork_vec.size(); i++){
			indeli << setprecision(2) << miswork_vec[i]->scores.indelscore[0] << endl;
		}
		indeli.close();
	}
	
	if(abstrandflag){
		ofstream strandi(strand_filename.c_str());
		if(!strandi.is_open()) cout << "outputfile failed create" << endl;
		for(uint64_t i=0; i<miswork_vec.size(); i++){
			strandi << setprecision(2) << miswork_vec[i]->scores.strandscore << endl;
		}
		strandi.close();
	}
	
	if(abisizeflag){
		ofstream inserti(insert_filename.c_str());
		if(!inserti.is_open()) cout << "outputfile failed create" << endl;
		for(uint64_t i=0; i<miswork_vec.size(); i++){
			inserti << setprecision(2) << miswork_vec[i]->scores.insertscore << endl;
		}
		inserti.close();
	}
	
	if (abmateflag)
	{
		ofstream matei(mate_filename.c_str());
		if(!matei.is_open()) cout << "outputfile failed create" << endl;
		for(uint64_t i=0; i<miswork_vec.size(); i++){
			matei << setprecision(2) << miswork_vec[i]->scores.matescore << endl;
		}
		matei.close();
	}
}




//void Evaluate::clusterfile(){
//	vector<string> files;
//	string directory = "./";
//	string resultDirectory = "./result";
//
//	string foldername;
//	foldername = outdir + resultDirectory;
//	if (access (foldername.c_str (), 0) != -1){   //fold is or not exist
//		mkdir(foldername.c_str(), S_IRWXU);
//	}
//
//	for (const auto& entry : fs::directory_iterator(directory)) {
//		 if (entry.is_regular_file() && entry.path().extension() == ".csv") {
//					files.push_back(entry.path().string());
//		}
//	}
//
//	 for (const auto& entry :fs::directory_iterator(currentPath)) {
//	        std::string filePath = entry.path().string();
//	 }
//	for (const auto& file : files) {
//		string cmd = "mlpack_kmeans -c 2 -i " + file + " -o ./result/" + file ;
//		cout << cmd << endl;
////		string cmd = "mlpack_kmeans -c 2 -i " + file + " -o ./result/";
//		system(cmd.c_str());
//	}
//
//}
void Evaluate::clusterfile(){
//	vector<string> files;
//	string directory = "./";
//	string resultDirectory = "./result";

    cout << "clustering the feature: " << endl;
	string cmd;
	if(covflag){
		cmd = "mlpack_kmeans -c 2 -i " + cov_filename + " -o " + cov_cluster_filename ;
//		cout << cmd << endl;
		system(cmd.c_str());
	}

	if(indelflag){
		cmd = "mlpack_kmeans -c 2 -i " + indel_filename + " -o " + indel_cluster_filename ;
//		cout << cmd << endl;
		system(cmd.c_str());
	}

	if(abstrandflag){
		cmd = "mlpack_kmeans -c 2 -i " + strand_filename + " -o " + strand_cluster_filename ;
//		cout << cmd << endl;
		system(cmd.c_str());
	}

	if(abisizeflag){
		cmd = "mlpack_kmeans -c 2 -i " + insert_filename + " -o " + insert_cluster_filename ;
//		cout << cmd << endl;
		system(cmd.c_str());
	}
	if (abmateflag){
		cmd = "mlpack_kmeans -c 2 -i " + mate_filename + " -o " + mate_cluster_filename ;
//		cout << cmd << endl;
		system(cmd.c_str());
	}


}

vector<vector<double>> Evaluate::formatfile(string& filename) {
    ifstream infile(filename);
    vector<vector<double>> fvec;
    bool flag = false;
    string line;

  //error
    vector<string> str_vec;
    vector<double> double_vec;
    string value;



    if (infile) {   //reline = infile.readline()
//        cout << line <<endl;

        while (getline(infile, line)) {
        	str_vec = split(line, ",");
        	for(size_t i = 0; i < str_vec.size(); i++){
        		double_vec.push_back(stod(str_vec.at(i)));
        	}

            fvec.push_back(double_vec);
            double_vec.clear();
        }
        infile.close();

//        for (size_t i = 0; i < fvec.size(); i++) {
//            for (size_t j = 0; j < fvec[i].size(); j++) {
//                cout << fvec[i][j] << " ";
//            }
//            cout << endl;
//        }

        for (const auto& j : fvec) {
            if ((j[0] > fvec[0][0] && j[1] < fvec[0][1]) || (j[0] < fvec[0][0] && j[1] > fvec[0][1])) {
                flag = true;
                break;
            }
        }

        if (flag) {
            for (auto& k : fvec) {
                k[1] = 1 - k[1];
            }
        }
    }
    else {
        cerr << "Error opening file: " << filename << endl;
    }
    return fvec;
}

void Evaluate::analysfile(){
//	vector<string> files;
	string line;
	string directory = "./";
	string result_path = outdir + "/"+ "result";

	vector<double> last ;
    vector<string> region;//save region data
    vector<string> regsplit;
    vector<string> tone;
    vector<vector<string>> regVec;

//	cout << "analysfile is start" << endl;

    ifstream reg(regionfile);
    if (!reg) {
        cerr << "Error opening file: " << regionfile << endl;
    }


    while (getline(reg, line)) {
    	region.emplace_back(line);
    }

    //delete the first element
    if (!region.empty()) {
    	region.erase(region.begin());
    }
    reg.close();

	ofstream clusterResult(cluster_stat_filename);
	if (!clusterResult) {
		cerr << "Error opening file: " << cluster_stat_filename << endl;
	}

	ofstream regionClassification(cluster_detail_filename);
	if (!regionClassification) {
		cerr << "Error opening file: " << final_result_filename << endl;
	}

	// write to the two files
//	cout << 1 << endl;
    string titleline = "#scaffold\tstartPos\tendPos";
	vector<vector<double>> covCluster,indelCluster,abstrandCluster,abisizeCluster,abmateCluster;  //restore
	if(covflag){
		covCluster = formatfile(cov_cluster_filename);
		titleline += "\tcov";
	}


	if(indelflag){
        indelCluster = formatfile(indel_cluster_filename);
		titleline += "\tindel";
	}

	if(abstrandflag){
        abstrandCluster = formatfile(strand_cluster_filename);
		titleline += "\tabstrand";
	}

	if(abisizeflag){
		abisizeCluster = formatfile(insert_cluster_filename);
		titleline += "\tabisize";
	}

	if (abmateflag){
		abmateCluster = formatfile(mate_cluster_filename);
		titleline += "\tabisize";
	}



	//title+mark

	if(covflag){
		titleline += "\tcovMark";
	}

	if(indelflag){
		titleline += "\tindelMark";
	}

	if(abstrandflag){
		titleline += "\tabstrandMark";
	}

	if(abisizeflag){
		titleline += "\tabisizeMark";
	}

	if (abmateflag){
		titleline += "\tabisizeMark";
	}
	titleline += "\tregionMark";
	regionClassification << titleline << endl ;

//	vector<int> last;   //final result
	double last_1 = 0;
	for(size_t i=0; i<region.size(); i++){
		string linei = region[i];
		if(covflag){
			last_1 += covCluster[i][1];
			linei = linei + '\t' + to_string(covCluster[i][0]);
		}

		if(indelflag){
			last_1 += indelCluster[i][1];
			linei = linei + '\t' + to_string(indelCluster[i][0]);
		}

		if(abstrandflag){
			last_1 += abstrandCluster[i][1];
			linei = linei + '\t' + to_string(abstrandCluster[i][0]);
		}

		if(abisizeflag){
			last_1 += abisizeCluster[i][1];
			linei = linei + '\t' + to_string(abisizeCluster[i][0]);
		}

		if (abmateflag){
			last_1 += abmateCluster[i][1];
			linei = linei + '\t' + to_string(abmateCluster[i][0]);
		}



		if(covflag){
			if(covCluster[i][1]!= 0)
				linei = linei + '\t' + 'P';
			else
				linei = linei + '\t' + 'N';
		}

		if(indelflag){
			if(indelCluster[i][1]!= 0)
				linei = linei + '\t' + 'P';
			else
				linei = linei + '\t' + 'N';
		}

		if(abstrandflag){
			if(abstrandCluster[i][1]!= 0)
				linei = linei + '\t' + 'P';
			else
				linei = linei + '\t' + 'N';
		}

		if(abisizeflag){
			if(abisizeCluster[i][1]!= 0)
				linei = linei + '\t' + 'P';
			else
				linei = linei + '\t' + 'N';
		}

		if (abmateflag){
			if(abmateCluster[i][1]!= 0)
				linei = linei + '\t' + 'P';
			else
				linei = linei + '\t' + 'N';
		}
		if (last_1 != 0) {
		    linei = linei + '\t' + 'P';
		} else {
		    linei = linei + '\t' + 'N';
		}
		regionClassification << linei << endl ;
		last.push_back(last_1);
		last_1 = 0;
	}
	regionClassification.close();

	for (size_t i = 0; i < last.size(); i++) {
	        std::string outline = std::to_string(i + 1) + ':' + std::to_string(last[i]) + '\n';
	        clusterResult << outline;
	    }

	clusterResult.close();

//	cout << "11111111111111111111111111111111111111111111" << endl;//


//        while (getline(infile, line)) {
//        	str_vec = split(line, ",");
//        	for(size_t i = 0; i < str_vec.size(); i++){
//        		double_vec.push_back(stod(str_vec.at(i)));
//        	}
//
//            fvec.push_back(double_vec);
//            double_vec.clear();
//        }
//        infile.close();

	ifstream fsplit(cluster_detail_filename);
    if (!fsplit.is_open()) {
        cerr << "Error opening file: " << final_result_filename << endl;
    }
    string regline;
    while (getline(fsplit, regline)) {
//		cout <<"line: "<< regline << endl;
    	tone = split(regline,"\t");
    	regVec.push_back(tone);
    	tone.clear();
    }


//    for (size_t i = 1; i < regsplit.size(); i++) {
//        line = regsplit[i];//
//        cout <<"line: "<< line << endl;
//        tone = split(line,"\t");
//        regVec.push_back(tone);
//
//    }

//    vector<int> regPostive;
//    vector<int> regNegative;

	ofstream final_result(final_result_filename);//positive
	final_result << "#scaffold\tstartPos\tendPos\tmark" << endl;
    for (size_t i = 1; i< regVec.size(); i++) {
    	final_result << regVec[i][0] << " " << regVec[i][1] << " " << regVec[i][2] << " ";
        if (regVec[i].back() == "P") {//return regVec[i]'s last element and jugde it
//            cout << "regPostive: " << i << endl;
            final_result << "\tmisassembly" << endl;
        } else {
//            cout << "regNegative: " << i << endl;
            final_result << "\tnormal" << endl;
        }
    }
    cout << "misClus work is completeded." << endl;
//
//	if (final_result.is_open()) {
//		final_result << "#scaffold\tstartPos\tendPos\tmark" << endl;
//		cout << "misassembly: " << endl;
//		for (size_t i= 0; i< regPostive.size(); i++){
//			final_result << region.at(regPostive.at(i)-1) << "\tmisassembly" << endl;
//		}
//		cout << "normal: " << endl;
//		for (size_t i= 0; i< regNegative.size(); i++){
//			final_result << region.at(regNegative.at(i)- 1) << "\tnormal" << endl;
//		}
//		final_result.close();
//	} else {
//		cout << "Failed to open final_result." << endl;
//	}



//	if(covflag){
//		 system("rm cov.csv");
//	}
//
//	if(indelflag){
//		system("rm indel.csv");
//	}
//
//	if(abstrandflag){
//		system("rm strand.csv");
//	}
//
//	if(abisizeflag){
//		system("rm insert.csv");
//	}
//
//	if (abmateflag){
//		system("rm mate.csv");
//	}
}









