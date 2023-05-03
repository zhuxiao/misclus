#include "Evaluate.h"

Evaluate::Evaluate(char** argv,faidx_t *fai,Paras *paras){
	this->paras = paras;
	this->fai = fai;
	regionfile = argv[3];
	bamFile = paras->inBamFile;
	string fainame = paras->refFile + ".fai";
	cout << "############# Parameters #############" << endl;
	if(!(strcmp(argv[4],"200")==0)){
		minRegsize = atoi(argv[4]);
	}
	cout << "min-regsize: " << minRegsize << endl;
	if(!(strcmp(argv[5],"2")==0)){
		exRegFold = atof(argv[5]);
	}
	cout << "extend-regsize-fold: " << exRegFold << endl;
	if((strcmp(argv[6],"off")==0)){
		covflag = false;
	}
	cout << "cov: " << covflag << endl;
	if(!(strcmp(argv[7],"0.5")==0) or !(strcmp(argv[8],"2")==0)){
		minCovFold = atof(argv[7]);
		maxCovFold = atof(argv[8]);
	}
	cout << "min-cov-fold: " << minCovFold << endl;
	cout << "max-cov-fold: " << maxCovFold << endl;
	if((strcmp(argv[9],"off")==0)){
		indelflag = false;
	}
	cout << "indels: " << indelflag << endl;
	if(!(strcmp(argv[10],"0.2")==0)){
		minLocRatio = atof(argv[10]);
	}	
	cout << "min-locs-ratio: " << minLocRatio << endl;
	if((strcmp(argv[11],"off")==0)){
		abstrandflag = false;
	}	
	cout << "abstrand: " << abstrandflag << endl;

	if(!(strcmp(argv[12],"0.3")==0)){
		minStrandRatio = atof(argv[12]);
	}
	cout << "min-abstrand-ratio: " << minStrandRatio << endl;
	if((strcmp(argv[13],"off")==0)){
		abisizeflag = false;
	}	
	cout << "abisize: " << abisizeflag << endl;
	if(!(strcmp(argv[14],"0.3")==0)){
		minisizeRatio = atof(argv[14]);
	}		
	cout << "min-abisize-ratio: " << minisizeRatio << endl;

	if(!(strcmp(argv[15],"3")==0)){
		IsizeSdevFold = atof(argv[15]);
	}	
	cout << "isize-sdev-fold: " << IsizeSdevFold << endl;

	if((strcmp(argv[16],"off")==0)){
		abmateflag = false;
	}		
	cout << "abmate: " << abmateflag << endl;
	
	if(!(strcmp(argv[17],"0.3")==0)){
		minMateRatio = atof(argv[17]);
	}	
	cout << "min-abmate-ratio: " << minMateRatio << endl;

	if(!(strcmp(argv[18],"misEval_out")==0)){
		outdir = argv[18];
	}

	if(!(strcmp(argv[19],"0.2")==0)){
		randomCoef = atof(argv[19]);
	}	
	cout << "NumNormalRegCoef: " << randomCoef << endl;	

	getinregion(regionfile);
	getchrregions(fainame);
	getranregion();
	loadAlnData();
}

Evaluate::~Evaluate(){
	if(!baseArr.empty()){
		vector<Base*>().swap(baseArr);
	}
	if(!alnDataVector.empty()){
		for(uint64_t i=0; i<alnDataVector.size(); i++){
			for(uint64_t j=0; j<alnDataVector[i].size(); j++){
				bam_destroy1(alnDataVector[i][j]);
			}
		}
		vector<vector<bam1_t*>>().swap(alnDataVector);
	}

}

void Evaluate::getinregion(string filename){
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

void Evaluate::splitstring(string origin, char flag, string &front, string &after){
	int pos = origin.find(flag);
	front = origin.substr(0, pos);
	after = origin.substr(pos+1, origin.size());
}

void Evaluate::getchrregions(string fai){
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
		tmp.startPos=atoi(startpos.c_str());
		tmp.endPos=atoi(endpos.c_str());
		if(tmp.endPos-tmp.startPos<minRegsize){
			tmp.startPos = (atoi(startpos.c_str()) + atoi(endpos.c_str()))/2 - minRegsize/2;
			tmp.endPos = (atoi(startpos.c_str()) + atoi(endpos.c_str()))/2 + minRegsize/2;
		}
		evaregions.push_back(tmp);
	}
	int count = 0;
	uint64_t j = 0;
	for(uint64_t i=0; i<evaregions.size(); i++){
		evaregions.at(i).num = count;
		for(; j<fairegion.size(); j++){
			if(evaregions.at(i).chrname==fairegion.at(j).chrname){
				if(i==0){
					chrregions.push_back(fairegion.at(j));
				}
				if(i>0 and evaregions.at(i).chrname!=evaregions.at(i-1).chrname){
					chrregions.push_back(fairegion.at(j));
					count++;
					evaregions.at(i).num = count;
				}
				break;
			}
		}
	}
	chrregions.shrink_to_fit();
	int len;
	for(uint64_t i=0; i<evaregions.size(); i++){
		len = evaregions.at(i).endPos - evaregions.at(i).startPos + 1;
		tmp.chrname = evaregions.at(i).chrname;
		if((evaregions.at(i).startPos - exRegFold*len)>1){
			tmp.startPos = (evaregions.at(i).startPos - exRegFold*len);
		}else{
			tmp.startPos = 1;
		}
		if((evaregions.at(i).startPos + exRegFold*len) < chrregions.at(evaregions.at(i).num).endPos){
			tmp.endPos = (evaregions.at(i).startPos + exRegFold*len);
		}else{
			tmp.endPos = chrregions.at(evaregions.at(i).num).endPos;
		}
		tmp.num = evaregions.at(i).num;
		regions.push_back(tmp);
	}
}

bool Evaluate::compareContigLen(region a,region b){
	return a.endPos > b.endPos;
}

bool Evaluate::isRegExist(region reg, vector <region> vec){
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

vector<Base*> Evaluate::initbaseinfo(){
	Base *chrbase;
	for(uint64_t i=0; i<chrregions.size(); i++){
		covLoader cov_loader(chrregions[i].chrname, chrregions[i].startPos,  chrregions[i].endPos, fai);
		chrbase = cov_loader.initBaseArray();
		baseArr.push_back(chrbase);
	}
	return baseArr;
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

void Evaluate::getreadsMarkregion(region r){

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
				sttmp.num = r.num;
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
				istmp.num = r.num;
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
				mttmp.num = r.num;
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

void Evaluate::getMultiReads(region r){
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

void Evaluate::preData(region r){
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
	vector<bam1_t *> ().swap(alnDataVec);
	vector<region>().swap(abisize);
	vector<region>().swap(abmate);
	vector<region>().swap(abstrand);
	chimeriCoef = 1;
}


double* Evaluate::getabreadsreatio(region r){
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
			sum = basearr[j - r.startPos].coverage.num_bases[5] + sum;
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

void Evaluate::getindicator(int s){
	preData(regions.at(s));
	indicators tmp;
	meancov=0;
	int blankregion=0, sum=0, count=0, nindel, basecov;
	double *abReadsRatio;
	double covscore=0, indelscore=0, minLocov;
	for(int64_t i=regions.at(s).startPos; i<=regions.at(s).endPos; i++){
		if(basearr[i-regions.at(s).startPos].coverage.idx_RefBase<4){
			sum = basearr[i-regions.at(s).startPos].coverage.num_bases[5] + sum;
			count++;
		}
	}
	meancov = sum/count;
	mincov = meancov*minCovFold;
	maxcov = meancov*maxCovFold;
	minLocov = meancov*minLocRatio;
	for(int64_t i = regions.at(s).startPos; i <= regions.at(s).endPos; i++){
		if(basearr[i-regions.at(s).startPos].coverage.idx_RefBase < 4){
			basecov = basearr[i-regions.at(s).startPos].coverage.num_bases[5] + basearr[i - regions.at(s).startPos].del_num_from_del_vec;
			nindel = basearr[i-regions.at(s).startPos].insVector.size() + basearr[i-regions.at(s).startPos].del_num_from_del_vec;
			if(basecov > minLocov){
				if(indelscore < double(nindel)/basecov){
					indelscore = double(nindel)/basecov;
				}
			}
		}
	}

	for(int64_t i=evaregions.at(s).startPos; i <= evaregions.at(s).endPos; i++){
		if(basearr[i-regions.at(s).startPos].coverage.idx_RefBase < 4){
			basecov = basearr[i-regions.at(s).startPos].coverage.num_bases[5];
			nindel = basearr[i-regions.at(s).startPos].insVector.size() + basearr[i-regions.at(s).startPos].del_num_from_del_vec;
			if(basecov < mincov or (basecov - 0.5*basearr[i-regions.at(s).startPos].coverage.num_multiReads)>maxcov){
				covscore += 1;
			}
		}else{
			blankregion++;
		}
	}
	abReadsRatio = getabreadsreatio(regions.at(s));
	if(blankregion < (evaregions.at(s).endPos - evaregions.at(s).startPos+1)){
		tmp.covscore[0] = chimeriCoef*(covscore)/(evaregions.at(s).endPos - evaregions.at(s).startPos-blankregion+1);
	}else{
		tmp.covscore[0] = 0;
	}
	tmp.indelscore[0] = indelscore;
	tmp.matescore = abReadsRatio[0];
	tmp.strandscore = abReadsRatio[1];
	tmp.insertscore = abReadsRatio[2];
	blankregion = 0;
	scores.push_back(tmp);
	clearData();

}

void Evaluate::extracfeature(){
	int tmp = 0;
	cout << endl << "############# Start region analysis #############" << endl;
	cout << "Total regions: " << evaregions.size() << endl;
	for(uint64_t i=0; i<evaregions.size(); i++){
		getindicator(i);
		if(100*i/evaregions.size()-tmp >= 1 ){
			cout << "Processed regions: " << i << '/' << evaregions.size();
			cout << " (" << setprecision(4) << 100*float(i)/evaregions.size() << '%' << ')' << endl;
			tmp = 100*i/evaregions.size();
		}
	}
	cout << "Processed regions: " << evaregions.size() << '/' << evaregions.size();
	cout << " (" << setprecision(4) << 100*float(evaregions.size())/evaregions.size() << '%' << ')' << endl;
	outputfile();
}

void Evaluate::outputfile(){
	string foldername;
	foldername = outdir;
	if (access (foldername.c_str (), 0) == -1){
		mkdir(foldername.c_str(), S_IRWXU);
	}

	if(covflag){
		string cov = "./" + foldername + "/cov.csv";
		ofstream covi(cov.c_str());
		if(!covi.is_open()) cout << "outputfile failed create" << endl;
		for(uint64_t i=0; i<scores.size(); i++){
			covi << setprecision(2) << scores[i].covscore[0] << endl;
		}
		covi.close();
	}

	if(indelflag){
		string indel = "./" + foldername + "/indel.csv";
		ofstream indeli(indel.c_str());
		if(!indeli.is_open()) cout << "outputfile failed create" << endl;
		for(uint64_t i=0; i<scores.size(); i++){
			indeli << setprecision(2) << scores[i].indelscore[0] << endl;
		}
		indeli.close();
	}
	
	if(abstrandflag){
		string strand = "./" + foldername + "/strand.csv";
		ofstream strandi(strand.c_str());
		if(!strandi.is_open()) cout << "outputfile failed create" << endl;
		for(uint64_t i=0; i<scores.size(); i++){
			strandi << setprecision(2) << scores[i].strandscore << endl;
		}
		strandi.close();
	}
	
	if(abisizeflag){
		string insert = "./" + foldername + "/insert.csv";
		ofstream inserti(insert.c_str());
		if(!inserti.is_open()) cout << "outputfile failed create" << endl;
		for(uint64_t i=0; i<scores.size(); i++){
			inserti << setprecision(2) << scores[i].insertscore << endl;
		}
		inserti.close();
	}
	
	if (abmateflag)
	{
		string mate = "./" + foldername + "/mate.csv";
		ofstream matei(mate.c_str());
		if(!matei.is_open()) cout << "outputfile failed create" << endl;
		for(uint64_t i=0; i<scores.size(); i++){
			matei << setprecision(2) << scores[i].matescore << endl;
		}
		matei.close();
	}
}
