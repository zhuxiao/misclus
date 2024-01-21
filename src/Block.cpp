#include "Block.h"

#include <sys/stat.h>
#include <sys/types.h>
#include <algorithm>

#include "covLoader.h"
#include "util.h"
using namespace std;

//pthread_mutex_t mutex_print = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t mutex_write_misAln = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t mutex_assem_work = PTHREAD_MUTEX_INITIALIZER;

// Constructor with parameters
Block::Block(string chrname, size_t startPos, size_t endPos, faidx_t *fai,  Paras *paras){
	string chrname_tmp;
	this->paras =  paras;
	this->chrname = chrname;
	this->chrlen = faidx_seq_len(fai, chrname.c_str()); // get the reference length
	this->startPos = startPos;
	this->endPos = endPos;
	this->fai = fai;
	chrname_tmp = preprocessPipeChar(chrname);
	workdir = chrname_tmp;
	outCovFile = chrname_tmp + "_" + to_string(startPos) + "-" + to_string(endPos) + ".bed";
	baseArr = NULL;
	winSize = 0;
//	winSize = paras->slideSize * 3;
	headIgnFlag = false;
	tailIgnFlag = false;
	process_flag = true;

	meanCov = 0;

	// detect output files
	snvDetectPrefix_local = "snv_cand";
	indelDetectPrefix_local = "indel_cand";
	clipRegDetectPrefix_local = "clipReg_cand";
	snvFilenameDetect = snvDetectPrefix_local + "_" + chrname_tmp + "_" + to_string(startPos) + "-" + to_string(endPos) + ".bed";
	indelFilenameDetect = indelDetectPrefix_local + "_" + chrname_tmp + "_" + to_string(startPos) + "-" + to_string(endPos) + ".bed";
	clipRegFilenameDetect = clipRegDetectPrefix_local + "_" + chrname_tmp + "_" + to_string(startPos) + "-" + to_string(endPos) + ".bed";

	// assemble output files
	indelAssemblePrefix_local = "contig";
	indelFilenameAssemble = indelAssemblePrefix_local + "_" + chrname_tmp + "_" + to_string(startPos) + "-" + to_string(endPos) + ".fa";

	var_cand_indel_file = NULL;
	misAln_reg_file = NULL;
	var_cand_clipReg_file = NULL;
}

// Destructor
Block::~Block(){
	if(baseArr) destroyBaseArray();
	if(!alnDataVector.empty()) destoryAlnData();
	if(!snvVector.empty()) destroySnvVector();
	if(!indelVector.empty()) destroyIndelVector();
	if(!clipRegVector.empty()) destroyClipRegVector();
}

// set the output directory
void Block::setOutputDir(string& out_dir_detect_prefix, string& out_dir_assemble_prefix, string& out_dir_call_prefix){
	out_dir_detect = out_dir_detect_prefix;
	out_dir_assemble = out_dir_assemble_prefix;
	out_dir_call = out_dir_call_prefix + "/" + chrname + "_" + to_string(startPos) + "-" + to_string(endPos);

	out_dir_call = preprocessPipeChar(out_dir_call);
}

// set limit regions
void Block::setLimitRegs(vector<simpleReg_t*> &sub_limit_reg_vec){
	for(size_t i=0; i<sub_limit_reg_vec.size(); i++) this->sub_limit_reg_vec.push_back(sub_limit_reg_vec.at(i));
	this->sub_limit_reg_vec.shrink_to_fit();
}

// set process flag
void Block::setProcessFlag(bool process_flag){
	this->process_flag = process_flag;
}

// initialize the base array of the block
Base *Block::initBaseArray(){
	covLoader cov_loader(chrname, startPos, endPos, fai);
	baseArr = cov_loader.initBaseArray();
	return baseArr;
}

// destroy the base array of the block
void Block::destroyBaseArray(){
	delete[] baseArr;
	baseArr = NULL;
}

// destroy the alignment data of the block
void Block::destoryAlnData(){
	vector<bam1_t*>::iterator aln;
	for(aln=alnDataVector.begin(); aln!=alnDataVector.end(); aln++)
		bam_destroy1(*aln);
	vector<bam1_t*>().swap(alnDataVector);
}

// destroy the SNV vector
void Block::destroySnvVector(){
	vector<size_t>().swap(snvVector);
}

// destroy the indel vector
void Block::destroyIndelVector(){
	vector<reg_t*>::iterator it;
	for(it=indelVector.begin(); it!=indelVector.end(); it++)
		delete (*it);
	vector<reg_t*>().swap(indelVector);
}

// destroy clip region vector
void Block::destroyClipRegVector(){
	vector<reg_t*>().swap(clipRegVector);
}

// set the region ingnore flag in the block
void Block::setRegIngFlag(bool headIgnFlag, bool tailIgnFlag){
	this->headIgnFlag = headIgnFlag;
	this->tailIgnFlag = tailIgnFlag;
}

// prepare the alignment data and fill the estimation data
void Block::blockFillDataEst(size_t op_est){
	// initialize the alignment data
	initBaseArray();
	loadAlnData();

	// compute the block base information, including coverage, insertions, deletions and clippings
	computeBlockBaseInfo();

	// fill the data for estimation
	fillDataEst(op_est);

	// compute mean read length
	if(op_est==SIZE_EST_OP){
		AddReadLenEstInfo();
	}
}

// fill the estimation data
void Block::fillDataEst(size_t op_est){
//	int64_t i;
//	size_t j, len;
//	Base *base;
	//vector<insEvent_t*>::iterator ins;
	//vector<delEvent_t*>::iterator del;
	//vector<clipEvent_t*>::iterator clip;

//	for(i=startPos; i<=endPos; i++){
//		base = baseArr + i - startPos;
//		if(base->coverage.idx_RefBase!=4){ // excluding 'Ns' gap region
//
//			if(op_est==SIZE_EST_OP){
//				// insertion
//				for(j=0; j<base->insVector.size(); j++){ // insSizeEstArr
//					len = base->insVector.at(j)->seq.size();
//					if(len>AUX_ARR_SIZE) len = AUX_ARR_SIZE;
//					paras->insSizeEstArr[len] ++;
//				}
//				// deletion
//				for(j=0; j<base->delVector.size(); j++){ // delSizeEstArr
//					len = base->delVector.at(j)->seq.size();
//					if(len>AUX_ARR_SIZE) len = AUX_ARR_SIZE;
//					paras->delSizeEstArr[len] ++;
//				}
//				// clipping
////				for(j=0; j<base->clipVector.size(); j++){ // clipSizeEstArr
////					len = stoi(base->clipVector.at(j)->seq);
////					if(len>AUX_ARR_SIZE){
////						len = AUX_ARR_SIZE;
////					}
////					paras->clipSizeEstArr[len] ++;
////				}
//			}else if(op_est==NUM_EST_OP){
//				// insertion
//				len = base->insVector.size();
//				if(len>AUX_ARR_SIZE) len = AUX_ARR_SIZE;
//				paras->insNumEstArr[len] ++; // insNumEstArr
//				// deletion
//				len = base->delVector.size();  // delNumEstArr
//				if(len>AUX_ARR_SIZE) len = AUX_ARR_SIZE;
//				paras->delNumEstArr[len] ++;
//				// clipping
//				len = base->clipVector.size();  // clipNumEstArr
//				if(len>AUX_ARR_SIZE){
//					//cout << "clipping, num=" << len << endl;
//					len = AUX_ARR_SIZE;
//				}
//				paras->clipNumEstArr[len] ++;
//			}else if(op_est==SNV_EST_OP){
//
//			}else{
//				cerr << __func__ << ", line=" << __LINE__ << ", invalid estimation op_flag: " << op_est << endl;
//				exit(1);
//			}
//		}
//	}
}

// add read length estimation information
void Block::AddReadLenEstInfo(){
//	bam1_t *b;
//	size_t total, num;
////
//	if(!alnDataVector.empty()){
//		total = num = 0;
//		for(size_t i=0; i<alnDataVector.size(); i++){
//			b = alnDataVector.at(i);
//			total += b->core.l_qseq;
//			num ++;
//		}
//		paras->mean_read_len += total;
//		paras->total_read_num_est += num;
//	}
}

// block process
void Block::blockDetect(){
//	pthread_mutex_lock(&mutex_print);
//	cout << chrname << ":" << startPos << "-" << endPos << endl;
//	pthread_mutex_unlock(&mutex_print);

	// initialize the alignment data
	initBaseArray();
	loadAlnData();

	// compute the block base information,
	// including coverage, insertions, deletions and clippings
	computeBlockBaseInfo();

	// mask misAln regions
	//if(paras->maskMisAlnRegFlag) maskMisAlnRegs();

	// compute abnormal signatures
	//computeAbSigs();

	//printRegVec(indelVector, "indelVector");

	removeFalseIndel();

	//printRegVec(indelVector, "indelVector");

	// remove false SNV
	removeFalseSNV();

	// sort the indels and clip regions
	sortRegVec(indelVector);
	sortRegVec(clipRegVector);

	// merge overlapped indels and clip regions
	mergeOverlappedReg(indelVector);
	mergeOverlappedReg(clipRegVector);

	// remove redundant items
	removeRedundantItems(indelVector);
	removeRedundantItems(clipRegVector);

	// save SV to file
	saveSV2File();

	// output base coverage
	//outputCovFile();

	// release memory
	if(baseArr) destroyBaseArray();
	if(!alnDataVector.empty()) destoryAlnData();
}

// load alignment data with specified region in the format like `chr2:100-200'
int Block::loadAlnData(){
	alnDataLoader data_loader(chrname, startPos, endPos, paras->inBamFile);
	data_loader.loadAlnData(alnDataVector);
	return 0;
}

// compute block base information
int Block::computeBlockBaseInfo(){

//	// generate the base coverage array
//	covLoader cov_loader(chrname, startPos, endPos, fai, paras->min_ins_size_filt, paras->min_del_size_filt);
//	cov_loader.generateBaseCoverage(baseArr, alnDataVector);
//
//	// update the coverage information for each base
//	computeBlockMeanCov();
	return 0;
}

int Block::outputCovFile(){
	string out_file_str = workdir + "/" + outCovFile;
	uint16_t *num_baseArr;
	ofstream ofile(out_file_str);

	mkdir(workdir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

	if(!ofile.is_open()){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file " << out_file_str << endl;
		exit(1);
	}

	int i, j, len;
	len = endPos - startPos + 1;
	for(i=0; i<len; i++){
		ofile << chrname << "\t" << startPos << "\t" << endPos;
		num_baseArr = baseArr[i].coverage.num_bases;
		for(j=0; j<5; j++)
			ofile << "\t" << num_baseArr[j];
		ofile << endl;
	}
	ofile.close();
	return 0;
}

// update the coverage information for each base in this block
void Block::computeBlockMeanCov(){
	int64_t pos, totalReadBeseNum = 0, totalRefBaseNum = 0;
	for(pos=startPos; pos<=endPos; pos++){
		// compute the meanCov excluding the gap regions
		if(baseArr[pos-startPos].coverage.idx_RefBase!=4){ // excluding 'N'
			totalReadBeseNum += baseArr[pos-startPos].coverage.num_bases[5];
			totalRefBaseNum ++;
		}
	}
	if(totalRefBaseNum) meanCov = (double)totalReadBeseNum/totalRefBaseNum;
	else meanCov = 0;
}

// extract abnormal signatures
//int Block::computeAbSigs(){
//	int64_t startPosWin, endPosWin, pos, tmp_endPos;
//
//	// process the head region
//	if(!headIgnFlag){
//		startPosWin = startPos;
//		endPosWin = startPosWin + 2 * paras->slideSize - 1;
//		if(endPosWin>endPos) endPosWin = endPos;
//		processSingleRegion(startPosWin, endPosWin, HEAD_REGION);
//	}
//
//	// process the inner regions
//	if(endPos>=2*(int64_t)paras->slideSize){
//		tmp_endPos = endPos - 2 * paras->slideSize;
//		for(pos=startPos; pos<=tmp_endPos; pos+=paras->slideSize){
//			startPosWin = pos;
//			endPosWin = startPosWin + winSize - 1;
//			if(endPosWin>endPos) endPosWin = endPos;
//			processSingleRegion(startPosWin, endPosWin, INNER_REGION);
//		}
//	}else pos = startPos;
//
//	// process the tail region
//	if(!tailIgnFlag){
//		startPosWin = pos;
//		endPosWin = startPosWin + winSize - 1;
//		if(endPosWin>endPos) endPosWin = endPos;
//		else if(endPosWin<endPos){
//			cerr << __func__ << "(), line=" << __LINE__ << ": endPosWin=" << endPosWin << ", endPos=" << endPos << ", error!" << endl; exit(1);
//		}
//		processSingleRegion(startPosWin, endPosWin, TAIL_REGION);
//	}
//
//	// merge adjacent zero coverage regions
//	mergeAdjacentReg(zeroCovRegVector, MIN_ADJACENT_REG_DIST);
//
//	// updateZeroCovReg
//	updateZeroCovRegUsingIndelReg(zeroCovRegVector, indelVector);
//
//	// update indel regions
//	updateSVRegUsingLongZeroCov();
//
//	return 0;
//}

// process single region
//int Block::processSingleRegion(int64_t startRpos, int64_t endRPos, int64_t regFlag){
//
//	vector<simpleReg_t*> sub_limit_reg_vec_tmp;
//	bool process_flag;
//
////	if(startRpos>4534233)
////		cout << "\t" << chrname << ":" << startRpos << "-" << endRPos << endl;
//
//	process_flag = true;
//	if(paras->limit_reg_process_flag){
//		if(sub_limit_reg_vec.size()){
//			//sub_limit_reg_vec_tmp = getSimpleRegs(chrname, startRpos, endRPos, sub_limit_reg_vec);
//			sub_limit_reg_vec_tmp = getOverlappedSimpleRegsExt(chrname, startRpos, endRPos, sub_limit_reg_vec, ASSEM_SLIDE_SIZE);
//			if(sub_limit_reg_vec_tmp.size()==0) process_flag = false;
//		}else process_flag = false;
//	}
//
//	if(process_flag){
//		// construct a region
//		Region tmp_reg(chrname, startRpos, endRPos, chrlen, startPos, endPos, baseArr+startRpos-startPos, regFlag, paras);
//		tmp_reg.setMeanBlockCov(meanCov);  // set the mean coverage
//
//		// compute the abnormal signatures in a region
//		if(tmp_reg.wholeRefGapFlag==false and isMisAlnReg(tmp_reg)==false){
//			tmp_reg.computeRegAbSigs();
//
//			// detect clip regions
//			tmp_reg.detectHighClipReg();
//
//			// detect indel candidate regions and SNV candidates
//			tmp_reg.detectIndelReg();
//			tmp_reg.detectSNV();
//
//			// copy the indels and SNVs
//			copySVEvents(tmp_reg);
//
//			// output the nums
//			//tmp_reg.printRegAbSigs();
//
//			//cout << "NO: " << tmp_reg.startMidPartPos << "-" << tmp_reg.endMidPartPos << endl;
//		}
//	}
//
//	return 0;
//}

// add SV events
//void Block::copySVEvents(Region &reg){
//	size_t i;
//	vector<reg_t*> regVec;
//	reg_t *reg_tmp;
//	vector<int64_t> posVec;
//	int64_t pos_tmp;
//
//	// copy clip region
//	regVec = reg.getClipRegVector();
//	for(i=0; i<regVec.size(); i++){
//		reg_tmp = regVec.at(i);
//		clipRegVector.push_back(reg_tmp);
//	}
//
//	// copy indel
//	regVec = reg.getIndelVector();
//	for(i=0; i<regVec.size(); i++){
//		reg_tmp = regVec.at(i);
//		indelVector.push_back(reg_tmp);
//	}
//
//	// copy SNV
//	posVec = reg.getSnvVector();
//	for(i=0; i<posVec.size(); i++){
//		pos_tmp = posVec.at(i);
//		snvVector.push_back(pos_tmp);
//	}
//
//	// compute zero coverage region
//	computeZeroCovReg(reg);
//}
//
//void Block::computeZeroCovReg(Region &reg){
//	size_t i, j;
//	reg_t *reg_tmp;
//	vector<int64_t> posVec;
//	int64_t pos1, pos2;
//	int64_t start_vec_idx, end_vec_idx;
//
//	posVec = reg.getZeroCovPosVector();
//	i = 0;
//	while(i<posVec.size()){
//		start_vec_idx = i;
//		end_vec_idx = -1;
//		for(j=i+1; j<posVec.size(); j++){
//			pos1 = posVec.at(j-1);
//			pos2 = posVec.at(j);
//			if(pos1+1==pos2){
//				end_vec_idx = j;
//			}else break;
//		}
//
//		if(start_vec_idx!=-1 and end_vec_idx!=-1){
//			reg_tmp = new reg_t();
//			reg_tmp->chrname = chrname;
//			reg_tmp->startRefPos = posVec.at(start_vec_idx);
//			reg_tmp->endRefPos = posVec.at(end_vec_idx);
//			reg_tmp->startLocalRefPos = reg_tmp->endLocalRefPos = 0;
//			reg_tmp->startQueryPos = reg_tmp->endQueryPos = 0;
//			reg_tmp->var_type = VAR_UNC;
//			reg_tmp->sv_len = 0;
//			reg_tmp->query_id = -1;
//			reg_tmp->blat_aln_id = -1;
//			reg_tmp->call_success_status = false;
//			reg_tmp->short_sv_flag = false;
//			reg_tmp->zero_cov_flag = false;
//
//			zeroCovRegVector.push_back(reg_tmp);
//			i = end_vec_idx + 1;
//		}else i ++;
//	}
//}

void Block::updateZeroCovRegUsingIndelReg(vector<reg_t*> &zeroCovRegVec, vector<reg_t*> &indelVec){
	size_t i, j;
	int32_t minValue, maxValue;
	reg_t *zero_cov_reg, *reg_tmp;
	int32_t regIdx;
	vector<int32_t> regIdx_vec;

	for(i=0; i<indelVec.size(); i++){
		reg_tmp = indelVec.at(i);
		regIdx = getOverlappedRegIdx(reg_tmp, zeroCovRegVec);
		regIdx_vec.push_back(regIdx);
	}

	// compute startRefPos and endRefPos
	for(i=0; i<zeroCovRegVec.size(); i++){
		zero_cov_reg = zeroCovRegVec.at(i);
		// compute startRefPos and endRefPos
		minValue = zero_cov_reg->startRefPos;
		maxValue = zero_cov_reg->endRefPos;
		for(j=0; j<indelVec.size(); j++){
			regIdx = regIdx_vec.at(j);
			if(regIdx==(int32_t)i){
				reg_tmp = indelVec.at(j);
				if(minValue>reg_tmp->startRefPos) minValue = reg_tmp->startRefPos;
				if(maxValue<reg_tmp->endRefPos) maxValue = reg_tmp->endRefPos;
			}
		}
		zero_cov_reg->startRefPos = minValue;
		zero_cov_reg->endRefPos = maxValue;
	}

	mergeOverlappedReg(zeroCovRegVec); // merge overlapped region
}

void Block::updateSVRegUsingLongZeroCov(){

//	cout << ">>>>>>>>>>> Before update clipRegVector: <<<<<<<<<<<<<<" << endl;
//	printRegVec(clipRegVector, "clipRegVector");  // print clipRegVector
//	printRegVec(zeroCovRegVector, "zeroCovRegVector");  // print zeroCovRegVector

	// update clipping regions
	updateClipRegUsingLongZeroCov(clipRegVector, zeroCovRegVector);

//	cout << ">>>>>>>>>>> After update clipRegVector: <<<<<<<<<<<<<<" << endl;
//	printRegVec(clipRegVector, "clipRegVector");  // print clipRegVector
//	printRegVec(zeroCovRegVector, "zeroCovRegVector");  // print zeroCovRegVector

//	cout << "<<<<<<<<<<<<<< Before update indelVector: >>>>>>>>>>>" << endl;
//	printRegVec(indelVector, "indelVector");  // print indelVector
//	printRegVec(zeroCovRegVector, "zeroCovRegVector");  // print zeroCovRegVector

	// update indel regions
	updateIndelRegUsingLongZeroCov(indelVector, zeroCovRegVector);

//	cout << "<<<<<<<<<<<<<< After update indelVector: >>>>>>>>>>>" << endl;
//	printRegVec(indelVector, "indelVector");  // print indelVector
//	printRegVec(zeroCovRegVector, "zeroCovRegVector");  // print zeroCovRegVector
}

void Block::updateClipRegUsingLongZeroCov(vector<reg_t*> &regVec, vector<reg_t*> &zero_cov_vec){
	size_t i;
	reg_t *reg_tmp;
	int32_t regIdx;

	for(i=0; i<regVec.size(); ){
		reg_tmp = regVec.at(i);
		regIdx = getOverlappedRegIdx(reg_tmp, zero_cov_vec);
		if(regIdx!=-1){
			delete reg_tmp;
			regVec.erase(regVec.begin()+i);
		}else i++;
	}
}

void Block::updateIndelRegUsingLongZeroCov(vector<reg_t*> &regVec, vector<reg_t*> &zero_cov_vec){
	size_t i;
	reg_t *zero_cov_reg, *reg_tmp;
	int32_t regIdx;

	for(i=0; i<regVec.size(); ){
		reg_tmp = regVec.at(i);
		regIdx = getOverlappedRegIdx(reg_tmp, zero_cov_vec);
		if(regIdx!=-1){ // remove item
			delete reg_tmp;
			regVec.erase(regVec.begin()+i);
		}else i++;
	}
	for(i=0; i<zero_cov_vec.size(); i++){
		zero_cov_reg = zero_cov_vec.at(i);
		zero_cov_reg->zero_cov_flag = true;
		regVec.push_back(zero_cov_reg);
	}
	zero_cov_vec.clear();
}

// remove false indels
void Block::removeFalseIndel(){
	bool flag = false;
	reg_t *reg, *foundReg;
	for(size_t i=0; i<indelVector.size(); ){
		foundReg = NULL;
		reg = indelVector.at(i);
		if(reg->zero_cov_flag==false)
			//foundReg = findVarvecItemExtSize(reg->startRefPos, reg->endRefPos, clipRegVector, CLIP_END_EXTEND_SIZE, CLIP_END_EXTEND_SIZE);
			foundReg = findVarvecItem(reg->startRefPos, reg->endRefPos, clipRegVector);
		if(foundReg) { flag = true; delete reg; indelVector.erase(indelVector.begin()+i); }
		else i++;
	}
	if(flag) indelVector.shrink_to_fit();
}

// remove false SNV
void Block::removeFalseSNV(){
	bool flag = false;
	vector<size_t>::iterator it;
	for(it=snvVector.begin(); it!=snvVector.end(); )
		if(isInReg(*it, indelVector)) { flag = true; it = snvVector.erase(it); }
		else it++;
	if(flag) snvVector.shrink_to_fit();
}

// sort the region vector
void Block::sortRegVec(vector<reg_t*> &regVector){
	size_t min_idx;
	reg_t *tmp;

	if(regVector.size()>=2){
		for(size_t i=0; i<regVector.size(); i++){
			min_idx = i;
			for(size_t j=i+1; j<regVector.size(); j++)
				if(regVector[min_idx]->startRefPos>regVector[j]->startRefPos) min_idx = j;
			if(min_idx!=i){
				tmp = regVector[i];
				regVector[i] = regVector[min_idx];
				regVector[min_idx] = tmp;
			}
		}
	}
}

// save SV to file
void Block::saveSV2File(){
	ofstream out_file;
	string out_file_str;
	size_t i;

	// creat the directory and open file for detect command
	mkdir(out_dir_detect.c_str(), S_IRWXU | S_IROTH);

	// save indel
	out_file_str = out_dir_detect + "/" + indelFilenameDetect;
	out_file.open(out_file_str);
	if(!out_file.is_open()){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file" << endl;
		exit(1);
	}
	reg_t* reg;
	for(i=0; i<indelVector.size(); i++){
		reg = indelVector.at(i);
		out_file << reg->chrname << "\t" << reg->startRefPos << "\t" <<  reg->endRefPos << endl;
	}
	out_file.close();

	// save snv
	out_file_str = out_dir_detect + "/" + snvFilenameDetect;
	out_file.open(out_file_str);
	if(!out_file.is_open()){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file" << endl;
		exit(1);
	}
	for(i=0; i<snvVector.size(); i++)
		out_file << chrname << "\t" << snvVector.at(i) << endl;
	out_file.close();

	// save clip region
	out_file_str = out_dir_detect + "/" + clipRegFilenameDetect;
	out_file.open(out_file_str);
	if(!out_file.is_open()){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file" << endl;
		exit(1);
	}
	for(i=0; i<clipRegVector.size(); i++){
		reg = clipRegVector.at(i);
		out_file << reg->chrname << "\t" << reg->startRefPos << "\t" <<  reg->endRefPos << endl;
	}
	out_file.close();
}

// set assembly information file
void Block::setVarCandFiles(ofstream *var_cand_indel_file, ofstream *var_cand_clipReg_file){
	this->var_cand_indel_file = var_cand_indel_file;
	this->var_cand_clipReg_file = var_cand_clipReg_file;

}
// reset assembly information file
void Block::resetVarCandFiles(){
	var_cand_indel_file = NULL;
	var_cand_clipReg_file = NULL;
}

// set misAln region file
void Block::setMisAlnRegFile(ofstream *misAln_reg_file){
	this->misAln_reg_file = misAln_reg_file;
}
// reset misAln region file
void Block::resetMisAlnRegFile(){
	misAln_reg_file = NULL;
}
void Block::blockAnalyse(){
	// initialize the alignment data
	initBaseArray();
	loadAlnData();

	// compute the block base information,
	// including coverage, insertions, deletions and clippings
	computeBlockBaseInfo();

}
void Block::checkpos(int64_t start,  int64_t end){
	int64_t i, j;
	for(i=start; i<=end; i++){
		for(j=0; j<6; j++){
			cout << baseArr[i-startPos].coverage.num_bases[j] << " ";
		}
		cout << baseArr[i-startPos].insVector.size() << " " << baseArr[i-startPos].del_num_from_del_vec << " " << baseArr[i-startPos].clipVector.size() << endl;
	}
}
