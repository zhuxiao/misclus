#include "misCluster.h"

misCluster::misCluster(Paras *paras){
	this->paras = paras;
	init();
	cout << "############# Parameters #############" << endl;

	getinregion(regionfile);
	getEvaRegions(fainame);
	getRandRegion();
	loadAlnData();
	getMisworkVec();
}

misCluster::~misCluster(){
	clearMisworkVec();
	fai_destroy(fai);
}

void misCluster::init(){
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
//	minCovFold = paras->minCovFold;
//	maxCovFold = paras->maxCovFold;
	indelflag = paras->indel_flag;
//	minLocRatio = paras->minLocRatio;
	abstrandflag = paras->abstrand_flag;
//	minStrandRatio = paras->minStrandRatio;
	abisizeflag = paras->abisize_flag;
//	minisizeRatio = paras->minisizeRatio;
	IsizeSdevFold = paras ->isizeSdevFold;
	abmateflag = paras->abmate_flag;
//	minMateRatio = paras->minmateRadio;

}
void misCluster::getinregion(string &filename){
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

void misCluster::splitstring(string &origin, char flag, string &front, string &after){
	int32_t pos = origin.find(flag);
	front = origin.substr(0, pos);
	after = origin.substr(pos+1, origin.size());
}

long int misCluster::getChrLen(string &chrname){
	for(size_t i=0; i<fairegion.size(); i++){
		if(chrname == fairegion.at(i).chrname){
			return fairegion.at(i).endPos;
		}
	}
	return 0;
}


void misCluster::getEvaRegions(string &fai){
	fstream in;
	string buf;
	region tmp;
	string nextchr, startpos, endpos, postmp;
	vector<string> r;
	int32_t len;

	in.open(fai, ios::in);
	if(!in.is_open()) cout << fai << " is not opened" << endl;

	while(getline(in,buf)){
		r = split(buf, "\t");
		tmp.chrname = r[0];
		tmp.startPos = 1;
		tmp.endPos = stoi(r[1]);
		fairegion.push_back(tmp);
	}
	fairegion.shrink_to_fit();

	for(size_t i=0; i<inregions.size(); i++){
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
	
	for(size_t i=0; i<regions.size(); i++){
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

bool misCluster::compareContigLen(region &a,region &b){
	return a.endPos > b.endPos;
}

bool misCluster::isRegExist(region &reg, vector <region> &vec){
	for(size_t i=0; i<vec.size(); i++){
		if (reg.chrname==vec.at(i).chrname and ((vec.at(i).startPos<reg.endPos and reg.endPos<vec.at(i).endPos) or (reg.startPos<vec.at(i).endPos and vec.at(i).endPos<reg.endPos))){
			return true;
		}
	}
	return false;
}

void misCluster::getRandRegion(){
	region tmp;
	int32_t len, randomFai, randomNum, end_skipsize;
	int64_t one, evaStart, evaEnd;
	vector <region> v1;
//	int32_t contig_len;
	len = regions.at(0).endPos - regions.at(0).startPos + 1;
	randomNum = int32_t(randomCoef*regions.size());
	one = 1;

	end_skipsize = 2000;
	v1 = fairegion;
	sort(v1.begin(), v1.end(), compareContigLen);
	//srand(time(NULL));

	for(int32_t i=0; i<randomNum; ){
		randomFai = (rand() % (int32_t(0.2*v1.size())));
		tmp.chrname = v1.at(randomFai).chrname;
		if(v1.at(randomFai).endPos <= 3*end_skipsize) break;
		tmp.startPos = (rand() % (v1.at(randomFai).endPos- 2*end_skipsize));//20400? or V1.at(V1.size() * percent).endPos
//		if(v1.at(randomFai).endPos-20400<0)
//			cout << v1.at(randomFai).endPos-20400 << endl;//contig_length < num
		tmp.startPos += end_skipsize;//startPos and endPos caculate this random num to stat away from paired ends of the scaffolds
		tmp.endPos = tmp.startPos + len;
		if(isRegExist(tmp, evaregions)) continue;
		regions.push_back(tmp);
		evaStart = tmp.startPos - 2*len;
		evaEnd = tmp.endPos + 2*len;

		tmp.startPos = max(one, evaStart);
		tmp.endPos = min(evaEnd, v1.at(randomFai).endPos);
		evaregions.push_back(tmp);
		i++;
	}
	
	regions.shrink_to_fit();
	evaregions.shrink_to_fit();
}

int32_t misCluster::loadAlnData(){
	double sdev = 0;
	InsertSizeEst isize(bamFile);

	sdev = 0;
	isize.estInsertSize(meaninsertsize, sdev);
	mininsertsize = meaninsertsize - IsizeSdevFold * sdev;
	maxinsertsize = meaninsertsize + IsizeSdevFold * sdev;

	cout << "minimal insert size: " << mininsertsize << ", maximal insert size: " << maxinsertsize << endl;
	return 0;
}

void misCluster::extracfeature(){
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
//		if(i<45000){
//			continue;
//		}
//		cout << i << endl;
		mis_work = miswork_vec.at(i);

		miswork_opt = new misWork_opt();
		miswork_opt->mis_work = mis_work;
		miswork_opt->work_id = i;
		miswork_opt->num_work = num_work;
		miswork_opt->num_work_percent = num_work_percent;
		miswork_opt->p_misworkDone_num = &(misworkDone_num);
		miswork_opt->p_mtx_misworkDone_num = &(mtx_misworkDone_num);

		hts_tpool_dispatch(p, q, processSingleMisWork, miswork_opt);
	}

	hts_tpool_process_flush(q);
	hts_tpool_process_destroy(q);
	hts_tpool_destroy(p);

	outputfile();  // test and then will restore
}

void misCluster::getMisworkVec(){
	miswork* tmp;
	for(size_t i = 0; i < evaregions.size(); i++){
		tmp = new miswork(regions.at(i), evaregions.at(i), fai, paras, mininsertsize, maxinsertsize);
		tmp->maxinsertsize = maxinsertsize;		//maxsize chushihua
		tmp->mininsertsize = mininsertsize;
		miswork_vec.push_back(tmp);
	}
}

void misCluster::clearMisworkVec(){
	for(size_t i = 0;i < miswork_vec.size();i++)
		delete miswork_vec.at(i);
	vector<miswork*>().swap(miswork_vec);
}

void misCluster::outputfile(){
	string foldername;
	ofstream out_file;
	size_t i;

	foldername = outdir;
	if(covflag){
		out_file.open(cov_filename);
		if(!out_file.is_open()) {
			cerr << "cannot create file '" << cov_filename << "', error!" << endl;
			exit(1);
		}
		for(i=0; i<miswork_vec.size(); i++)  //score[]
			out_file << setprecision(2) << miswork_vec[i]->scores.covscore[0] << endl;
		out_file.close();
	}
	if(indelflag){
		out_file.open(indel_filename);
		if(!out_file.is_open()) {
			cerr << "cannot create file '" << indel_filename << "', error!" << endl;
			exit(1);
		}
		for(i=0; i<miswork_vec.size(); i++)
			out_file << setprecision(2) << miswork_vec[i]->scores.indelscore[0] << endl;
		out_file.close();
	}
	
	if(abstrandflag){
		out_file.open(strand_filename);
		if(!out_file.is_open()) {
			cerr << "cannot create file '" << strand_filename << "', error!" << endl;
			exit(1);
		}
		for(i=0; i<miswork_vec.size(); i++)
			out_file << setprecision(2) << miswork_vec[i]->scores.strandscore << endl;
		out_file.close();
	}
	
	if(abisizeflag){
		out_file.open(insert_filename);
		if(!out_file.is_open()) {
			cerr << "cannot create file '" << insert_filename << "', error!" << endl;
			exit(1);
		}
		for(i=0; i<miswork_vec.size(); i++)
			out_file << setprecision(2) << miswork_vec[i]->scores.insertscore << endl;
		out_file.close();
	}
	
	if (abmateflag)	{
		out_file.open(mate_filename);
		if(!out_file.is_open()) {
			cerr << "cannot create file '" << mate_filename << "', error!" << endl;
			exit(1);
		}
		for(i=0; i<miswork_vec.size(); i++)
			out_file << setprecision(2) << miswork_vec[i]->scores.matescore << endl;
		out_file.close();
	}
}

void misCluster::clusterfile(){
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

vector<vector<double>> misCluster::formatfile(string& filename) {
	ifstream infile(filename);
	vector<vector<double>> fvec;
	bool flag;
	string line, value;
	vector<string> str_vec;
	vector<double> double_vec;

	flag = false;
	if (infile.is_open()) {
		while (getline(infile, line)) {
			str_vec = split(line, ",");
			for(size_t i = 0; i < str_vec.size(); i++)
				double_vec.push_back(stod(str_vec.at(i)));
			fvec.push_back(double_vec);
			double_vec.clear();
		}
		infile.close();

		for (const auto& j : fvec) {
			if ((j[0] > fvec[0][0] && j[1] < fvec[0][1]) || (j[0] < fvec[0][0] && j[1] > fvec[0][1])) {
				flag = true;
				break;
			}
		}

		if (flag) {
			for (auto& k : fvec)
				k[1] = 1 - k[1];
		}
	}
	else {
		cerr << "Error opening file: " << filename << endl;
		exit(1);
	}

	return fvec;
}

void misCluster::analysfile(){
	string line, linei, titleline, outline;
	vector<double> last;
	vector<string> region, regsplit, tone;
	vector<vector<string>> regVec;
	vector<vector<double>> covCluster, indelCluster, abstrandCluster, abisizeCluster, abmateCluster;  //restore
	double last_1;
	size_t i;

	ifstream reg(regionfile);
	if (!reg) {
		cerr << "Error opening file: " << regionfile << endl;
		exit(1);
	}

	while (getline(reg, line))
		region.emplace_back(line);

	//delete the first element
	if (!region.empty())
		region.erase(region.begin());
	reg.close();

	ofstream clusterResult(cluster_stat_filename);
	if (!clusterResult) {
		cerr << "Error opening file: " << cluster_stat_filename << endl;
		exit(1);
	}
	ofstream regionClusterDetail(cluster_detail_filename);
	if (!regionClusterDetail) {
		cerr << "Error opening file: " << cluster_detail_filename << endl;
		exit(1);
	}

	// write to the two files
	titleline = "#scaffold\tstartPos\tendPos";
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
		titleline += "\tabmate";
	}

	if(covflag) titleline += "\tcovMark";
	if(indelflag) titleline += "\tindelMark";
	if(abstrandflag) titleline += "\tabstrandMark";
	if(abisizeflag) titleline += "\tabisizeMark";
	if (abmateflag) titleline += "\tabmateMark";

	titleline += "\tregionMark";
	regionClusterDetail << titleline << endl ;

	last_1 = 0;
	for(i=0; i<region.size(); i++){
		linei = region.at(i);
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
			if(covCluster[i][1]!= 0) linei = linei + '\t' + 'P';
			else linei = linei + '\t' + 'N';
		}

		if(indelflag){
			if(indelCluster[i][1]!= 0) linei = linei + '\t' + 'P';
			else linei = linei + '\t' + 'N';
		}

		if(abstrandflag){
			if(abstrandCluster[i][1]!= 0) linei = linei + '\t' + 'P';
			else linei = linei + '\t' + 'N';
		}

		if(abisizeflag){
			if(abisizeCluster[i][1]!= 0) linei = linei + '\t' + 'P';
			else linei = linei + '\t' + 'N';
		}

		if (abmateflag){
			if(abmateCluster[i][1]!= 0) linei = linei + '\t' + 'P';
			else linei = linei + '\t' + 'N';
		}
		if (last_1 != 0) linei = linei + '\t' + 'P';
		else linei = linei + '\t' + 'N';

		regionClusterDetail << linei << endl ;
		last.push_back(last_1);
		last_1 = 0;
	}
	regionClusterDetail.close();

	for (i = 0; i<last.size(); i++) {
		outline = to_string(i + 1) + ':' + to_string(last.at(i)) + '\n';
		clusterResult << outline;
	}
	clusterResult.close();

	ifstream fsplit(cluster_detail_filename);
	if (!fsplit.is_open()) {
		cerr << "Error opening file: " << final_result_filename << endl;
		exit(1);
	}

	while (getline(fsplit, line)) {
		tone = split(line, "\t");
		regVec.push_back(tone);
		tone.clear();
	}
	fsplit.close();

	ofstream final_result(final_result_filename); //positive
	final_result << "#scaffold\tstartPos\tendPos\tmark" << endl;
	for (i = 1; i< regVec.size(); i++) {
		final_result << regVec[i][0] << " " << regVec[i][1] << " " << regVec[i][2] << " ";
		if (regVec[i].back() == "P") {//return regVec[i]'s last element and jugde it
//			cout << "regPostive: " << i << endl;
			final_result << "\tmisassembly" << endl;
		} else {
//			cout << "regNegative: " << i << endl;
			final_result << "\tnormal" << endl;
		}
	}
}
