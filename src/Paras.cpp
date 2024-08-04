#include "Paras.h"

#include <string.h>
#include <vector>
#include <htslib/hts.h>
#include <cmath>

#include "util.h"

// Constructor with parameters
Paras::Paras(){
	init();

}

// Constructor with parameters
Paras::Paras(int argc, char **argv){
	init();
	if(parseParas(argc, argv)!=0) exit(1);

	// check bam file
	if(checkBamFile()!=0) exit(1);
}

//Destructor
Paras::~Paras(){
	if(!limit_reg_vec.empty())
		destroyLimitRegVector(limit_reg_vec);
}

// initialization
void Paras::init(){
	command = "";
	refFile = "";
	inBamFile = "";
	outFilePrefix = "";
	regionFile = "";
	outDir = OUT_DIR;
	num_threads = 0;

	minRegsize = 0;
	exRegFold = 0;
	minCovFold = 0;
	maxCovFold = 0;
	minLocRatio = 0;
	minStrandRatio = 0;
	minisizeRatio = 0;
	isizeSdevFold = 0;
	random_norm_reg_percent = 0;

	cov_flag = abstrand_flag = abisize_flag = abmate_flag = IDC_flag = true;
	num_parts_progress = NUM_PARTS_PROGRESS;
}

// get the Canu program version
//string Paras::getCanuVersion(){
//	string canu_version_str, canu_version_cmd, canu_version_filename, line, canu_prog_name;
//	ifstream infile;
//	vector<string> str_vec;
//
//	canu_version_filename = "canu_version";
//	canu_version_cmd = "canu -version > " + canu_version_filename;
//	system(canu_version_cmd.c_str());
//
//	infile.open(canu_version_filename);
//	if(!infile.is_open()){
//		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file " << canu_version_filename << endl;
//		exit(1);
//	}
//
//	canu_version_str = "";
//	getline(infile, line);
//	if(line.size()){
//		str_vec = split(line, " ");
//		canu_prog_name = str_vec.at(0);
//
//		if(canu_prog_name.compare("Canu")==0) canu_version_str = str_vec.at(1);
//		else{
//			cout << "Canu may be not correctly installed, please check its installation before running this program." << endl;
//			exit(1);
//		}
//	}
//
//	infile.close();
//
//	return canu_version_str;
//}

void Paras::show_version(){
	cout << PROG_VERSION << endl;
}


// check Bam file, and generate the BAM index if it is unavailable
int Paras::checkBamFile(){
	samFile *in = 0;
	string idx_filename;
	int ret = 1;

	if ((in = sam_open(inBamFile.c_str(), "r")) == 0) {
		cerr << __func__ << ": failed to open " << inBamFile.c_str() << " for reading" << endl;
		exit(1);
	}

	hts_idx_t *idx = sam_index_load(in, inBamFile.c_str()); // load index
	if (idx == NULL) { // index is unavailable, then generate it
		if(in) sam_close(in);

		cout << __func__ << ": BAM index is unavailable, now generate it, please wait ...\n" << endl;

		// construct the index file name
		if(inBamFile.substr(inBamFile.size()-4).compare(".bam")==0)
			idx_filename = inBamFile.substr(0, inBamFile.size()-4) + ".bai";
		else
			idx_filename = inBamFile + ".bai";

		ret = sam_index_build3(inBamFile.c_str(), idx_filename.c_str(), 0, num_threads);
		switch (ret) {
			case 0:
				break;
			case -2:
				cerr << __func__ << ", failed to open " << inBamFile << endl;
				break;
			case -3:
				cerr << __func__ << inBamFile << " is in a format that cannot be usefully indexed" << endl;
				break;
			case -4:
				cerr << __func__ << ", failed to create or write index" << endl;
				break;
			default:
				cerr << __func__ << ", failed to create index for " << inBamFile << endl;
				break;
		}
	}else{
		ret = 0;
		hts_idx_destroy(idx); // destroy the BAM index
		if(in) sam_close(in);
	}

	return ret;
}

// parse the parameters
int Paras::parseParas(int argc, char **argv){
	if (argc < 2) { showUsage(); return 1; }

	if (strcmp(argv[1], "-h") == 0 or strcmp(argv[1], "--help") == 0) {
		if (argc == 2) {
			showUsage(); exit(0);
		}
	}else if (strcmp(argv[1], "-v") == 0 or strcmp(argv[1], "--version") == 0){
		show_version(); exit(0);
	}

//	int opt, mask_val;
	int opt;
	int option_index;
	int32_t threadNum_tmp;
	outDir = OUT_DIR;
	simpleReg_t *simple_reg;
	string simple_reg_str;

	minRegsize = MIN_REG_SIZE;//
	exRegFold = EXREG_FOLD; //
	minCovFold = MIN_COV_FOLD;//m
	maxCovFold = MAX_COV_FOLD;//M
	minLocRatio = MIN_LOCAL_RATIO;//l
	minStrandRatio = MIN_STRAND_RATIO;//s
	minisizeRatio = MIN_SIZE_RATIO;//
	isizeSdevFold = ISIZE_DEV_FOLD;//
	minmateRadio = MIN_MATE_RATIO;//
	random_norm_reg_percent = RAND_NORM_REG_PERCENT;//
	threadNum_tmp = NUM_THREADS_PER_WORK;

	static struct option lopts[] = {
		{ "cov-off", no_argument, NULL, 0 },
//		{ "indel-off", no_argument, NULL, 0 },
		{ "IDC-off", no_argument, NULL, 0 },
		{ "ab-strand-off", no_argument, NULL, 0 },
		{ "ab-isize-off", no_argument, NULL, 0 },
		{ "ab-mate-off", no_argument, NULL, 0 },

		{ "min-regsize", required_argument, NULL, 0 },
		{ "extend-regsize-fold", required_argument, NULL, 0 },
		{ "min-cov-fold", required_argument, NULL, 0 },
		{ "max-cov-fold", required_argument, NULL, 0 },
		{ "min-locs-ratio", required_argument, NULL, 0 },
		{ "min-abstrand-ratio", required_argument, NULL, 0 },
		{ "min-abisize-ratio", required_argument, NULL, 0 },
		{ "min-abmate-ratio", required_argument, NULL, 0 },
		{ "norm-reg-percent", required_argument, NULL, 0 },
		{ "isize-sdev-fold", required_argument, NULL, 0 },
//		{ "max_proc_running_minutes_call", required_argument, NULL, 0 },
//		{ "technology", required_argument, NULL, 0 },
//		{ "include-decoy", required_argument, NULL, 0 },
//		{ "gt_min_sig_size", required_argument, NULL, 0 },
//		{ "gt_size_ratio_match", required_argument, NULL, 0 },
//		{ "gt_min_alle_ratio", required_argument, NULL, 0 },
//		{ "gt_max_alle_ratio", required_argument, NULL, 0 },
//		{ "version", no_argument, NULL, 'v' },
		{ "help", no_argument, NULL, 'h' },
		{ NULL, 0, NULL, 0 }
	};

	while( (opt = getopt_long(argc, argv, ":o:t:vh", lopts, &option_index)) != -1 ){
		switch(opt){
//			case 'a': ;refFile = stoi(optarg); break;
//			case 'b': inBamFile = stoi(optarg); break;
//
//			case 'r': regionFile = stoi(optarg); break;
			case 'o': outDir = optarg; break;
			case 't': threadNum_tmp = stoi(optarg);break;
//			case 'r': max_seg_size_ratio_usr = stod(optarg); break;
//			case 'e': minClipEndSize = stoi(optarg); break;
//			case 'x': expected_cov_assemble = stod(optarg); break;
//			case 'o': outDir = optarg; break;
//			case 'p': outFilePrefix = optarg; break;
//			case 't': threadNum_tmp = stoi(optarg); break;
			case 'v': show_version(); exit(0);
			case 'h': showUsage(); exit(0);
			case '?':
				if(optopt) cout << "unknown option '-" << (char)optopt << "'" << endl;
				else{ // Bad long option
					if (optind > 0 && strncmp(argv[optind-1], "--", 2) == 0) {
						cout << "unknown option '" << argv[optind-1] << "'" << endl;
					}
				}
				exit(1);
			case ':':
				if(optopt) cout << "the option '-" << (char)optopt << "' needs a value" << endl;
				else{
					if (optind > 0 && strncmp(argv[optind-1], "--", 2) == 0)
						cout << "the option '" << argv[optind-1] << "' needs a value" << endl;
				}
				exit(1);
			case 0: // long options
				if(parse_long_opt(option_index, optarg, lopts)!=0){
					showUsage();
					return 1;
				}
				break;
			default:
				cout << "Error: please specify the correct options for 'all' command" << endl;
				showUsage();
				return 1;
		}
	}

	load_from_file_flag = true;
	if(threadNum_tmp==0) num_threads = sysconf(_SC_NPROCESSORS_ONLN);
	else num_threads = (threadNum_tmp>=sysconf(_SC_NPROCESSORS_ONLN)) ? sysconf(_SC_NPROCESSORS_ONLN) : threadNum_tmp;

	outDir = deleteTailPathChar(outDir);

//	if(mask_val==1) maskMisAlnRegFlag = true;
//	else if(mask_val==0) maskMisAlnRegFlag = false;
//	else{
//		cout << "Error: Please specify the correct mask flag for mis-aligned regions using -M option." << endl << endl;
//		showUsage();
//		return 1;
//	}


	//non-option-input regionFile refFile inBamFile
	opt = argc - optind; // the reference and BAM file on the command line
	if(opt>=3) {
		regionFile = argv[optind];
		refFile = argv[optind+1];
		inBamFile = argv[optind+2];

		for(int i=optind+2; i<argc; i++){
			simple_reg_str = argv[i];
			simple_reg = allocateSimpleReg(simple_reg_str);
			if(simple_reg) limit_reg_vec.push_back(simple_reg);
		}
		if(limit_reg_vec.size()) limit_reg_process_flag = true;
	}else{
		cout << "Error: Please specify the reference file and coordinate-sorted BAM file." << endl << endl;
		showUsage();
		return 1;
	}

	return 0;
}

int Paras::parse_long_opt(int32_t option_index, const char *optarg, const struct option *lopts){
	int ret = 0;
//	size_t find_pos;
	string lower_str;

	string opt_name_str = lopts[option_index].name;
	if(opt_name_str.compare("sample")==0){ // "sample"
		if(optarg) sample = optarg;
		else{
			cout << "Error: Please specify the correct sample name using --sample option." << endl << endl;
			ret = 1;
		}
	}else if(opt_name_str.compare("cov-off")==0){ // cov-true
		cov_flag = false;
//	}else if(opt_name_str.compare("indel-off")==0){ //
//		indel_flag = false;
	}else if(opt_name_str.compare("ab-strand-off")==0){ //
		abstrand_flag = false;
	}else if(opt_name_str.compare("IDC-off")==0){ //
		IDC_flag = false;
	}else if(opt_name_str.compare("ab-isize-off")==0){ //
		abisize_flag = false;
	}else if(opt_name_str.compare("ab-mate-off")==0){ //
		abmate_flag = false;
		//
	}else if(opt_name_str.compare("min-regsize")==0){ // "min-regsize"
		minRegsize = stoi(optarg);
	}else if(opt_name_str.compare("extend-regsize-fold")==0){ // "mask-noisy-region"
		exRegFold = stod(optarg);
	}else if(opt_name_str.compare("min-cov-fold")==0){ // min-cov-fold
		minCovFold = stod(optarg);
	}else if(opt_name_str.compare("max-cov-fold")==0){ // max-cov-fold
		maxCovFold = stod(optarg);
	}else if(opt_name_str.compare("min-locs-ratio")==0){ //
		minLocRatio = stod(optarg);
	}else if(opt_name_str.compare("min-abstrand-ratio")==0){ //
		minStrandRatio = stod(optarg);
	}else if(opt_name_str.compare("min-abisize-ratio")==0){ //
		minisizeRatio = stod(optarg);
	}else if(opt_name_str.compare("min-abmate-ratio")==0){ //min-abmate-ratio
		minmateRadio = stoi(optarg);
	}else if(opt_name_str.compare("norm-reg-percent")==0){ //norm-reg-percent
		random_norm_reg_percent = stod(optarg);
	}else if(opt_name_str.compare("isize-sdev-fold")==0){ //
		isizeSdevFold = stod(optarg);
//	}else if(opt_name_str.compare("max_proc_running_minutes_assemble")==0){ // max_proc_running_minutes_assemble
//		max_proc_running_minutes_assemble = stoi(optarg);
//		if(max_proc_running_minutes_assemble<ULTRA_LOW_PROC_RUNNING_MINUTES){
//			cout << "Error: The specified maximum process running minutes is too small '" << max_proc_running_minutes_assemble << "', please specify a larger one at least " << ULTRA_LOW_PROC_RUNNING_MINUTES << "." << endl << endl;
//			ret = 1;
//		}
//	}else if(opt_name_str.compare("max_proc_running_minutes_call")==0){ // max_proc_running_minutes_call
//		max_proc_running_minutes_call = stoi(optarg);
//		if(max_proc_running_minutes_call<ULTRA_LOW_PROC_RUNNING_MINUTES){
//			cout << "Error: The specified maximum process running minutes is too small '" << max_proc_running_minutes_call << "', please specify a larger one at least " << ULTRA_LOW_PROC_RUNNING_MINUTES << "." << endl << endl;
//			ret = 1;
//		}
//	}else if(opt_name_str.compare("technology")==0){ // technology
//		technology = optarg;
//		lower_str.resize(technology.size());
//		transform(technology.begin(), technology.end(), lower_str.begin(), ::tolower);
//		if(lower_str.compare(PACBIO_CLR_TECH_STR)==0 or lower_str.compare(PACBIO_CCS_TECH_STR)==0 or lower_str.compare(NANOPORE_TECH_STR)==0){
//			technology = lower_str;
//		}else{
//			cout << "Error: Please specify the correct sequencing technology using '--technology' option." << endl << endl;
//			ret = 1;
//		}
//	}else if(opt_name_str.compare("include-decoy")==0){ // include-decoy
//		include_decoy = true;
//	}else if(opt_name_str.compare("gt_min_sig_size")==0){ // "gt_min_sig_size"
//		gt_min_sig_size = stoi(optarg);
//	}else if(opt_name_str.compare("gt_size_ratio_match")==0){ // "gt_size_ratio_match"
//		gt_size_ratio_match = stof(optarg);
//	}else if(opt_name_str.compare("gt_min_alle_ratio")==0){ // "gt_min_alle_ratio"
//		gt_min_alle_ratio = stof(optarg);
//	}else if(opt_name_str.compare("gt_max_alle_ratio")==0){ // "gt_max_alle_ratio"
//		gt_max_alle_ratio = stof(optarg);
	}


	return ret;
}

void Paras::showUsage(){
	cout << "Program: " << PROG_NAME << " (" << PROG_DESC << ")" << endl;
	cout << "Version: " << PROG_VERSION << " (using htslib " << hts_version() << ")" << endl;
	cout << "Usage:  misclus [options] <REG_FILE> <REF_FILE> <BAM_FILE>" << endl << endl;

	cout << "Description:" << endl;
	cout << "   REG_FILE   Region file that need to be processed (required)." << endl;
	cout << "              Each region in REG_FILE should be specified in the format" << endl;
	cout << "              per line: CHR|CHR:START-END." << endl;
	cout << "   REF_FILE   Reference file (required)" << endl;
	cout << "   BAM_FILE   Coordinate-sorted BAM file (required)" << endl;

	cout << "Options: " << endl;
	cout << "   -t INT     Number of threads [" << NUM_THREADS_PER_WORK << "]. 0 for the maximal number" << endl;
	cout << "              of threads in machine" << endl;
	cout << "   -o DIR     Output directory [" << OUT_DIR << "]" << endl;

	cout << "   --norm-reg-percent DOUBLE" << endl;
	cout << "              The percentage of total regions to generate random normal regions" << "[" << RAND_NORM_REG_PERCENT << "]" << endl;
	cout << "              These randomly generated normal regions are used for comparison with" << endl;
	cout << "              candidate misassembly regions during clustering." << endl;
	cout << "   --min-regsize  INT" << endl;
	cout << "              The minimal region size for feature extraction. The region will be extended to" << endl;
	cout << "              INT bp on both sides if it is shorter than INT bp." << "[" << MIN_REG_SIZE << "]" << endl;
	cout << "   --extend-regsize-fold DOUBLE" << endl;
	cout << "              Auxiliary region is used for an candidate misassembly region." << endl;
	cout << "              The auxiliary region is determined by extending the candidate region by" << endl;
	cout << "              DOUBLE fold of the size of candidate region on both sides." << "[" << EXREG_FOLD << "]" << endl;
	cout << "   --cov-off  disable the clustering on coverage feature" << endl;
	cout << "      --min-cov-fold DOUBLE" << endl;
	cout << "              Coverage is anomalous when it is lower than DOUBLE*meancov for an auxiliary region," << endl;
	cout << "              where meancov is the average coverage of the auxiliary region." << "[" << MIN_COV_FOLD << "]" << endl;
	cout << "              This option takes effect only if '--cov-off' option is not specified." << endl;
	cout << "      --max-cov-fold DOUBLE" << endl;
	cout << "              Coverage is anomalous when it is higher than DOUBLE*meancov for an auxiliary region." << "[" << MAX_COV_FOLD << "]" << endl;
	cout << "              This option takes effect only if '--cov-off' option is not specified." << endl;
	cout << "   --IDC-off  disable the clustering on IDC feature" << endl;
	cout << "      --min-locs-ratio DOUBLE" << endl;
	cout << "              Mask the ultra-low coverage locations whose coverage is lower" << endl;
	cout << "              than DOUBLE*meancov. When indels and clips appear at ultra-low coverage location," << endl;
	cout << "              the performance on IDC is unreliable." << "[" << MIN_LOCAL_RATIO << "]" << endl;
	cout << "              This option takes effect only if '--IDC-off' option is not specified." << endl;
	cout << "   --ab-strand-off  disable the clustering on anomalous orientation feature" << endl;
	cout << "      --min-abstrand-ratio DOUBLE" << endl;
	cout << "              Search the anomalous orientation clumps in which the anomalous" << endl;
	cout << "              orientation reads take a higher ratio than DOUBLE at each location." << "[" << MIN_STRAND_RATIO << "]" << endl;
	cout << "              This option takes effect only if '--ab-strand-off' option is not specified." << endl;
	cout << "   --ab-isize-off  disable the clustering on anomalous insert size feature" << endl;
	cout << "      --min-abisize-ratio DOUBLE" << endl;
	cout << "              Search the anomalous insert size clumps in which the anomalous" << endl;
	cout << "              insert size reads take a higher ratio than DOUBLE at each location." << "[" << MIN_SIZE_RATIO << "]" << endl;
	cout << "              This option takes effect only if '--ab-isize-off' option is not specified." << endl;
	cout << "      --isize-sdev-fold DOUBLE" << endl;
	cout << "              Normal insert size is usually in the range of " << endl;
	cout << "              [isize-DOUBLE*sdev, isize+DOUBLE*sdev] bp, otherwise it is anomalous." << "[" << ISIZE_DEV_FOLD << "]" << endl;
	cout << "              This option takes effect only if '--ab-isize-off' option is not specified." << endl;
	cout << "   --ab-mate-off  disable the clustering on abnormal mate read pair feature"  << endl;
	cout << "      --min-abmate-ratio DOUBLE" << endl;
	cout << "              Search the mate read pair clumps in which the abnormal mate read pairs" << endl;
	cout << "              take a higher ratio than DOUBLE at each location." << "[" << MIN_MATE_RATIO << "]" << endl;
	cout << "              This option takes effect only if '--ab-mate-off' option is not specified." << endl;

	cout << "   -v         show version information" << endl;
	cout << "   -h         show this help message and exit" << endl;
}

void Paras::outputParas(){

	cout << "Program: " << PROG_NAME << " (" << PROG_DESC << ")" << endl;
	cout << "Version: " << PROG_VERSION << " (using htslib " << hts_version() << ")" << endl << endl;

	if(regionFile.size()) cout << "Region file: " << regionFile << endl;
	if(refFile.size()) cout << "Reference file: " << refFile << endl;
	if(inBamFile.size()) cout << "Bam file: " << inBamFile << endl;

	cout << "Number of threads: " << num_threads << endl;
	if(outDir.size()) cout << "Output directory: " << outDir << endl;

	cout << "Percentage of regions to randomly generate normal regions: " << random_norm_reg_percent << endl;
	cout << "Minimal region size for feature extraction: " << minRegsize << endl;
	cout << "Auxiliary region is extended by the the folds of candidate misassembly region: " << exRegFold << endl;
	if(cov_flag){
		cout << "Coverage clustering: on" << endl;
		cout << "\tMinimal coverage fold to determine normal coverage region: " << minCovFold << endl;
		cout << "\tMaximal coverage fold to determine normal coverage region: " << maxCovFold << endl;
	}
//	if(indel_flag){
//		cout << "Indels clustering: on" << endl;
//		cout << "\tUltra-low coverage ratio threshold: " << minLocRatio << endl;
//	}
	if(abstrand_flag){
		cout << "Abnormal strand orientation clustering: on" << endl;
		cout << "\tMinimal ratio threshold for abnormal strand orientation reads: " << minStrandRatio << endl;
	}
	if(abisize_flag){
		cout << "Abnormal insert size clustering: on" << endl;
		cout << "Maximal ratio threshold for anomalous insert size clumps: " << minisizeRatio << endl;
		cout << "\tFolds of standard deviation for anomalous insert size clumps: " << isizeSdevFold << endl;// need to polish
	}
	if(abmate_flag){
		cout << "Abnormal mate read pair clustering: on" << endl;
		cout << "\tMaximal ratio for abnormal mate read pair clumps: " << minmateRadio << endl;
	}
	if(IDC_flag){
		cout << "IDC clustering: on" << endl;
		cout << "\tUltra-low coverage ratio threshold: " << minLocRatio << endl;
	}
}

// estimate single parameter
size_t Paras::estimateSinglePara(size_t *arr, size_t n, double threshold, size_t min_val){
	size_t i, total, val = 1;
	double total1;

	for(i=0, total=0; i<AUX_ARR_SIZE; i++) total += arr[i];
	if(total>0)
//		cout << "total=" << total << endl;
		for(i=0, total1=0; i<AUX_ARR_SIZE; i++){
			total1 += arr[i];
//			cout << "\t" << i << ": " << arr[i] << ", " << total1/total << endl;
			if(total1/total>=threshold){
				val = i + 1;
				break;
			}
		}
	if(val<min_val) val = min_val;

	return val;
}
