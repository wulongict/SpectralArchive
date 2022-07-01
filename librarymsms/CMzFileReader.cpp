//
// Created by wulong on 8/1/18.
//
#include "spdlog/spdlog.h"
#include <iostream>
#include "Util.h"
#include "CMzFileReader.h"
#include "ProteomicsDataTypes.h"
#include "CThreadsPool.h"
#include "PeakList.h"
#include <memory>
#include "XMLFileParser.h"
#include "ICQuery.h"


CMzFileReader::CMzFileReader(string mzxml_filelist, bool overwrite, bool islist, bool rmParentIon,
                             double localMaxHalfWidth, int minPeakNum, bool verbose, string mzFileName) {
    m_mzFileName = mzFileName; // to be used.
    m_mzs = nullptr;
    m_minPeakNum = minPeakNum;
    m_localMaxHalfWidth = localMaxHalfWidth;
    m_PeakNumPerSpectrum = 50;
    m_mzxml_filelist = mzxml_filelist;
    m_norm2_ptr=nullptr;

    init(overwrite, islist, rmParentIon);
}

CMzFileReader::CMzFileReader(DataFile &df, bool overwrite, bool rmParentIon, SpectraFilter *sf,
                             double localMaxHalfWidth, int minPeakNum, bool verbose) {
                                m_mzs = nullptr;

    m_norm2_ptr=nullptr;
    m_minPeakNum = minPeakNum;
    m_localMaxHalfWidth = localMaxHalfWidth;
	cout << "sf = " << sf << endl;
	m_PeakNumPerSpectrum = 50;
	m_mzxml_filelist = df.getSourceFileName();

	if(not File::isExist(getMzFilename()) or overwrite)    {
	    ofstream fout(getMzFilename(), ios::out | ios::binary);
//        m_scanInfo.appendFile(df.getSourceFileName());
        toMzScanFilePair(fout, df, rmParentIon, sf, false);
        m_scanInfo.clear(); // remove all the information.
//        m_scanInfo.getScanFileReader(m_scanInfo.getFileNum() - 1).write_specinfo(); // todo together
    }

    loadScanMzToMemory(false, false);
	spdlog::get("A")->info("mz file is created!\n");
}

// dump scaninfo to disk.
CMzFileReader::~CMzFileReader() {
    m_scanInfo.exportToCombinedFile();
    cout << "Releasing space of all spectra from memory" << endl;
    delete[] m_mzs;
}

void CMzFileReader::init(bool redo, bool is_file_list, bool rmParentIon) {
    //cout << "Starting loading mzfile " << endl;
//    spdlog::get("A")->info("loading file ...");
    if (not File::isExist(getMzFilename(),true) or redo) {
        vector<string> filelist;
        if (is_file_list) {
            string mzXMLList = getListFile();
            filelist = readlines(mzXMLList);
        } else {
            filelist = {m_mzxml_filelist};
        }
        createWithList(rmParentIon, filelist);
    }
// reload the mz file again.
    loadScanMzToMemory(is_file_list, false);
    initPeakNum(false);
    spdlog::get("A")->info("mz/scan file is created!\n");
}


void CMzFileReader::toMzScanFilePair(ofstream &fout_mzFile, DataFile &df, bool rmParentIon, SpectraFilter *sf,
                                     bool verbose) {
    m_scanInfo.appendFile(df.getSourceFileName());

    // note: the SpectraFilter is not used;
    //  sf = nullptr;
    const int PeakNum = getPeakNumPerSpec();
	int totalspecnum = df.getSpectrumNum();
	int ms2specnum = df.getMS2SpectrumNum();
	int ms2count = 0;

	vector<vector<uint16_t >> batchMS(ms2specnum, vector<uint16_t>(PeakNum));
	//Progress ps(totalspecnum, "Generating mz format");

	for (int spec_id_k = 0; spec_id_k < totalspecnum; spec_id_k++) {
		//ps.increase();
		CSpectrum *spec = df.getSpectrum(spec_id_k);

		if(nullptr!=sf and sf->skip(spec)){
			continue;
		}

		if (spec == nullptr) {
			cout << "empty spec : " << spec << endl;
			continue;
		}
		if (spec->getMSLevel() != 2) {
			continue;
		}
		vector<double> mz, intensity;
		bool removeLowIntensePeaks = true;
		bool rmIsotopicPeaks = true;
        spec->getAllPeaks(mz, intensity, removeLowIntensePeaks, rmParentIon, rmIsotopicPeaks, m_localMaxHalfWidth);

		if (mz.empty()) {
			//cout << "Bad mz" << endl;
		}

		PeakList pl;
		pl.setM_intensityList(intensity);
		pl.setM_mzList(mz);

		pl.KeepTopN(PeakNum);
		pl.rankingAsIntensity(PeakNum);

		vector<double> mzs, intens;
		mzs = pl.getM_mzList();
		intens = pl.getM_intensityList();
		int _peaknum = mzs.size();
		if (mzs.size() < PeakNum and _peaknum > m_minPeakNum) {
			mzs.resize(PeakNum, 0);
			intens.resize(PeakNum, 0);
		} else if(_peaknum <= m_minPeakNum)   {
		    // put dummy peak there. so every spectra is on the hyper sphere now
		    mzs.assign(PeakNum, 0);
		    intens.assign(PeakNum,0);
            _peaknum = 1;
            mzs[0] = 0.5;
            intens[0] = 1000;
        }
		const double *mzptr = mzs.data();
		const double *intensityptr = intens.data();
		getCompactForm(mzptr, intensityptr, batchMS[ms2count]);
		m_scanInfo.appendSpecInfo(_peaknum, spec->m_rt_in_sec, spec->m_parentMz, spec->m_precursor_charge, spec->m_scanNum, ms2count);
		ms2count++;
	}
	batchMS.resize(ms2count);
	if(verbose)cout << "Writing mz file to disk" << endl;
	for (auto & k : batchMS) {
		fout_mzFile.write((char *)k.data(), sizeof(uint16_t) * PeakNum);
	}

    m_scanInfo.getScanFileReader(m_scanInfo.getFileNum() - 1).write_specinfo();
    m_specnum += ms2count; // ms2count updated.
}

int CMzFileReader::getPeakNumPerSpec() const {
    return m_PeakNumPerSpectrum;
}

long CMzFileReader::getSpecNum() { return m_specnum; }

// control all of the scan file put them into one single file.
void CMzFileReader::loadScanMzToMemory(bool isList, bool verbose) {
      //CScanInfo

    if(verbose) cout << "start loading scan mz to memory" << endl;
    loadMzFile(false);

    m_scanInfo.load(m_mzxml_filelist, isList, verbose);

    if(m_scanInfo.getTotalScanNum() != m_specnum)    {
        cout << "spectrum num in mz file : " << m_specnum << " while in meta info " << m_scanInfo.getTotalScanNum() << endl;
        cout << "[Error] Scan Header and mz file are inconsistent" << endl;
        throw runtime_error("Error: Scan header file and mz file are inconsistent!") ;
    }
}

void CMzFileReader::loadMzFile(bool verbose) {
    string mzfile = getMzFilename();
    long filebytes = 0;
    File::getfilesize(mzfile, filebytes);
    // cout << "mz file and filebytes: " << mzfile << " " << filebytes << endl; 
    long peaknum = filebytes / 2;
    m_specnum = peaknum / m_PeakNumPerSpectrum;

    // if(verbose)
    // cout << "filebytes: " << filebytes << " peaknum " << peaknum << " specnum " << m_specnum << endl;
    if(m_mzs != nullptr) {
        // the mz file object is not empty.
        // cout << "will release the pointer first!" << " " << m_mzs << " pointer" << endl;
        delete[] m_mzs;
        m_mzs = nullptr;
    }
    m_mzs = new uint16_t[peaknum];
    // cout << "creating a pointer for spectra num: " << m_specnum << endl;
    ifstream fin(mzfile.c_str(), ios::in | ios::binary);
    fin.read((char *) m_mzs, filebytes);
    fin.close();
}

uint16_t *CMzFileReader::getSpecBy(long queryindex)  {
    return m_mzs + queryindex * m_PeakNumPerSpectrum;
}


void CMzFileReader::append(DataFile &df, bool rmParentIon, SpectraFilter *sf){
    ofstream fout(getMzFilename(), ios::app | ios::binary);
//    m_scanInfo.appendFile(df.getSourceFileName());
    toMzScanFilePair(fout, df, rmParentIon, sf, false);
//    m_scanInfo.getScanFileReader(m_scanInfo.getFileNum() - 1).write_specinfo();
    // write the  whole file
    // do we need to export it here?
    // m_scanInfo.exportToCombinedFile();  // export in destructor.
    // you have to close the file before you open it.
    fout.close();
    

// // why reload?
//     cout << "Start to reload" << endl;


//     // update the object.
//     m_scanInfo.clear();
//     loadScanMzToMemory(true, false);
//     initPeakNum(false);
//     spdlog::get("A")->info("mz/scan object is refreshed!\n");

}

void CMzFileReader::createWithList(bool rmParentIon, const vector<string> &filelist) {
    ofstream fout(getMzFilename(), ios::out | ios::binary);
    int total_file_num = filelist.size();

    for (int j = 0; j < total_file_num; j++) {
        cout << "Processing file " << j+1 << " / " << total_file_num << " : " << filelist[j] << endl;
        string filename = filelist[j];

        DataFile df(filename,0,-1);

        toMzScanFilePair(fout, df, rmParentIon, nullptr, false);

    }
    m_scanInfo.clear();
}

string CMzFileReader::getListFile() {
    return m_mzxml_filelist;
}

string CMzFileReader::getMzFilename() {
    return m_mzxml_filelist + ".mz";
}

// to calculate the score for all
void CMzFileReader::calcDotProduct(int TopNPeak, int tol, uint16_t *queryspec, int blockSize,
                                   vector<long> &indexlist, vector<int> &scores) {

    if(indexlist.empty()){
        cout << "Error: empty index list" << endl;
        return;
    } else{

    }
    scores.assign(indexlist.size(),0);
    vector<int> vecformY;
    get_vector_form(queryspec, tol, vecformY);

    bool debug = false;
    for(int i = 0; i < indexlist.size(); i ++)    {
            if (indexlist.at(i) == -1) cout << "Warning: index list contains -1, which is invalid index" << endl;
            scores[i]= calculate_dot_product_with_vecfrom(indexlist.at(i), vecformY, TopNPeak, debug, 0);
    }
}

class CDistTask:public ICThreadTask {
    vector<long> &retIdx;
    double tolerance;
    uint16_t *m_queryptr;
    vector<float> &accDist;
    CMzFileReader *m_p;
    vector<int> &dpscore;
public:
    CDistTask(CMzFileReader *p, vector<float> &acc_dist, uint16_t *query_ptr, double tol, vector<long> &retIdx,
              vector<int>&dp_score): accDist(acc_dist), retIdx(retIdx), dpscore(dp_score){
        tolerance = tol;
        m_queryptr = query_ptr;
        m_p = p;
    }
    void run() override{
        if(retIdx.empty())   {
//            cout << "Empty index list in " << __FUNCTION__ << endl;
                return;
        }
         m_p->distOnSpec(m_queryptr, retIdx, tolerance, accDist,dpscore);
    }
};


void CMzFileReader::dpscore(double tolerance, vector<vector<long>> &allRetIdx, int threadnum,
                            vector<vector<float>> &accDist, ICQuery &query,vector<vector<int>> &dpscores) {
    vector<CDistTask> vtask;
    for(int i = 0; i < query.size(); i ++)    {
        vtask.emplace_back(CDistTask(this, accDist[i],query.getPtrUint16(i),tolerance, allRetIdx[i],dpscores[i]));
    }

    vector<ICThreadTask*> vTasks;
    for(auto & i : vtask)    {
        vTasks.push_back(&i);
    }
    CTaskPool runDist(vTasks,false,false);
    runDist.start(threadnum, "L2 distance");
}

void CMzFileReader::initPeakNum(bool verbose) {
//    Progress ps(getSpecNum(),"initialize peaknum vector");
    m_peaknum.assign(getSpecNum(),0);
    for(long i = 0; i < getSpecNum(); i ++)  {
//        ps.increase();
        m_peaknum[i] = getMzSpec(i).getPeakNum();
    }

    vector<long> histogram(getPeakNumPerSpec()+1,0);
    for(auto n: m_peaknum) histogram[n] ++;

    if(verbose)    {
        cout << "PeakNum\t#Spectra\tFrequency" << endl;
        for(int i = 0; i < histogram.size(); i ++)   {
            cout << i << "\t" << histogram[i] << "\t" << histogram[i] * 1.0 / getSpecNum() << endl;
        }
    }
}

float CMzFileReader::getSquaredNorm(long queryindex) {
    if(m_norm2_ptr==nullptr)   {
//        cout << "First call--- initializing..." << endl;
        m_norm2_ptr=make_shared<vector<float>>(getSpecNum(),0);

        Progress ps(getSpecNum(), "squared norm");
        for(long i = 0; i < getSpecNum(); i ++)  {
            ps.increase();
            uint16_t  *p=getSpecBy(i);
            (*m_norm2_ptr)[i]=ICMzFile::getSquaredNorm(p, getPeakNumPerSpec());
        }
//        cout << "initialization done" << endl;
    }
    return (*m_norm2_ptr)[queryindex];
}

int CMzFileReader::getPeakNum(long queryindex) {
    if(m_peaknum.empty())    {
        initPeakNum(false);
    }
    if(queryindex <0 or queryindex >= getSpecNum()){
        cout << "Error: query index " << queryindex  <<" out of range: [" << 0 << ":"  << getSpecNum()-1 << "]" << endl;
    }
    return m_peaknum[queryindex];
}

MzSpecInfo &CMzFileReader::getSpecMetaInfo(long queryidx, string &filename) {
    return m_scanInfo.getMzSpecInfo(queryidx, filename);
}

string CMzFileReader::getClassName() {
    return "CMzFileReader";
}

float CMzFileReader::getSquaredNorm(uint16_t *p) {
    return ICMzFile::getSquaredNorm(p);
}

// read mzXML/mzML file
// add scan info to a vector. MzSpecInfo object
// write scan info to file
void CScanFile::init() {
    long cnts = m_resIdxOffset;
    cout << "Creating scan header file" << endl;
    DataFile df(m_datafilename,0,-1);
    int totalspecnum = df.getSpectrumNum();
    SpectraFilter *sf = nullptr;
    int ms2count = 0;

    for (int i = 0; i < totalspecnum; i++) {
        CSpectrum *spec = df.getSpectrum(i);
        if((nullptr!=sf and sf->skip(spec)) or spec==nullptr or spec->getMSLevel() !=2) {
            continue;
        }

        cnts ++;
        infos.emplace_back(cnts, -1, spec->m_rt_in_sec, spec->m_parentMz, spec->m_precursor_charge, spec->m_scanNum, ms2count);
        ms2count++;
    }
    write_specinfo();
}

void CScanFile::readLines(istream &fin, long len) {
    if(isCreated()) {
        throw logic_error("Error: should be empty");
    }
    copy_n(istream_iterator<MzSpecInfo>(fin),len,back_inserter(infos));
}

// can we avoid to do this when db is available?
void CScanFile::create(long &cnts, bool verbose) {
    m_resIdxOffset = cnts;
    if(isCreated()) {
        return;
    }

    if(not load()){
        init();
    }
    cnts = getLastQueryIDX() ;
    if(verbose)cout <<"total scan number: " <<getLastQueryIDX() << "\nscans added: " <<  getScanNum() << " from file" << m_datafilename << endl;
}

void CScanFile::fixResIdx() {
    // it assumes that the cnt value in the infos are correct/ this could be wrong.
//    cout << "[Info] fixing the universal idx from current scan file" << endl;
    int fixed_residx_num = 0;
    for(auto &info : infos){
        // residx fixed
        if(info.residx != m_resIdxOffset + 1) {
            info.residx = m_resIdxOffset + 1;
            fixed_residx_num ++;
        }
        m_resIdxOffset ++;
    }
//    cout << "[Info] fixed universal idx in current scan file " << fixed_residx_num << endl;
}

MzSpecInfo &CScanFile::binarySearch(long queryidx) {
    MzSpecInfo t(queryidx,-1,-1,-1,-1,-1,-1);
    auto it = lower_bound(infos.begin(), infos.end(), t,[](const MzSpecInfo &x, const MzSpecInfo &y){return x.residx<y.residx;});
    if(it == infos.end())    {
        cout << "Error: queryidx " << queryidx << " was not found " << endl;
        throw logic_error("Error: queryidx not was found");
    } else{
        return *it;
    }
}

void CScanFile::write_specinfo(ofstream &fout) {
    // output all of the scan info into a single file.
    // specific file
    for(int i = 0; i < infos.size() ; i ++){
        fout << infos[i] << endl;
    }
}

void CScanFile::write_specinfo() {
    ofstream scanout(m_scanFilename, ios::out);
    write_specinfo(scanout);
    spdlog::get("A")->info("scan file saved as {}", m_scanFilename);
}

bool CScanFile::isCreated() {
    return getScanNum()>0;
}

int CScanFile::getScanNum() const {
    return infos.size();
}

string CScanFile::getFileName() {
    return m_datafilename;
}

long CScanFile::getLastQueryIDX() {
    if(infos.empty()){
        return m_resIdxOffset;
    }
    return infos.back().residx;
    cout << "Error: the index is not valid:" << endl;
    throw logic_error("Error! Invalid index") ;
}

CScanFile::CScanFile(string _filename) {
    m_datafilename = std::move(_filename);
    m_scanFilename = m_datafilename + ".scan";
    m_resIdxOffset = -1;
}

CScanFile::CScanFile(const CScanFile &other) {
    m_datafilename = other.m_datafilename;
    m_scanFilename = other.m_scanFilename;
    infos = other.infos;
    m_resIdxOffset = other.m_resIdxOffset;
}

bool CScanFile::load() {
    bool ret = false;
    if(File::isExist(m_scanFilename))    {
        ifstream fin(m_scanFilename, ios::in);
        cout << "loading file: " << m_scanFilename << endl;
        copy(istream_iterator<MzSpecInfo>(fin), istream_iterator<MzSpecInfo>(),
             back_inserter(infos));
        fixResIdx();
        ret = true;
    }
    return ret;
}


bool SpectraFilter::skip(CSpectrum *spec) {
    if(this) {
        cout << "its not nullptr " << this << endl;
    }else{
        cout << "its null " << this << endl;
    }
    if(!this or  this==nullptr) return false;
    cout << "this filter is " << (void *)this << " " << " is this the difference? " << (nullptr==this) << " not  equal :"  << (nullptr!=this)<< endl;
    spec->print();
    bool ret = false;
    if(spec == nullptr)    {
        cout << "Empty" << endl;
        ret = true;
    }  else if(spec->getMSLevel() !=2)   {
        cout << "MS1" << endl;
        ret = true;
    }  else   {
        cout << "[Warning] This part is not ready!!" << __FUNCTION__ << "--> line: " << __LINE__  << endl;
        int scan = spec->m_scanNum;
        cout << "Start scan searching size = " << m_psms.size() << endl;
        auto x = find_if(m_psms.begin(), m_psms.end(), [&scan](const PSMInfo &psminfo){return psminfo.start_scan==scan;});
        cout << "End of scan searching" << endl;
        if(x==m_psms.end())   {
            ret = true;
        }  else {
            m_scans.push_back(scan);
        }
        cout << "[testing here]" << endl;
    }
    return ret;
}

SpectraFilter::~SpectraFilter() {
    ofstream fout("outfile");
    for(int m_scan : m_scans) {
        fout<< m_scan << endl;
    }
}

SpectraFilter::SpectraFilter(string pepxml) {
    m_pepxml = pepxml;
    m_ppp = make_shared<PeptideProphetParser>(pepxml);
    m_ppp->filter_with_FDR(0.01, m_psms);
}

// save the combined scan file.
void CScanInfo::exportToCombinedFile() {
    ofstream fout(getCombinedScanFilename(), ios::out);
    cout << "[Info] exporting " << getFileNum() << " .scan files to " << getCombinedScanFilename() << endl;

    long lastResidx = -1;

    for (auto &eachFile: m_readers) {
        // cout << "processing " << eachFile.getFileName() << endl;
        lastResidx += eachFile.getScanNum();
        fout << lastResidx << "\t" << eachFile.getFileName() << endl;
        eachFile.write_specinfo(fout);
    }
    // cout << "scan file saved as: " << getCombinedScanFilename() << endl;
    // todo: No need to keep these scan files...
}

void CScanInfo::appendSpecInfo(int _peaknum, float rt, float mz, int charge, int scan, int ms2count) {
    m_lastResIdxOffset ++;
    m_readers.back().infos.emplace_back(MzSpecInfo(m_lastResIdxOffset, _peaknum, rt, mz, charge, scan, ms2count));
}

string CScanInfo::getFileName(int fileid) {
    return getScanFileReader(fileid).getFileName();
}

int CScanInfo::getFileId(long queryidx) {
    for(int i = 0; i < getFileNum(); i ++)    {
        long last_Idx = getScanFileReader(i).getLastQueryIDX();
        if(queryidx <= last_Idx) return i;
    }
    cout << "Error: queryidx out of range" << endl;
    throw "Error";
}

MzSpecInfo &CScanInfo::getMzSpecInfo(long queryidx, string &filename) {
    int fileid = getFileId(queryidx);
    filename = getFileName(fileid);
    return getScanFileReader(fileid).binarySearch(queryidx);
}

CScanFile &CScanInfo::getScanFileReader(int fileid) {
    if( fileid >= 0 and  getFileNum() > fileid )    {
        return m_readers[fileid];
    } else {
        cout << "Error: invalid file index: " << fileid << " out of range " << 0 << "~" << getFileNum() << endl;
        throw "Error: OOR";
    }
}

int CScanInfo::getFileNum() {
    return m_readers.size();
}


long CScanInfo::getTotalScanNum() const {
    return m_lastResIdxOffset + 1;
}

// append a new data file to get its scan info.
void CScanInfo::appendFile(string file) {
    m_readers.emplace_back(CScanFile(std::move(file)));
}

// create a list of readers.
void CScanInfo::loadFiles(const vector<string>& files) {
    for(const auto& file: files){
        appendFile(file);
    }
    for(auto &eachScanFile: m_readers){
        eachScanFile.create(m_lastResIdxOffset, false);  // Do not use varialbe m_cnt like this.
    }
}

CScanInfo::~CScanInfo() = default;

CScanInfo::CScanInfo() {
    m_lastResIdxOffset = -1;
}

// loading the combiend file.
// if it does not exist [fail to init]
// then try load individual scan file, then combine them.
void CScanInfo::loadCombinedFile() {
    if(not initWithCombinedFile()){
        vector<string> files= readlines(m_mzXMLListfile);
        loadFiles(files);
        exportToCombinedFile();
    }
}

bool CScanInfo::initWithCombinedFile() {
    if(File::isExist(getCombinedScanFilename(),true))    {
        string combined_filename = getCombinedScanFilename();
        long filesize;
        File::getfilesize(combined_filename, filesize);
        char *buf = new char[filesize+1];
        buf[filesize]='\0';
        {
        ifstream fin(getCombinedScanFilename(), ios::in);
        fin.read(buf,filesize);
        fin.close();
    }

        istringstream  iss(buf);

        cout << "reading file " << getCombinedScanFilename() << endl;
        string filename;
        long lastval = -1;
        iss >> m_lastResIdxOffset >> filename;
//        CountProgress pc(10,"loading scan file");
        long totalCounts = 0;
        while(!iss.eof())  {
//            pc.increase();
            long len = m_lastResIdxOffset - lastval;
            lastval = m_lastResIdxOffset;
            appendFile(filename);
            totalCounts += len;
    //            cout << "[Error] loading " <<getFileNum() <<"-ith file " << filename << " len " << len  << " total len: " << totalCounts << endl;
            if(len > 0) {
                m_readers.back().readLines(iss, len);
            }
            iss >> m_lastResIdxOffset >> filename;
        }
        delete buf;

        return true;
    }else{
        return false;
    }
}

// void CScanInfo::setCombinedFilename(string mzxmlfiles) {
//     m_mzXMLListfile = std::move(mzxmlfiles);
// }

void CScanInfo::load(string m_mzxml_filelist, bool isList, bool verbose) {
    if(verbose)cout << "loading meta info for mz files : isList: " << isList << endl;
    if(isList)    {
        if(verbose)cout << "loading list of mzXML files" << m_mzxml_filelist << endl;
        m_mzXMLListfile = m_mzxml_filelist;
        
        loadCombinedFile();

    } else {
        vector<string> files;
        if(verbose)cout << "loading single raw data file : " << m_mzxml_filelist << endl;
        files.push_back(m_mzxml_filelist);
        loadFiles(files);
    }
}

bool CScanInfo::isEmpty() {
    return m_readers.empty();
}

void CScanInfo::clear() {
    m_lastResIdxOffset = -1;
    vector<CScanFile>().swap(m_readers);
//    m_readers.swap(vector<CScanFile>());
}

shared_ptr<ICMzFile> CMzFileReaderFactory::create() {
    cout << "[Info] Creating CPU scorer Object" << endl;
    return getMzFileObject();
}

double MzSpecInfo::getPrecursorNeutralMass() const {
    const double PROTON_MASS = 1.007276;
    return (precursormz - PROTON_MASS) * charge ;
}

long MzSpecInfo::getResIdx() const {
    return residx;
}

istream &operator>>(istream &ism, MzSpecInfo &sf) {
    ism >> sf.residx >> sf.ms2idx >> sf.scan >> sf.precursormz >> sf.charge >> sf.rt;
    return ism;
}

ostream &operator<<(ostream &osm, const MzSpecInfo &sf) {
    osm << sf.residx << "\t" << sf.ms2idx << "\t" << sf.scan << "\t" << sf.precursormz << "\t" << sf.charge << "\t" << sf.rt;
    return osm;
}

MzSpecInfo::MzSpecInfo() :MzSpecInfo(-1,-1,-1,-1,-1,-1,-1)  {}

MzSpecInfo::MzSpecInfo(long _residx, int _peaknum, float _rt, float _parentMz, int _chg, int _scan, int _ms2idx) {
    residx = _residx;
    peaknum = _peaknum;
    rt = _rt;
    precursormz = _parentMz;
    charge = _chg;
    scan = _scan;
    ms2idx = _ms2idx;
}
