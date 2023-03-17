//
// Created by wulong on 12/2/19.
//

#include "ICGtInfoUpdate.h"
#include "ProteomicsDataTypes.h"
#include "Util.h"
#include "CPSMAnnotation.h"
#include <cassert>


bool CsvAnnotation::updateGtInfo(SPsmAnnotation &gtinfo) {
    bool ret = false;
    string key=gtinfo.mzxml_filename+ "."+to_string(gtinfo.ms2_scan);
    int rownum = m_psm.getRowByKey(key,m_psm.getColByHeader("spectrumname"));
    if(rownum!=-1)    {
        gtinfo.rfscore = atof(m_psm.getEntry(rownum, m_psm.getColByHeader("rf+")).c_str());
        ret = true;
    }


    return ret;

    // update value
//    m_psm.getRowByKey()

    return ret;
}

CsvAnnotation::CsvAnnotation(CTable &psm) :m_psm(psm){
    // read table

    // make new filename
    //map<string, int> col2index={{"filename",0}, {"ms2scan",2}};
    m_psm.appendHeader({"spectrumname"});
    for(int i = 0; i < m_psm.m_row; i ++)
    {
        // for each row , change the first column
        string filename = m_psm.getEntry(i,m_psm.getColByHeader("filename"));
        string scannum = m_psm.getEntry(i, m_psm.getColByHeader("scan"));
        string spectrumname = File::CFile(filename).filename +"."+scannum;
        m_psm.appendEntry(i,spectrumname);

    }
    //
    m_psm.build_table_index(m_psm.getColByHeader("spectrumname"));
    cout << "index created " << endl;
    m_psm.printRow(0);
    m_psm.printRow(1);

}

// Significance will be 1 for mgf/sptxt file.
bool AnnotationDataFile::updateGtInfo(SPsmAnnotation &gtinfo) {
    bool ret = false;
    CSpectrum *spec = m_df.getSpectrum(gtinfo.ms2idx);
    if(spec!= nullptr and spec->getScanNum()==gtinfo.ms2_scan)    {
        assert(spec->getScanNum()==gtinfo.ms2_scan);
        string peptidestr = spec->getSpectrumName();
        getmodificationfrompeptidestring(spec->getSpectrumName(),gtinfo);
        // added hcd
        gtinfo.m_collision_energy = spec->m_collision_energy;
        string filename=m_df.getSourceFileName();
        string ext = File::CFile(filename).ext;
        transform(ext.begin(), ext.end(), ext.begin(), ::tolower);

        gtinfo.significance=1;


        if(ext=="sptxt"){
            gtinfo.protein = spec->getProtein();// update protein
            if (string::npos == gtinfo.protein.find("DECOY")){
                // not found DECOY label, it is a target protein.
            }
            else{
                gtinfo.isDecoy = 1;
                gtinfo.significance=0;

            }
        }

        ret = true;
    } else  {
        cout << "spec -> " << spec << endl;
        cout << "spec->getScanNum()==gtinfo.ms2_scan " << (spec->getScanNum()==gtinfo.ms2_scan) << endl;
    }
    return ret;
}

AnnotationDataFile::AnnotationDataFile(DataFile &df) : m_df(df){}

ICFileParser::~ICFileParser() = default;

RawDataFile::RawDataFile(DataFile &df) : m_df(df){}

bool RawDataFile::updateGtInfo(SPsmAnnotation &gtinfo) {
    return false;
}

