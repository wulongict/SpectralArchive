//
// Created by wulong on 10/20/19.
//

#include "randomSpecParser.h"
#include <cmath>
#include<string>
#include "MSReader.h"
#include "Spectrum.h"
//#include "MSObject.h"
#include <iostream>

#include "../../External/SpectraST/SpectraSTLib.hpp"
#include "../../librarymsms/Util.h"

// random-access spectrum reads not allowed with MGF format.
// it will retrieve the first scan in the mgf file.
bool getPeakListMzxmlMzml(string filename, int scannum, vector<double> &mz, vector<double> &intensity){
    MSToolkit::MSReader mstreader;
    vector<MSToolkit::MSSpectrumType> m{MSToolkit::MS1, MSToolkit::MS2};
    mstreader.setFilter(m);
    MSToolkit::Spectrum s;
//    cout << "Start reading file with mstoolkit, mstreader, readfile "  << endl;
    cout << filename << "\t" << scannum << endl;
    string datafilename;
    if(File::CFile(filename).ext=="scanlist"){
        ifstream fin(filename, ios::in);
        fin >> datafilename;
        cout << "datafilename is :" << datafilename << endl;
    }else{
        datafilename = filename;
    }
    bool status = mstreader.readFile(datafilename.c_str(), s, scannum);
//    cout << "after finished reading the file, status " << status << endl;
    if (status) {
//        cout << "Got spectra s  " << endl;
        mz.assign(s.size(), 0);
        intensity.assign(s.size(), 0);
        for (int i = 0; i < s.size(); i++) {
            mz[i] = s.at(i).mz;
            intensity[i] = s.at(i).intensity;
        }
    } else {
//        cout << "Fail" << endl;
    }
    return status;
}



bool getPeakListFromSplib(string filename, int scannum, vector<double> &mz, vector<double> &intensity) {
    bool status = false;
    SpectraSTSearchParams s;
    SpectraSTLib splib(filename,&s);
    SpectraSTMzLibIndex *p = splib.getMzLibIndexPtr();
    long entry_num = p->getEntryCount();
    cout << "entry num: " << entry_num << endl;
    p->sortEntriesByOffsets();

    fstream::off_type  offset = 0;
    int k = 1;
    while(p->nextSortedFileOffset(offset) )
    {
//        cout << k << endl;
        if (k < scannum)   {   k++;   continue;   }

        // Remeber to release the pointers

        SpectraSTLibEntry * entry = p->thisSortedEntry();
        unsigned int lib_id = entry->getLibId();
        string entry_name = entry->getName();

        SpectraSTPeakList *pklist = entry->getPeakList(); // no need to delete
        cout << "SpectraST lib entry: " << lib_id << "\t" << entry_name << endl;
//        cout << "SpectraST lib entry(same?): " << entry->getLibId() << "\t" << entry->getName() << endl;
//        pklist->printPeaks();
        unsigned  int pknum=pklist->getNumPeaks();
        mz.assign(pknum, 0);
        intensity.assign(pknum,0);
        for(unsigned int i = 0; i < pknum ; i ++)
        {
            Peak peak;
            pklist->getPeak(i,peak);
            mz[i]=peak.mz;
            intensity[i] = peak.intensity;
        }
        status = true;
        delete entry;
        // break
        break;
    }

    return status;
}

void printSpectrum(MSToolkit::Spectrum &s) {
    for (int i = 0; i < s.size(); i++) {
        cout << s.at(i).mz << "\t" << s.at(i).intensity << endl;
    }
}
bool getPeakList(string filename, int scannum, vector<double> &mz, vector<double> &intensity) {
    cout << "[Info] Looking for real peaks " << endl;
    bool status = false;
    if(File::CFile(filename).ext=="sptxt")
    {
        cout << "looking for splib file, load spectra from splib " << endl;
        filename = filename.substr(0,filename.size()-5)+"splib";
        status = getPeakListFromSplib(filename,scannum, mz,intensity);
    }
     else{
        cout << "loading mzXML/mzML" << endl;
        status = getPeakListMzxmlMzml(filename, scannum, mz,intensity);
    }
    return status;

}