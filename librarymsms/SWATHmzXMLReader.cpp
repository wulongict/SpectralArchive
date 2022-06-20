//
// Created by wulong on 7/12/15.
//

#include "SWATHmzXMLReader.h"
#include "SpectraST_cramp.hpp"
#include "PeakList.h"
using SpectraST_msms::cRamp;
using SpectraST_msms::rampRunInfo;
using SpectraST_msms::rampScanInfo;
using SpectraST_msms::rampPeakList;


vector<PeakList *> SWATHmzXMLReader::ReadSWATHmzXMLToPeakLists(mzXMLFilename f, int currentcyclenum, int cycle) {
    cout << "[Info] Loading SWATH " << currentcyclenum << "/" << cycle << " from file " << f << ". " << flush;
    int emptyScan = 0;

    cRamp *cramp = new cRamp(f.c_str());
    if (!cramp->OK()) {
        cout << "Cannot open file \"" << f << "\". File skipped." << endl;
        delete (cramp);
        throw "can not open file";
    }
    rampRunInfo *runInfo = cramp->getRunInfo();

    if (!runInfo) {
        cout << "Cannot open file \"" << f << "\". File skipped." << endl;

        throw "can not open file";

    }
    //rampInstrumentInfo* instr = cramp->getInstrumentInfo();
    vector<PeakList *> vpl;
    for (int k = 1; k <= cramp->getLastScan(); k++) {
        if ((k + cycle - currentcyclenum) % cycle != 0) continue;
        rampScanInfo *scanInfo = cramp->getScanHeaderInfo(k);



        // If the scan is NULL, skip this scan
        if (!scanInfo) {
            emptyScan++;
            //cout << "[Info] Number of empty scan: " << emptyScan  << endl;

            vpl.push_back(NULL);


            continue;
        }
        PeakList *p = new PeakList();
        p->setRTinSeconds(scanInfo->getRetentionTimeSeconds());


        // If the msLevel is not equal to MSLEVEL, skip this scan

        //cout << k << endl;

        rampPeakList *peaks = cramp->getPeakList(scanInfo->m_data.acquisitionNum);
        delete scanInfo;
        //cout << scanInfo->m_data.precursorScanNum<< " " << scanInfo->m_data.msLevel << endl;
        if (!peaks) {

            emptyScan++;
            //
            vpl.push_back(p);

            continue;

        }
        //double sumint = 0;

        for (int j = 0; j < peaks->getPeakCount(); j++) {
            double mz = peaks->getPeak(j)->mz;
            double intensity = (double) (peaks->getPeak(j)->intensity);
            p->InsertPeak(mz, intensity);

        }
        delete peaks;

        vpl.push_back(p);


    }
    if (emptyScan > 0)
        cout << "#(Empty scan): " << emptyScan << endl;
    else
        cout << endl;
    delete cramp;

    return vpl;
}

int SWATHmzXMLReader::GetSwathCycleNum(const mzXMLFilename& f) {
    cRamp *cramp = new cRamp(f.c_str());
    if (!cramp->OK()) {
        cout << "Cannot open file \"" << f << "\". File skipped." << endl;
        delete cramp;
        exit(0);
    }
    rampRunInfo *runInfo = cramp->getRunInfo();

    if (!runInfo) {
        cout << "Cannot open file \"" << f << "\". File skipped." << endl;
        exit(0);

    }
    int cycle = 0;
    rampScanInfo *scanInfo = nullptr;
    for (int k = 1; k <= cramp->getLastScan(); k++) {
        if (scanInfo) {
            delete scanInfo;
            scanInfo = nullptr;
        }
        scanInfo = cramp->getScanHeaderInfo(k);
        if (!scanInfo) { continue; }
        if (scanInfo->m_data.msLevel == 1 && cycle != 0) break;
        else cycle++;

        // If the scan is NULL, skip this scan


        // If the msLevel is not equal to MSLEVEL, skip this scan
        if (scanInfo->m_data.msLevel != 1) continue;
    }

    if (scanInfo) {
        delete scanInfo;
        scanInfo = nullptr;
    }

    delete cramp;
    cout << "[Resu] " << f << ": Cycle = " << cycle << endl;
    return cycle;

}

SWATHmzXMLReader::SWATHmzXMLReader() : mzXMLReader("1") {}

SWATHmzXMLReader::~SWATHmzXMLReader() = default;