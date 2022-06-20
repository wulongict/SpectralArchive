//
// Created by wulong on 7/12/15.
//

#ifndef PROJECT_SWATHMZXMLREADER_H
#define PROJECT_SWATHMZXMLREADER_H

#include "mzXMLReader.h"

class SWATHmzXMLReader : public mzXMLReader {
public:
    SWATHmzXMLReader();
    ~SWATHmzXMLReader();
    int GetSwathCycleNum(const mzXMLFilename& f);
    vector<PeakList *> ReadSWATHmzXMLToPeakLists(mzXMLFilename f, int currentcyclenum, int cycle);
};

#endif //PROJECT_SWATHMZXMLREADER_H
