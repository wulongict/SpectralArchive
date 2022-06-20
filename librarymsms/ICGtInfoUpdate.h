//
// Created by wulong on 12/2/19.
//

#ifndef MYTOOL_ICGTINFOUPDATE_H
#define MYTOOL_ICGTINFOUPDATE_H

class DataFile;
class SPsmAnnotation;
class CTable;



class ICGtInfoUpdate {
public:
    virtual ~ICGtInfoUpdate()= default;
    virtual bool updateGtInfo(SPsmAnnotation &gtinfo)=0;
};

class AnnotationDataFile : public ICGtInfoUpdate
{
    DataFile &m_df;
public:
    explicit AnnotationDataFile(DataFile &df);
    bool updateGtInfo(SPsmAnnotation &gtinfo) override;
};

class CsvAnnotation: public ICGtInfoUpdate
{
    CTable &m_psm;
public:
    explicit CsvAnnotation(CTable &psm);
    bool updateGtInfo(SPsmAnnotation &gtinfo) override;
};


class RawDataFile : public ICGtInfoUpdate
{
    DataFile &m_df;
public:
    explicit RawDataFile(DataFile &df);
    bool updateGtInfo(SPsmAnnotation &gtinfo) override;
};


#endif //MYTOOL_ICGTINFOUPDATE_H
