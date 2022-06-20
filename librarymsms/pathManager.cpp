//
// Created by wulong on 1/21/19.
//

#include "pathManager.h"
#include "Util.h"
#include <sys/stat.h>

// todo: to be merged with FILE name space
std::string GetCurrentWorkingDir(void) {
    char buff[FILENAME_MAX];
    GetCurrentDir(buff, FILENAME_MAX);
    std::string current_working_dir(buff);
    return current_working_dir;
}

bool createDirectory(string directoryname) {
    File::CFile fileObj(directoryname);
    struct stat info;
    // create parent folder, recursively?
    if(stat(fileObj.path.c_str(), &info)!=0){
        bool x = createDirectory(fileObj.path);
    }
    // Creating a directory
    bool ret = true;
    if (mkdir(directoryname.c_str(), 0777) == -1) {
        //cerr << "Warning in creating directory "<< directoryname << ":  " << strerror(errno) << endl;
        ret = false;
    }
    return ret;
}

