#include "dependentlibs/msbasic/commonfun.h"
#include "dependentlibs/msbasic/CDebugMode.h"
#include "../librarymsms/CMzFileReader.h"
#include "FastCgiInterface.h"
#include "dependentlibs/spectralIndex/CKMeans.h"
#include "CSpectralArchive.h"
#include <thread>
#include "CTimerSummary.h"
#include "Util.h"
#include <filesystem>
namespace fs = std::filesystem;

std::vector<fs::path> get_mz_file_list(const fs::path& dir_path) {
    std::vector<fs::path> file_list;
    for (const auto& entry : fs::recursive_directory_iterator(dir_path)) {
        if (entry.is_regular_file() && (entry.path().extension() == ".mzXML" || entry.path().extension() == ".mzML")) {
            file_list.push_back(entry.path());
        }
    }
    return file_list;
}

// check if one string is end as another string as suffix. 
// return true if the string is end with the suffix
bool isEndWith(const std::string& str, const std::string& suffix)
{
    return str.size() >= suffix.size() &&
           str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0;
}

// search for specific file extension, including
// .mzXML, .pep.xml, .mzML, .mgf, .sptxt, .pepXML, .MGF etc. specifiec in lower case
// we store the extension in a vector, and then search for the file
// with the extension in the vector.
std::vector<fs::path> get_file_list(const fs::path& dir_path, std::vector<std::string> ext_list) {
    std::vector<fs::path> file_list;
    for (const auto& entry : fs::recursive_directory_iterator(dir_path)) {
        if (entry.is_regular_file()) {
            for (auto ext : ext_list){
                string fullpath = entry.path().string();
                transform(fullpath.begin(), fullpath.end(), fullpath.begin(), ::tolower);
                if (isEndWith(fullpath, ext)){
                    file_list.push_back(entry.path());
                } else{
                    // cout << "skipping unsupported file format in " << entry.path() << endl;
                }
            }
        }
    }
    return file_list;
}