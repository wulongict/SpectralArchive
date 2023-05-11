#ifndef __FILEUTILITY_H__
#define __FILEUTILITY_H__
// #include "dependentlibs/msbasic/commonfun.h"
// #include "dependentlibs/msbasic/CDebugMode.h"
// #include "../librarymsms/CMzFileReader.h"
// #include "FastCgiInterface.h"
// #include "dependentlibs/spectralIndex/CKMeans.h"
// #include "CSpectralArchive.h"
// #include <thread>
// #include "CTimerSummary.h"
// #include "Util.h"
#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;
// #include <filesystem>
// namespace fs = std::filesystem;


std::vector<fs::path> get_mz_file_list(const fs::path& dir_path);
std::vector<fs::path> get_file_list(const fs::path& dir_path, std::vector<std::string> ext_list);
#endif /* __FILEUTILITY_H__ */
