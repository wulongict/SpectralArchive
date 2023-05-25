#ifndef __FILEUTILITY_H__
#define __FILEUTILITY_H__

#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;


std::vector<fs::path> get_mz_file_list(const fs::path& dir_path);
std::vector<fs::path> get_file_list(const fs::path& dir_path, std::vector<std::string> ext_list);
#endif /* __FILEUTILITY_H__ */
