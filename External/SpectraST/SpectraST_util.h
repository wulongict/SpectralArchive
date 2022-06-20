/*
 *
 * declarations for TPP's util.c
 *
 */

#ifndef SPECTRAST_UTIL_H_INCLUDED
#define SPECTRAST_UTIL_H_INCLUDED

#include <string>

#ifdef _MSC_VER
#define strcasecmp stricmp
#define strncasecmp strnicmp
#endif

void getword(char *word, char *line, char stop);

char *makeword(char *line, char stop);

char *fmakeword(FILE *f, char stop, int *cl);

char x2c(char *what);

void unescape_url(char *url);

void plustospace(char *str);

int rind(char *s, char c);

// copy one file to another, return 0 on success
int copy_file(const char *fromName, const char *toName);

void send_fd(FILE *f, FILE *fd);


int get_line(char *s, int n, FILE *f); // was "getline", which conflicts with stdio.h

/** Basic routine to escape spaces in a string with
 *  the "\" character.  This can be used to create
 *  a POSIX path/filename from a string.  The escaped
 *  space path/filename can then be used in a "system()"
 *  call.
 */
char *escape_spaces(char *inStr);

#define fixpath(in) (in)

//
// return nonzero if fname appears to be absolute
// does not check for actual existence of file
//
int isAbsolutePath(const char *fname);

// accept either .xml or .pepXML as valid pepXML filename extensions
// return NULL if invalid, else return pointer to .ext
const char *hasValidPepXMLFilenameExt(const char *fname);

// accept either .xml or .protXML as valid protXML filename extensions
// return NULL if invalid, else return pointer to .ext
const char *hasValidProtXMLFilenameExt(const char *fname);

//
// get the WEBSERVER_ROOT env var, with path separators cleaned up
//
const char *getWebserverRoot();

//
// in-place replacement of path with newpath
// so /blarg/foo becomes newpath/foo
// 
void replace_path(char *path, int pathmaxlen, const char *newpath);

void replace_path(std::string &path, const char *newpath);

//
// return a copy of the input filename with an absolute path in the webserver_root tree
// caller must free() the result
//
char *translate_relative_webserver_root_path_to_absolute_filesystem_path(const char *path);

//
// remove the webserver root portion of the input path, if any
// so /inetpub/wwwroot/foo/blah becomes /foo/blah
//
void translate_absolute_filesystem_path_to_relative_webserver_root_path(std::string &path);

void translate_absolute_filesystem_path_to_relative_webserver_root_path(char *path);


//
// return a copy of the input filename with the filesystem webserver root prepended
// caller must free() the result
//
char *prepend_webserver_root(const char *path);

// fix up path's root dir if needed
void resolve_root_dir(char *path, int buflen);

// do we need to fix up the path at all?
// returns:
//   malloc'd path of file that can be opened,
//   or just a copy of input path
// in any case caller must free()
char *resolve_root(const char *path);

//
// get the WEBSERVER_URL env var, with path separators cleaned up
//
const char *getWebserverUrl();

// get the canonical ppepxml filename extension, including period
const char *get_pepxml_dot_ext();

// get the canonical protxml filename extension, including period
const char *get_protxml_dot_ext();

//
// case insensitive strstr
//
char *strstri(const char *str, const char *strCharSet);

char *strstrir(const char *str, const char *strCharSet); // finds rightmost match

//
// getcwd with path separators cleaned up
//
char *safepath_getcwd(char *buf, int buflen);

std::string safepath_getcwd(); // std::string version

//
// is c a path seperator for linux or windows?
//
int isPathSeparator(char c); // return nonzero if c is a path seperator / or \ .

//
// does path end with a path seperator?
//
bool endsWithPathSeparator(const char *path);

int isAbsolutePath(const char *fname); // return nonzero if fname appears to a full path

char *makeFullPath(const char *fname); // return a full path version of fname - caller must free()
void makeFullPath(std::string &fname); // make a full path version of fname, in place

std::string XMLEscape(const std::string &s);

char *findRightmostPathSeparator(char *path); // return pointer to rightmost / or \ , or NULL.
const char *findRightmostPathSeparator_const(const char *path); // return pointer to rightmost / or \ , or NULL.
int findRightmostPathSeparator(const std::string &str); // return position to rightmost / or \, or npos.


#endif // UTIL_H_INCLUDED
