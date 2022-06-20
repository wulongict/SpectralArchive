/*
 *

 * Revision 1.13  2006/05/25 22:25:43  brendanx
 * Win32 build for VC++7.1 includes:
 * - All relative paths (to new module win_lib)
 * - Compile time warnings fixed
 * - Build instructions
 * - All char[line_width_] allocations moved off stack to avoid stack overflow
 * - Better use of PATH for XML_ONLY build
 *
 * Also in this check-in (wish they could have been separate):
 * - PeptideProphet support for X! Tandem Native scoring (v0.01)
 *
 * Revision 1.11  2006/03/17 23:27:57  pcbrefugee
 * added read_dta_or_out_from_tgz_file() to assist in Mascot2XML regression tests
 *
 * Revision 1.10  2006/01/06 20:05:23  pcbrefugee
 * There were a few places in the TPP code with the comment "//Convert to same case first for lack of strstri" wherein filenames being checked for inclusion of the webserver_root path got trashed from the *nix case sensitive filename point of view.  The solution was just to implement strstri() (in util.c).
 *
 * Revision 1.9  2005/12/16 20:03:45  pcbrefugee
 * added a bit more sophistication to regression test handling: more crossplatform tolerant, more informative error messages
 *
 * Revision 1.8  2005/12/08 22:05:10  pcbrefugee
 * need to include <errno.h> for linux build
 *
 * Revision 1.7  2005/12/08 21:58:32  pcbrefugee
 * add a utility routine to copy one file by name to another
 *
 * Revision 1.6  2005/11/21 23:48:30  pcbrefugee
 * portability tweak (stricmp vs strcasecmp)
 *
 * Revision 1.5  2005/11/21 23:38:02  pcbrefugee
 * portability tweak (stricmp vs strcasecmp)
 *
 * Revision 1.4  2005/11/21 23:11:36  pcbrefugee
 * 1) Added integrated regression test support for major C++ implemented TPP components called by xinteract.  See TESTING doc for details.
 * 2) Moving toward use of .pepXML instead of .xml as standard pepXML filename extension (much as we do with .mzXML), but still actually using .xml at this time - the actual choice has been virtualized, though.  Still need to modify perl CGIs to handle the switch (or rather, to handle any extension thrown at them, as continued support for .xml extension for pepXML files is required).
 * 3) Some cleanups of unused variables etc that were causing noisy compiles, and portability tweaks.
 *
 * Revision 1.3  2005/10/20 00:27:59  pcbrefugee
 * rename util.c's getline() to get_line(), which won't conflict with stdio.h
 *
 * Revision 1.2  2005/10/19 23:49:12  pcbrefugee
 * TPP components now display release version info and build numbers.  This text from the Makefile explains:
 *
 * #
 * # A note about TPP Version and Build Numbers
 * #
 * # The file src/common/TPPVersion.h controls the version info displayed by TPP components.
 * # This info is furnished to the various exe and cgi files by linking with an object file
 * # that also contains a build number (timestamp, actually).  The source code for this object
 * # file is generated automatically, see target TPPVersionInfo.cxx: below.  For perl scripts
 * # the info is generated as an include file TPPVersionInfo.pl.
 * #
 * # The build number (timestamp) is refreshed when src/common/TPPVersion.h is changed,
 * # or during a "make clean all".  The intent is that the build number be the same across
 * # all TPP components delivered to the end user.  It could of course be refreshed with
 * # every "make all" but that would mean every app would relink every time, which would bog
 * # down the development cycle.  If you wish to force a build number update without doing
 * # a "make clean", use "make new_buildnum".
 * #
 * #
 *
 * Revision 1.1  2005/06/20 19:21:20  dshteyn
 * First commit after reorg
 *
 * Revision 1.2  2005/02/09 20:13:36  adkeller
 * update before TPP branch
 *
 * Revision 1.1  2003/03/11 00:46:43  rhubley
 *   Several changes:
 * 	- Moved the location of util.c since it is now shared by
 * 	  both bin and cgi applications.
 *         - Made sure programs would compile appropriately under
 *           windows and unix.
 *         - Fixed problem with missing type definitions in Jimmy's
 *           recent changes.
 *         - Final checkin before release 6.0
 *
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>
#include <sys/stat.h>
#include <iostream>
#include <string>

#ifdef _MSC_VER // if under Visual Studio

#include <direct.h> // for getcwd 
#else

#include <unistd.h>

#endif

#include "SpectraST_util.h"
#include "SpectraST_constants.h"

#define LF 10
#define CR 13

using namespace std;

void getword(char *word, char *line, char stop) {
    int x = 0, y;

    for (x = 0; ((line[x]) && (line[x] != stop)); x++)
        word[x] = line[x];

    word[x] = '\0';
    if (line[x]) ++x;
    y = 0;

    while (line[y++] = line[x++]);
}

char *makeword(char *line, char stop) {
    int x = 0, y;
    char *word = (char *) malloc(sizeof(char) * (strlen(line) + 1));

    for (x = 0; ((line[x]) && (line[x] != stop)); x++)
        word[x] = line[x];

    word[x] = '\0';
    if (line[x]) ++x;
    y = 0;

    while (line[y++] = line[x++]);
    return word;
}

char *fmakeword(FILE *f, char stop, int *cl) {
    int wsize;
    char *word;
    int ll;

    wsize = 102400;
    ll = 0;
    word = (char *) malloc(sizeof(char) * (wsize + 1));

    while (1) {
        word[ll] = (char) fgetc(f);
        if (ll == wsize) {
            word[ll + 1] = '\0';
            wsize += 102400;
            word = (char *) realloc(word, sizeof(char) * (wsize + 1));
        }
        --(*cl);
        if ((word[ll] == stop) || (feof(f)) || (!(*cl))) {
            if (word[ll] != stop) ll++;
            word[ll] = '\0';
            return word;
        }
        ++ll;
    }
}

char x2c(char *what) {
    register char digit;

    digit = (what[0] >= 'A' ? ((what[0] & 0xdf) - 'A') + 10 : (what[0] - '0'));
    digit *= 16;
    digit += (what[1] >= 'A' ? ((what[1] & 0xdf) - 'A') + 10 : (what[1] - '0'));
    return (digit);
}

void unescape_url(char *url) {
    register int x, y;

    for (x = 0, y = 0; url[y]; ++x, ++y) {
        if ((url[x] = url[y]) == '%') {
            url[x] = x2c(&url[y + 1]);
            y += 2;
        }
    }
    url[x] = '\0';
}

void plustospace(char *str) {
    register int x;

    for (x = 0; str[x]; x++) if (str[x] == '+') str[x] = ' ';
}

int rind(char *s, char c) {
    register int x;
    for (x = (int) strlen(s) - 1; x != -1; x--)
        if (s[x] == c) return x;
    return -1;
}

int get_line(char *s, int n, FILE *f) { // "getline" conflicts with stdio.h
    register int i = 0;

    while (1) {
        s[i] = (char) fgetc(f);

        if (s[i] == CR)
            s[i] = fgetc(f);

        if ((s[i] == 0x4) || (s[i] == LF) || (i == (n - 1))) {
            s[i] = '\0';
            return (feof(f) ? 1 : 0);
        }
        ++i;
    }
}

//
// case insensitive strstr
//
char *strstri(const char *str, const char *strCharSet) {
    char *lstr = strdup(str);
    char *lstrCharSet = strdup(strCharSet);
    char *result;
    char *cp;
    // strlwr
    for (cp = lstr; *cp; cp++) {
        *cp = tolower(*cp);
    }
    for (cp = lstrCharSet; *cp; cp++) {
        *cp = tolower(*cp);
    }
    // search
    result = strstr(lstr, lstrCharSet);
    if (result) {
        // result pointer must be to original data
        result = (char *) str + (result - lstr);
    }
    // no leaks!
    free(lstr);
    free(lstrCharSet);
    return result;
}

//
// case insensitive strstr, finds rightmost match
//
char *strstrir(const char *str, const char *strCharSet) {
    char *search = strstri(str, strCharSet);
    char *result = search;
    while (search && *search) {
        search = strstri(result + 1, strCharSet);
        if (search) {
            result = search;
        }
    }
    return result;
}

// copy one file to another, return 0 on success
int copy_file(const char *fromName, const char *toName) {
    int result = 0;
    FILE *from = fopen(fromName, "rb");
    if (from) {
        FILE *to = fopen(toName, "wb");
        if (to) {
            send_fd(from, to);
            fclose(to);
        } else {
            result = errno;
        }
        fclose(from);
    } else {
        result = errno;
    }
    return result;
}

void send_fd(FILE *f, FILE *fd) {
    char c[0xfff];
    int n;
    while (n = (int) fread(c, 1, sizeof(c), f)) {
        fwrite(c, 1, n, fd);
    }
}

string XMLEscape(const string &s) {
    string ret;
    for (size_t i = 0; i < s.length(); i++) {
        char ch = s.at(i);
        switch (ch) {
            case '<':
                ret.append("&lt;");
                break;
            case '>':
                ret.append("&gt;");
                break;
            case '&':
                ret.append("&amp;");
                break;
            case '"':
                ret.append("&quot;");
                break;
            case '\'':
                ret.append("&apos;");
                break;
            case '\n':
                ret.append(" ");
                break;
            case '\r':
                ret.append(" ");
                if (i + 1 < s.length() && s.at(i + 1) == '\n')
                    i++;
                break;
            default:
                ret.append(1, ch);
                break;
        }
    }
    return ret;
}

/** Basic routine to escape spaces in a string with
 *  the "\" character.  This can be used to create
 *  a POSIX path/filename from a string.  The escaped
 *  space path/filename can then be used in a "system()"
 *  call.
 */
char *escape_spaces(char *inStr) {
    int i;
    int space_count = 0;
    int inStrLen = 0;
    char *pinStr = inStr;
    char *outStr;
    char *poutStr;

    while (pinStr = strchr(pinStr, ' ')) {
        pinStr++;
        space_count++;
    }

    if (space_count > 0) {
        poutStr = outStr =
                (char *) malloc(sizeof(char) * (strlen(inStr) + 1 + space_count));

        if (outStr == NULL) {
            printf("Error mallocing memory in escape_spaces()!\n");
            return ('\0');
        }

        for (i = 0; inStr[i]; i++) {
            if (inStr[i] == ' ') {
                *poutStr++ = '\\';
            }
            *poutStr++ = inStr[i];
        }
        *poutStr = '\0';

        return (outStr);

    }

    return (inStr);

}

// accept either .xml or .pepXML as valid pepXML filename extensions
// return NULL if invalid, else return pointer to .ext
const char *hasValidPepXMLFilenameExt(const char *fname) {
    const char *result = fname ? strstrir(fname, get_pepxml_dot_ext()) : NULL; // usual ext?
    if (!result) { // no, see what's there, though
        result = fname ? strrchr(fname, '.') : NULL;
    }
    if (result) {
        if (strcasecmp(result, ".xml") && // pre Jan 2008
            strcasecmp(result, ".pepXML") && // never really took off...
            strcasecmp(result, get_pepxml_dot_ext())) { // user config, default .pep.xml
            result = NULL;
        }
    }
    return result;
}

// accept either .xml or .protXML as valid protXML filename extensions
// return NULL if invalid, else return pointer to .ext
const char *hasValidProtXMLFilenameExt(const char *fname) {
    const char *result = fname ? strstrir(fname, get_protxml_dot_ext()) : NULL; // usual ext?
    if (!result) { // no, see what's there, though
        result = fname ? strrchr(fname, '.') : NULL;
    }
    if (result) {
        if (strcasecmp(result, ".xml") && // pre Jan 2008 (actually -prot.xml, usually)
            strcasecmp(result, ".protXML") && // never really took off
            strcasecmp(result, get_protxml_dot_ext())) { // user config, default .prot.xml
            result = NULL;
        }
    }
    return result;
}

//
//
// is c a path separator for linux or windows?
//
int isPathSeparator(char c) {
    return (c == '/');
}

// convert rel path to abs path - caller must free() result
char *makeFullPath(const char *fname) {
    char *result;
    if (isAbsolutePath(fname)) {
        result = strdup(fname);
    } else {
        char szFile[SIZE_FILE];
        int len;
        safepath_getcwd(szFile, SIZE_FILE);
        len = (int) strlen(szFile);
        strncat(szFile, "/", SIZE_FILE - len++);
        strncat(szFile, fname, SIZE_FILE - len);
        result = strdup(szFile);
    }
    return result;
}

void makeFullPath(std::string &fname) {
    char *str = makeFullPath(fname.c_str());
    fname = str;
    free(str);
}


//
// local helper func for cached env read, with path separators cleaned up
//
static char *check_env_var(const char *env, char *buf, int buflen, bool bIsPath, const char *defaultval = NULL) {
    if (!buf[0]) { // first access?
        const char *e = getenv(env);
        if (e) {
            strncpy(buf, e, buflen);
        } else {
            buf[0] = -1; // mark as failed
        }
        if ((buf[0] <= 0) && defaultval) { // no env var value, use default if provided
            strncpy(buf, defaultval, buflen);
        }
        // tidy up path seps if needed
        if (bIsPath) {
            char *cp;
            for (cp = buf; *cp; cp++) {
                if ('\\' == *cp) {
                    *cp = '/';
                }
            }
            // add trailing path sep if needed
            if ((cp > buf) && ('/' != *(cp - 1))) {
                if ((cp + 1 - buf) < buflen) {
                    *cp++ = '/';
                    *cp = 0;
                }
            }
        }
    }
    return (buf[0] > 0) ? buf : NULL; // buf[0]<0 means env var does not exist
}

//
// get the WEBSERVER_ROOT env var, with path separators cleaned up
//
static char wsroot[1024] = {0};

const char *getWebserverRoot() {
    return check_env_var("WEBSERVER_ROOT", wsroot, sizeof(wsroot), true);
}

//
// return a copy of the input filename with a webserver root relative path
// caller must free() the result
//
char *translate_relative_webserver_root_path_to_absolute_filesystem_path(const char *path) {
    char *result;
    char *pathcopy = strdup(path);
    const char *szWebserverRoot = getWebserverRoot();
    if (!szWebserverRoot) {
        return (strdup(path));
    }
    char *szWebserverRootDup = strdup(szWebserverRoot);
    const char *pStr;
    if (!findRightmostPathSeparator(pathcopy)) { // no path info
        char *bbuf = (char *) malloc(1025 + strlen(pathcopy));
        safepath_getcwd(bbuf, 1023);
        strcat(bbuf, "/");
        strcat(bbuf, pathcopy);
        free(pathcopy);
        pathcopy = bbuf;
    }
    pStr = strstri(pathcopy, szWebserverRootDup);
    if (pStr == NULL) {
        free(szWebserverRootDup);
        free(pathcopy);
        return strdup(path); // presumably set up already
    }
    result = (char *) malloc(strlen(PEPXML_STD_XSL) + strlen(pathcopy) + 1);
    pStr = pathcopy + strlen(szWebserverRootDup);
    sprintf(result, "%s%s", PEPXML_STD_XSL, pStr);
    free(szWebserverRootDup);
    free(pathcopy);
    return result;
}

//
// in-place replacement of path with newpath
// so /blarg/foo becomes newpath/foo
//
void replace_path(char *path, int pathmaxlen, const char *newpath) {
    char *fname = findRightmostPathSeparator(path);
    fname = strdup(fname ? (fname + 1) : path);
    strncpy(path, newpath, pathmaxlen);
    strncat(path, fname, pathmaxlen - strlen(path));
    free(fname);
}

void replace_path(std::string &path, const char *newpath) {
    int len;
    char *tmp = (char *) malloc(len = (int) (path.length() + strlen(newpath) + 2));
    replace_path(tmp, len, newpath);
    path = tmp;
    free(tmp);
}


//
// remove the webserver root portion of the input path, if any
// so /inetpub/wwwroot/foo/blah becomes /foo/blah
//
void translate_absolute_filesystem_path_to_relative_webserver_root_path(std::string &path) {
    char *p = strdup(path.c_str());
    translate_absolute_filesystem_path_to_relative_webserver_root_path(p);
    path = p;
    free(p);
}

void translate_absolute_filesystem_path_to_relative_webserver_root_path(char *path) {
    const char *szWebserverRoot = getWebserverRoot();
    if (!szWebserverRoot) {
        return;
    }
    char *tmpStr = strstri(path, szWebserverRoot);
    if (tmpStr != NULL) {
        tmpStr += strlen(szWebserverRoot);
        char *relpath = strdup(tmpStr);
        if ('/' != *relpath) {
            *path++ = '/';
        }
        strcpy(path, relpath);
        free(relpath);
    }
}

//
// get the WEBSERVER_ROOT env var, with path separators cleaned up
//
static char wsurl[1024] = {0};

const char *getWebserverUrl() {
    if (!wsurl[0]) {
        const char *e = getenv("WEBSERVER_URL");
        if (e) {
            int len;
            strncpy(wsurl, e, sizeof(wsurl));
            len = (int) strlen(wsurl);
            if (wsurl[len - 1] != '/') {
                wsurl[len] = '/';
                wsurl[len + 1] = '\0';
            }
        }
    }
    return wsurl[0] ? wsurl : NULL;
}

//
// getcwd with path separators cleaned up
//
char *safepath_getcwd(char *buf, int buflen) {
    char *cp;
    char *result;
    result = getcwd(buf, buflen);
    for (cp = buf; *cp; cp++) {
        if ('\\' == *cp) {
            *cp = '/';
        }
    }
    return result;
}

std::string safepath_getcwd() {
    char buf[1024];
    std::string result = safepath_getcwd(buf, sizeof(buf));
    return result;
}

//
// return nonzero if fname appears to be absolute
// does not check for actual existence of file
//
int isAbsolutePath(const char *fname) {
    return (isPathSeparator(*fname));
}

//
// does path end with a path separator?
//
bool endsWithPathSeparator(const char *path) {
    return path && *path && isPathSeparator(path[strlen(path) - 1]);
}

// get the canonical pepxml filename extension, including period
static char pepxmldotext[1024] = {0};

const char *get_pepxml_dot_ext() {
    if (!pepxmldotext[0]) { // first pass
        check_env_var("PEPXML_EXT", pepxmldotext, sizeof(pepxmldotext), false, DEFAULT_PEPXML_FILENAME_DOTEXT);
        if (!strchr(pepxmldotext, '.')) { // user didn't put a dot in there?
            memmove(pepxmldotext + 1, pepxmldotext, strlen(pepxmldotext) + 1);
            pepxmldotext[0] = '.';
        }
    }
    return pepxmldotext;
}

// get the canonical protxml filename extension, including period
static char protxmldotext[1024] = {0};

const char *get_protxml_dot_ext() {
    if (!protxmldotext[0]) { // first pass
        check_env_var("PROTXML_EXT", protxmldotext, sizeof(protxmldotext), false, DEFAULT_PROTXML_FILENAME_DOTEXT);
        if (!strchr(protxmldotext, '.')) { // user didn't put a dot in there?
            memmove(protxmldotext + 1, protxmldotext, strlen(protxmldotext) + 1);
            protxmldotext[0] = '.';
        }
    }
    return protxmldotext;
}

const char *findRightmostPathSeparator_const(const char *path) { // return pointer to rightmost / or \ .
    const char *result = path + strlen(path);
    while (result-- > path) {
        if (isPathSeparator(*result)) {
            return result;
        }
    }
    return NULL; // no match
}

int findRightmostPathSeparator(const std::string &str) {
    int slashPos = (int) str.find_last_of('/');
    if (slashPos == string::npos) {
        slashPos = (int) str.find_last_of('\\');
    }
    return slashPos;
}

char *findRightmostPathSeparator(char *path) { // return pointer to rightmost / or \ .
    return (char *) findRightmostPathSeparator_const(path);
}
