#ifndef SPECTRAST_CONSTANTS_H
#define SPECTRAST_CONSTANTS_H

#define SIZE_BUF             8192
#define SIZE_FILE            1024
#define SIZE_PEP             128
#define INIT_INTERACT_LINES  1000
#define PERCENTAGE           0.75

#define USE_LOCAL_TIME 1

#define SIZE_BUF 8192

#ifndef Boolean
typedef unsigned int Boolean;
#endif

//
// CONFIGURATION DEFINITIONS
//

//
// Linux Installation
//

// Web Paths
#define WEBSERVER_HOME_URL "/tpp/"
#define CGI_BIN "/tpp/cgi-bin/"
#define HELP_DIR "/tpp/html/"

// Filesystem Paths
#define CGI_FULL_BIN "/local/TPP/cgi-bin/"
#define LOCAL_BIN "/local/TPP/bin/"
#define COMETLINKSDIR "/local/TPP/etc/"
#define LOCAL_HTML "/local/TPP/html/"

// schema
#define PEPXML_STD_XSL "/local/TPP/schema/"  // a webserver reference, to the directory, filename is hard coded

// schema namespace constants
#define PEPXML_NAMESPACE "http://regis-web.systemsbiology.net/pepXML"

// pepXML stuff that's platform independent
//#define PEPXML_FILENAME_DOTEXT ".pepXML"
#define DEFAULT_PEPXML_FILENAME_DOTEXT ".pep.xml"
#define PEPXML_NAMESPACE_PX "pepx"
#define PEPXML_SCHEMA "pepXML_v18.xsd"

// protXML stuff that's platform independent
//#define PROTXML_FILENAME_DOTEXT ".protXML"
#define DEFAULT_PROTXML_FILENAME_DOTEXT ".prot.xml"


#endif
