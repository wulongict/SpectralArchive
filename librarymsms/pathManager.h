//
// Created by wulong on 1/21/19.
//

#ifndef MYTOOL_PATHMANAGER_H
#define MYTOOL_PATHMANAGER_H


#include <stdio.h>  /* defines FILENAME_MAX */
// #define WINDOWS  /* uncomment this line to use it for windows.*/
#ifdef WINDOWS
#include <direct.h>
const auto GetCUrrentDir = _getcwd;
//#define GetCurrentDir _getcwd
#else
#include <unistd.h>
// this is how we rename a function in c++11
const auto GetCurrentDir = getcwd;
//#define GetCurrentDir getcwd
#endif
#include <iostream>
#include <bits/stdc++.h>
#include <iostream>
#include <sys/stat.h>
#include <sys/types.h>
using namespace std;

std::string GetCurrentWorkingDir( void );
bool createDirectory(string directoryname );

#endif //MYTOOL_PATHMANAGER_H
