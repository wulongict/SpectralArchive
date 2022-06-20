#ifndef SPECTRASTSIMSCORES_HPP_
#define SPECTRASTSIMSCORES_HPP_

#include <fstream>



/*

Program       : Spectrast
Author        : Henry Lam <hlam@systemsbiology.org>                                                       
Date          : 03.06.06 


Copyright (C) 2006 Henry Lam

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307, USA

Henry Lam
Institute for Systems Biology
401 Terry Avenue North 
Seattle, WA  98109  USA
hlam@systemsbiology.org

*/

/* Class: SpectraSTSimScores
 * 
 * Manages the similarity scores. 
 * Needs to be modified for other scoring functions and for different output types.
 */


using namespace std;

class SpectraSTSimScores {

public:
    SpectraSTSimScores();

    SpectraSTSimScores(SpectraSTSimScores &other);

    ~SpectraSTSimScores();

    SpectraSTSimScores &operator=(SpectraSTSimScores &other);

    // printing methods for different output types
    void printFixedWidth(ofstream &fout);

    void printTabDelimited(ofstream &fout);

    void printPepXML(ofstream &fout);

    void printHtml(ofstream &fout);

    static void printHeaderTabDelimited(ofstream &fout);

    static void printHeaderFixedWidth(ofstream &fout);

    static void printHeaderHtml(ofstream &fout);

    // the various scores
    double dot;
    double delta;
    double dotBias;
    double precursorMzDiff;
    double pValue;
    double KSScore;
//  double pms;
//  double fftDot;
    unsigned int hitsNum;
    double hitsMean;
    double hitsStDev;
    double fval;
    pair<double, string> openMod;

    double calcFval(double fractionDelta, bool useDotBias);

    double calcOldFval();

    double calcSP5Fval(unsigned int sampleSize);

    double imposeDotBiasPenalty();

    unsigned int firstNonHomolog;

};

#endif /*SPECTRASTSIMSCORES_HPP_*/
