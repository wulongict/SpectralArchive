#include "SpectraSTXHunterLibImporter.hpp"
#include "SpectraSTLog.hpp"
#include "SpectraSTConstants.hpp"
#include "FileUtils.hpp"
#include "Peptide.hpp"
#include "ProgressCount.hpp"
#include <iostream>
#include <sstream>
#include <string.h>

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

/* Class: SpectraSTXHunterLibImporter
 * 
 * Implements a library importer for the X!Hunter file format (ref http://www.thegpm.org/hunter/format.html).
 * Note that there is no guarantee that X!Hunter libraries will work well with SpectraST!
 */


extern bool g_verbose;
extern bool g_quiet;
extern SpectraSTLog *g_log;

// constructor
SpectraSTXHunterLibImporter::SpectraSTXHunterLibImporter(vector<string> &impFileNames, SpectraSTLib *lib,
                                                         SpectraSTCreateParams &params) :
        SpectraSTLibImporter(impFileNames, lib, params),
        m_numSkipped(0) {}

// destructor
SpectraSTXHunterLibImporter::~SpectraSTXHunterLibImporter() {

}

// import - prints the preamble, then loops over all files and import them one by one
void SpectraSTXHunterLibImporter::import() {

    for (vector<string>::iterator i = m_impFileNames.begin(); i != m_impFileNames.end(); i++) {
        string fullName(*i);
        makeFullPath(fullName);
        string quoted("\"" + fullName + "\"");
        string desc = m_params.constructDescrStr(quoted, ".hlf");
        m_preamble.push_back(desc);
    }

    m_lib->writePreamble(m_preamble);

    for (vector<string>::iterator i = m_impFileNames.begin(); i != m_impFileNames.end(); i++) {
        readFromFile(*i);
    }
}

// readFromCurFile - reads one .hlf file
void SpectraSTXHunterLibImporter::readFromFile(string &impFileName) {

    ifstream fin;
    if (!myFileOpen(fin, impFileName)) {
        g_log->error("CREATE", "Cannot open .hlf (X!Hunter) file \"" + impFileName + "\" for reading. File skipped.");
        return;
    }

    g_log->log("HLF IMPORT", "Importing .hlf file \"" + impFileName + "\".");

    if (g_verbose) {
        cout << "\nImporting spectra from .hlf (X!Hunter)library file..." << endl;
    }

    char buffer[MAX_LINE];

    // header
    int dummy = 0;
    unsigned int numSpectra = 0;
    fin.read((char *) (&dummy), sizeof(int));
    fin.read((char *) (&numSpectra), sizeof(unsigned int));
    fin.read(buffer, 248);

    // start the progress count
    ProgressCount pc(!g_quiet && !g_verbose, 1, (int) numSpectra);
    pc.start("\nImporting spectra from .hlf (X!Hunter) library file");

    for (unsigned int s = 0; s < numSpectra && !fin.eof(); s++) {

        pc.increment();

        double precursorMH = 0.0;
        int precursorCharge = 0;
        float spectralMagnitude = 0.0;
        float medianExpectation = 0.0;
        int peptideLen = 0;
        string peptideSeq("");
        int numPeaks = 0;
        vector<int> intensities;
        vector<float> mzs;

        fin.read((char *) (&precursorMH), sizeof(double));
        fin.read((char *) (&precursorCharge), sizeof(int));
        fin.read((char *) (&spectralMagnitude), sizeof(float));
        fin.read((char *) (&medianExpectation), sizeof(float));

        double precursorMz = precursorMH / (double) precursorCharge;

        fin.read((char *) (&peptideLen), sizeof(int));
        memset(buffer, '\0', MAX_LINE);
        fin.read(buffer, peptideLen);
        peptideSeq = buffer;

        fin.read((char *) (&numPeaks), sizeof(int));

        int i = 0;
        for (i = 0; i < numPeaks; i++) {
            unsigned char c;
            fin.read((char *) (&c), sizeof(unsigned char));
            intensities.push_back((int) c);
        }
        for (i = 0; i < numPeaks; i++) {
            float mz;
            fin.read((char *) (&mz), sizeof(float));
            mzs.push_back(mz);
        }

        SpectraSTPeakList *peakList = new SpectraSTPeakList(precursorMz, precursorCharge, numPeaks);
        peakList->setNoiseFilterThreshold(m_params.rawSpectraNoiseThreshold);

        for (i = 0; i < numPeaks; i++) {
            peakList->insert((double) (mzs[i]), (float) intensities[i], "", "");
        }

        vector<double> modMasses(peptideSeq.length(), 0.0);

        int numMods = 0;
        fin.read((char *) (&numMods), sizeof(int));
        int m = 0;
        for (m = 0; m < numMods; m++) {
            int pos;
            double modMass;
            fin.read((char *) (&pos), sizeof(int));
            fin.read((char *) (&modMass), sizeof(double));

            if (pos >= 1 && pos <= peptideSeq.length()) {
                modMasses[pos - 1] += modMass;
            }
        }

        stringstream modifiedPeptide;

        for (int aa = 0; aa < peptideSeq.length(); aa++) {

            if (modMasses[aa] > 0.001 || modMasses[aa] < 0.001) {

                int d = 0;
                // round the mod mass to nearest integer
                if (modMasses[aa] >= 0) {
                    d = (int) (modMasses[aa] + 0.5);
                } else {
                    d = (int) (modMasses[aa] - 0.5);
                }

                if (aa == 0 && (d == 42 || d == 144 || d == 1)) {
                    // assumed to be N-terminal mods
                    modifiedPeptide << "n[" << d + 1 << "]";
                } else {
                    double AAplusModMass = modMasses[aa] + (*Peptide::AAMonoisotopicMassTable)[peptideSeq[aa]];
                    modifiedPeptide << peptideSeq[aa] << '[' << (int) (AAplusModMass + 0.5) << ']';
                }

            } else {
                modifiedPeptide << peptideSeq[aa];
            }

        }

        int numProteins = 0;
        vector<string> proteins;
        fin.read((char *) (&numProteins), sizeof(int));
        int p = 0;
        for (p = 0; p < numProteins; p++) {
            int proteinLen;
            int startPos;
            memset(buffer, '\0', MAX_LINE);
            fin.read((char *) (&proteinLen), sizeof(int));
            fin.read(buffer, proteinLen);
            fin.read((char *) (&(startPos)), sizeof(int));
            proteins.push_back(buffer);
        }

        Peptide *pep = createPeptide(modifiedPeptide.str(), precursorCharge, "", "", "HLF");

        if (pep) {

            if (m_params.setDeamidatedNXST) {
                setDeamidatedNXST(pep);
            }

            stringstream comments;
            comments << "Spec=Consensus ";
            comments << "Fullname=" << pep->interactStyleFullWithCharge() << ' ';
            comments << "Mods=" << pep->mspMods() << ' ';
            comments << "MedianExpectation=" << medianExpectation << ' ';
            comments << "Prob=1.0000 ";
            comments << "Protein=\"" << numProteins;
            for (p = 0; p < numProteins; p++) {
                comments << '/' << proteins[p];
            }
            comments << "\" ";

            SpectraSTLibEntry *entry = new SpectraSTLibEntry(pep, comments.str(), "Normal", peakList);

            if (insertOneEntry(entry, "HLF")) {
                if (g_verbose) {
                    cout << "Importing record: " << pep->interactStyleWithCharge() << endl;
                }
                m_count++;
            }
            delete (entry);

        } else {
            m_numSkipped++;
        }

    }

    pc.done();
}





