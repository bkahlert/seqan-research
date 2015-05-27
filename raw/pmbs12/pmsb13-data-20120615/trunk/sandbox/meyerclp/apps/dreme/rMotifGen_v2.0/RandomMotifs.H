
//============================================================================
//--------------------------- BEGIN OF RandomMotifs.H ------------------------
//============================================================================

//==========================================================================
// rMotifGen
// Eric C. Rouchka
// C. Timothy Hardin
// (c) 2004-2007, University of Louisville
//
// FILE:  RandomMotifs.H
//
//
//    This file is part of rMotifGen.
//
//    rMotifGen is free software; you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation; either version 2 of the License, or
//    (at your option) any later version.
//
//    rMotifGen is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with rMotifGen; if not, write to the Free Software
//    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
//
//==========================================================================



#ifndef MY_RANDOMMOTIFS_H
#define MY_RANDOMMOTIFS_H

#include "SequenceClass.H"
#include "CommonRoutines.H"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>
using namespace std;

//------------------- LOCAL FUNCTION PROTOTYPES --------------------
void printWelcome();
void printOverallHeader();
void printMotifHeader(int currMotif);

void allocateMotifArrays(bool *&motif_isRandomArr,
                           int *&motif_lenArr,
                           string *&motif_seqArr,
                           double *&motif_pctSeqConsArr,
                           PAM **&motif_PAMArr,
                           double **&motifPctArr,
                           double *&motif_pctConsensusConsArr,
                           int numMotifs,
                           int numResidues,
                           bool *&motif_isPalindromicArr);

void getMotifData(bool * &motif_isPalindromicArr,
                  bool * &motif_isRandomArr,
                  int  * &motif_lenArr,
                  string * &motif_seqArr,
                  int    * &motif_pctSeqConsArr,
                  int    * &motif_pctConsensusConsArr,
                  int    * &motif_pctAArr,
                  int    * &motif_pctCArr,
                  int    * &motif_pctGArr,
                  int    * &motif_pctTArr,
                  int    numMotifs);

void getUserDefinedMotifInfo(ifstream &infile, bool *&motif_isRandomArr,
                           string *&motif_seqArr, double *&motif_pctSeqConsArr,
                           double *&motif_pctConsensusConsArr,
                           PAM *&motif_PAMArr, char seqType, int i);

string setMotifLabel(int i);

void getRandomMotifInfo(ifstream &infile, bool *&motif_isRandomArr,
                        string *&motif_seqArr, double *&motif_pctSeqConsArr,
                        double *&motif_pctConsensusConsArr,
                        PAM **&motif_PAMArr, bool *&motif_isPalindromicArr,
                        int *&motif_lenArr, char seqType, int i,
                        int numResidues, char *ResidueLabels, double **&motifPctArr);

void checkValidPercentages(int pctA, int pctC, int pctG, int pctT, int errorCode);

void printToFile(SequenceClass *&S);

void fileMode(string FN, char seqType);

void getResidueFrequencies(ifstream &inFile, double **&motifPctArr,
                           char *ResidueLabels, string motifLabel, 
                           int i, int numResidues);

void getPAMValue(PAM **&motif_PAMArr, ifstream &infile,
                 int i, string motifLabel);

void getPctConserved(double *&motif_pctConsensusConsArr, ifstream &infile,
                 int i, string motifLabel);

#endif

//============================================================================
//------------------------- END OF RandomMotifs.H  ---------------------------
//============================================================================



