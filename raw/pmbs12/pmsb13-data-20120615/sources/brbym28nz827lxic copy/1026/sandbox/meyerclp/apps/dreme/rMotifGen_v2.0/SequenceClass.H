
//============================================================================
//------------------------ BEGIN OF SequenceClass.H --------------------------
//============================================================================

//==========================================================================
// rMotifGen
// Eric C. Rouchka
// C. Timothy Hardin
// (c) 2004-2007, University of Louisville
//
// FILE: SequenceClass.H
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

//==========================================================================

#ifndef MY__SEQCLASS_H
#define MY__SEQCLASS_H

#include "SubstitutionMatrices.H"
#include <string>
#include <cstdlib>
#include <iostream>
using namespace std;

string itos(int i);

class SequenceClass
{
   public:
      SequenceClass(int s, int nm, int ns, char 
                      *labelArr, double *tmpFreqArr, char st, int nr);
      ~SequenceClass();
      int getSeqLens();
      double getOverallPct(char residue);
      int getNumMotifs();
      int getNumSeqs();
      void MotifsToInclude();
      string randomSeq(int len, char *labelArr, long double *FreqArr, bool isPalindrome);
      void createRandomSequences();
      void createConsensusMotifs();
      bool isValid(int pos, bool *used, int len);
      void generateSequences();
      string formatSeq(string s, int width, int blockWidth, bool numbersOn);
      int calculateDigits(int len);
      void setValues(double **&FreqArr, int *&mLenArr, 
                     PAM **&mPAMArr, bool *&useRandomConsArr, 
                     string *&userConsensusArr, double *&nMotifsConsArr,
                     double *&PctConsensusConserverdArr, char seqType,
                     bool *&isPalindromicArr);
      void printMotifConsensus(ostream &os);
      int getIndex(char aminoAcid);
      friend ostream &operator << (ostream &os, SequenceClass s);
      int seqLens;
      int numMotifs;
      int numSeqs;
      string *seqs;
      string *seqHdrs;
      string *consMotifArr;
      int **motifPosArr;
      bool **includeMotifArr;

   private:
      char *LabelArr;
      long double *overallFreqArr;
      long double **FreqArr;
      int *motifLenArr;
      double *pctConservedArr;   // Percent conservation of nts inside of single motif
      double *pctMotifsConsArr;  // Percentage of motifs conserved per sequence
      bool *palindromeArr;
      bool *randomConsArr;
      bool *userConsArr;
      string *userConsSeqArr;
      PAM *PAMArr;
      char seqType;
      int numResidues;
      bool destructorCalled;
};

#endif

//============================================================================
//------------------------ END OF SequenceClass.H --------------------------
//============================================================================


