

//============================================================================
//------------------------ BEGIN OF SequenceClass.cpp ------------------------
//============================================================================

//==========================================================================
// rMotifGen
// Eric C. Rouchka
// C. Timothy Hardin
// (c) 2004-2007, University of Louisville
//
// FILE:  SequenceClass.cpp
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

#include "SubstitutionMatrices.H"
#include "SequenceClass.H"
#include <iostream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <time.h>
using namespace std;

//-----------------------------------------------------------------------------
ostream &operator << (ostream &os, SequenceClass s) {

   //**************************************************************************
   // PRE:  os is a defined output stream (console, error, or file)
   //       s  is a valid instance of the sequence class
   // POST: the returned value will print out the header and sequence
   // DESCRIPTION: This is the overloaded output operator for the 
   //              SequenceClass
   //**************************************************************************

   for(int i = 0; i < s.numSeqs; i++) {
      os << s.seqHdrs[i] << endl;
      string currSeq = s.formatSeq(s.seqs[i], 60, 60, false);
      os << currSeq;
   }
   return(os);
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
void SequenceClass::setValues(double **&mFreqArr, int *&mLenArr,
                                PAM **&mPAMArr, bool *&useRandomConsArr,
                                string *&userConsensusArr, double *&nMotifsConsArr, double *&pctConsConsArr, char seqType, bool *&isPalindromicArr) {

   //**************************************************************************
   // PRE:  NONE
   // POST: NONE
   // DESCRIPTION: Sets the values for the data members of the SequenceClass
   //**************************************************************************
   for(int i = 0; i < numMotifs; i++) 
   {
      for(int j = 0; j < numResidues; j++) {
         FreqArr[i][j] = (long double)mFreqArr[i][j];
      }
      motifLenArr[i] = mLenArr[i];
      if(seqType == 'P') {
      PAMArr[i]      = *(mPAMArr[i]);
      }
      palindromeArr[i] = isPalindromicArr[i]; 

      randomConsArr[i] = useRandomConsArr[i];
      userConsArr[i] = !(useRandomConsArr[i]);
      userConsSeqArr[i] = userConsensusArr[i];
      pctMotifsConsArr[i] = nMotifsConsArr[i];
      pctConservedArr[i] = pctConsConsArr[i];
   }
   generateSequences();
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
void SequenceClass::printMotifConsensus(ostream &os) {

   //**************************************************************************
   // PRE:  NONE
   // POST: NONE
   // DESCRIPTION: Prints the current motif sequence for each motif
   //              To the specified output stream
   //**************************************************************************

   os << endl;
   os << "NUMMOTIFS: " << numMotifs << endl;
   for(int i = 0; i < numMotifs; i++) {
      int currMotif = i + 1;
      os << "Motif " << currMotif << " consensus: " << consMotifArr[i] << endl;
   }
   os << endl;
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
string SequenceClass::randomSeq(int len, char *labelArr, long double *FreqArr, bool isPalindrome)
{
   //**************************************************************************
   // PRE:  NONE
   // POST: NONE
   // DESCRIPTION:
   //**************************************************************************
   string s = "";
   int i = 0;

   if(isPalindrome) {
      int halflen = len / 2;
      string post = "";
      while(i < halflen) {
         int val = lrand48() % 101;
         long double cumSum = 0;
         for(int j = 0; j < numResidues; j++) {
            cumSum += FreqArr[j];
            if(val < cumSum) {
               s += labelArr[j];
               if(labelArr[j] == 'A') { post = "T" + post; }
               if(labelArr[j] == 'C') { post = "G" + post; }
               if(labelArr[j] == 'G') { post = "C" + post; }
               if(labelArr[j] == 'T') { post = "A" + post; }
               j = numResidues;
            }
         }
         i++;
      }
      if((len % 2) != 0) {
         int val = lrand48() % 101;
         long double cumSum = 0;
         for(int j = 0; j < numResidues; j++) {
            cumSum += FreqArr[j];
            if(val < cumSum) {
               s += labelArr[j];
               j = numResidues;
            }
         }
      }
      s += post;
   }
   else {
      while (i < len)
      {
         int val = lrand48() % 101;
         long double cumSum = 0;
         for(int j = 0; j < numResidues; j++) {
            cumSum += FreqArr[j];
            if(val < cumSum) {
               s += labelArr[j];
               j = numResidues;
            }
         }
         i++;
      }
   }
   return(s);
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
void SequenceClass::createRandomSequences() 
{
   //**************************************************************************
   // PRE:  NONE
   // POST: NONE
   // DESCRIPTION: creates random sequences
   //**************************************************************************
   string currSeq;
   for(int i = 0; i < numSeqs; i++) 
   {
      currSeq = "";
      currSeq = randomSeq(seqLens, LabelArr, overallFreqArr, false);
      seqs[i] = currSeq;
   }
}
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
void SequenceClass::createConsensusMotifs() 
{
   //**************************************************************************
   // PRE:  NONE
   // POST: NONE
   // DESCRIPTION: Creates consensus motif patterns
   //**************************************************************************

   // NOW CREATE RANDOM MOTIFS
   string currSeq;
	
   for(int i = 0; i < numMotifs; i++) 
   {
      if(randomConsArr[i]) {
         currSeq = "";
         currSeq = randomSeq(motifLenArr[i], LabelArr, FreqArr[i], palindromeArr[i]);
         consMotifArr[i] = currSeq;
      }
      else {
         consMotifArr[i] = userConsSeqArr[i];
         motifLenArr[i] = consMotifArr[i].length();
      }
   }
}

//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
bool SequenceClass::isValid(int pos, bool *used, int len)
{
   //**************************************************************************
   // PRE:  NONE
   // POST: NONE
   // DESCRIPTION: This is the default constructor for the sequence class
   //**************************************************************************
   for (int i = pos; i < (pos + len); i++)
   {
       if (used[i])
          return (false);   // Position is occupied
   }
       return (true);            // Position is unoccupied
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
void SequenceClass::generateSequences() 
{
   //**************************************************************************
   // PRE:  NONE
   // POST: NONE
   // DESCRIPTION: 
   //**************************************************************************

// FIRST, CREATE RANDOM SEQUENCES
   MotifsToInclude();        // Set which motifs are included 
   createRandomSequences();  
   createConsensusMotifs();
   bool *usedArr;

   for(int i = 0; i < numSeqs; i++) 
   {
      // Initialize all positions to be valid motif start positions
      usedArr = new bool[seqLens];
      for (int j = 0; j < seqLens; j++)
      {
         usedArr[j] = false;
      }

      for(int j = 0; j < numMotifs; j++) 
      {
         if (includeMotifArr[i][j])
         {
            int newLen = seqLens - motifLenArr[j];
            int currLoc = (int)lrand48() % newLen;
            while (!isValid(currLoc, usedArr, motifLenArr[j]))
            {
               currLoc = (int)lrand48() % newLen;
            }
            motifPosArr[i][j] = currLoc;
            for (int k = currLoc; k < (currLoc + motifLenArr[j]); k++)
            {
               usedArr[k] = true;
            }
         }
         else
         {
            motifPosArr[i][j] = -1;
         }
      }

      for(int j = 0; j < numMotifs; j++) 
      {
         if (motifPosArr[i][j] != -1)
         {
            string currMotif;
            if(seqType == 'P') {
               currMotif = mutate(consMotifArr[j], PAMArr[j]);
            }
            else {
               currMotif = mutate(consMotifArr[j], pctConservedArr[j]);
            }
            seqs[i] = seqs[i].replace(motifPosArr[i][j], motifLenArr[j], currMotif);
         }
      }
   }

   // ---------------------------------
   // Create the random sequence header
   // ---------------------------------

   for(int i = 0; i < numSeqs; i++) 
   {
      string currhdr = "";
      int currSeq = i + 1;
      string mName = "rMotifGen_RandSeq_" + itos(currSeq);
      currhdr = ">" + mName + " " + itos(numMotifs) + " ";
      for(int j = 0; j < numMotifs; j++)
         currhdr += itos(motifPosArr[i][j]) + " ";
      seqHdrs[i] = currhdr;
   }

   delete [] usedArr;
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
string itos(int i) {
   stringstream s;
   s << i;
   return(s.str());
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
string SequenceClass::formatSeq(string s, int width, int blockWidth, 
                        bool numbersOn) 
{
   //**************************************************************************
   // PRE:  NONE
   // POST: NONE
   // DESCRIPTION: 
   //**************************************************************************
   if(width % blockWidth != 0) 
      return(s);

   string formatSeq = "";
   int len = (int)s.length();
   int numLines = (int)(len/width);
   int remainder = len % width;
   string currSegment = "";
   int begNum = 0;
   int endNum = 0;

   int numDigits = calculateDigits(len);

   for(int currLine = 0; currLine < numLines; currLine++) 
   {
      currSegment = "";
      if(numbersOn) 
      {
         begNum = 1 + width * currLine;
         endNum = begNum + width - 1;
         int currDigits = calculateDigits(begNum);
         for(int i = 0; i < (numDigits - currDigits); i++) 
         {
            currSegment += " ";
         }
         currSegment += begNum + " ";
      }
      int numBlocks = (width / blockWidth);
      for(int currBlock = 0; currBlock < numBlocks; currBlock++) 
      {
         currSegment += s.substr(currLine * width + currBlock * blockWidth, blockWidth);
         currSegment += " ";
      }
      if(numbersOn) 
      {
         currSegment += " ";
         currSegment += endNum;
      }
      formatSeq += currSegment + "\n";

   }
   if(remainder > 0) 
   {
      int rem = remainder;
      string remainingSeq = s.substr(numLines * width, remainder);
      len = remainingSeq.length();
      int numBlocks = (int)(len / blockWidth);
      remainder = len % blockWidth;
   
      if(numbersOn) 
      {
         begNum = 1 + width * numLines;
         endNum = (int)s.length();
         int currDigits = calculateDigits(begNum);
         for(int i = 0; i < numDigits - currDigits; i++) 
         {
            formatSeq += " ";
         }
         formatSeq += begNum + " ";
      }
      for(int currBlock = 0; currBlock < numBlocks; currBlock++) 
      {
         formatSeq += remainingSeq.substr(currBlock * blockWidth, blockWidth);
         formatSeq += " ";
      }
      if(remainder != 0) 
      {
         formatSeq += remainingSeq.substr(numBlocks * blockWidth, remainder);
      }
      if(numbersOn) 
      {
         int blocks = (int)(width / blockWidth) - numBlocks;
         int numSpaces = width - rem + blocks + 1;
         for(int i = 0; i < numSpaces; i++) 
         {
            formatSeq += " ";
         }
         formatSeq += endNum;
      }
      formatSeq += "\n";
   }
   return(formatSeq);
}

//-----------------------------------------------------------------------------
		
//-----------------------------------------------------------------------------
int SequenceClass::calculateDigits(int len) 
{
   //**************************************************************************
   // PRE:  NONE
   // POST: NONE
   // DESCRIPTION: returns the number of digits in the number len
   //**************************************************************************
   int num = len;
   int n = 0;
   while((int)(num/10) != 0) 
   {
      n++;
      num /= 10;
   }
   return(n);
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
void SequenceClass::MotifsToInclude() {

   //**************************************************************************
   // PRE: numSeqs has been set (>0) as has numMotifs
   // POST: includeMotifArr set appropriately
   // DESCRIPTION: This function is used to determine which motifs will be
   //              included in the sequences based upon the percentage of
   //              motifs that are conserved
   //**************************************************************************

   for(int i = 0; i < numSeqs; i++) {
       for(int j = 0; j < numMotifs; j++) {
           int randVal = lrand48() % 101;
           if(randVal > pctMotifsConsArr[j]) {
               includeMotifArr[i][j] = false;
           }
           else {
               includeMotifArr[i][j] = true;
           }
       }
   }
   return;
}
//-----------------------------------------------------------------------------
        
//-----------------------------------------------------------------------------
int SequenceClass::getSeqLens()     { return(seqLens); }
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
double SequenceClass::getOverallPct(char residueToFind) {
   int index = getIndex(residueToFind);
   if(index == -1) 
      return(-1);
   else
      return(overallFreqArr[index]);
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
int SequenceClass::getIndex(char residueToFind) {
   for(int i = 0; i < numResidues; i++) {
      if(LabelArr[i] == residueToFind) {
         return(i);
      }
   }
   return(-1);
}
//-----------------------------------------------------------------------------

int SequenceClass::getNumMotifs()   { return(numMotifs); }
int SequenceClass::getNumSeqs()     { return(numSeqs); }
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
SequenceClass::~SequenceClass() {

   //**************************************************************************
   //PRE:         NONE
   //POST:        Memory dynamically allocated for Sequence class is free
   //DESCRIPTION: Destructor for SequenceClass
   //**************************************************************************

   if(!destructorCalled) {
      delete [] LabelArr;
      delete [] overallFreqArr;
      delete [] motifLenArr;
      if(PAMArr)
         delete [] PAMArr;
      delete [] pctMotifsConsArr;
      delete [] randomConsArr;
      delete [] userConsArr;
      delete [] userConsSeqArr;
      delete [] palindromeArr;
      for(int i = 0; i < numMotifs; i++) {
         delete[] FreqArr[i];
      }
      delete[] FreqArr;
      for(int i = 0; i < numSeqs; i++) {
         delete [] motifPosArr[i];
         delete [] includeMotifArr[i];
      } 
      delete [] motifPosArr;
      delete [] includeMotifArr;
      destructorCalled = true;
   }
}
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
SequenceClass::SequenceClass(int s, int nm, int ns, char *labelArr,
                                 double *tmpFreqArr, char st, int nr)
{
   //**************************************************************************
   //PRE: s is the length of the sequences to produce (>10)
   //     tmpFreqArr is the amino acid frequency array which must sum to 100
   //     nm is the number of motifs (>0)
   //     ns is the number of sequences (>0))
   //POST:
   //DESCRIPTION: This is a constructor for the sequence class
   //**************************************************************************

   seqLens = s;
   numMotifs = nm;
   numSeqs = ns;
   seqType = st;
   numResidues = nr;
   destructorCalled = false;

   //------------------------------------------------------------
   // ALLOCATE SPACE FOR THE PRIVATE DATA MEMBERS FOR THE CLASS
   // Values should be set in the setValues function
   //------------------------------------------------------------

   PAMArr           = new PAM[numMotifs];    // Substitution matrix array
   randomConsArr    = new bool[numMotifs];   // Random consensus motif?
   userConsArr      = new bool[numMotifs];   // User consensus motif?
   userConsSeqArr   = new string[numMotifs]; // User consensus motif sequence
   motifLenArr      = new int[numMotifs];    // length of motifs
   pctMotifsConsArr = new double[numMotifs];    // Pct motifs conserved/sequence
   LabelArr       = new char[numResidues]; // Residue labels
   overallFreqArr = new long double[numResidues];   // Overall residue frequencies
   pctConservedArr  = new double[numMotifs]; // Pct conservation inside motif
   palindromeArr    = new bool[numMotifs];  // Palindromic motif?

   for(int i = 0; i < numResidues; i++) {
      LabelArr[i] = labelArr[i];
      overallFreqArr[i] = (long double)tmpFreqArr[i];
   }

   //------------------------------------------------------------
   // ALLOCATE SPACE FOR THE PUBLIC DATA MEMBERS FOR THE CLASS
   // Values will be assigned in the generateSequences Function
   //------------------------------------------------------------
   seqs             = new string[numSeqs];   // Overall Sequences
   seqHdrs          = new string[numSeqs];   // Sequence Headers
   consMotifArr     = new string[numMotifs]; // Consensus Motif
   motifPosArr      = new int *[numSeqs];    // Motif Positions
   includeMotifArr  = new bool *[numSeqs];   // Is motif included in sequence

   for (int i = 0; i < numSeqs; i++)
   {
      motifPosArr[i] = new int[numMotifs];
      includeMotifArr[i] = new bool[numMotifs];
   }
   FreqArr   = new long double *[numMotifs];  
   for(int i = 0; i < numMotifs; i++) {
      FreqArr[i] = new long double[numResidues];
   }
   srand48(time(NULL));
}

//============================================================================
//------------------------ END OF SequenceClass.cpp ------------------------
//============================================================================



