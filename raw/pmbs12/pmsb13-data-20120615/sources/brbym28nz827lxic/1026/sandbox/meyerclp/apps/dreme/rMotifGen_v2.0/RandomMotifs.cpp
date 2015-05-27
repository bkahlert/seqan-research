
//============================================================================
//------------------------- BEGIN OF RandomMotifs.cpp ------------------------
//============================================================================

//==========================================================================
// rMotifGen
// Eric C. Rouchka
// C. Timothy Hardin
// (c) 2004-2007, University of Louisville
//
// FILE: RandomMotifs.cpp
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

#include "RandomMotifs.H"
#include "SequenceClass.H"
#include "CommonRoutines.H"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>
using namespace std;

//---------------------------------------------------------------------------
void printWelcome() {

   //**********************************************************************
   // PRE  : NONE
   // POST : Prints out a welcome screen information about rMotifGen
   //**********************************************************************

   system("clear");
   cout << "          ===========================================" << endl;
   cout << "          | rMotifGen (v2.0)                        |" << endl;
   cout << "          |                                         |" << endl;
   cout << "          | Eric C. Rouchka                         |" << endl;
   cout << "          | C. Timothy Hardin                       |" << endl;
   cout << "          |                                         |" << endl;
   cout << "          | (c) 2004-2007, Eric C. Rouchka and      |" << endl;
   cout << "          |     C. Timothy Hardin                   |" << endl;
   cout << "          |     University of Louisville            |" << endl;
   cout << "          |                                         |" << endl;
   cout << "          | rMotifGen comes with ABSOLUTELY NO      |" << endl;
   cout << "          | WARRANTY.  It is distruited under the   |" << endl;
   cout << "          | GNU General Public License.  Please see |" << endl;
   cout << "          | the LICENSE file for more details.      |" << endl;
   cout << "          |                                         |" << endl;
   cout << "          | Please cite:                            |" << endl;
   cout << "          |   Rouchka EC, Hardin CT.  (2006)        |" << endl;
   cout << "          |   rMotifGen: random motif generator for |" << endl;
   cout << "          |   genomic sequences.  (Under review).   |" << endl;
   cout << "          ===========================================" << endl;
   cout << endl << endl;
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
void printOverallHeader() {

   //**********************************************************************
   // PRE  :
   // POST :
   //**********************************************************************

   for(int j = 0; j < 76; j++) 
       cout << "-";
   cout << endl;
   for(int j = 0; j < 24; j++)
       cout << "-";
   cout << " OVERALL SEQUENCE SETTINGS ";
   for(int j = 0; j < 25; j++)
       cout << "-";
   cout << endl;
   for(int j = 0; j < 76; j++) 
       cout << "-";
   cout << endl;
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
void printMotifHeader(int currMotif) {

   //**********************************************************************
   // PRE  :
   // POST :
   //**********************************************************************

   cout << endl << endl;
   for(int j = 0; j < 76; j++) 
       cout << "=";
   cout << endl;
   for(int j = 0; j < 33; j++)
       cout << "=";
   cout << " MOTIF " << currMotif << "  ";
   for(int j = 0; j < 33; j++)
       cout << "=";
   cout << endl;
   for(int j = 0; j < 76; j++) 
       cout << "=";
   cout << endl;
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
void getOverallData(int &numSeqs, int &seqLen, double *&pctArr, 
                    char *LabelArr, int numResidues, int &numMotifs) {

   //**********************************************************************
   // PRE  :
   // POST :
   //**********************************************************************

   pctArr = new double[numResidues];
   numSeqs = getIntInput("Enter in the number of sequences to generate (10): ",
                         10, 0, 1001);
   seqLen  = getIntInput("Enter in the sequence length (500) : ", 500, 9, 10001);
   double CumSum = 0;
   for(int i = 0; i < numResidues; i++) {
      char *prompt;
      prompt = new char[500];
      sprintf(prompt, "Enter in the overall percent %c : ", LabelArr[i]);
      pctArr[i] = getDoubleInput(prompt, 25, -1, 101);
      CumSum += pctArr[i];
   }
   CumSum /= 100.0;
   for(int i = 0; i < numResidues; i++) {
      pctArr[i] /=  CumSum;
   }
   numMotifs = getIntInput("Enter in the number of motifs per sequence (3) : ",
                           3, -1, 11);
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
void allocateMotifArrays(bool *&motif_isRandomArr, int *&motif_lenArr, 
                         string *&motif_seqArr, double *&motif_pctSeqConsArr, 
                         PAM **&motif_PAMArr, double **&motifPctArr,
                         double *&motif_pctConsensusConsArr, int numMotifs,
                         int numResidues, bool *&motif_isPalindromicArr) {

   //**********************************************************************
   // PRE  :
   // POST :
   //**********************************************************************

   motif_isRandomArr         = new bool[numMotifs];
   motif_lenArr              = new int[numMotifs];
   motif_seqArr              = new string[numMotifs];
   motif_pctSeqConsArr       = new double[numMotifs];
   motif_PAMArr              = new PAM *[numMotifs];
   motif_pctConsensusConsArr = new double[numMotifs];
   motifPctArr               = new double *[numMotifs];
   motif_isPalindromicArr    = new bool[numMotifs];
   for(int i = 0; i < numMotifs; i++) {
      motifPctArr[i] = new double[numResidues];
   }
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
void getMotifData(bool *&motif_isPalindromicArr, 
                  bool *&motif_isRandomArr, 
                  int  *&motif_lenArr,
                  string *&motif_seqArr, 
                  double *&motif_pctSeqConsArr, 
                  double *&motif_pctConsensusConsArr, 
                  double **&motif_pctArr,
                  PAM **&motif_PAMArr, 
                  int numMotifs,
                  char seqType,
                  int numResidues,
                  char *LabelArr) {

   //**********************************************************************
   // PRE  : space has been allocated for the motif information
   // POST : valid data has been stored for each motif
   //**********************************************************************

   for(int i = 0; i < numMotifs; i++) {
      int currMotif = i + 1;
      printMotifHeader(currMotif);
      string tmp = " ";
      tmp = getCharInput("Is this motif (U)ser defined, or (Random)? (Please enter (U) or (R)) : ", "U", "R");
      
      if(tmp == "R") {
         //--------------------------------------------------------------------
         // RANDOM MOTIF -- must get the motif length and the motif composition
         //                 as well as if the motif is palindromic
         //--------------------------------------------------------------------

         motif_isRandomArr[i] = true;
         motif_lenArr[i] = getIntInput("Enter in the motif length (10) : ",
                                       10, 0, 101);
         if(seqType == 'N') {
            string tmp2 = " ";
            tmp2 = getCharInput("Is this motif palindromic (Y) or (N)? ", "Y", "N");
            if(tmp2 == "Y")  {  motif_isPalindromicArr[i] = true; }
            else             {  motif_isPalindromicArr[i] = false; }
         }
         double sum = 0.0;
         for(int j = 0; j < numResidues; j++) {
            char *prompt;
            prompt = new char[500];
            sprintf(prompt, "Enter in the motif percent %c : ", LabelArr[j]);
            double tmpdouble = getDoubleInput(string(prompt), 100.0 / numResidues, -1, 101);
            sum += tmpdouble;
            motif_pctArr[i][j] = tmpdouble;
         }
         for(int j = 0; j < numResidues; j++) {
            motif_pctArr[i][j] /= sum;
            motif_pctArr[i][j] *= 100.0;
         }
      }
      else {
         // USER INPUT MOTIF -- RETRIEVE THE CONSENSUS FROM THE USER
         motif_isRandomArr[i] = false;
         motif_seqArr[i] = getSequence("Enter in the motif sequence : ", seqType);
      }
      string prompt = "Enter in the percentage of sequences containing the motif : ";
      motif_pctSeqConsArr[i] = getIntInput(prompt, 100, -1, 101);
      if(seqType == 'N') {
         prompt = "Enter in the percent conservation of the motif consensus sequence : ";
         motif_pctConsensusConsArr[i] = getIntInput(prompt, 100, -1, 101);
      }
      else { // MUST BE AMINO ACID SEQUENCE
         prompt = "Enter in the number for the PAM matrix to be used (eg 250 for PAM250) : ";
         int PAMnum = getIntInput(prompt, 250, -1, 501);
         motif_PAMArr[i] = new PAM(PAMnum);
      }
   }
}

//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
void printToFile(SequenceClass *&S) {

   //**********************************************************************
   // PRE  : S is the sequence class that has been created
   // POST : retrieves a file name from the user and then 
   //        prints out the sequence class to that file
   //**********************************************************************

   ofstream outfile;
   string outfn;

   bool done = false;
   while(!done) {
      cout << "Enter the file name to print to: ";
      cin >> outfn;

      // CHECK TO SEE IF FILE EXISTS
      if(fileExists(outfn)) {
         string tmp = getCharInput("File exists: Overwrite? (Y or N)", "Y", "N");
         if(tmp == "Y") {
            outfile.open(outfn.c_str(), ios::trunc);
            if(outfile) {
               outfile << (*S) << endl;
               done = true;
               outfile.close();
            }
            else {
               cerr << "ERROR ENCOUNTERED: CANNOT WRITE TO " << outfn << endl;
            } 
         }
      }
      else {     
         outfile.open(outfn.c_str(), ios::out);

         if(outfile) {
             outfile << (*S) << endl;
             done = true;
             outfile.close();
         } 
         else {
            cerr << "ERROR ENCOUNTERED: CANNOT WRITE TO " << outfn << endl;
         }
      }
   }
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
int setResidues(char *&ResidueLabels, double *&ResiduePcts,
                char seqType) {

   //**********************************************************************
   // PRE  : seqType is 'N' for DNA, 'P' for proteins
   // POST : allocates the residue labels and percentages arrays
   //**********************************************************************

   if(seqType == 'N') {
      ResidueLabels = new char[4];
      ResiduePcts   = new double[4];
      ResidueLabels[0] = 'A';   // ADENINE
      ResidueLabels[1] = 'C';   // CYTOSINE
      ResidueLabels[2] = 'G';   // GUANINE
      ResidueLabels[3] = 'T';   // THYMINE
      return(4);
   }
   else {
      ResidueLabels = new char[20];
      ResiduePcts   = new double[20];
      ResidueLabels[0] = 'A';  // ALANINE
      ResidueLabels[1] = 'R';  // ARGININE
      ResidueLabels[2] = 'N';  // ASN
      ResidueLabels[3] = 'D';  // ASPARTAMINE
      ResidueLabels[4] = 'C'; 
      ResidueLabels[5] = 'Q';
      ResidueLabels[6] = 'E';
      ResidueLabels[7] = 'G';
      ResidueLabels[8] = 'H';
      ResidueLabels[9] = 'I';
      ResidueLabels[10] = 'L';
      ResidueLabels[11] = 'K';
      ResidueLabels[12] = 'M';
      ResidueLabels[13] = 'F';
      ResidueLabels[14] = 'P';
      ResidueLabels[15] = 'S';
      ResidueLabels[16] = 'T';
      ResidueLabels[17] = 'W';
      ResidueLabels[18] = 'Y';
      ResidueLabels[19] = 'V';
      return(20);
   }
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
void getOverallBackgroundResiduePct(ifstream &infile, char *ResidueLabels, 
                                    double *&ResiduePcts, int numResidues) {

   //**********************************************************************
   // PRE  : infile is an open input file stream 
   // POST : reads in the overall background residue frequencies from infile
   //**********************************************************************
   if(infile.eof()) { exit(26); }
   double cumSum = 0.0;
   for(int i = 0; i < numResidues; i++) {
      char *tmpch;
      tmpch = new char[1];
      sprintf(tmpch, "%c", ResidueLabels[i]);
      string tmpst(tmpch); 
      string ResidueLabel = "OAPct" + tmpst;
      double tmpdouble;
      if(infile.eof()) { exit(26); }
      getValidDoubleTokens(tmpdouble, infile, -1.0, 100.0, ResidueLabel, 11);
      ResiduePcts[i] = tmpdouble;
      cumSum += ResiduePcts[i];
      delete [] (tmpch);
   }
  
   // Normalize Residue Percentages to sum to 100
   cumSum /= 100.0;
   for(int i = 0; i < numResidues; i++) {
      ResiduePcts[i] /= cumSum;
   }
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
void getUserDefinedMotifInfo(ifstream &infile, bool *&motif_isRandomArr, 
                           string *&motif_seqArr, double *&motif_pctSeqConsArr,
                           double *&motif_pctConsensusConsArr, 
                           PAM **&motif_PAMArr, char seqType, int i) {

   //**********************************************************************
   // PRE  : infile is an open input file stream 
   // POST : reads in the overall background residue frequencies from infile
   //**********************************************************************

   // Need to update this for differences between protein and dna
   // sequences
   
   motif_isRandomArr[i] = false;
   string motifLabel = setMotifLabel(i);
   string currLabel = motifLabel + "Seq";
   string tmpstr;
   if(infile.eof()) { exit(26); }
   if(seqType == 'P') { // PROTEIN SEQUENCE
      getValidAASequenceTokens(tmpstr, infile, currLabel, 15);
   }
   else { // DNA SEQUENCE
      getValidDNASequenceTokens(tmpstr, infile, currLabel, 15);
   }
   motif_seqArr[i] = tmpstr;
   currLabel = motifLabel + "PctSeqsContain";
   double tmpdouble;
   if(infile.eof()) { exit(26); }
   getValidDoubleTokens(tmpdouble, infile, -1, 101, currLabel, 16);
   motif_pctSeqConsArr[i] = tmpdouble;
   if(infile.eof()) { exit(26); }
   if(seqType == 'P') { // PROTEIN SEQUENCE
      getPAMValue(motif_PAMArr, infile, i, motifLabel);
   }
   else { // DNA SEQUENCE
      getPctConserved(motif_pctConsensusConsArr, infile, i, motifLabel);
   }
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
void getPctConserved(double *&motif_pctConsensusConsArr, ifstream &infile, 
                 int i, string motifLabel) {

   //**********************************************************************
   // PRE  : infile is open and motif_PAMArr has been allocated
   // POST : reads in the PAM value from infile, allocates the array
   //      : and sets it accordingly
   //**********************************************************************

   string currLabel = motifLabel + "PctCons";
   double tmpdouble;
   getValidDoubleTokens(tmpdouble, infile, -1, 101, currLabel, 17);
   motif_pctConsensusConsArr[i] = tmpdouble;
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
string setMotifLabel(int i) {

   //**********************************************************************
   // PRE  : i is the motif number
   // POST :  returns a string giving motif followed by i
   //**********************************************************************

   char *tmpNum;
   tmpNum = new char[2];
   sprintf(tmpNum, "%i", i+1);
   string tmpNum2(tmpNum);
   delete [] tmpNum;
   string motifLabel = "motif" + tmpNum2;
   return(motifLabel);
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
void getRandomMotifInfo(ifstream &infile, bool *&motif_isRandomArr, 
                        string *&motif_seqArr, double *&motif_pctSeqConsArr,
                        double *&motif_pctConsensusConsArr,
                        PAM **&motif_PAMArr, bool *&motif_isPalindromicArr,
                        int *&motif_lenArr, char seqType, int i, int numResidues,
 char *ResidueLabels, double **&motifPctArr) {


    motif_isRandomArr[i] = true;
    string motifLabel = setMotifLabel(i);
    string currLabel = motifLabel + "Len";

    int tmpint;
    if(infile.eof()) { exit(26); }
    getValidIntTokens(tmpint, infile, 0, 101, currLabel, 18);
    motif_lenArr[i] = tmpint;
    
    if(seqType == 'N') { // DNA SEQUENCE
       currLabel = motifLabel + "IsPalindrome";
       vector <string>validChVector;
       validChVector.clear();
       validChVector.push_back("Y");
       validChVector.push_back("N");
       string isPalindrome;
       if(infile.eof()) { exit(26); }
       getValidCharTokens(isPalindrome, infile, validChVector, currLabel, 19);
       if(isPalindrome == "Y") {
          motif_isPalindromicArr[i] = true;
       }
       else {
          motif_isPalindromicArr[i] = false;
       }
    }
    getResidueFrequencies(infile, motifPctArr, ResidueLabels, motifLabel, i, numResidues);

    currLabel = motifLabel + "PctSeqsContain";
    if(infile.eof()) { exit(26); }
    getValidIntTokens(tmpint, infile, -1, 101, currLabel, 16);
    motif_pctSeqConsArr[i] = tmpint;
  
    if(infile.eof()) { exit(26); }
    if(seqType == 'P') { // PROTEIN SEQUENCE
        getPAMValue(motif_PAMArr, infile, i, motifLabel);
    }
    else { // NUCLEOTIDE SEQUENCE
      getPctConserved(motif_pctConsensusConsArr, infile, i, motifLabel);
    }
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
void getResidueFrequencies(ifstream &inFile, double **&motifPctArr, 
                           char *ResidueLabels, string motifLabel, int i,
                           int numResidues) {

   //**********************************************************************
   // PRE  : infile is open and motif_PctArr has been allocated
   // POST : reads in the residue frequencies from the infile, and sets
   //      : the values accordingly 
   //**********************************************************************
    double cumSum = 0.0;
    for(int j = 0; j < numResidues; j++) {
       string Label = motifLabel + "Pct" + ResidueLabels[j];
       double tmpdouble;
       if(inFile.eof()) { exit(26); }
       getValidDoubleTokens(tmpdouble, inFile, -1.0, 100.0, Label, 11);
       motifPctArr[i][j] = tmpdouble;
       cumSum += motifPctArr[i][j];
    }
    cumSum /= 100.0;
    for(int j = 0; j < numResidues; j++) {
       motifPctArr[i][j] /= cumSum;
    }
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
void getPAMValue(PAM **&motif_PAMArr, ifstream &infile, 
                 int i, string motifLabel) {

   //**********************************************************************
   // PRE  : infile is open and motif_PAMArr has been allocated
   // POST : reads in the PAM value from infile, allocates the array
   //      : and sets it accordingly
   //**********************************************************************

    string currLabel = motifLabel + "PAM";
    int pamVal;
    if(infile.eof()) { exit(26); }
    getValidIntTokens(pamVal, infile, -1, 501, currLabel, 17);
    motif_PAMArr[i] = new PAM(pamVal);
}
//---------------------------------------------------------------------------

 
//---------------------------------------------------------------------------
void fileMode(string FN, char seqType){

   //**********************************************************************
   // PRE  : FN is a file name that exists
   // POST : reads in the motif generation information from a file with
   //      : a format of label=value
   //      : checking to make sure each is appropriately listed
   //**********************************************************************


   double *ResiduePcts;
   char   *ResidueLabels;
   int seqLen;
   int numMotifs;
   int numSeqs;
   bool finished;
   bool   *motif_isRandomArr;
   bool   *motif_isPalindromicArr;
   int    *motif_lenArr;
   string *motif_seqArr;
   double *motif_pctSeqConsArr;
   double *motif_pctConsensusConsArr;
   PAM    **motif_PAMArr;
   double **motifPctArr;
   vector <string>validChVector;
   string currLabel;
   
   int numResidues = setResidues(ResidueLabels, ResiduePcts, seqType);
   if(!fileExists(FN)) {
      exit(7);
   }
   ifstream infile;
   infile.open(FN.c_str(), ios::in);
   if(!infile) {
      exit(7);
   }
   string currLine;
   char *lhTok, *rhTok;
   getValidIntTokens(numSeqs, infile, 0, 1001, "numSeq", 8);
   getValidIntTokens(seqLen, infile, 9, 10001, "seqLen", 10);
   getOverallBackgroundResiduePct(infile, ResidueLabels, 
                                  ResiduePcts, numResidues);
   getValidIntTokens(numMotifs, infile, -1, 11, "numMotifs", 13);
   allocateMotifArrays(motif_isRandomArr, motif_lenArr, motif_seqArr,
                       motif_pctSeqConsArr, motif_PAMArr, motifPctArr,
                       motif_pctConsensusConsArr, numMotifs, numResidues,
                       motif_isPalindromicArr);

   for(int i = 0; i < numMotifs; i++) {
      // RETRIEVE THE MOTIF-SPECIFIC INFORMATION
      string motifLabel = setMotifLabel(i);
      currLabel = motifLabel +  "Type";
      string motiftype;
      validChVector.clear();
      validChVector.push_back("U");   
      validChVector.push_back("R");
 
      getValidCharTokens(motiftype, infile, validChVector, currLabel, 14);
      if(motiftype == "U") {
          getUserDefinedMotifInfo(infile, motif_isRandomArr, motif_seqArr, 
                              motif_pctSeqConsArr, motif_pctConsensusConsArr,
                              motif_PAMArr, seqType, i);
      }
      else {
          getRandomMotifInfo(infile, motif_isRandomArr, motif_seqArr,
                             motif_pctSeqConsArr, motif_pctConsensusConsArr,
                             motif_PAMArr, motif_isPalindromicArr, 
                             motif_lenArr, seqType, i, numResidues, 
                             ResidueLabels, motifPctArr);
      }
   }
   validChVector.clear();
   validChVector.push_back("F");
   validChVector.push_back("S");
   string outputType;
   currLabel="output";
   if(infile.eof()) { exit(26); }
   getValidCharTokens(outputType, infile, validChVector, currLabel, 22);
   ofstream outfile;
   if(outputType == "F") {
       string outfn;
       currLabel = "outFile";
       if(infile.eof()) { exit(26); }
       getStringTokens(outfn, infile, currLabel, 23);
       outfile.open(outfn.c_str(), ios::out);
       if(!outfile) {
          exit(24);
       }
   }
   infile.close();
   SequenceClass *S;
   S = new SequenceClass(seqLen, numMotifs, numSeqs, ResidueLabels, 
                           ResiduePcts, seqType, numResidues);
   S->setValues(motifPctArr, motif_lenArr, motif_PAMArr, motif_isRandomArr,
                motif_seqArr, motif_pctSeqConsArr, motif_pctConsensusConsArr, seqType, motif_isPalindromicArr);

   if(outputType == "F") {
      S->printMotifConsensus(outfile);
      outfile << *S << endl;
      outfile.close();
   }
   else {
      S->printMotifConsensus(cout);
      cout << *S << endl;
   }
   delete [] ResiduePcts;
   delete [] ResidueLabels;
   delete [] motif_isRandomArr; 
   delete [] motif_lenArr;
   delete [] motif_seqArr;
   delete [] motif_pctSeqConsArr;
   for(int i = 0; i < numMotifs; i++) {
      if(seqType == 'P') 
         delete [] motif_PAMArr[i];
      if(seqType == 'N')
         delete [] motifPctArr[i];
   }
   if(motif_PAMArr) {
      delete [] motif_PAMArr;
   }
   delete [] motifPctArr;
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
void interactiveMode(){

   //**********************************************************************
   // PRE  : NONE
   // POST : reads in the motif generation information from the keyboard
   //**********************************************************************

   int    seqLen;
   int    numMotifs;
   int    numSeqs;
   bool   finished;
   double *ResiduePcts;
   char   *ResidueLabels;
   PAM    **motif_PAMArr;
   double **motifPctArr;
   bool   *motif_isPalindromicArr;
   bool   *motif_isRandomArr;
   int    *motif_lenArr;
   string *motif_seqArr;
   double *motif_pctSeqConsArr;
   double *motif_pctConsensusConsArr;
   double *residueBackgroundFrequencies;
   char seqType;

   printWelcome();
   printOverallHeader();
   string tmp = getCharInput("Type of sequence to generate (N)ucleotide or (P)rotein : ", 
                             "N", "P");
   if(tmp == "P") { seqType = 'P'; }
   if(tmp == "N") { seqType = 'N'; }

   int numResidues = setResidues(ResidueLabels, ResiduePcts, seqType);

   getOverallData(numSeqs, seqLen, ResiduePcts, ResidueLabels, numResidues, numMotifs);

   allocateMotifArrays(motif_isRandomArr, motif_lenArr, motif_seqArr,
                       motif_pctSeqConsArr, motif_PAMArr, motifPctArr,
                       motif_pctConsensusConsArr, numMotifs, numResidues,
                       motif_isPalindromicArr);

   getMotifData(motif_isPalindromicArr, motif_isRandomArr, motif_lenArr,
                motif_seqArr, motif_pctSeqConsArr, 
                motif_pctConsensusConsArr, 
                motifPctArr,
                motif_PAMArr, numMotifs, seqType, numResidues, ResidueLabels);

   SequenceClass *S;
   S = new SequenceClass(seqLen, numMotifs, numSeqs, ResidueLabels, 
                           ResiduePcts, seqType, numResidues);
   S->setValues(motifPctArr, motif_lenArr, motif_PAMArr, motif_isRandomArr,
                motif_seqArr, motif_pctSeqConsArr, motif_pctConsensusConsArr, 
                seqType, motif_isPalindromicArr);

   S->printMotifConsensus(cout);
   tmp = getCharInput("Output to (S)creen or (F)ile : ", "S", "F");
   if(tmp == "S") 
      cout << *S << endl;
   else {
      printToFile(S);
   }
   delete [] ResiduePcts;
   delete [] ResidueLabels;
   delete [] motif_isPalindromicArr;
   delete [] motif_isRandomArr;
   delete [] motif_lenArr;
   delete [] motif_seqArr;
   delete [] motif_pctSeqConsArr;
   delete [] motif_pctConsensusConsArr;
   delete [] residueBackgroundFrequencies;
   delete [] motif_PAMArr;
   for(int i = 0; i < numMotifs; i++) 
      delete [] motifPctArr[i];
   delete [] motifPctArr;
   delete S;
}
//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
int main(int argc, char *argv[]) {

   if(argc == 1) {
      interactiveMode();
   }
   else if(argc != 4) {
      exit(5);
   }
   else {
      char residueType='N';
      string arg1(argv[1]);
      if(!((arg1 == "N") || (arg1 == "P"))) 
         exit(5);
      if(arg1 == "P")
         residueType = 'P';

      if(argc == 4) {
         string arg2(argv[2]);
         if(arg2 != "-f") {
            exit(5);
         }
         string inFN(argv[3]);
         fileMode(inFN, residueType);
      }
   }
   return(EXIT_SUCCESS);
}
//---------------------------------------------------------------------------

//============================================================================
//--------------------------- END OF RandomMotifs.cpp ------------------------
//============================================================================



