//===========================================================================
//=====================  BEGIN OF CommonRoutines.cpp ========================
//===========================================================================


//==========================================================================
// rMotifGen
// Eric C. Rouchka
// C. Timothy Hardin
// (c) 2004-2007, University of Louisville
//
// FILE: CommonRoutines.cpp
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

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include "SequenceClass.H"
#include "CommonRoutines.H"
using namespace std;


//---------------------------------------------------------------------------
bool isNum(char ch) {

   //**********************************************************************
   // PRE : ch is a character
   // POST: returns true if ch is a numeral; false otherwise
   //**********************************************************************

   if((ch == '1') || (ch == '2') || (ch == '3') || (ch == '4') ||
      (ch == '5') || (ch == '6') || (ch == '7') || (ch == '8') ||
      (ch == '9') || (ch == '0')) {
      return(true);
   }
   return(false);
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
bool isNumOrSign(char ch) {

   //**********************************************************************
   // PRE : ch is a character
   // POST: returns true if ch is a numeral, +, or -; false otherwise
   //**********************************************************************

   if(isNum(ch) || (ch == '+') || (ch == '-')) {
      return(true);
   }
   return(false);
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
bool getDoubleFromString(double &retVal, string s, 
                         double lowBounds, double hiBounds) {

   //**********************************************************************
   // PRE  : defaultValue has been set
   // POST : returns true if a double has been read in -- the
   //        resulting value is stored in the parameter retVal
   //**********************************************************************
   char *lhTok, *rhTok;

   lhTok = strtok((char *)s.c_str(), ".");
   rhTok = strtok(NULL, "\n");
   
   int lhs, rhs;
   if(rhTok == NULL) {
      if(getIntFromString(lhs, lhTok, -1, 101)) {
         retVal = (double) lhs;
         return(true);
      }
      else {
         return(false);
      }
   }

   if(getIntFromString(lhs, lhTok, -1, 101) && 
      getIntFromString(rhs, rhTok, -1, 10001)) {

      // rhs should be divided accordingly
      int divisor = 1;
      int tmp = rhs;
      while(tmp > 0) {
         divisor *= 10;
         tmp = tmp / 10;
      }
      double rhs2 = (double)rhs / divisor;

      retVal = lhs + rhs2;
      return(true);
   }
   else {
      return(false);
   }
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
bool getIntFromString(int &retVal, string s, int lowBounds, int hiBounds) {

   //**********************************************************************
   // PRE  : defaultValue has been set
   // POST : returns true if an integer has been read in -- the
   //        resulting value is stored in the parameter retVal
   //**********************************************************************

   int len;
   retVal = 0;
   len = s.length();

   if(s == "") {
      return(true);
   }

   for(int i = 0; i < len; i++) {
      char currCh = s.at(i);
      if(!isNum(currCh))
         return(false);
      else {
         int currMult = (int)pow(10,len-i-1);
         retVal += currMult * ((int)currCh - 48);
      }
   }
   if((retVal > lowBounds)  && (retVal < hiBounds)) 
      return(true);
   else
      return(false);
   // check to make sure it is a valid number
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
bool getDouble(double &retVal, double defaultValue) {

   //**********************************************************************
   // PRE  : defaultValue has been set
   // POST : returns true if a double has been read in -- the
   //        resulting value is stored in the parameter retVal
   //**********************************************************************

   string s;
   int len;

   retVal = 0;
   cin >> s;
   len = s.length();
   if(s == "") {
      retVal = defaultValue;
      return(true);
   }

   int numDecimals=0;
   int numDecimalPlaces = 0;

   for(int i = 0; i < len; i++) {
      char currCh = s.at(i);
      if(!isNum(currCh)) {
         if((currCh == '.') && ((numDecimals == 0) && (i > 0))) {
            numDecimals++;
         }
         else {
            return(false);
         }
      }
   }  

   if(getDoubleFromString(retVal, s, -1000000.0, 1000000.0)) {
       return(true);
    }
    else {
       return(false);
    }
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
bool getInt(int &retVal, int defaultValue) {

   //**********************************************************************
   // PRE  : defaultValue has been set
   // POST : returns true if an integer has been read in -- the
   //        resulting value is stored in the parameter retVal
   //**********************************************************************

   string s;
   int len;

   retVal = 0;
   cin >> s;
   len = s.length();
   if(s == "") {
      retVal = defaultValue;
      return(true);
   }

   for(int i = 0; i < len; i++) {
      char currCh = s.at(i);
      if(!isNum(currCh))
         return(false);
      else {
         int currMult = (int)pow(10,len-i-1);
         retVal += currMult * ((int)currCh - 48);
      }
   }

   return(true);
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
void getValidAASequenceTokens (string &retVal, ifstream &ifs, 
                             string validLabel, int errorCode) {

   //**********************************************************************
   // PRE  : validLabel contains a label for the next valid field to be
   //        read in from the file
   // POST : Reads in the next line from the file, storing the value in
   //        the retVal parameter if the label matches the validLabel,
   //        otherwise exiting with the defined errorCode
   //**********************************************************************

   string currLine;

   ifs >> currLine;
   char *lhTok, *rhTok;

   lhTok = strtok((char *)currLine.c_str(), "=");
   rhTok = strtok(NULL, "\n");
   string label = string(lhTok);
   if(label != validLabel)
      exit(errorCode);
   retVal = toUpper(string(rhTok));
   if(!validAASeq(retVal))
      exit(errorCode);
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
void getValidDNASequenceTokens (string &retVal, ifstream &ifs, 
                             string validLabel, int errorCode) {

   //**********************************************************************
   // PRE  : validLabel contains a label for the next valid field to be
   //        read in from the file
   // POST : Reads in the next line from the file, storing the value in
   //        the retVal parameter if the label matches the validLabel,
   //        otherwise exiting with the defined errorCode
   //**********************************************************************

   string currLine;

   ifs >> currLine;
   char *lhTok, *rhTok;

   lhTok = strtok((char *)currLine.c_str(), "=");
   rhTok = strtok(NULL, "\n");
   string label = string(lhTok);
   if(label != validLabel)
      exit(errorCode);
   retVal = toUpper(string(rhTok));
   if(!validDNASeq(retVal))
      exit(errorCode);
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
void getStringTokens (string &retVal, ifstream &ifs, string validLabel, 
                      int errorCode) {

   //**********************************************************************
   // PRE  : validLabel contains a label for the next valid field to be
   //        read in from the file
   // POST : Reads in the next line from the file, storing the value in
   //        the retVal parameter if the label matches the validLabel,
   //        otherwise exiting with the defined errorCode
   //**********************************************************************

   string currLine;
   ifs >> currLine;
   char *lhTok, *rhTok;

   lhTok = strtok((char *)currLine.c_str(), "=");
   rhTok = strtok(NULL, "\n");
   string label = string(lhTok);
   string val = toUpper(string(rhTok));
   if(label != validLabel)
      exit(errorCode);
   retVal = val;
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
void getValidCharTokens (string &retVal, ifstream &ifs, 
                         vector<string> validVals, string validLabel, 
                         int errorCode) {

   //**********************************************************************
   // PRE  : validLabel contains a label for the next valid field to be
   //        read in from the file
   // POST : Reads in the next line from the file, storing the value in
   //        the retVal parameter if the label matches the validLabel,
   //        and the input is listed in the validVals vector,
   //        otherwise exiting with the defined errorCode
   //**********************************************************************

   string currLine;
   ifs >> currLine;
   char *lhTok, *rhTok;

   lhTok = strtok((char *)currLine.c_str(), "=");
   rhTok = strtok(NULL, "\n");
   string label = string(lhTok);
   string val = toUpper(string(rhTok));
   if(label != validLabel)
      exit(errorCode);
   bool isFound = false;
   for(int i = 0; i < validVals.size(); i++) {
      if(val == validVals.at(i)) {
         isFound = true;
         i = validVals.size();
      }
   }  
   if(!isFound) 
      exit(errorCode);
   retVal = val;
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
void getValidDoubleTokens (double &retVal, ifstream &ifs, double low,
                           double hi, string validLabel, int errorCode) {

   //**********************************************************************
   // PRE  : validLabel contains a label for the next valid field to be
   //        read in from the file
   // POST : Reads in the next line from the file, storing the value in
   //        the retVal parameter if the label matches the validLabel,
   //        and is a double in the range from low to hi,
   //        otherwise exiting with the defined errorCode
   //**********************************************************************

   string currLine;

   ifs >> currLine;
   char *lhTok, *rhTok;

   lhTok = strtok((char *)currLine.c_str(), "=");
   rhTok = strtok(NULL, "\n");
   string label = string(lhTok);

   if((label != validLabel) || !getDoubleFromString(retVal, rhTok, low, hi))
      exit(errorCode);
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
void getValidIntTokens (int &retVal, ifstream &ifs, int low,
                        int hi, string validLabel, int errorCode) {

   //**********************************************************************
   // PRE  : validLabel contains a label for the next valid field to be
   //        read in from the file
   // POST : Reads in the next line from the file, storing the value in
   //        the retVal parameter if the label matches the validLabel,
   //        and is an integer in the range from low to hi,
   //        otherwise exiting with the defined errorCode
   //**********************************************************************

   string currLine;

   ifs >> currLine;
   char *lhTok, *rhTok;

   lhTok = strtok((char *)currLine.c_str(), "=");
   rhTok = strtok(NULL, "\n");
   string label = string(lhTok);

   if((label != validLabel) || !getIntFromString(retVal, rhTok, low, hi))
      exit(errorCode);
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
double getDoubleInput(string prompt, double defaultVal, 
                      double lowBounds, double hiBounds) {

   //**********************************************************************
   // PRE  : prompt is a string that is printed out to the user
   //        defaultVal is the default value if enter by itself is used
   //        lowBounds is the first illegal value
   //        highBounds is the last illegal value
   // POST : returns the integer value entered, once it has been
   //        validated according to the parameters
   //**********************************************************************

   bool finished = false;
   double retVal;

   while(!finished) {
      cout << prompt;
      if(getDouble(retVal, defaultVal)) {
         if((retVal > lowBounds)  && (retVal < hiBounds)) {
            finished = true;
         }
      }
      if(!finished) {
         cout << endl << endl;
         cout << "ERROR: NUMBER BETWEEN " << lowBounds + 1 << " AND " <<
              hiBounds - 1 << " EXPECTED!!" << endl << endl;
      }
   }
   return(retVal);
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
int getIntInput(string prompt, int defaultVal, int lowBounds, int hiBounds) {

   //**********************************************************************
   // PRE  : prompt is a string that is printed out to the user
   //        defaultVal is the default value if enter by itself is used
   //        lowBounds is the first illegal value
   //        highBounds is the last illegal value
   // POST : returns the integer value entered, once it has been
   //        validated according to the parameters
   //**********************************************************************

   bool finished = false;
   int retVal;

   while(!finished) {
      cout << prompt;
      if(getInt(retVal, defaultVal)) {
         if((retVal > lowBounds)  && (retVal < hiBounds)) {
            finished = true;
         }
      }
      if(!finished) {
         cout << endl << endl;
         cout << "ERROR: NUMBER BETWEEN " << lowBounds + 1 << " AND " <<
              hiBounds - 1 << " EXPECTED!!" << endl << endl;
      }
   }
   return(retVal);
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
char toUpper(char ch) {

   //**********************************************************************
   // PRE  : ch is the character to change to uppercase
   // POST : returns the uppercase equivalent of ch by using ASCII values
   //**********************************************************************

   int val;
   val = (int) ch;
   if((val > 96) && (val < 123)) {
      val -= 32;
   }
   ch = (char)val;
   return(ch);
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
string toUpper(string s) {

   //**********************************************************************
   // PRE  : s is the string to change to uppercase
   // POST : returns the uppercase equivalent of s
   //**********************************************************************

   int len;
   int val;
   string retVal = "";
   len = s.length();

   for(int i = 0; i < len; i++) {
      char ch = s.at(i);
      ch = toUpper(ch);
      retVal.append(1,ch);
   }
   return(retVal);
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
bool validAASeq(string s) {

   //**********************************************************************
   // PRE  : s is a sequence to be tested
   // POST : returns true if s is a valid Amino Acid sequence; 
   //        false if it is not
   //**********************************************************************

   int len = s.length();

   for(int i = 0; i < len; i++) {
      char ch = s.at(i);
      if(!((ch == 'A') || (ch == 'R') || (ch == 'N') || (ch == 'D') || 
           (ch == 'C') || (ch == 'Q') || (ch == 'E') || (ch == 'G') ||
           (ch == 'H') || (ch == 'I') || (ch == 'L') || (ch == 'K') ||
           (ch == 'M') || (ch == 'F') || (ch == 'P') || (ch == 'S') ||
           (ch == 'T') || (ch == 'W') || (ch == 'Y') || (ch == 'V'))) {
         return(false);
      }
   }
   return(true);
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
bool validDNASeq(string s) {

   //**********************************************************************
   // PRE  : s is a sequence to be tested
   // POST : returns true if s is a valid DNA sequence; false if it is not
   //**********************************************************************

   int len = s.length();

   for(int i = 0; i < len; i++) {
      char ch = s.at(i);
      if(!((ch == 'A') || (ch == 'C') || (ch == 'G') || (ch == 'T'))) {
         return(false);
      }
   }
   return(true);
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
string getSequence(string prompt, char seqType) {

   //**********************************************************************
   // PRE  : prompt is the message to be printed to the user
   // POST : reads in a valid DNA sequence and returns that value
   //**********************************************************************

   string tmp = "1";
   bool done = false;
   while(!done) {
      cout << prompt;
      cin >> tmp;
      tmp = toUpper(tmp);
      if(seqType == 'N') {
         if(validDNASeq(tmp)) {
            done = true;
         }
         else {
            cout << endl << endl << "ERROR: ONLY A, C, G, T ARE VALID" << endl
            << endl;
         }
      }
      else {
         if(validAASeq(tmp)) {
            done = true;
         }
         else {
            cout << endl << endl << "ERROR: INVALID AMINO ACID CODE" << endl 
            << endl;
         }
      }
   }
   return(tmp);
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
string getCharInput(string prompt, string legalVal1, string legalVal2) {

   //**********************************************************************
   // PRE  : reads in the prompt and two possible legal values
   // POST : loops unitl a valid value has been entered -- it is then 
   //        returned as the function's return value
   //**********************************************************************

   string tmp = " ";
   while((tmp != legalVal1) && (tmp != legalVal2)) {
      cout << prompt;
      cin >> tmp;
      tmp = toUpper(tmp);
      if(!((tmp == legalVal1) || (tmp == legalVal2))) {
          cout << endl << endl;
          cout << "*** ERROR: MUST BE " << legalVal1 << " OR ";
          cout << legalVal2 << " ***" << endl << endl;
      }
   }
   return(tmp);
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
bool fileExists(string fn) {

   //**********************************************************************
   // PRE  : fn is a filename
   // POST : returns true if the file fn exists; false if it does not
   //**********************************************************************

   ifstream fin;
   fin.open(fn.c_str());
   if(fin.fail())
      return(false);
   fin.close();
   return(true);
}
//---------------------------------------------------------------------------

//===========================================================================
//=====================  END OF CommonRoutines.cpp ==========================
//===========================================================================


