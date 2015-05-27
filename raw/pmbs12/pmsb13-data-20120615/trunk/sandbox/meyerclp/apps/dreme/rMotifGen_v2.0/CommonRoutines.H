//===========================================================================
//=====================  BEGIN OF CommonRoutines.H   ========================
//===========================================================================


//==========================================================================
// rMotifGen
// Eric C. Rouchka
// C. Timothy Hardin
// (c) 2004-2007, University of Louisville
//
// FILE: CommonRoutines.H
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

#ifndef MY_COMMONROUTINES_H
#define MY_COMMONROUTINES_H

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include "SequenceClass.H"
#include "CommonRoutines.H"
using namespace std;


bool isNum(char ch);

bool isNumOrSign(char ch);

bool getDoubleFromString(double &retVal, string s,
                         double lowBounds, double hiBounds);

bool getDouble(double &retVal, double defaultValue);

bool getInt(int &retVal, int defaultValue);

bool getIntFromString(int &retVal, string s, int lowBounds, int hiBounds);

int getIntInput(string prompt, int defaultVal, int lowBounds, int hiBounds);

double getDoubleInput(string prompt, double defaultVal, 
                      double lowBounds, double hiBounds);

char toUpper(char ch);

string getCharInput(string prompt, string legalVal1, string legalVal2);

string toUpper(string s);

bool validDNASeq(string s);
bool validAASeq(string s);

string getSequence(string prompt, char seqType);

string getCharInput(string prompt, string legalVal1, string legalVal2);

bool fileExists(string fn);

void getValidDoubleTokens (double &retVal, ifstream &ifs, double low,
                           double hi, string validLabel, int errorCode);

void getValidIntTokens (int &retVal, ifstream &ifs, int low,
                        int hi, string validLabel, int errorCode);

void getValidCharTokens (string &retVal, ifstream &ifs, 
                         vector<string> validVector, string validLabel, 
                         int errorCode);

void getStringTokens (string &retVal, ifstream &ifs, string validLabel, 
                      int errorCode);

void getValidDNASequenceTokens (string &retVal, ifstream &ifs, 
                             string validLabel, int errorCode);

void getValidAASequenceTokens (string &retVal, ifstream &ifs, 
                             string validLabel, int errorCode);
#endif

//===========================================================================
//=======================  END OF CommonRoutines.H   ========================
//===========================================================================


