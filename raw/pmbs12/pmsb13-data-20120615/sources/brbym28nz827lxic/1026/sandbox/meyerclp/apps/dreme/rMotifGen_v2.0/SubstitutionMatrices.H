
//==========================================================================
// -----------------   BEGIN OF  SubstitutionMatricies.H -------------------
//==========================================================================

//==========================================================================
// rMotifGen
// Eric C. Rouchka
// C. Timothy Hardin
// (c) 2004-2007, University of Louisville
//
// FILE: SubstitutionMatrices.H
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
#ifndef MY_PAM_H
#define MY_PAM_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <time.h>
using namespace std;

class PAM {
   public:
      PAM();
      PAM(int val);
      int findAAIndex(char ch); 
      void setInitialMatrix();     // Initializes the PAM1 matrix
      void printMatrix();
      void normalizeMatrix();
      char substitute(char c);     // subsitutes a single character using PAM probabilities
      long double rowSum(int r);   // Finds the sum of the row probabilities
      long double colSum(int c);   // Finds the sum of the column probabilities

   private:
      long double matrix[20][20];  // PAM Matrix (in terms of probabilities)
      char labels[20];             // Amino Acid Labels for row/column header
};

string mutate(string s, PAM p);
string mutate(string s, double pctConserved);

#endif

//==========================================================================
// -------------------   END OF  SubstitutionMatricies.H -------------------
//==========================================================================


