
//==========================================================================
// -----------------   BEGIN OF  SubstitutionMatricies.cpp -----------------
//==========================================================================

//==========================================================================
// rMotifGen
// Eric C. Rouchka
// C. Timothy Hardin
// (c) 2004-2007, University of Louisville
//
// FILE: SubstitutionMatrices.cpp
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

#include "SubstitutionMatrices.H"

//==========================================================================
long double PAM::rowSum(int r) {

   //***********************************************************************
   // PRE:  r is the row number (between 0 and 19) to find the sum for
   // POST: returns the sum of the probabilities for row r
   // DESCRIPTION: 
   //***********************************************************************

   if((r < 0) || (r > 19)) {
      exit(0);
   }

   long double rS = 0;
   for(int c = 0; c < 20; c++) {
      rS += matrix[r][c];
   }   
   return(rS);
}
//==========================================================================

//==========================================================================
long double PAM::colSum(int c) {

   //***********************************************************************
   // PRE:  c is the column number (between 0 and 19) to find the sum for
   // POST: returns the sum of the probabilities for column c
   // DESCRIPTION: 
   //***********************************************************************

   if((c < 0) || (c > 19)) {
      exit(0);
   }

   long double cS = 0;
   for(int r = 0; r < 20; r++) {
      cS += matrix[r][c];
   }
   return(cS);
}
//==========================================================================

//==========================================================================
char PAM::substitute(char c) {

   //***********************************************************************
   // PRE:  
   // POST:
   // DESCRIPTION: 
   //***********************************************************************

   int index = findAAIndex(c);
   long double rS = rowSum(index);
   long double runningSum = 0;
   double rVal = drand48();
   for(int col = 0; col < 20; col++) {
      runningSum += matrix[index][col] / rS;
      if(rVal < runningSum) {
         return(labels[col]);
      }
   }
   return(c);
}
//==========================================================================

//==========================================================================
PAM::PAM() {
   //***********************************************************************
   // PRE:  
   // POST:
   // DESCRIPTION: 
   //***********************************************************************
   setInitialMatrix();
}
//==========================================================================

//==========================================================================
PAM::PAM(int val) {

   //***********************************************************************
   // PRE:   val is the number of the PAM matrix, between 1 and 500
   // POST:  sets the matrix data member
   // DESCRIPTION: note that PAM matrices are calculated by multiplying
   //              the PAM1 n times to get PAMn -- if PAM is < 1 or > 500,
   //              them PAM1 is used
   //***********************************************************************

   setInitialMatrix();              // calculates PAM1
   if((val > 1) && (val < 501)) {

      // FIRST CALCULATE PAM1
      long double PAM1[20][20];
      for(int r = 0; r < 20; r++) {
         for(int c = 0; c < 20; c++) {
            PAM1[r][c] = matrix[r][c];
         }
      }
 
      // THEN FIND PAMn BY MULTIPLYING BY PAM1 n TIMES
      long double tmpMatrix[20][20];
      for(int i = 2; i <= val; i++) {
         for(int r = 0; r < 20; r++) {
            for(int c = 0; c < 20; c++) {
               long double cVal = 0;
               for(int j = 0; j < 20; j++) {
                  cVal += PAM1[r][j] * matrix[j][c];
               }
               tmpMatrix[r][c] = cVal;
            }
         }
         for(int r = 0; r < 20; r++) {
            for(int c = 0; c < 20; c++) {
               matrix[r][c] = tmpMatrix[r][c];
            }
         }
      }
  }
  else {
    if(val == 0) {
      for(int r = 0; r < 20; r++) {
         for(int c = 0; c < 20; c++) {
            if(r == c) {
               matrix[r][c] = 1;
            }
            else {
               matrix[r][c] = 0;
            }
         }
      }
    }
     // ELSE DO NOTHING -- PAM1 BY DEFAULT
  } 
}
//==========================================================================

//==========================================================================
void PAM::normalizeMatrix() {

   //***********************************************************************
   // PRE:   NONE
   // POST:  Normalizes a row by dividing by the rowSum
   // DESCRIPTION:  NOT A GOOD IDEA TO USE, SINCE COLS ARE NO LONGER VALID
   //***********************************************************************

   for(int r = 0; r < 20; r++) {
      long double rS = rowSum(r);
      for(int c = r; c < 20; c++) {
         matrix[r][c] /= rS;
         matrix[c][r] /= rS;
      }
  }
}
//==========================================================================

//==========================================================================
void PAM::printMatrix() {

   //***********************************************************************
   // PRE:   the pam matrix has been allocated and set
   // POST:  prints out row by row the values in the PAM matrix
   // DESCRIPTION: 
   //***********************************************************************

   cout << setw(6);
   for(int r = 0; r < 20; r++) {
      for(int c = 0; c < 20; c++) {
         cout << matrix[r][c] << " ";
      }
      cout << endl;
   }
}

//==========================================================================

//==========================================================================
int PAM::findAAIndex(char ch) {

   //**************************************************************************
   // PRE:   ch is a valid (uppercase) amino acid character
   // POST:  returns the index of the character in the labels array
   //        or -1 if it is not found
   // DESCRIPTION: 
   //**************************************************************************

   for(int i = 0; i < 20; i++) {
      if(labels[i] == ch) {
         return(i);
      }
   }
   return(-1);
}

//==========================================================================

//==========================================================================
void PAM::setInitialMatrix() {

   //**************************************************************************
   // PRE:   the file PAM1.prob exists, and contains the PAM1 matrix in 
   //        terms of probabilities multiplied by 10000
   // POST:
   // DESCRIPTION: 
   //**************************************************************************

   ifstream infile;
   infile.open("PAM1.prob", ios::in);
   if(!infile) {
      exit(7);
   }
   for(int c = 0; c < 20; c++) {
      char ch; 
      infile >> ch;
      labels[c] = ch;
   }
   for(int r = 0; r < 20; r++) {
      char ch;
      infile >> ch;
      for(int c = 0; c < 20; c++) {
         int tmpVal;
         infile >> tmpVal;
         matrix[r][c] = tmpVal / 10000.0; 
      }
   } 
   for(int r = 0; r < 20; r++) {
      for(int c = r+1; c < 20; c++) {
         long double avg = (matrix[r][c] + matrix[c][r]) / 2.0;
         matrix[r][c] = matrix[c][r] = avg;
      }
   }
   infile.close();
}
//==========================================================================

//==========================================================================
string mutate(string s, PAM p) {

   //**************************************************************************
   // PRE:  s is a valid amino acid sequence; p is a PAM matrix
   // POST: returns the string s mutated according to the PAM matrix
   // DESCRIPTION:  This function takes in a string and a mutation matrix
   //               and mutates each position in the string independently
   //               returning the mutated string
   //**************************************************************************

   int len = s.length();
   string newStr;
   for(int i = 0; i < len; i++) {
      char ch = s.at(i);
      char newCh = p.substitute(ch);  // substitute according to the PAM
      newStr += newCh;                // append character onto the end
   }
   return(newStr);
}
//==========================================================================

//==========================================================================
string mutate(string s, double pctConserved)
{
   //**************************************************************************
   // PRE:  NONE
   // POST: NONE
   // DESCRIPTION: Mutates a DNA string according to the percent conservation
   //**************************************************************************

   string retVal;
   retVal = s;
   int len = (int)s.length();

   for(int i = 0; i < len; i++)
   {
      int randVal = lrand48() % 101;
      if(randVal > pctConserved)
      {
         string oldCh = retVal.substr(i, 1);
         string newCh = oldCh;
         while(oldCh == newCh)
         {
            int randVal2 = lrand48() % 101;
            if(randVal2 < 25)      newCh = "A";
            else if(randVal2 < 50) newCh = "C";
            else if(randVal2 < 75) newCh = "G";
            else                   newCh = "T";
         }
         retVal.replace(i, 1, newCh);
      }
   }
   return(retVal);
}
//==========================================================================
    

//==========================================================================
// -----------------   END OF  SubstitutionMatricies.cpp -------------------
//==========================================================================



