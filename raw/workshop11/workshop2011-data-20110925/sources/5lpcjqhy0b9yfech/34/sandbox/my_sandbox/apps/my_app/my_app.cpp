// ==========================================================================
//                                   my_app
// ==========================================================================
// Copyright (c) 2006-2011, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Your Name <your.email@example.net>
// ==========================================================================

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include <seqan/misc/misc_cmdparser.h>

#include "my_app.h"

using namespace seqan;

// Program entry point
int main(int argc, char const ** argv)
{
    // Setup command line parser.
    CommandLineParser parser;
    Options options;
    setupCommandLineParser(parser, options);
    
    // Then, parser the command line and handle the cases where help display
    // is requested or errornoeous parameters were given.
    int ret = parseCommandLineAndCheck(options, parser, argc, argv);
    if (ret != 0)
        return ret;
    if (options.showHelp || options.showVersion)
        return 0;
    
    // Finally, launch the program.
    ret = mainWithOptions(options);
    return ret;
}
/*
// Hello World

#include <iostream>
#include <seqan/sequence.h>  // CharString, ...
#include <seqan/file.h>      // to stream a CharString into cout

int main(int, char **) {
    std::cout << "Hello World!" << std::endl;
    seqan::CharString mySeqanString = "Hello SeqAn!";
    std::cout << mySeqanString << std::endl;
    return 1;
}

// letterCount

#include <seqan/sequence.h>
#include <seqan/basic.h>
#include <seqan/file.h>
#include <iostream>
using namespace seqan;

void countOneMers(String<char> str) {

	std::cout << "word: " << str << std::endl;
	String<int> letterNum;
	resize(letterNum, 256, 0);
	
	for(int i=0; i<length(str);i++) {
		char thisLetter = str[i];
		letterNum[thisLetter]++;
	}

	for(int i=0; i<length(letterNum);i++) {
		int thisNum = letterNum[i];
		if(thisNum>0){
			std::cout << "letter " << (char)i << ": " << thisNum << std::endl;
		}
	}
	

}

int main(int, char **) {

	String<char> str1 = "hello world";
	String<char> str2 = "banana";
	String<char> str3 = "mississippi";

	countOneMers(str1);
	countOneMers(str2);
	countOneMers(str3);

    return 1;
}

// generic letter count 

#include <seqan/sequence.h>
#include <seqan/basic.h>
#include <seqan/file.h>
#include <iostream>
using namespace seqan;

template <typename T> 
void countOneMers(T & str) {

	std::cout << "word: " << str << std::endl;
	String<int> letterNum;
	resize(letterNum, 256, 0);
	
	for(int i=0; i<length(str);i++) {
		typename Value<T>::Type thisLetter = str[i];
		letterNum[thisLetter]++;
	}

	for(int i=0; i<length(letterNum);i++) {
		int thisNum = letterNum[i];
		if(thisNum>0){
			std::cout << "letter " << (typename Value<T>::Type)i << ": " << thisNum << std::endl;
		}
	}
	

}

int main(int, char **) {

	String<char> str1 = "hello world";
	String<Dna> str2 = "TATACGCTA";
	String<AminoAcid> str3 = "MQDRVKRPMNAFIVWSRDQRRKMALEN";

	countOneMers(str1);
	countOneMers(str2);
	countOneMers(str3);

    return 1;
}

//crazy alphabet

#include <seqan/sequence.h>
#include <seqan/basic.h>
#include <seqan/file.h>
#include <iostream>
using namespace seqan;


void permutations2(int len, String<char> str){
	if(len==2){
		for(char i = 'a'; i<= 'z';i++){
			String<char> newString;
			append(newString,str);
			append(newString,i);
			permutations2(1,newString);
		}
	}else{
		// now want third letter
		for(char i = 'a'; i<= 'z';i++){
			std::cout << str << i << ","; 
		}
		std::cout << std::endl;
	}
}

void printPermutations(int len) {

	for(char i = 'a'; i<= 'z';i++){
		
		String<char> myStr;
		append(myStr,i);
		permutations2(2,myStr);

	}

}


int main(int, char **) {

	printPermutations(3);

    return 1;
}

//replace letter

#include <seqan/sequence.h>
#include <seqan/basic.h>
#include <seqan/file.h>
#include <iostream>
using namespace seqan;


void replaceAs(String<char> str){

	std::cout << str << std::endl;

	typedef Iterator<String<char> >::Type TIterator;
	for (TIterator it = begin(str); it != end(str); ++it){
	
		if(*it == 'a'){
			assignValue(it,'X');
		}
		
	}

	std::cout << str << std::endl;

}


int main(int, char **) {

	String<char> str1 = "abcdefghijklmnopqrstuvxyz";
	String<char> str2 = "Hello SeqAn!";
	String<char> str3 = "Hello Seqan!";

	replaceAs(str1);
	replaceAs(str2);
	replaceAs(str3);

    return 1;
}

// alignment

#include <seqan/sequence.h>
#include <seqan/basic.h>
#include <seqan/file.h>
#include <seqan/align.h>
#include <iostream>
using namespace seqan;


int main(int, char **) {

	typedef String<Dna> TSequence;  // sequence type
    typedef Align<TSequence,ArrayGaps>  TAlign;     // align type
	typedef Row<TAlign>::Type TRow;  

	String<Dna> str1 = "acgtacgtact";
	String<Dna> str2 = "actactacgt";

	TAlign align;
    resize(rows(align), 2);
    assignSource(row(align,0),str1);
    assignSource(row(align,1),str2);

	int score = globalAlignment(align,Score<int>(1,-1,-1,-1));
	
	TRow &row1 = row(align,0);
    TRow &row2 = row(align,1);

	::std::cout << align;

	::std::cout << ::std::endl << "ViewPos1: ";
    for(unsigned i = 0; i < length(source(row1)); ++i){
		if(isGap(row1,i)){
			::std::cout << toViewPosition(row1, i) << ",";
		}
	}
    ::std::cout << ::std::endl << "ViewPos2: ";
    for(unsigned i = 0; i < length(source(row2)); ++i){
		if(isGap(row2,i)){
			::std::cout << toViewPosition(row2, i) << ",";
		}
	}
    ::std::cout << ::std::endl;

    return 1;
}



*/