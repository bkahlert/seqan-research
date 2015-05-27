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
#include <seqan/file.h>

#include <seqan/misc/misc_cmdparser.h>

#include "my_app.h"

using namespace seqan;

template <typename TAlphabet>
void showAllLettersOfMyAlphabet (TAlphabet const &) {
	typedef typename Size<TAlphabet>::Type TSize;
	TSize alphSize = ValueSize<TAlphabet>::VALUE;
	for (TSize i=0; i<alphSize; ++i)
		std::cout << /*i << "," << */TAlphabet(i) << "  ";
	std::cout << std::endl;
}

template <typename T>
void countOneMers (String<T> const& str) {
	typedef typename Size<T>::Type TSize;
	TSize alphSize = ValueSize<T>::VALUE;
	String<int> counter;
	resize(counter,alphSize,0);
	typedef typename Iterator<String<T> >::Type TIterator;
	for (TIterator it=begin(str); it!=end(str); ++it) {
		counter[(*it)]++;
	}
	for (int i=0; i<alphSize; i++)
		if (counter[i]>0) std::cout << T(i) << "," << counter[i] << "  ";
	std::cout << std::endl;
	
}

bool incPerm (CharString& str, const int pos) {
	if (pos==0) return true;
	if (str[pos-1]=='z') {
		str[pos-1]='a';
		return incPerm (str, pos-1);
	}
	else str[pos-1]++;
	return false;
}

void printPermutations(const int len) {
	CharString str;
	resize(str,len);
	for (int i=0; i<len; i++) str[i]='a';
	bool finished=false;
	while (!finished) {
		std::cout << str << " ";
		finished = incPerm(str,len);
	}
}

// Program entry point
int main(int argc, char const ** argv)
{
    // Setup command line parser.
    CommandLineParser parser;
    Options options;
    setupCommandLineParser(parser, options);
    
    // Then, parse the command line and handle the cases where help display
    // is requested or erroneous parameters were given.
    int ret = parseCommandLineAndCheck(options, parser, argc, argv);
    if (ret != 0)
        return ret;
    if (options.showHelp || options.showVersion)
        return 0;
    
		
    // Finally, launch the program.
//     ret = mainWithOptions(options);
// 		tut Basics - 1
		showAllLettersOfMyAlphabet(AminoAcid());
		showAllLettersOfMyAlphabet(Dna());
		showAllLettersOfMyAlphabet(Dna5());
		showAllLettersOfMyAlphabet(Iupac());
		
// 		task Basics - 1a
		Peptide str_b_1="MQDRVKRPMNAFIVWSRDQRRKMALEN";
		std::cout << str_b_1 << std::endl;
		
// 		task Basics - 1b
		typedef Iterator<Peptide>::Type TIterator;
		for (TIterator it=begin(str_b_1); it != end(str_b_1); ++it)
			if (value(it)=='R') value(it)='A';
		std::cout << str_b_1 << std::endl;
		
// 		task Basics - 1c
		String<int> counter;
		typedef Size<AminoAcid>::Type TSize;
		TSize alphSize = ValueSize<AminoAcid>::VALUE;
		resize(counter,alphSize,0);
		for (TIterator it=begin(str_b_1); it != end(str_b_1); ++it)
			counter[(*it)]++;
		
// 		task Basics - 1d
		showAllLettersOfMyAlphabet(AminoAcid());
		for (int i=0; i<alphSize; i++) std::cout << counter[i] << "  ";
		std::cout << std::endl;
		
// 		tut Basics - 2
		Allocator<SimpleAlloc<> > mySimpleAlloc;
		Allocator<MultiPool<> > myMultiPoolAlloc;
		
		double Tstart, Tend;
		int memsize=1;
		std::cout << "                ";
		for (int i=1; i<=3; i++) {
			memsize=memsize*10;
			std::cout << memsize << "alloc   " << memsize << "clear   ";
		}
		int runs=1000000;
		memsize=1;
		std::cout << "\nSimple Alloc    ";
		char* myString_arr;
		for (int i=1; i<=3; i++) {
			memsize=memsize*10;
			double t_a=0, t_c=0;
				Tstart=sysTime();
			for (int r=0; r<runs; r++) {
				allocate (mySimpleAlloc, myString_arr, memsize,TagAllocateTemp());
			}
				Tend=sysTime();
				t_a+=Tend-Tstart;
				Tstart=sysTime();
				clear (mySimpleAlloc);
				Tend=sysTime();
				t_c+=Tend-Tstart;
			std::cout << t_a << "  ";
			std::cout << t_c << "  ";
		}
		memsize=1;
		std::cout << "\nMultiPool Alloc ";
		for (int i=1; i<=3; i++) {
			memsize=memsize*10;
			double t_a=0, t_c=0;
				Tstart=sysTime();
			for (int r=0; r<runs; r++) {
				allocate (myMultiPoolAlloc, myString_arr, memsize,TagAllocateTemp());
			}
				Tend=sysTime();
				t_a+=Tend-Tstart;
				Tstart=sysTime();
				clear (myMultiPoolAlloc);
				Tend=sysTime();
				t_c+=Tend-Tstart;
			std::cout << t_a << "  ";
			std::cout << t_c << "  ";
		}
		std::cout << std::endl;
		
// 		task Sequences 1
		CharString input_list[3]={"hello world", "banana", "mississippi"};
		for (int i=0; i<3; i++) countOneMers(input_list[i]);
		
// 		task Sequences 1b
		CharString str1="Hello world!";
		countOneMers(str1);
		DnaString str2="TATACGCTA";
		countOneMers(str2);
		Peptide str3="MQDRVKRPMNAFIVWSRDQRRKMALEN";
		countOneMers(str3);
		
// 		task Sequences 2
		printPermutations (3);
    return 0;
}
