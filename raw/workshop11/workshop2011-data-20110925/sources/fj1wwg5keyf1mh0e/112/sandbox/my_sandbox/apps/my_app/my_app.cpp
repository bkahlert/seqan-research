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
#include <seqan/align.h>
#include <iostream>

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
	std::cout << std::endl;
}

void replaceAs (CharString& str, const char from, const char to) {
	std::cout << str << " => ";
	typedef Iterator<CharString,Rooted>::Type TIterator;
	TIterator it=begin(str);
	while (!atEnd(it)) {
		if (value(it)==from) value(it)=to;
		++it;
	}
	std::cout << str << std::endl;
}

struct MyFunctor : public std::unary_function<char, char> {
	inline char operator()(char x) const {
		if (('a'<=x) && (x<='z')) return ('z'-x+'a');
		if (('A'<=x) && (x<='Z')) return ('Z'-x+'A');
		return x;
	}
};

template <typename T>
void globalAlign (String<T> str1, String<T> str2) {
		typedef String<T> TSequence;
		typedef typename Align<TSequence, ArrayGaps> TAlign;
 		typedef Row<TAlign>::Type TRow;
		typedef Iterator<TRow>::Type TIterator;
		TAlign align;
		resize(rows(align),2);
		assignSource(row(align,0),str1);
		assignSource(row(align,1),str2);
		std::cout << align;
// 		TRow &row1 = row(align,0), &row2=row(align,1);
// 		insertGap(row1,2);
// 		insertGap(row1,5);
		int score = globalAlignment(align,Score<int>(1,-1,-1,-1));
		std::cout << "Score = " << score << std::endl << align;
		
		for (unsigned int i=0; i<length(rows(align)); i++) {
			TIterator it=begin(row(align,i));
			TIterator itEnd=end(row(align,i));
			int pos=0;
			std::cout << "Row " << i << ": ";
			while (it!=itEnd) {
				if (isGap(it)) std::cout << pos << " ";
				++it, ++pos;
			}
			std::cout << std::endl;
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
		int runs=100000;
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
// 		printPermutations (3);

// 		tut Sequences
		CharString str_s_1="start_middle_end";
		String<char, CStyle> input_s_1;
		infix(str_s_1,6,12)=input_s_1;
		std::cout << str_s_1 << std::endl;

// 		task Sequences 6
		CharString str_s_6="abcdefghijklmnopqrstuvxyz"; 
		replaceAs (str_s_6,'a','X');
		str_s_6 = "Hello SeqAn!";
		replaceAs (str_s_6,'a','X');
		str_s_6 = "Hello Seqan!";
		replaceAs (str_s_6,'a','X');
		
// 		tut Modifiers
		CharString str_m_1="A man, a plan, a canal-Panama";
		ModifiedString<CharString, ModReverse> mod_m_1(str_m_1);
		ModifiedString<CharString, ModView<MyFunctor> > mod_m_2(str_m_1);
		std::cout << str_m_1 << "; " << mod_m_1 << "; " << mod_m_2 << std::endl;
		infix(str_m_1, 9, 9) ="master ";
		std::cout << str_m_1 << "; " << mod_m_1 << "; " << mod_m_2 << std::endl;
		
		DnaString str_m_2 = "attacgg";
		DnaStringReverseComplement mod_m_3(str_m_2);
		std::cout << str_m_2 << std::endl;
		std::cout << mod_m_3 << std::endl;
		infix (str_m_2, 1,1)="cgt";
		std::cout << str_m_2 << std::endl;
		std::cout << mod_m_3 << std::endl;
		
// 		tut Alignments
		typedef String<char> TSequence;
		typedef Align<TSequence, ArrayGaps> TAlign;
		typedef Row<TAlign>::Type TRow;
		TSequence seq1="CDFGHC", seq2="CDEFGAHC";
		TAlign align;
		resize(rows(align),2);
		assignSource(row(align,0),seq1);
		assignSource(row(align,1),seq2);
		std::cout << align;
		TRow &row1 = row(align,0), &row2=row(align,1);
// 		insertGap(row1,2);
// 		insertGap(row1,5);
		int score = globalAlignment(align,Score<int>(1,-1,-1,-1));
		std::cout << "Score = " << score << std::endl << align;
		
		seq1 = str_m_1;
		seq2 = mod_m_1;
		assignSource(row(align,0),seq1);
		assignSource(row(align,1),seq2);
		score = globalAlignment(align,Score<int>(1,-1,-1,-1), Hirschberg());
		std::cout << "Score = " << score << std::endl << align;
		
// 		task Alignments 1
		DnaString dna_a_1="acgtacgtact", dna_a_2="actactacgt";
		globalAlign(dna_a_1, dna_a_2);
    return 0;
}
