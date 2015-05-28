// ==========================================================================
//                                 alignment2
// ==========================================================================
// Copyright (c) 2006-2012, Knut Reinert, FU Berlin
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

#include <iostream>
#include <seqan/align.h>
#include <seqan/graph_msa.h>


using namespace seqan;

void msaDemo()
{
  typedef String<AminoAcid> TSequence;
  typedef StringSet<TSequence> TStringSet;

  typedef StringSet<TSequence, Dependent<> > TDepStringSet;
  typedef Graph<Alignment<TDepStringSet> > TAlignGraph;

  TStringSet seq;
  appendValue(seq,"DPKKPRGKMSSYAFFVQTSREEHKKKHPDASVNFSEFSKKCSERWKTMSAKEKGKFEDMAKADKARYEREMKTYIPPKGE");
  appendValue(seq,"RVKRPMNAFIVWSRDQRRKMALENPRMRNSEISKQLGYQWKMLTEAEKWPFFQEAQKLQAMHREKYPNYKYRPRRKAKMLPK");
  appendValue(seq,"FPKKPLTPYFRFFMEKRAKYAKLHPEMSNLDLTKILSKKYKELPEKKKMKYIQDFQREKQEFERNLARFREDHPDLIQNAKK");
  appendValue(seq,"HIKKPLNAFMLYMKEMRANVVAESTLKESAAINQILGRRWHALSREEQAKYYELARKERQLHMQLYPGWSARDNYGKKKKRKREK");

  TAlignGraph alignG(seq);

  int score = globalMsaAlignment(alignG, Blosum62(-1, -11));
  std::cout << "score= " << score << std::endl;
  std::cout << alignG << ::std::endl;
}


void localAlignmentDemo()
{
  typedef String<Dna> TString;
  typedef Align<TString, ArrayGaps> TAlign;

  std::cout << std::endl;

  TString seq1 = "ataagcgtctcg";
  TString seq2 = "tcatagagttgc";

  

  TAlign align;
  resize(rows(align), 2);
  assignSource(row(align, 0), seq1);
  assignSource(row(align, 1), seq2);

  Score<int> scoring(2, -1, -2, 0);
  int score = localAlignment(align, scoring);


  std::cout << "local: " << std::endl;
  std::cout << score << std::endl;
  std::cout << align << std::endl;
  
  

}


void globalAlignmentDemo()
{
  typedef String<char> TString;
  typedef Align<TString, ArrayGaps> TAlign;

  std::cout << std::endl;

  TString seq1 = "blablubalu";
  TString seq2 = "abba";

  TAlign align;
  resize(rows(align), 2);
  assignSource(row(align, 0), seq1);
  assignSource(row(align, 1), seq2);

  Score<int, Simple> scoring(1, -1, -1);
  int score = globalAlignment(align, scoring, AlignConfig<true, true, true, true>());

  std::cout << score << std::endl;
  std::cout << align << std::endl;
  
  

}




int main()
{

  typedef String<char> TString;
  typedef Align<TString, ArrayGaps> TAlign;

  TString seq1 = "CDFGDC";
  TString seq2 = "CDEFGAHGC";

  TAlign align;
  resize(rows(align), 2);
  assignSource(row(align, 0), seq1);
  assignSource(row(align, 1), seq2);

  insertGap(row(align, 0), 2);
  insertGaps(row(align, 0), 5, 2);


  std::cout << align << std::endl;

  std::cout << "View: " << std::endl;

  for (unsigned i = 0; i < length(row(align, 0)); ++i)
    std::cout << toSourcePosition(row(align, 0), i) << " ";
  std::cout << std::endl;

  for (unsigned i = 0; i < length(row(align, 1)); ++i)
    std::cout << toSourcePosition(row(align, 1), i) << " ";
  std::cout << std::endl;


  std::cout << "View2: " << std::endl;

  for (unsigned i = 0; i < length(source(row(align, 0))); ++i)
    std::cout << toViewPosition(row(align, 0), i) << " ";
  std::cout << std::endl;

  for (unsigned i = 0; i < length(source(row(align, 1))); ++i)
    std::cout << toViewPosition(row(align, 1), i) << " ";
  std::cout << std::endl;



  globalAlignmentDemo();
  localAlignmentDemo();


}



