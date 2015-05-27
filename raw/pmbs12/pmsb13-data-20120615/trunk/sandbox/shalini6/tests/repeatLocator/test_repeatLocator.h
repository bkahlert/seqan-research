// ==========================================================================
//                               repeatLocator
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
// Author: Vishalini Vimalakanthan <vishalini.vimalakanthan@fu-berlin.de>
// ==========================================================================

// ################################ DESCRIPTION #####################################
// Given a genome and a set of reads, this program will determine the locations 
// of possible repeatsegments and save them in two separate files, one comprising
// the repeatpositions of the genome and the other those of the readset. In addition
// this program will softmask the genome and save the masked sequence in a third file.
// Designed for FastA -files containing DNA-sequences.
// ##################################################################################

#ifndef SANDBOX_SHALINI6_TESTS_REPEATLOCATOR_TEST_REPEATLOCATOR_H_
#define SANDBOX_SHALINI6_TESTS_REPEATLOCATOR_TEST_REPEATLOCATOR_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <seqan/file.h>
#include <seqan/index.h>
#include <algorithm>
#include <seqan/stream.h>

using namespace seqan;

    //TEST parts of main function
	SEQAN_DEFINE_TEST(test_repeatLocator_reading_files)
{
		//Read genome and store its ID and sequence in Strings
		std::fstream fstrm;
	    fstrm.open("BacillusGenom.fasta", ::std::ios_base::in | ::std::ios_base::binary);
	    RecordReader<std::fstream, SinglePass<> > fstrmReader(fstrm);
	    CharString id; //genomeID
	    String<Dna> seq; // genomesequence
	 
		
	   if (readRecord(id, seq, fstrmReader, Fasta()) != 0){
	        std::cerr << "ERROR: Could not read FASTA!\n";
	       
        }
	    
	    fstrm.close();


		//Read Reads and its IDs into StringSets
		MultiSeqFile multiSeqFile;
		open(multiSeqFile.concat, "CandidatusReads.fasta", OPEN_RDONLY);
        
		AutoSeqFormat format;
		guessFormat(multiSeqFile.concat, format);
		split(multiSeqFile, format);

		unsigned seqCount = length(multiSeqFile);

		StringSet<String<Dna> > readseqs;
		StringSet<CharString> readseqIDs;
		reserve(readseqs, seqCount, Exact());
		reserve(readseqIDs, seqCount, Exact());
		String<Dna> readseq;
		CharString readid;

		for (unsigned i = 0; i < seqCount; ++i){
			assignSeq(readseq, multiSeqFile[i], format);   
			assignSeqId(readid, multiSeqFile[i], format);   

			appendValue(readseqs, readseq, Generous());
			appendValue(readseqIDs, readid, Generous());
		}

	
		//========================================================================================

		//Test reading genomefile, EXAMPLE: first and last line of BacillusGenom.fasta

		String<Dna> genomstr1=""; //first line of input genome
	    for(unsigned i=0; i<70;++i)
		appendValue(genomstr1,seq[i]);

		CharString strGenom1 = "ATCTTTTTCGGCTTTTTTTAGTATCCACAGAGGTTATCGACAACATTTTCACATTACCAACCCCTGTGGA";
		SEQAN_ASSERT_EQ(genomstr1, strGenom1);


		String<Dna> genomstr2="";
		for(unsigned j=(length(seq)-66); j<length(seq); ++j)
		appendValue(genomstr2,seq[j]);


		CharString strGenom2 = "GATTCCTTAATTTTACGGAAAAAAGACAAATTCAAACAATTTGCCCCTAAAATCACGCATGTGGAT"; //last line
		SEQAN_ASSERT_EQ(genomstr2, strGenom2);
		

		//Test reading readfile, EXAMPLE: first and last read of CandidatusReads.fasta

		CharString strRead1 =  "ATGAAAATGTTTTATGAAAAAGATGCAGATGTAGATTTAATTAAAAGTAAAAAAATAGCAATTTTTGGTTATGGAAGCCAAGGTCATGCACATGCACTAAATTTAAAAGATAGTGGAGCCAAAGAAGTTGTTGTAGCGCT";
		append(strRead1,"TAGAGATGGGTCAGCAAGTAAAGCAAAAGCTGAGTCAAAAGGATTAAGAGTTATGAATATGTCAGATGCTGCTGAATGGGCAGAAGTTGCAATGATCTTAACACCTGATGAATTGCAAGCTTCAATTTATAAAAATCATA");
		append(strRead1,"TTGAGCAAAGAATCAAACAAGGAACAAGTTTAGCTTTCGCTCATGGATTAAATATTCATTACAAGCTGATTGATGCTAGAAAAGATTTAGATGTATTTATGGTTGCTCCAAAAGGACCGGGTCACTTGGTTAGAAGTGAG");
		append(strRead1,"TTTGAAAGAGGTGGTGGAGTTCCATGTTTATTTGCAGTTCATCAAGATGGAACAGGTAAAGCAAGAGACCTTGCATTATCATATGCCTCAGCAATAGGTGGGGGGAAATCTGGTATTATTGAAACTACCTTTAAAGATGA");
		append(strRead1,"GTGTGAAACAGATTTATTTGGTGAACAATCAGTTTTATGTGGTGGATTAGTTGAGCTTATTAAAAATGGCTTTGAAACCCTAACTGAAGCTGGTTATGAACCAGAAATGGCATACTTTGAATGTCTACATGAAGTAAAGT");
		append(strRead1,"TGATAGTAGACTTAATTTACGAAGGTGGAATTGCTAATATGAATTACTCAATTTCCAATACAGCAGAGTATGGAGAGTATGTATCTGGTAAAAAGGTTGTTGATAGTGAGAGCAAGAAAAGAATGAAGGAAGTTCTAGCT");
		append(strRead1,"GATATTCAATCTGGTAAGTTCACAAAAGACTGGATGAAAGAGTGTGAAGGTGGCCAAAAGAACTTTTTAAAAATGAGAAAAGATCTTGCGGATCATCCAATTGAAAAGGTTGGTGCTGAACTTAGAGCCATGATGCCTTG");
		append(strRead1,"GATTGGTAAGAAAAAATTAATAGATAGTGATAAGAGTTAA"); //first read 
		SEQAN_ASSERT_EQ(readseqs[0], strRead1);

		
		CharString strRead2 = "ATGTGTAAGATAGTTGAGACAGCAATCAGCAGTGGAGCCAAAACAATAAATATTCCTGACACTGTTGGTTATACAATTTATCAGCCTTGA"; //last read
		SEQAN_ASSERT_EQ(readseqs[length(readseqs)-1], strRead2);
	    
}
	//Function to identify repeatpositions in a given readset
	//Repeatsections are written into a txt-file
	//Functionparameters: reads -Stringset, String containing positions of all overrepresented qGrams, readsIDs-StringSet and qGramlength q
	//Outputfile: "RepeatlocationsInReads.txt"
	void getOccReads(StringSet<String<Dna> > &readset, String<int> &posOfOverrepQgrams, StringSet<CharString> &readIDs, int const q){
	

	//Create QGramIndex for reads
	typedef Index<StringSet<String<Dna> >, IndexQGram<UngappedShape<10> > > ReadIndex;
	typedef Fibre<ReadIndex, QGramDir>::Type ReadQgramDir;
	typedef Iterator<ReadQgramDir, Rooted>::Type ReadQDirIter;
	ReadIndex rQIndex(readset);
	indexRequire(rQIndex, QGramDir());


	ReadQDirIter readqdirStart = begin(indexDir(rQIndex),Rooted());
	ReadQDirIter readqdirEnd = end(indexDir(rQIndex), Rooted());
	
	String<Dna> readdirqgrams; //unhashed qgrams of the Read-QGramDir are temporarily stored in this String 
	
	String< SAValue<ReadIndex>::Type > repeatInRead; //repeatpositions in reads are stored in here  
	
	
	
	//iterate over posOverrepQGrams String which is passed on by countQGrams-Function, use unhash to 
	//obtain associated qGram and hash that to acquire the positions in Reads

		
		for(unsigned a = 0; a<length(posOfOverrepQgrams); ++a){
			

				unhash(readdirqgrams, posOfOverrepQgrams[a],q);
				hash(indexShape(rQIndex),begin(readdirqgrams)); 
		       
				for (unsigned i = 0; i < length(getOccurrences(rQIndex, indexShape(rQIndex))); ++i){
					
					appendValue(repeatInRead, getOccurrences(rQIndex, indexShape(rQIndex))[i]);
			    
				}
		
	
		}
		

	


	   //writing repeatpositions and corresponding readIDs into a .txt file
	   
	   orderOccurrences(repeatInRead); 
	   std::ofstream readRepeatsToFile; 
	   readRepeatsToFile.open ("RepeatlocationsInReads.txt");
	   readRepeatsToFile<<"Repeatpositions in Reads"<<"\n";
	   readRepeatsToFile<<"\n";
	   for(unsigned w = 0; w<length(repeatInRead); ++w)
	   readRepeatsToFile<<readIDs[repeatInRead[w].i1]<<"\t"<<"Repeatrange:"<<repeatInRead[w].i2<<"-"<<repeatInRead[w].i2+(q-1)<<"\n";
	   readRepeatsToFile.close();
	  
		
}

	//TEST function 'getOccReads'
	SEQAN_DEFINE_TEST(test_repeatLocator_function_getOccurrencesReads)
{

	StringSet<String<Dna> > reads; //1st function parameter
	appendValue(reads, "ATGAAAATGTTTTATGAAAAAGATGCAGATGTAGATTTAATTAAAAGTAAAAAAATAGCAATTTTTGGTTATGGAAGCCAAGGTCATGCACATGCACTAAATTTAAAAGATAGTGGAGCCAAAGAAGTTGTTGTAGCGCT");
	appendValue(reads, "ATGGTAAAATCTAAATCCGCATATAGCACACCAAGTAAATCTGGAAAAATTGATACACACATTATTGTTGTATGGGTTGATAACGAAGCAAGTGTTTTATCCAGAGTGGTAGGGCTATTTTCAGGTAGAGGATATAATAT");
	appendValue(reads, "ATGCCAAAATTATATTCAGGAGCAGAGATAGTATTCAAGTGTCTTGAAGACCAAAAAGTTGAACATATTTTTGGATATCCAGGTGGAGCAGTTCTTCCTATCTATGATGAACTAAAAAATCATCCAACTATAAAACATAT");
	appendValue(reads, "ATGGACAAGCAGTCCAAAATTATCTTAATTTCAGGACCGACTGCATCAGGTAAATCAAATTTTGCTGTTAAGATTGCAAAAAAAATTGAAGGAGAAATCATCAATGCCGACAGTATGCAGGTTTACAAAAAATTAAAAAT");
	appendValue(reads, "ATGATAAACTCTGAATATTTAAACAACCTGAATAATGCTCAAAAAGAAGCTGTTTTATATCTTGATGGACCGTTACTAATTGTTGCAGGAGCTGGATCTGGAAAAACGAAAGTTCTCACTTCAAGAATTGCTCATATAAT");

	String<int> positions = 0; //2nd function parameter
	appendValue(positions, 5);
	appendValue(positions, 10);
	appendValue(positions, 20);
	appendValue(positions, 50);

	StringSet<CharString> readID; //3rd function parameter
	appendValue(readID, "id1");
	appendValue(readID, "id2");
	appendValue(readID, "id3");
	appendValue(readID, "id3");
	appendValue(readID, "id4");
	
	
	int const q = 10; //4th function parameter
	
	getOccReads(reads, positions, readID, 10);


}

//Function to softmask the given genomesequence and save softmasked Genome in File
//Furthermore Repeatpositions are written into a .txt file
//Functionparameters: toMask = genomesequence, repeatposInGenome = repeatpositions in genome, readIDs-StringSet and genomeID
//outputfiles : "RepeatlocationsInGenome.txt" and "SoftmaskedGenome.txt"
template <typename TString>
void masking(TString &toMask, String<int> &repeatposInGenome,int q, StringSet<CharString> &readIDs, String<char> &genID){ 
	

	std::ofstream locationsGenome;
	locationsGenome.open("RepeatlocationsInGenome.txt");
	locationsGenome<<"Repeatpositions in Genome:"<<genID<<"\n";
	locationsGenome<<"\n";

	std::ofstream maskieren;
	maskieren.open ("SoftmaskedGenome.txt"); //softmasked genome

	orderOccurrences(repeatposInGenome);
	
	for(unsigned j = 0; j < length(repeatposInGenome); ++j){
		locationsGenome<< repeatposInGenome[j]<<"-"<<repeatposInGenome[j]+(q-1)<<"\n";

		String<char> repeat; // stores q bases of pos[j] to pos[j]-(q-1), in order to replace capitals in small letters
		
		for(unsigned x = 0; x<q; ++x){ 
			appendValue(repeat, value(toMask, repeatposInGenome[j]+x));
		}
		
		for(unsigned l = 0; l<length(repeat); ++l){
		
			switch(repeat[l]){

				case 'A' : repeat[l] = 'a'; break;

				case 'C' : repeat[l] = 'c'; break;

				case 'G' : repeat[l] = 'g'; break;

				case 'T' : repeat[l] = 't'; break;

				case 'a' : repeat[l] = 'a'; break;

				case 'c' : repeat[l] = 'c'; break;

				case 'g' : repeat[l] = 'g'; break;

				case 't' : repeat[l] = 't'; break;

		
			}
		
		
		 }
		replace(toMask, repeatposInGenome[j],repeatposInGenome[j]+q, repeat); // actual replacement of repeatsegments by small letters using 'repeat'

		clear(repeat);
	  }
	   
	   maskieren<<"Softmasking of Genome :"<< genID<<"\n";
	   maskieren<<"\n";
	   maskieren <<toMask<<"\n"; //writing the replaced genomesequence into a txt-file
	   
	   maskieren.close();
	   
	   locationsGenome.close();
	  

}

    //TEST function 'masking'
	SEQAN_DEFINE_TEST(test_repeatLocator_function_masking)
{
	//1st function parameter
	CharString genomToMask =  "ATGAAAATGTTTTATGAAAAAGATGCAGATGTAGATTTAATTAAAAGTAAAAAAATAGCAATTTTTGGTTATGGAAGCCAAGGTCATGCACATGCACTAAATTTAAAAGATAGTGGAGCCAAAGAAGTTGTTGTAGCGCT";
				append(genomToMask,"TAGAGATGGGTCAGCAAGTAAAGCAAAAGCTGAGTCAAAAGGATTAAGAGTTATGAATATGTCAGATGCTGCTGAATGGGCAGAAGTTGCAATGATCTTAACACCTGATGAATTGCAAGCTTCAATTTATAAAAATCATA");
				append(genomToMask,"TTGAGCAAAGAATCAAACAAGGAACAAGTTTAGCTTTCGCTCATGGATTAAATATTCATTACAAGCTGATTGATGCTAGAAAAGATTTAGATGTATTTATGGTTGCTCCAAAAGGACCGGGTCACTTGGTTAGAAGTGAG");
				append(genomToMask,"TTTGAAAGAGGTGGTGGAGTTCCATGTTTATTTGCAGTTCATCAAGATGGAACAGGTAAAGCAAGAGACCTTGCATTATCATATGCCTCAGCAATAGGTGGGGGGAAATCTGGTATTATTGAAACTACCTTTAAAGATGA");
				append(genomToMask,"GTGTGAAACAGATTTATTTGGTGAACAATCAGTTTTATGTGGTGGATTAGTTGAGCTTATTAAAAATGGCTTTGAAACCCTAACTGAAGCTGGTTATGAACCAGAAATGGCATACTTTGAATGTCTACATGAAGTAAAGT");
				append(genomToMask,"TGATAGTAGACTTAATTTACGAAGGTGGAATTGCTAATATGAATTACTCAATTTCCAATACAGCAGAGTATGGAGAGTATGTATCTGGTAAAAAGGTTGTTGATAGTGAGAGCAAGAAAAGAATGAAGGAAGTTCTAGCT");
				append(genomToMask,"GATATTCAATCTGGTAAGTTCACAAAAGACTGGATGAAAGAGTGTGAAGGTGGCCAAAAGAACTTTTTAAAAATGAGAAAAGATCTTGCGGATCATCCAATTGAAAAGGTTGGTGCTGAACTTAGAGCCATGATGCCTTG");
				append(genomToMask,"GATTGGTAAGAAAAAATTAATAGATAGTGATAAGAGTTAA");

	String<int> positions = 10; //2nd function parameter
				appendValue(positions, 25);
				appendValue(positions, 150);
				appendValue(positions, 300);
				appendValue(positions, 350);

	int q = 10; //3rd finction parameter

	StringSet<CharString> readID; //4th function parameter
	appendValue(readID, "id1");
	appendValue(readID, "id2");
	appendValue(readID, "id3");
	appendValue(readID, "id3");
	appendValue(readID, "id4");

	CharString genID = "GenomID"; //5th function parameter

	masking(genomToMask, positions, 10, readID, genID);



}


//Function to count qGrams and identify overrepresented repeats
//calls functions to retrieve repeatpositions in reads, softmask the given genome and save the results in files
//Functionparameters: genome-String, QGramIndex, reads-StringSet, genome-CharString (for masking purpose), readIDs-StringSet, genomeID and qlength q

template <typename TString, typename TIndexSpec>
void countQGrams(TString &genomestr, TIndexSpec, StringSet<String<Dna> > &reads, String<char> &maskgenom, StringSet<CharString> &readID,String<char> &genId, int const q){

	//Create QGramIndex for genome
	typedef Index<TString, TIndexSpec > TIndex;
	typedef Fibre<TIndex, QGramDir>::Type TQgramDir;
	typedef Iterator<TQgramDir, Rooted>::Type TIterQDir;
	

	TIndex qindex(genomestr);
	indexRequire(qindex, QGramDir());
	
	TIterQDir qiterstart = begin(indexDir(qindex),Rooted());
	TIterQDir qiterend = end(indexDir(qindex), Rooted());
	

	String<Dna> qgram; //unhashed qgrams of the genome-QGramDir are temporarily stored in this String
	String<int> counts; //stores number of hits for each qgram
	
		
		//iterate over QGramDir and count number of occurrences for each of its entry
		for(TIterQDir k= qiterstart; k!=qiterend; ++k){
			
			unhash(qgram,position(k),q);
				
			hash(indexShape(qindex), begin(qgram));
				
			int ocr = countOccurrences(qindex, indexShape(qindex));
			appendValue(counts, ocr);
				
		}

		
	

	String<int> posOfOverrepQgrams; //stores positions of overrepresented qgrams to be passed on to 'getOccReads'

	
	int max = counts[0];       

     for(int i = 1; i<length(counts); ++i)
     {
          if(counts[i] > max)
             max = counts[i];
     }
	 int threshold = max*0.25; //computing threshold
	
	  for(unsigned i= 0; i<length(counts);++i){

			if(counts[i]>= threshold){   //threshold fixed after analysing various histograms and testing coverage of masked genome

				appendValue(posOfOverrepQgrams, i);

			}
	  }



	  String<Dna> tmp; //unhashed qgrams of the genome-QGramDir are temporarily stored in this String
	  String<int> posOfOverrepInGenome; //stores positions of overrepresented q-Grams -> to be passed on to 'masking'  
	 
	  

		  //retrieve only positions of overrepresented qgrams
		  for(unsigned o = 0; o<length(posOfOverrepQgrams); ++o){
			
			
			  
			  unhash(tmp,posOfOverrepQgrams[o],q);
			  
			  hash(indexShape(qindex), begin(tmp));  
			  
			  String<int> positions = getOccurrences(qindex, indexShape(qindex)); 
			   for(unsigned i = 0; i < length(positions); ++i){
			  
				 appendValue(posOfOverrepInGenome, positions[i]); 
			   }
			  
		     }
		  
	   

	  
	   
	  

}	


//TEST function 'countQGrams'
SEQAN_DEFINE_TEST(test_repeatLocator_function_countQGrams)
{
	//1st function parameter
	String<Dna> genom = "ATGAAAATGTTTTATGAAAAAGATGCAGATGTAGATTTAATTAAAAGTAAAAAAATAGCAATTTTTGGTTATGGAAGCCAAGGTCATGCACATGCACTAAATTTAAAAGATAGTGGAGCCAAAGAAGTTGTTGTAGCGCT";
	append(genom,"TAGAGATGGGTCAGCAAGTAAAGCAAAAGCTGAGTCAAAAGGATTAAGAGTTATGAATATGTCAGATGCTGCTGAATGGGCAGAAGTTGCAATGATCTTAACACCTGATGAATTGCAAGCTTCAATTTATAAAAATCATA");
	append(genom,"TTGAGCAAAGAATCAAACAAGGAACAAGTTTAGCTTTCGCTCATGGATTAAATATTCATTACAAGCTGATTGATGCTAGAAAAGATTTAGATGTATTTATGGTTGCTCCAAAAGGACCGGGTCACTTGGTTAGAAGTGAG");
	append(genom,"TTTGAAAGAGGTGGTGGAGTTCCATGTTTATTTGCAGTTCATCAAGATGGAACAGGTAAAGCAAGAGACCTTGCATTATCATATGCCTCAGCAATAGGTGGGGGGAAATCTGGTATTATTGAAACTACCTTTAAAGATGA");
	append(genom,"GTGTGAAACAGATTTATTTGGTGAACAATCAGTTTTATGTGGTGGATTAGTTGAGCTTATTAAAAATGGCTTTGAAACCCTAACTGAAGCTGGTTATGAACCAGAAATGGCATACTTTGAATGTCTACATGAAGTAAAGT");
	append(genom,"TGATAGTAGACTTAATTTACGAAGGTGGAATTGCTAATATGAATTACTCAATTTCCAATACAGCAGAGTATGGAGAGTATGTATCTGGTAAAAAGGTTGTTGATAGTGAGAGCAAGAAAAGAATGAAGGAAGTTCTAGCT");
	append(genom,"GATATTCAATCTGGTAAGTTCACAAAAGACTGGATGAAAGAGTGTGAAGGTGGCCAAAAGAACTTTTTAAAAATGAGAAAAGATCTTGCGGATCATCCAATTGAAAAGGTTGGTGCTGAACTTAGAGCCATGATGCCTTG");
	append(genom,"GATTGGTAAGAAAAAATTAATAGATAGTGATAAGAGTTAA");


	StringSet<String<Dna> > reads; //3st function parameter
	appendValue(reads, "ATGAAAATGTTTTATGAAAAAGATGCAGATGTAGATTTAATTAAAAGTAAAAAAATAGCAATTTTTGGTTATGGAAGCCAAGGTCATGCACATGCACTAAATTTAAAAGATAGTGGAGCCAAAGAAGTTGTTGTAGCGCT");
	appendValue(reads, "ATGGTAAAATCTAAATCCGCATATAGCACACCAAGTAAATCTGGAAAAATTGATACACACATTATTGTTGTATGGGTTGATAACGAAGCAAGTGTTTTATCCAGAGTGGTAGGGCTATTTTCAGGTAGAGGATATAATAT");
	appendValue(reads, "ATGCCAAAATTATATTCAGGAGCAGAGATAGTATTCAAGTGTCTTGAAGACCAAAAAGTTGAACATATTTTTGGATATCCAGGTGGAGCAGTTCTTCCTATCTATGATGAACTAAAAAATCATCCAACTATAAAACATAT");
	appendValue(reads, "ATGGACAAGCAGTCCAAAATTATCTTAATTTCAGGACCGACTGCATCAGGTAAATCAAATTTTGCTGTTAAGATTGCAAAAAAAATTGAAGGAGAAATCATCAATGCCGACAGTATGCAGGTTTACAAAAAATTAAAAAT");
	appendValue(reads, "ATGATAAACTCTGAATATTTAAACAACCTGAATAATGCTCAAAAAGAAGCTGTTTTATATCTTGATGGACCGTTACTAATTGTTGCAGGAGCTGGATCTGGAAAAACGAAAGTTCTCACTTCAAGAATTGCTCATATAAT");


	CharString maskGenom = genom; //4th function parameter

	StringSet<CharString> readID; //5th function parameter
	appendValue(readID, "id1");
	appendValue(readID, "id2");
	appendValue(readID, "id3");
	appendValue(readID, "id3");
	appendValue(readID, "id4");

	CharString genID = "GenomID"; //6th function parameter

	int const q = 11; //7th function parameter

	countQGrams(genom, IndexQGram<UngappedShape<q> >(), reads, maskGenom, readID, genID, q);



}


#endif  // SANDBOX_SHALINI6_TESTS_REPEATLOCATOR_TEST_REPEATLOCATOR_H_
