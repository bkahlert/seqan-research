// ################################ DESCRIPTION #######################################
// Given a genome and a set of reads, this program will determine the locations 
// of possible repeatsegments and save them in two separate files, one comprising
// the repeatpositions of the genome and the other those of the readset. In addition
// this program will softmask the genome and save the masked sequence in a third file.
// Designed for FastA -files containing DNA-sequences.
// ####################################################################################
#include <iostream>
#include <fstream>
#include <cstdio>
#include <seqan/sequence.h>
#include <seqan/file.h>
#include <seqan/index.h>
#include <algorithm>
#include <seqan/stream.h>
	
using namespace seqan;


//Function to identify repeatpositions in a given readset
//Repeatsections are written into a txt-file
//Functionparameters: reads -Stringset, String containing positions of all overrepresented qGrams, readsIDs-StringSet and qGramlength q
//Outputfile: "RepeatlocationsInReads.txt"

void getOccReads(StringSet<String<Dna> > &readset, String<int> &posOfOverrepQgrams, StringSet<CharString> &readIDs, int const q){
	
	std::cout<<"Obtaining Repeatpositions in Reads"<<std::endl;

	//Create QGramIndex for reads
	typedef Index<StringSet<String<Dna> >, IndexQGram<UngappedShape<11> > > ReadIndex;
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
	   std::cout<<"Writing Repeatpositions of Reads in File 'RepeatlocationsInReads.txt '"<<std::endl;
	   orderOccurrences(repeatInRead); 
	   std::ofstream readRepeatsToFile;
	   readRepeatsToFile.open ("RepeatlocationsInReads.txt");
	   readRepeatsToFile<<"Repeatpositions in Reads"<<"\n";
	   readRepeatsToFile<<"\n";
	   for(unsigned w = 0; w<length(repeatInRead); ++w)
	   readRepeatsToFile<<readIDs[repeatInRead[w].i1]<<"\t"<<"Repeatrange:"<<repeatInRead[w].i2<<"-"<<repeatInRead[w].i2+(q-1)<<"\n";
	   readRepeatsToFile.close();
	   std::cout<<"Writing COMPLETED"<<std::endl;
	   std::cout<<"\n"<<std::endl;
		
}

//Function to softmask the given genomesequence and save softmasked Genome in File
//Furthermore Repeatpositions are written into a .txt file
//Functionparameters: toMask = genomesequence, repeatposInGenome = repeatpositions in genome, readIDs-StringSet and genomeID
//Outputfiles : "RepeatlocationsInGenome.txt" and "SoftmaskedGenome.txt"
template <typename TString>
void masking(TString &toMask, String<int> &repeatposInGenome,int q, StringSet<CharString> &readIDs, String<char> &genID){ 
	std::cout<<"Softmasking Genome"<<std::endl;

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
	   std::cout<<"Softmasking Genome COMPLETED"<<std::endl;
	   std::cout<<"Softmasked Genome saved as 'SoftmaskedGenome.txt'"<<std::endl;
	   locationsGenome.close();
	   std::cout<<"'Repeatlocations of Genome to be found in File 'RepeatlocationsInGenome.txt'"<<std::endl;

}




//Function to count qGrams and identify overrepresented qgrams
//calls the other functions to retrieve repeatpositions in reads, softmask the given genome and save the results in files
//Functionparameters: genome-String, QGramIndex, reads-StringSet, genome-CharString (for masking purpose), readIDs-StringSet, genomeID and qlength q
//returns String containing positions of overrepresented qGrams, which is used for the summery output
template <typename TString, typename TIndexSpec>
String<int> countQGrams(TString &genomestr, TIndexSpec, StringSet<String<Dna> > &reads, String<char> &maskgenom, StringSet<CharString> &readID,String<char> &genId, int const q){

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
	 int threshold = max*0.25;  //computing threshold
	
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
	

	  
	   
	   //call other functions
	   masking(maskgenom, posOfOverrepInGenome, q, readID,genId);

	   
	   getOccReads(reads, posOfOverrepQgrams, readID, q);

	   return posOfOverrepInGenome;

}	




		

int main(int argc, char const ** argv){
		
	//Check whether number of input parameters correct.
	if(argc!=3){

			std::cerr << "Error: Too many arguments"<<std::endl
					  << "USAGE: " << argv[0] << " FILE" << std::endl;

		    return 1;
	}
		std::cout<<"Start RepeatLocator Program"<<std::endl;
		std::cout<<""<<std::endl;
		std::cout<<"Reading Genomefile"<<"\t"<<argv[1]<<std::endl;

		//Read genome and store its ID and sequence in Strings
		std::fstream fstrm;
	    fstrm.open(argv[1], ::std::ios_base::in | ::std::ios_base::binary);
	    RecordReader<std::fstream, SinglePass<> > fstrmReader(fstrm);
	    CharString id; //genomeID
	    String<Dna> seq; // genomesequence
	 
	    if (readRecord(id, seq, fstrmReader, Fasta()) != 0){
	        std::cerr << "ERROR: Could not read FASTA!\n";
	        return 1;
        }
	    
	    fstrm.close();


		String<char> chargenom = seq;  //converting DNAString to CharString for softmasking-function
		
		std::cout<<"Reading Genomefile COMPLETED"<<std::endl;
		 
		//Read readset and store its IDs and sequences in StringSets
		std::cout<<"Reading Readsfile"<<"\t"<<argv[2]<<std::endl;
		
		MultiSeqFile multiSeqFile;
		open(multiSeqFile.concat, argv[2], OPEN_RDONLY);
        
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

		std::cout<<"Reading Readsfile COMPLETED"<<std::endl;
		
		int const q = 11;
		
		String<int> posOfOverrepGenome = countQGrams(seq,IndexQGram<UngappedShape<q> >(), readseqs, chargenom, readseqIDs, id, q ); 
		

	   std::ofstream results;
	   results.open ("SoftmaskedGenome.txt", ::std::ios::in | ::std::ios::out | ::std::ios::ate);
	   results<<"\n";
	   results<<"\n";
	   results<<"Summary of Repeatmasking"<<"\n";
	   results<<"\n";
	   results<<"Genome"<<"\t"<<"Reads"<<"\t"<<"Genomelength"<<"\t"<<"Complete Repeatlength of Genome"<<"\t"<<"Softmasked Genome in %"<<"\n";
	   results<<"\n";
	   results<<argv[1]<<"\t"<<argv[2]<<"\t"<<length(seq)<<"\t"<<length(posOfOverrepGenome)*q<<"\t"<<(float)((length(posOfOverrepGenome)*q)*100)/length(seq) <<"\n";
	   results.close();
	  
	     
		std::cout<<"Summary of Repeatmasking:"<<std::endl;
		std::cout<<"======================================="<<std::endl;
		std::cout<<"Genome:"<<"\t"<<argv[1]<<std::endl;
		std::cout<<"Reads:"<<"\t"<<argv[2]<<std::endl;
		std::cout<<"Genomelength:"<<"\t"<<length(seq)<<std::endl;
		std::cout<<"Complete Repeatlength of Genome:"<<"\t"<<length(posOfOverrepGenome)*q<<std::endl;
		std::cout<<"Softmasked Genome in %:"<<"\t"<<(float)((length(posOfOverrepGenome)*q)*100)/length(seq)<<std::endl;
		std::cout<<"======================================="<<std::endl;
	
		return 0;

	}