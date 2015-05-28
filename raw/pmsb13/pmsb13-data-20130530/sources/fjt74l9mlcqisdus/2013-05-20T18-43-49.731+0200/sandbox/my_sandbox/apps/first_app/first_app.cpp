// ==========================================================================
//                                 createData
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
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
#include <fstream>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/file.h>      // For printing SeqAn Strings.
#include <seqan/random.h>

using namespace seqan;

// von euch geklaut
int hash(int x,int y,int z){
	/**@brief hash gibt fuer jedes Codon einen eindeutigen integer wieder.
	*@param x Zahlenwert der Aminosäure an 1. Position 
	*@param y Zahlenwert der Aminosäure an 2. Position
	*@param z Zahlenwert der Aminosäure an 3. Position
	*@return erg gibt die Position in einem Array an , welcher die Aminosaeure gespeichert hat 
	*/
	int erg = (x*16)+(y*4)+z;
	if (erg>=0 && erg<=63) return erg;
	else {
        std::cerr << "invalid input in hash-function" <<std::endl;
		return -1;
	}
}

int get_Amino_Acid_Pos(int pos){
	/**@brief get_Amino_Acid_Pos beinhaltet ein Array, das die Adresse der jeweiligen Aminosaeure beinhaltet.
	*@param pos Position der gesuchten Aminosaeure
	*@return x[pos] Gibt die Adresse der Aminosaeure zurueck
	*/
	if (pos>=0 && pos<=63){
		int x [64] = {11,2,11,2,16,16,16,16,1,15,1,15,9,9,12,9,5,8,5,8,14,14,14,14,1,1,1,1,10,10,10,10,6,3,6,3,0,0,0,0,7,7,7,7,19,19,19,19,20,18,20,18,15,15,15,15,20,4,17,4,10,13,10,13};
		return x[pos];
	}
	else{
        std::cerr << "invalid input in get_Amino_Acid_Pos"<<std::endl;
		return -1;
	}
}

// von euch geklaut
template <typename TText>
int get_translate_from_codon(TText actual_triplet){
	/**@brief get_translate_from_codon gibt die Gruppennummer des reduzierten Alphabet zurueck der jeweiligen Amnosauere
	*@param actual_triplet Triplet, welches fuer die Aminosaeure codiert
	*@param alphabet Reduziertes Alphabet
	*@param reduced Parameter 1 wenn reduziertes Alphabet gewuenscht ist, 0 wenn nicht
	@return amino_pos Gruppennummer der Aminosaeure
	*/
	int hash_value = hash((int) actual_triplet[0],(int) actual_triplet[1],(int) actual_triplet[2]);
	if (hash_value!=-1){
		int amino_pos = get_Amino_Acid_Pos(hash_value);
		if (amino_pos!=-1){
			return amino_pos;
		}
		else return -1;
	}
	else return -1;
}

int main(int argc, char const ** argv)
{
    // statistik
	StringSet<unsigned> read_id;
	StringSet<unsigned> protein_id;
	StringSet<unsigned> start_protein;
	StringSet<unsigned> end_protein;
	
	// typedefs
    typedef StringSet<String<AminoAcid> > TDB;
    typedef StringSet<String<Dna> > TReads;

    // number of different proteins and their max length
    unsigned numProtDB = 10000;
    unsigned maxProtLength = 1000;

    // init of random number generator
    Rng<MersenneTwister> rng(42);

    // creating the data base
    TReads tempDb;
    resize(tempDb, numProtDB);
    for (unsigned i = 0; i < length(tempDb); ++i)
    {
        unsigned localProtLength = (pickRandomNumber(rng) % maxProtLength) * 3 + 150;
        for (unsigned j = 0; j < localProtLength; ++j)
            appendValue(tempDb[i], (Dna)(pickRandomNumber(rng) % 4));
    }

    // number and length of reads
    unsigned readLen = 100;
    unsigned numReads = 1000;

    // creating the reads without errors
    TReads reads;
    resize(reads, numReads);
    for (unsigned i = 0; i < length(reads); ++i)
    {
        unsigned protId = pickRandomNumber(rng) % numProtDB;
      	unsigned readStart = pickRandomNumber(rng) % (length(tempDb[protId]) - readLen + 1);
        reads[i] = infix(tempDb[protId], readStart, readStart + readLen);
		appendValue(read_id,i);
		appendValue(protein_id,protId);
		appendValue(start_protein,readStart/3);
		appendValue(end_protein,(readStart + readLen)/3);
    }

    std::ofstream readFile("reads.fasta", std::ios_base::binary);
    for (unsigned i = 0; i < length(reads); ++i)
    {
        readFile << ">";
        readFile << i;
        readFile << '\n';
        readFile << reads[i];
        readFile << '\n';
    }
    readFile.close();
    
    std::ofstream dbFile("db.fasta", std::ios_base::binary);
    for (unsigned i = 0; i < length(tempDb); ++i)
    {   
        dbFile << ">";
        dbFile << i;
        dbFile << '\n';
        
        for (unsigned j = 0; j < length(tempDb[i]) - 2; j=j+3)
        {
//             std::cerr << infix(tempDb[i], j, j+3) << " " << length(tempDb[i]) << " " << i << std::endl;
            dbFile << (AminoAcid)get_translate_from_codon(infix(tempDb[i], j, j+3));
        }
            
        dbFile << '\n';
    }
    dbFile.close();

	std::ofstream outfile;
	outfile.open("statistik.txt");
    if (outfile.is_open()){
		outfile << "read_id\tprot_id\tprot_begin\tprot_end"<<endl;
		for (unsigned i=0;i<length(read_id);++i){
			outfile << read_id[i]<<"\t"<<protein_id<<"\t"<<start_protein[i]<<"\t"<<end_protein[i]<<endl;
		}
		outfile.close();
	}
	


    return 0;
}
