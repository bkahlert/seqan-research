// ==========================================================================
//                           fragmentStore_tutorial
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

#include <fstream>
#include <iostream>
#include <seqan/sequence.h>
#include <seqan/file.h>
#include <seqan/store.h>

using namespace seqan;
int main()
{
    FragmentStore<> store;
    std::ifstream file("/home/franzi/SeqAn/assignment_annotations.gtf", std::ios_base::in | std::ios_base::binary);
    read(file, store, Gtf());
    // Create AnnotationTree iterator
    Iterator<FragmentStore<>, AnnotationTree<> >::Type it;
    it = begin(store, AnnotationTree<>());
    // Move iterator one node down
	
	int counterAllExon = 0;
	int counterRna = 0;
	int counterLengthAll = 0;
	int counterGene = 0;
	

     do{
	int counter = 0;
	
	while(goDown(it)){
		if(getType(it) != "gene"){
			++counterGene;
		}	
	}

	++counter;

	while(goRight(it)){
		++counter;
		counterLengthAll += getAnnotation(it).endPos -getAnnotation(it).beginPos;	
	}
	
	while(goUp(it) && !goRight(it));
	
	std::cout << "Counted: " << counter << '\n';
	++counterRna;
	counterAllExon += counter;

     }while(!isRoot(it));
 	
    	Iterator<FragmentStore<>, AnnotationTree<> >::Type it2;
    	it2 = begin(store, AnnotationTree<>());

	while(getType(it2) != "exon"){
		if(!goDown(it2))
		{std::cerr << "No exon!" << '\n';}
	}

	std::cout << "Type: " << getType(it2) << '\n' << "Begin: " << getAnnotation(it2).beginPos << '\n';
	std::cout << "End: " << getAnnotation(it2).endPos << '\n' << "ID: " << value(it2) << '\n';
	goUp(it2);
	std::cout << "ParentID: " << value(it2) << '\n' << getName(it2) << '\n';

	std::cout << "Average exonNum per mRna: " << (counterAllExon/counterRna) << '\n';
	std::cout << "Average exonLength: " << (counterLengthAll/counterAllExon) << '\n';
	std::cout << "Average number of rnas per gene: " << (counterRna/counterGene) << '\n';

    return 0;
}