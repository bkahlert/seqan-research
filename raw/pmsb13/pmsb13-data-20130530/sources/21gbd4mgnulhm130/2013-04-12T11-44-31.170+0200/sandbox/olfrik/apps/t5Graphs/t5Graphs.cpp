// ==========================================================================
//                                  t5Graphs
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
// ====
// ======================================================================
// Author: Your Name <your.email@example.net>
// ==========================================================================
#include <seqan/basic.h>
#include <seqan/sequence.h>

#include <seqan/arg_parse.h>

#include <iostream>
#include <seqan/graph_types.h>
#include <seqan/graph_algorithms.h>
using namespace seqan;

int digraph ()
{
    typedef unsigned int TCargo;
    typedef Graph<Directed<TCargo> > TGraph;
    typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
    TGraph g;
    TVertexDescriptor vertices[] = {1,0, 0,4, 2,1, 4,1, 5,1, 6,2, 3,2, 2,3, 7,3, 5,4, 6,5, 5,6, 7,6, 7,7};
    int size =  (sizeof( vertices)/sizeof(TVertexDescriptor))/2; 
    std::cout << "sizeof " << size << std::endl;
    addEdges(g, vertices, size); 



    FILE* strmWrite = fopen("graph.dot", "w");
    write(strmWrite, g, DotDrawing());
    fclose(strmWrite);

    std::cout << g << ::std::endl;

    return 0;
}

template <typename TAlphabet>
int hmmGraph(TAlphabet const &){
    typedef unsigned int TCargo;
    typedef Graph<Hmm<Dna,TCargo> > TGraph;
    typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
    typedef double TProb; 
    typedef typename Size<TAlphabet>::Type TSize;
    
    TGraph g;
    TVertexDescriptor vertices[]  = {exonState, spliceState, intronState, startState, endState};
   
/*
    addVertex(g, exonState);
    addVertex(g, spliceState);
    addVertex(g, intronState);
    addVertex(g, startState);
    addVertex(g, endState);
*/
    /*
	                 A	     C	     G	     T
    exon state	    0.25	 0.25	 0.25	 0.25
    splice state	0.05	 0.0	 0.95	 0.0
    intron state	0.4	     0.1	 0.1	 0.4
     *
     *
     */

    TProb emmission_table[3][ValueSize<TAlphabet>::VALUE] = {
        {.25,.25,.25,.25},
        {.05,.0,.95,.0},
        {.4,.1,.1,.1}
    };

    std::cout <<  << std::endl;
   return 0; 
    //assignEmissionProbability(g, exonState, );

    assignTransitionProbability(g,startState, exonState, 1.0);

    assignTransitionProbability(g,exonState,exonState, 0.9);
    assignTransitionProbability(g,exonState,spliceState,0.1);
    
    assignTransitionProbability(g,spliceState, intronState, 1.0);

    assignTransitionProbability(g,intronState, intronState, 0.9);
    assignTransitionProbability(g,intronState, endState, 0.1);



    assignBeginState(g, startState);
//    assignEndState(g, intronState);

    return 0;
}

int main(){
    return hmmGraph(Dna());
}
