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

int hmmGraph(){
    typedef unsigned int TCargo;
    typedef Graph<Hmm<Dna,TCargo> > TGraph;
    typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;

    TGraph g;
    TVertexDescriptor exonState, spliceState, intronState;
    
    addVertex(g, exonState);
    addVertex(g, spliceState);
    addVertex(g, intronState);

    assignBeginState(g, exonState);
//    assignEndState(g, intronState);

    return 0;
}

int main(){
    return hmmGraph();
}
