#include <iostream>
#include <seqan/graph_types.h>
#include <seqan/graph_algorithms.h>
#include <seqan/basic.h>

using namespace seqan;

int main ()
{
	typedef LogProb<> TCargo;
	typedef Graph<Hmm<Dna, TCargo> > TGraph;
	typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef String<char> TStateName;
    typedef String<TStateName> TProperties;

	TGraph g;

	TVertexDescriptor start = addVertex(g);
	TVertexDescriptor exon = addVertex(g);
	TVertexDescriptor intron = addVertex(g);
	TVertexDescriptor splice = addVertex(g);
	TVertexDescriptor end = addVertex(g);

	TProperties stateNames;
	resizeVertexMap(g, stateNames);

	assignProperty(stateNames, exon, "Exon");
	assignProperty(stateNames, intron, "Intron" );
	assignProperty(stateNames, splice, "Splice" );
	assignProperty(stateNames, start, "Start" );
	assignProperty(stateNames, end, "End" );

	
	addEdge(g, start, exon, (TCargo)1.0);
	addEdge(g, exon, exon, (TCargo)0.9);
	addEdge(g, exon, splice, (TCargo)0.1);
	addEdge(g, splice, intron, (TCargo)1.0);
	addEdge(g, intron, intron, (TCargo)0.9);
	addEdge(g, intron, end, (TCargo)0.1);

	assignBeginState(g, start);
	assignEndState(g, end);

	emissionProbability(g, exon, (Dna) 'A') = 0.25;
	emissionProbability(g, exon, (Dna) 'C') = 0.25;
	emissionProbability(g, exon, (Dna) 'G') = 0.25;
	emissionProbability(g, exon, (Dna) 'T') = 0.25;
	emissionProbability(g, intron, (Dna) 'A') = 0.4;
	emissionProbability(g, intron, (Dna) 'C') = 0.1;
	emissionProbability(g, intron, (Dna) 'G') = 0.1;
	emissionProbability(g, intron, (Dna) 'T') = 0.4;
	emissionProbability(g, splice, (Dna) 'A') = 0.05;
	emissionProbability(g, splice, (Dna) 'C') = 0.0;
	emissionProbability(g, splice, (Dna) 'G') = 0.95;
	emissionProbability(g, splice, (Dna) 'T') = 0.0;
	
	typedef Iterator<TGraph, VertexIterator>::Type TVertexIterator;
    TVertexIterator it(g);
    for(;!atEnd(it);goNext(it)) {
        std::cout << value(it) << ": " << getProperty(stateNames, value(it)) << std::endl;
    }

	::std::cout << g << ::std::endl;

    
    return 0;
}