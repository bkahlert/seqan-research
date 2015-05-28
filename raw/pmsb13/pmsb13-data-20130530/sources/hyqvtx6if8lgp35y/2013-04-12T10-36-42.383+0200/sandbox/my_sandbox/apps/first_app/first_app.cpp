#include <iostream>
#include <seqan/graph_types.h>
#include <seqan/graph_algorithms.h>

using namespace seqan;

typedef Graph<Hmm<Default, LogProb, Default> > THmm;
typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
typedef Alphabet<Dna> TDna;

int main()
{
	TGraph hmm;

	TVertexDescriptor beginState = addVertex(hmm);	
	TVertexDescriptor exonState = addVertex(hmm);
	TVertexDescriptor spliceState = addVertex(hmm);
	TVertexDescriptor intronState = addVertex(hmm);
	TVertexDescriptor endState = addVertex(hmm);

	Tdna dnaA = 'A';
	Tdna dnaC = 'C';
	Tdna dnaG = 'G';
	Tdna dnaT = 'T';

	assignBeginState(hmm, beginState);
	assignEndState(hmm, endState);

	assignEmissionProbability(hmm, exonState, dnaA, 0.25);
	assignEmissionProbability(hmm, exonState, dnaC, 0.25);
	assignEmissionProbability(hmm, exonState, dnaG, 0.25);
	assignEmissionProbability(hmm, exonState, dnaT, 0.25);

	assignEmissionProbability(hmm, spliceState, dnaA, 0.05);
	assignEmissionProbability(hmm, spliceState, dnaC, 0.00);
	assignEmissionProbability(hmm, spliceState, dnaG, 0.95);
	assignEmissionProbability(hmm, spliceState, dnaT, 0.00);

	assignEmissionProbability(hmm, intronState, dnaA, 0.4);
	assignEmissionProbability(hmm, intronState, dnaC, 0.1);
	assignEmissionProbability(hmm, intronState, dnaG, 0.1);
	assignEmissionProbability(hmm, intronState, dnaT, 0.4);

	addEdge(hmm, beginState, exonState, 1.0);	
	addEdge(hmm, exonState, exonState, 0.9);
	addEdge(hmm, exonState, spliceState, 0.1);
	addEdge(hmm, spliceState, intronState, 1.0);
	addEdge(hmm, intronState, intronState, 0.9);
	addEdge(hmm, intronState, endState, 0.1);

	std::cout << hmm << ::std::endl;
	
	return 0;
}
