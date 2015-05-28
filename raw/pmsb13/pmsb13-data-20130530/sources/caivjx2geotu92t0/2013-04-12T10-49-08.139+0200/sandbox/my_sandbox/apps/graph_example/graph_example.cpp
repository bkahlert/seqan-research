#include <iostream>
#include <seqan/graph_types.h>
#include <seqan/graph_algorithms.h>
#include <array>
#include <seqan/basic.h>

using namespace seqan;
using namespace std;

int main() {
	typedef unsigned int TCargo;
    typedef Graph<Directed<TCargo> > TGraph;
    typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;

	TGraph graph;
	
	TVertexDescriptor vertZero = addVertex(graph);
	TVertexDescriptor vertOne = addVertex(graph);
    TVertexDescriptor vertTwo = addVertex(graph);
    TVertexDescriptor vertThree = addVertex(graph);
    TVertexDescriptor vertFour = addVertex(graph);
    TVertexDescriptor vertFive = addVertex(graph);
	TVertexDescriptor vertSix = addVertex(graph);
	TVertexDescriptor vertSeven = addVertex(graph);

	//(1,0), (0,4), (2,1), (4,1), (5,1), (6,2), (3,2), (2,3), (7,3), (5,4), (6,5), (5,6), (7,6), (7,7)
	TVertexDescriptor con []  = {1,0,0,4,2,1,4,1,5,1,6,2,3,2,2,3,7,3,5,4,6,5,5,6,7,6, 7,7};
	int i = 14;
	addEdges(graph, con, i);
	
	std::cout << graph << ::std::endl;

	typedef Graph<Hmm<> > THmmGraph;
	typedef VertexDescriptor<THmmGraph>::Type THmmVertexDescriptor;
	Dna a = 'A';
	Dna c = 'C';
	Dna g = 'G';
	Dna t = 'T';
	THmmGraph h;
	THmmVertexDescriptor exon =  addVertex(h);
	assignBeginState(h, exon);
	assignEmissionProbability(h, exon, a, 0.25);
	assignEmissionProbability(h, exon, c, 0.25);
	assignEmissionProbability(h, exon, g, 0.25);
	assignEmissionProbability(h, exon, t, 0.25);
	THmmVertexDescriptor splice = addVertex(h);
	THmmVertexDescriptor intron = addVertex(h);

	std::cout << h << ::std::endl;


	return 0;
}