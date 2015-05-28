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
	String<char> pm;
	char prop[] = {'a','b','c','d','e','f','g','h'};
	assignVertexMap(graph, pm, prop);
	
	std::cout << graph << ::std::endl;

    typedef Iterator<TGraph, VertexIterator>::Type TVertexIterator;
    TVertexIterator itV(graph);
    for(;!atEnd(itV);goNext(itV)) {
        std::cout << value(itV) << ':' << getProperty(pm, value(itV)) << std::endl;
    }


	// Assignment 2
	/*
	typedef Graph<Hmm<> > THmmGraph;
	typedef VertexDescriptor<THmmGraph>::Type THmmVertexDescriptor;
	Dna a = 'A';
	Dna c = 'C';
	Dna g = 'G';
	Dna t = 'T';
	THmmGraph h;
	THmmVertexDescriptor exon =  addVertex(h);
	emissionProbability(h, exon, a) = 0.25;
	emissionProbability(h, exon, c) = 0.25;
	emissionProbability(h, exon, g) = 0.25;
	emissionProbability(h, exon, t) = 0.25;
	THmmVertexDescriptor splice = addVertex(h);
	emissionProbability(h, splice, a) = 0.05;
	emissionProbability(h, splice, c) = 0.0;
	emissionProbability(h, splice, g) = 0.95;
	emissionProbability(h, splice, t) = 0.0;
	THmmVertexDescriptor intron = addVertex(h);
	emissionProbability(h, intron, a) = 0.4; 	 	 	
	emissionProbability(h, intron, c) = 0.1;
	emissionProbability(h, intron, g) = 0.1;
	emissionProbability(h, intron, t) = 0.4;
	THmmVertexDescriptor start = addVertex(h);
	assignBeginState(h, start);
	THmmVertexDescriptor end = addVertex(h);
	assignEndState(h, end);
	addEdge(h, start, exon, 1);
	addEdge(h, exon, exon, 0.9);
	addEdge(h, intron, intron, 0.9);
	addEdge(h, exon, splice, 0.1);
	addEdge(h, splice, intron, 0.1);
	addEdge(h, intron, end, 1);
	
	std::cout << h << ::std::endl;
	*/
	// Assignment 3


	return 0;
}