#include <iostream>
#include <seqan/graph_types.h>
#include <seqan/graph_algorithms.h>
#include <seqan/find_motif.h>


using namespace seqan;

int main ()
{
	typedef unsigned int TCargo;
    typedef Graph<Undirected<TCargo> > TGraph;
    typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
	TGraph g;
	TVertexDescriptor edges[]={1,0,0,4,2,1,4,1,5,1,6,2,3,2,2,3,7,3,5,4,6,5,5,6,7,6,7,7};
	addEdges(g,edges,14);
	
	char letters[]={'a','b','c','d','e','f','g','h'};
	String<char>map;
	assignVertexMap(g, map, letters);
	 TVertexDescriptor start = 0;
    typedef Iterator<TGraph, DfsPreorder>::Type TDfsIterator;
    TDfsIterator dfsIt(g, start);

    std::cout << "Iterate from '" << getProperty(nameMap, start) << "' in depth-first-search ordering: ";
    while(!atEnd(dfsIt)) {
        std::cout << getProperty(nameMap, getValue(dfsIt)) << ", ";
        goNext(dfsIt);
    }
    std::cout << std::endl;
	
	::std::cout << g << ::std::endl;
	return 0;
}