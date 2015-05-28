#include <iostream>
#include <seqan/graph_types.h>
#include <seqan/graph_algorithms.h>

using namespace seqan;

typedef unsigned int TCargo;
typedef Graph<Directed<TCargo> > TGraph;
typedef Size<TGraph>::Type TSize;
typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
typedef Iterator<TGraph, DfsPreorder>::Type TIterator;

int main()
{
	TGraph g;
	TSize size = 7;
	TVertexDescriptor edges [] = {1,0,0,4,2,1,4,1,5,1,6,2,3,2,2,3,7,3,5,4, 6,5,5,6,7,6,7,7};
	addEdges(g, edges, size);
	
	CharString vertexMap;
	char names [] = {'a','b','c','d','e','f','g'};
	assignVertexMap(g, vertexMap, names);

	TVertexDescriptor start = 0;
    TIterator it (g, start);
    for(; !atEnd(it); goNext(it)) 
	{
        std::cout << value(it) << ':' << getProperty(names, value(it)) << std::endl;
    }

	return 0;
}

