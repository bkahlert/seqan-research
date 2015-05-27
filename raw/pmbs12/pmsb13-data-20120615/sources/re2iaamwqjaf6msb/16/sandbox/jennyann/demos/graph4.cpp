#include <iostream>
#include <seqan/graph_types.h>
#include <seqan/graph_algorithms.h>

using namespace seqan;

int main ()
{
    typedef Graph<Directed<>> TGraph;
    typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;

	TVertexDescriptor edges[] = {1,0, 0,4, 2,1, 4,1, 5,1, 6,2, 3,2, 2,3, 7,3, 5,4, 6,5, 5,6, 7,6, 7,7};
	TGraph g;

	addEdges(g, edges, 14);

	String<char> nameMap;
    char names[] = {'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h'};
    assignVertexMap(g,nameMap, names);

	/*TVertexDescriptor start = 0;
    typedef Iterator<TGraph, DfsPreorder>::Type TDfsIterator;
    TDfsIterator dfsIt(g, start);
    ::std::cout << "Iterate from '" << getProperty(nameMap, start) << "' in depth-first-search ordering: ";
    while(!atEnd(dfsIt)) {
        ::std::cout << getProperty(nameMap, getValue(dfsIt)) << ", ";
        goNext(dfsIt);
    }*/

	String<unsigned int> component;
    stronglyConnectedComponents(g, component);

	typedef Iterator<TGraph, VertexIterator>::Type VIterator;
    VIterator VIt(g);
    while(!atEnd(VIt)) {
		::std::cout << "Vertex " << getProperty(nameMap, getValue(VIt)) << " : ";
        ::std::cout << getProperty(component, getValue(VIt)) << ::std::endl;
        goNext(VIt);
    }

    ::std::cout << ::std::endl;

	return 0;
}