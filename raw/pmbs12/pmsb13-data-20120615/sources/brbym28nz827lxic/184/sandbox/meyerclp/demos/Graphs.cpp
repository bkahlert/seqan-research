#include <iostream>
#include <seqan/graph_types.h>
#include <seqan/graph_algorithms.h>
using namespace seqan;

int main (){

	typedef unsigned int TCargo;
    typedef Graph<Directed<TCargo> > TGraph;
    typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
	TGraph g;
	TVertexDescriptor edges[]={1,0,0,4,2,1,4,1,5,1,6,2,3,2,2,3,7,3,5,4,6,5,5,6,7,6,7,7};
	addEdges(g,edges,14);
	//FILE* strmWrite = fopen("graph.dot", "w");
 //   write(strmWrite, g, DotDrawing());
 //   fclose(strmWrite);
	//::std::cout << g << ::std::endl;

	typedef String<char> TCityName;
	typedef String<TCityName> TProperties;
	TProperties cityNames;
	resizeVertexMap(g,cityNames);
	assignProperty(cityNames,0,"a");
	assignProperty(cityNames,1,"b");
	assignProperty(cityNames,2,"c");
	assignProperty(cityNames,3,"d");
	assignProperty(cityNames,4,"e");
	assignProperty(cityNames,5,"f");
	assignProperty(cityNames,6,"g");
	assignProperty(cityNames,7,"h");

	::std::cout << g << ::std::endl;
	::std::cout << cityNames << ::std::endl;

	typedef Iterator<TGraph, DfsPreorder>::Type TVertexIterator;
    TVertexIterator itV(g);
    for(;!atEnd(itV);goNext(itV)) {
        ::std::cout << value(itV) << ':' << getProperty(cityNames, value(itV)) << ::std::endl;
    }
	return 0;

}