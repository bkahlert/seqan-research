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
	assignValue(cityNames,0,"a");
	assignValue(cityNames,1,"b");
	assignValue(cityNames,2,"c");
	assignValue(cityNames,3,"d");
	assignValue(cityNames,4,"e");
	assignValue(cityNames,5,"f");
	assignValue(cityNames,6,"g");
	assignValue(cityNames,7,"h");
	/*
	String<char> nameMap;
    char names[] = {'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h'};
    assignVertexMap(g,nameMap, names);

	*/

	::std::cout << g << ::std::endl;
	::std::cout << cityNames << ::std::endl;

	typedef Iterator<TGraph, DfsPreorder>::Type TVertexIterator;
	TVertexDescriptor start;
	for(int i =0;i<length(cityNames);++i){
	start=i;
    TVertexIterator itV(g,start);
    for(;!atEnd(itV);goNext(itV)) {
        ::std::cout << value(itV) << ':' << getProperty(cityNames, value(itV)) << ::std::endl;
    }
	::std::cout<<::std::endl;
	}
	String<TCargo> component;



	stronglyConnectedComponents(g,component);
	
	TVertexIterator TvI(g,0);
	while(!atEnd(TvI)){

		::std::cout<<getProperty(cityNames,value(TvI))<<" ";
		::std::cout<<getProperty(component,value(TvI))<<::std::endl;

		++TvI;
	}


	return 0;

}