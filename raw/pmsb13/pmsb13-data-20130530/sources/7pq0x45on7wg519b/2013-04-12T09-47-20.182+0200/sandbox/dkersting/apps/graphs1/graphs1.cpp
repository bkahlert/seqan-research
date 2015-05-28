#include <iostream>
#include <seqan/graph_types.h>
#include <seqan/graph_algorithms.h>
using namespace seqan;
using namespace std;
int main ()
{
    typedef unsigned int TCargo;
    typedef Graph<Directed<> > > TGraph;
    typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
    TVertexDescriptor vertBerlin = addVertex(g);
    TGraph g;

  	TVertexDescriptor 0 = addVertex(g);
	TVertexDescriptor 1 = addVertex(g);
	TVertexDescriptor 2 = addVertex(g);
	TVertexDescriptor 3 = addVertex(g);
	TVertexDescriptor 4 = addVertex(g);
	TVertexDescriptor 5 = addVertex(g);
	TVertexDescriptor 6 = addVertex(g);
	TVertexDescriptor 7 = addVertex(g);

    addEdge(g, 1, 0);
    addEdge(g, 0, 4);
    addEdge(g, 2, 1);
    addEdge(g, 4, 1);
    addEdge(g, 5, 1);
    addEdge(g, 6, 2);
    addEdge(g, 3, 2);
    addEdge(g, 2, 3);
    addEdge(g, 7, 3);
    addEdge(g, 5, 4);
    addEdge(g, 6, 5);
    addEdge(g, 5, 6);
    addEdge(g, 7, 6);
    addEdge(g, 7, 7);


    cout << g <<endl;
