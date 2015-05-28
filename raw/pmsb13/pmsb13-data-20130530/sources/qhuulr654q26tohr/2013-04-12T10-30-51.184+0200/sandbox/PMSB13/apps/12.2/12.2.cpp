#include <iostream>
#include <seqan/graph_types.h>
#include <seqan/graph_algorithms.h>
using namespace seqan;
using namespace std;

int main ()
{
    typedef Graph<Hmm<> > TGraph;
    typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
    
    TGraph g;
    
    TVertexDescriptor vIntron = addVertex(g);
    TVertexDescriptor vExon = addVertex(g);
    TVertexDescriptor vSplice = addVertex(g);
    
    addEdge(g, vIntron, vIntron, 0.9);
    addEdge(g, vExon, vExon, 0.9);
    addEdge(g, vSplice, vSplice, 0.1);
    addEdge(g, vExon, vSplice, 0.1);
    addEdge(g, vIntron, vSplice, 0.9);
    addEdge(g, vIntron, vIntron, 0.9);
    addEdge(g, vIntron, vIntron, 0.9);
    addEdge(g, vIntron, vIntron, 0.9);
    addEdge(g, vIntron, vIntron, 0.9);
    addEdge(g, vIntron, vIntron, 0.9);
    addEdge(g, vIntron, vIntron, 0.9);
    
    FILE* strmWrite = fopen("graph.dot", "w");
    write(strmWrite, g, DotDrawing());
    fclose(strmWrite);
}