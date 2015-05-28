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
    TVertexDescriptor vEnd = addVertex(g);
    
    assignBeginState(g, vExon);
    addEdge(g, vExon, vExon, 0.9);
    addEdge(g, vExon, vSplice, 0.1);
    addEdge(g, vSplice, vIntron, 1);
    addEdge(g, vIntron, vIntron, 0.9);
    addEdge(g, vIntron, vEnd, 0.1);
    assignEndState(g, vEnd)
    
    assignEmissionProbability(g, vExon, 'A', 0.25);
    assignEmissionProbability(g, vExon, 'C', 0.25);
    assignEmissionProbability(g, vExon, 'G', 0.25);
    assignEmissionProbability(g, vExon, 'T', 0.25);
    assignEmissionProbability(g, vSplice, 'A', 0.05);
    assignEmissionProbability(g, vSplice, 'G', 0.95);
    assignEmissionProbability(g, vIntron, 'A', 0.4);
    assignEmissionProbability(g, vIntron, 'C', 0.1);
    assignEmissionProbability(g, vIntron, 'G', 0.1);
    assignEmissionProbability(g, vIntron, 'T', 0.4);
    
    int numberOfTestruns = 3;
    StringSet<CharString> sequences;
    resize(sequences,numberOfTestruns);
    StringSet<String<TVertexDescriptor> > states;
    resize(states,numberOfTestruns);
    
    generateSequence(g, sequences, states, numberOfTestruns, 1000)
    for (int i=0;i<numberOfTestruns;++i){
	cout << states << endl;
	cout << sequences << endl << endl;
    }
    FILE* strmWrite = fopen("graph.dot", "w");
    write(strmWrite, g, DotDrawing());
    fclose(strmWrite);
}