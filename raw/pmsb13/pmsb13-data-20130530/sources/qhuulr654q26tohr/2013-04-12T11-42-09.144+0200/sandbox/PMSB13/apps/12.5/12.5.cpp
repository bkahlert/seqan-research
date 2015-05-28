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
    
    TVertexDescriptor vStart = addVertex(g);
    TVertexDescriptor vExon = addVertex(g);
    TVertexDescriptor vSplice = addVertex(g);
    TVertexDescriptor vIntron = addVertex(g);
    TVertexDescriptor vEnd = addVertex(g);
    
    assignBeginState(g, vStart);
    addEdge(g, vStart, vExon, 1);
    addEdge(g, vExon, vExon, 0.9);
    addEdge(g, vExon, vSplice, 0.1);
    addEdge(g, vSplice, vIntron, 1);
    addEdge(g, vIntron, vIntron, 0.9);
    addEdge(g, vIntron, vEnd, 0.1);
    assignEndState(g, vEnd);
    
    assignEmissionProbability(g, vExon, Dna('A'), 0.25);
    assignEmissionProbability(g, vExon, Dna('C'), 0.25);
    assignEmissionProbability(g, vExon, Dna('G'), 0.25);
    assignEmissionProbability(g, vExon, Dna('T'), 0.25);
    assignEmissionProbability(g, vSplice, Dna('A'), 0.05);
    assignEmissionProbability(g, vSplice, Dna('G'), 0.95);
    assignEmissionProbability(g, vIntron, Dna('A'), 0.4);
    assignEmissionProbability(g, vIntron, Dna('C'), 0.1);
    assignEmissionProbability(g, vIntron, Dna('G'), 0.1);
    assignEmissionProbability(g, vIntron, Dna('T'), 0.4);
    
    /*
    int numberOfTestruns = 10;
    StringSet<CharString> sequences;
    resize(sequences,numberOfTestruns);
    StringSet<String<TVertexDescriptor> > states;
    resize(states,numberOfTestruns);
    
    generateSequence(g, sequences, states, numberOfTestruns, 1000);
    for (int i=0;i<numberOfTestruns;++i){
	cout << states[i] << endl;
	cout << " " << sequences[i] << endl << endl;*/
    
    DnaString s="CTTCATGTGAAAGCAGACGTAAGTCA";
    String<TVertexDescriptor> statePath;
    viterbiAlgorithm(g, s, statePath);
    cout <<" "<< s << endl << statePath << endl;
    
}