#include <iostream>
#include <fstream>
#include <seqan/basic.h>
#include <seqan/graph_algorithms.h>

using namespace seqan;

typedef LogProb<> TProbability;
typedef Dna TAlphabet;    
typedef Size<TAlphabet>::Type TSize;
typedef Graph<Hmm<TAlphabet, TProbability, Default> > THmm;
typedef VertexDescriptor<THmm>::Type TVertexDescriptor;
typedef EdgeDescriptor<THmm>::Type TEdgeDescriptor;  



int main(){
  
  Dna dnaA = Dna('A');
  Dna dnaC = Dna('C');
  Dna dnaG = Dna('G');
  Dna dnaT = Dna('T');

  THmm hmm;

  TVertexDescriptor begState = addVertex(hmm);
  assignBeginState(hmm, begState);

  TVertexDescriptor exonState = addVertex(hmm);
    emissionProbability(hmm, exonState, dnaA) = 0.25;
    emissionProbability(hmm, exonState, dnaC) = 0.25;
    emissionProbability(hmm, exonState, dnaG) = 0.25;
    emissionProbability(hmm, exonState, dnaT) = 0.25;

    TVertexDescriptor spliceState = addVertex(hmm);
    emissionProbability(hmm, spliceState, dnaA) = 0.05;
    emissionProbability(hmm, spliceState, dnaC) = 0.0;
    emissionProbability(hmm, spliceState, dnaG) = 0.95;
    emissionProbability(hmm, spliceState, dnaT) = 0.0;

    TVertexDescriptor intronState = addVertex(hmm);
    emissionProbability(hmm, intronState, dnaA) = 0.4;
    emissionProbability(hmm, intronState, dnaC) = 0.1;
    emissionProbability(hmm, intronState, dnaG) = 0.1;
    emissionProbability(hmm, intronState, dnaT) = 0.4;
  
  return 0;
}  