#include <iostream>
#include <vector>
#include <seqan/graph_types.h>
#include <seqan/graph_algorithms.h>
#include <seqan/stream.h>
#include <cstdio>
#include <fstream>
/*
#include <seqan/basic.h>
#include <seqan/index.h>
/*
#include <seqan/sequence.h>
#include <seqan/file.h>
#include <seqan/modifier.h>
#include <seqan/align.h>
#include <seqan/graph_msa.h>
//*/
using namespace seqan;

int nl() {
	std::cout << "\n";
	return 0;
}

int pP(int len) {
	if(len <= 0) return 1;

	std::vector<int> vec1 (len, 0);
	bool a = true;
	String<char> string1;
	int nextup = 0;

	while(a) {
		string1 = "";
		//a = false;
		for(int i=0;i<=(len-1);++i) {
			string1 += 'a'+vec1[i];

			if(nextup == i) {
				vec1[i]++;
			}
			if(vec1[i] >= 26) {
				a = false;
			}
		}
		std::cout << string1;
		nl();

		nextup++;
		if(nextup >= len) {
			nextup = 0;
		}
	}

	return 0;
}

int replaceAB(String<char> &string1, char a, char b) {
	Iterator< String<char> >::Type it = begin(string1); 
	Iterator< String<char> >::Type itEnd = end(string1); 
	
	while(it != itEnd) {
		if(value(it) == a) {
			value(it) = b;
		}
		++it;
	}


	return 0;
}

int count1mers(String<char> string1) {
	toLower(string1);
	std::vector<int> vec1 (ordValue( ValueSize<char>::VALUE ),0);
	Iterator< String<char> >::Type it = begin(string1); 
	Iterator< String<char> >::Type itEnd = end(string1); 

	while(it != itEnd) {
		vec1[ordValue(value(it))]++;
		++it;
	}

	for(int i=0;i!=ordValue( ValueSize<char>::VALUE );++i) {
		if(vec1[i] > 0) std::cout << char(i) << ": " << vec1[i] << "\n";
	}

	return 0;
}

int stream_task(char* type, char rw, char* filename) {
	if(rw == 'w') {
		if(type == "file" || true) { 
			FILE* strmWrite = fopen(filename, "w");
			fprintf(strmWrite, "Hello World\n");
			fclose(strmWrite);
		}
	}
	else {
		std::ifstream file(filename, std::ios_base::in | std::ios_base::binary);
		do {
			char c = '\0';
			int res = streamReadChar(c, file);
			if (streamEof(file))
				break;  // It's over! But no error.
			if (res)
				return res;  // Pass error code to caller.
			std::cout << c;
		} while (true);
		nl();
	}

	return 0;
}

CharString my_parser(char* filename, CharString mode) {
	CharString outs = "";
	std::ifstream file(filename, std::ios_base::in | std::ios_base::binary);
	do {
		char c = '\0';
		int res = streamReadChar(c, file);
		if (streamEof(file))
			break;  // It's over! But no error.
		if (res)
			return "-1";  // Pass error code to caller.
		if(mode == ".") {
			if( c>=33 && c<=47 ) outs += c;
		}
		else { //hex
			if(c>=97 && c<=102) c-=32;
			if( (c>=48 && c<=57) || (c>=65 && c<=70) ) outs += c;
		}
	} while (true);

	return outs;
}

int main() { 
	/*
	typedef String<AminoAcid> AString;
	AString string1 = "MQDRVKRPMNAFIVWSRDQRRKMALEN";
	Iterator<AString>::Type it = begin(string1); 
	Iterator<AString>::Type itEnd = end(string1); 

	std::vector<int> vec1 (ordValue( ValueSize<AminoAcid>::VALUE ),0);

	while(it != itEnd) {
		vec1[ordValue(value(it))]++;
		//std::cout << value(it) << " - " << vec1[ordValue(value(it))] << "\n";

		if(value(it) == 'R') {
			//value(it) = 'A';
		}
		std::cout << value(it);
		++it;
	}

	nl();
	//std::cout << "\n";

	for(int i=0;i!=ordValue( ValueSize<AminoAcid>::VALUE );++i) {
		std::cout << AminoAcid(i) << ": " << vec1[i] << "\n";
	}

	std::cout << "\n";
	*/

	//pP(3);

	//String<char> string1 = "abcde";
	//replaceAB(string1, 'a', 'X');
	//std::cout << string1;
	//nl();

	//count1mers("hallo lAchse!");

	/*
    typedef String<AminoAcid> TSequence;
    StringSet<TSequence> seq;

    appendValue(seq,"DPKKPRGKMSSYAFFVQTSREEHKKKHPDASVNFSEFSKKCSERWKTMSAKEKGKFEDMAKADKARYEREMKTYIPPKGE");
    appendValue(seq,"RVKRPMNAFIVWSRDQRRKMALENPRMRNSEISKQLGYQWKMLTEAEKWPFFQEAQKLQAMHREKYPNYKYRPRRKAKMLPK");
    appendValue(seq,"FPKKPLTPYFRFFMEKRAKYAKLHPEMSNLDLTKILSKKYKELPEKKKMKYIQDFQREKQEFERNLARFREDHPDLIQNAKK");
    appendValue(seq,"HIKKPLNAFMLYMKEMRANVVAESTLKESAAINQILGRRWHALSREEQAKYYELARKERQLHMQLYPGWSARDNYGKKKKRKREK");

	Graph<Alignment<StringSet<TSequence, Dependent<> > > > aliG(seq);
    globalMsaAlignment(aliG, Blosum80(-1, -11));
    std::cout << aliG << "\n";
	//*/

	/*
	Align< String<AminoAcid> > align;
	appendValue(rows(align),"DPKKPRGKMSSYAFFVQTSREEHKKKHPDASVNFSEFSKKCSERWKTMSAKEKGKFEDMAKADKARYEREMKTYIPPKGE");
    appendValue(rows(align),"RVKRPMNAFIVWSRDQRRKMALENPRMRNSEISKQLGYQWKMLTEAEKWPFFQEAQKLQAMHREKYPNYKYRPRRKAKMLPK");
    appendValue(rows(align),"FPKKPLTPYFRFFMEKRAKYAKLHPEMSNLDLTKILSKKYKELPEKKKMKYIQDFQREKQEFERNLARFREDHPDLIQNAKK");
    appendValue(rows(align),"HIKKPLNAFMLYMKEMRANVVAESTLKESAAINQILGRRWHALSREEQAKYYELARKERQLHMQLYPGWSARDNYGKKKKRKREK");
	globalMsaAlignment(align, Blosum80(-1, -11));
    std::cout << align << "\n";
	*/

	//StringSet< String<char> > strset1;
	/*appendValue(strset1, "tobeornottobe");
	appendValue(strset1, "thebeeonthecomb");
	appendValue(strset1, "beingjohnmalkovich");
	typedef Index< StringSet< String<char> > > INDEX_T;
	INDEX_T index1(strset1);
	Iterator< INDEX_T, TopDown< ParentLinks<Postorder> > >::Type it1(index1); //*/
	/*std::cout << "test\n";
	appendValue(strset1, "CDFGHC");
	appendValue(strset1, "CDEFGAHC");
	typedef Index< StringSet< String<char> > > INDEX_T;
	INDEX_T index1(strset1);
	Iterator< INDEX_T, Mums >::Type it1(index1);
	goBegin(it1);
	while (!atEnd(it1)) {
        std::cout << representative(it1) << "\n";
        ++it1;
    }*/

	/*
	Index< DnaString, IndexQGram<GenericShape> > index1("CATGATTACATA");

	stringToShape(indexShape(index1), "1101");
	hash(indexShape(index1), "ATAA");

    for (unsigned i = 0; i < length(getOccurrences(index1, indexShape(index1))); ++i)
        std::cout << getOccurrences(index1, indexShape(index1))[i] << "\n";
	//*/
	/*
	typedef Graph< Directed<> > GRAPH_DIR;
	typedef Graph<Hmm<Dna, LogProb<>, Default> > GRAPH_HMM;
	LogProb<> p;

	GRAPH_DIR g;
	VertexDescriptor<GRAPH_DIR>::Type es[] = {1,0, 0,4, 2,1, 4,1, 5,1, 6,2, 3,2, 2,3, 7,3, 5,4, 6,5, 5,6, 7,6, 7,7};
	int j = 14;
	addEdges(g, es, j);
	String< CharString > map;

	CharString names[] = {"a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l"};
	assignVertexMap(g, map, names);
	property(map, 0) = "z"; // erstes element umbenennen
    std::cout << g << "\n\n";//*/

	/*Iterator<GRAPH_DIR, DfsPreorder>::Type it(g, 0);
	while(!atEnd(it)) {
		std::cout << getProperty(map, getValue(it)) << "\n";
		++it;
	}//*/
	/*
	String<unsigned int> component;
    stronglyConnectedComponents(g, component); // stronglycc: erreichbar über ger. graph

	Iterator<GRAPH_DIR, VertexIterator>::Type it(g);
	while(!atEnd(it)) {
		std::cout << "Knoten: " << getProperty(map, getValue(it)) << " -> " << getProperty(component, getValue(it)) << "\n";
		++it;
	}//*/

	/*//Size<Dna>::Type
	GRAPH_HMM h;
	VertexDescriptor<GRAPH_HMM>::Type begState = addVertex(h);
	assignBeginState(h, begState);
	VertexDescriptor<GRAPH_HMM>::Type s1 = addVertex(h);
	VertexDescriptor<GRAPH_HMM>::Type s2 = addVertex(h);
	VertexDescriptor<GRAPH_HMM>::Type s3 = addVertex(h);
	VertexDescriptor<GRAPH_HMM>::Type endState = addVertex(h);
    assignEndState(h, endState);

	emissionProbability(h, s1, Dna('A')) = 0.25;
	emissionProbability(h, s1, Dna('C')) = 0.25;
	emissionProbability(h, s1, Dna('G')) = 0.25;
	emissionProbability(h, s1, Dna('T')) = 0.25;
	emissionProbability(h, s2, Dna('A')) = 0.05;
	emissionProbability(h, s2, Dna('C')) = 0.00;
	emissionProbability(h, s2, Dna('G')) = 0.95;
	emissionProbability(h, s2, Dna('T')) = 0.00;
	emissionProbability(h, s3, Dna('A')) = 0.4;
	emissionProbability(h, s3, Dna('C')) = 0.1;
	emissionProbability(h, s3, Dna('G')) = 0.1;
	emissionProbability(h, s3, Dna('T')) = 0.4;

	addEdge(h, begState, s1, 1.0);
    addEdge(h, s1, s1, 0.9);
    addEdge(h, s1, s2, 0.1);
    addEdge(h, s2, s3, 1.0);
    addEdge(h, s3, s3, 0.9);
    addEdge(h, s3, endState, 0.1);

	std::cout << h << "\n\n";

	//viterbi, fw, bw
	String<Dna> sequence = "CTTCATGTGAAAGCAGACGTAAGTCA";
    String<VertexDescriptor<GRAPH_HMM>::Type> path;

	p = viterbiAlgorithm(h, sequence, path);
	std::cout << p << " <- viterbi\n";
	for(int i=0;i<length(sequence);++i) std::cout << sequence[i] << ",";
	nl();
	for(int i=0;i<length(path);++i) std::cout << path[i] << ",";
	nl();

	nl();
	p = forwardAlgorithm(h, sequence);
	std::cout << p << " <- fw\n";
	p = backwardAlgorithm(h, sequence);
	std::cout << p << " <- bw\n";//*/

	//I/O task 1:
	//char* filename = "test123123123.txt";
	//stream_task("file", 'w', filename);

	//parser lachs:
	//CharString lachse = my_parser(filename, ".");
	std::cout << "test" << "\n";

	return 0;
}
