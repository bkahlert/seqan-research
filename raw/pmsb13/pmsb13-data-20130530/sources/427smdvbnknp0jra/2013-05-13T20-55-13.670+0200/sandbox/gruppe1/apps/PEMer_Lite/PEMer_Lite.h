/*
Name: saminput
Author: Lars Zerbe <larszerbe@live.de>
Author: Stephan Peter <s.peter@fu-berlin.de>
Maintainer: Lars Zerbe <larszerbe@live.de>
License: GPL v3
Copyright: 2008-2012, FU Berlin
Status: under development
*/

#ifndef PEMER_LITE_H_
#define PEMER_LITE_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include <algorithm>
#include <iterator>
#include <iostream>
#include <seqan/bam_io.h>
#include <seqan/file.h>
#include <seqan/arg_parse.h>
#include <math.h>
#include <vector>
#include <fstream>

struct candidate
{
	enum form{u,d,i};	// form definiert den Type der Kandidaten ob u:undefiniert(0); d:deletion(1); i:insertion(2)

	std::vector<unsigned> pos;	// speichert die Positon
	std::vector<unsigned> ref;	// speichert die Referenz
	std::vector<unsigned> len;	// speichert die Laenge
	std::vector<form> type;		// speichert den Type
	unsigned e;					// speichert den Durchschnitt MEDIAN/AVERRAGE wahlweise
	unsigned sd;				// speichert die Standardabweichung

	unsigned length()	// gibt die larnge aus *.length()
	{
		return len.size();
	}

	void clear()	// leert den Kandidaten
	{
		len.clear();
		pos.clear();
		ref.clear();
		type.clear();
	}
};

struct modifier
{
    unsigned x,y,sd,e;
    bool p;
    seqan::CharString inputFile;
    modifier() :
        x(1), y(1),e(0),sd(0),p(false)
    {}
};

int median(candidate &input);

unsigned average(candidate &input);

seqan::ArgumentParser::ParseResult commands(modifier &options, int argc, char const ** argv);

int analyse(candidate &input);

int saminput(candidate &save, char *file);

int devide(candidate &input,candidate &result,modifier &options);

int tsv(std::vector<std::vector<unsigned>> &input,candidate &ref, char *out);

int findCluster(candidate &input, std::vector<std::vector<unsigned> > &dest, candidate::form indel);

#endif  // PEMER_LITE_HEADER_H_
