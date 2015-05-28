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
#include <string>


struct medianTree
{
private:
   int leftSize, rightSize, wert;

    medianTree() :
        leftSize(0), rightSize(0),wert(NULL)
    {}
public:
   void add(int x)
   {
	   if(x<=wert){leftSize+=1;}
	   else{rightSize+=1;}
	   if((leftSize-rightSize)>1){wert=x;leftSize-=1;rightSize+=1;}
	   if((leftSize-rightSize)<0){wert=x;rightSize-=1;leftSize+=1;}
   }

   int getMedian(){return wert;}
};


struct candidate
{
	enum form{undefine,deletion,insertion};	// form definiert den Type der Kandidaten ob u:undefiniert(0); d:deletion(1); i:insertion(2)

	unsigned pos;	// speichert die Positon
	unsigned ref;	// speichert die Referenz
	unsigned len;	// speichert die Laenge
	form type;	// speichert den Type

	void addvalues(unsigned r,unsigned p,unsigned l,form t)	
	{
		ref=r;
		pos=p;
		len=l;
		type=t;
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

struct statisticals
{
	unsigned e_value;
	unsigned standardDegression;
	statisticals() :
        standardDegression(0)
    {}
};

typedef std::vector<candidate> candidates;

bool compareCandidatesLength(candidates x,candidate y);

unsigned median(candidate &input);

seqan::ArgumentParser::ParseResult commands(modifier &options, int argc, char const ** argv);

int analyse(candidate &input);

int saminput(candidate &save, char *file);

int devide(candidate &input,candidate &result,modifier &options);

int tsv(std::vector<std::vector<unsigned>> &input,candidate &ref, char *out);

int findCluster(candidate &input, std::vector<std::vector<unsigned> > &dest, candidate::form indel);

#endif  // PEMER_LITE_HEADER_H_
