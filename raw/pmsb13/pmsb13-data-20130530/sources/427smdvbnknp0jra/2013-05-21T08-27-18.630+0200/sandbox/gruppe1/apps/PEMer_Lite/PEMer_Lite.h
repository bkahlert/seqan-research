/*
Name: PEMer_Lite
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


struct medianTree	// struktur welche den Median Speichert, abh�ngig von allen in eingegebenen Werten mit konstantem Speicher.
{
private:
   int leftSize, rightSize, wert;	//Initialisierung: wert speichert den aktuellen Median, leftSize und rightSize die gr��e der elemente die kleiner oder g��er des Medians sind.
public:
    medianTree() :	//konstruktor der beim erzeugen  'leftSize, rightSize, wert' mit Anfangswerten definiert.
        leftSize(0), rightSize(0),wert(NULL)
    {}

   void add(int x)	// add f�gt eine neues Element x ein und �berpr�ft ob es ein Median ist.
   {
	   if(x<=wert){leftSize+=1;}	//ist x kleiner/gleich als der aktuelle Median wird leftSize erh�ht.
	   else{rightSize+=1;}			//ist x gr��er als der aktuelle Median wird rightSize erh�ht.
	   // Jetzt wir dgepr�ft ob der aktuelle wert immer noch der Median ist oder ob x an seine Stelle kommt.
	   // Dazu wird gepr�ft ob die differenz der Anzahl der Elemente die gr��er oder kleiner sind gr��er als 1 ist.
	   if((leftSize-rightSize)>1){wert=x;leftSize-=1;rightSize+=1;}	// x ist neuer Median, da die kleineren Elemente in der Anzahl gr��er sind als die gro�en Elemente. Es findet eine Verschiebung statt und der alte Median wird ersetzt, sowie die gr��en von kleiner und gr��eren Elementen angepasst
	   if((leftSize-rightSize)<0){wert=x;rightSize-=1;leftSize+=1;}	// x ist neuer Median, da die gr��eren Elemente in der Anzahl gr��er sind als die kleinen Elemente. Es findet eine Verschiebung statt und der alte Median wird ersetzt, sowie die gr��en von kleiner und gr��eren Elementen angepasst
   }

   int getMedian(){return wert;}	// gibt den Median aus.
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
    unsigned x,y,sd,e,s,i,d;
    bool p;
    seqan::CharString inputFile;
    modifier() :
        x(1), y(1),e(0),sd(0),p(false),s(false)
    {}
};

struct statisticals
{
	unsigned e_value;
	unsigned standardDegression;
	unsigned coverage;
	statisticals() :
        standardDegression(0)
    {}
};

typedef std::vector<candidate> candidates;

bool candidateSortFunction(candidate x, candidate y);

int samSort(candidates &save);

seqan::ArgumentParser::ParseResult commands(modifier &options, int argc, char const ** argv);

unsigned generateStandardDegression(candidates &input,int e_value);

int saminput(candidates &save, modifier &options, statisticals &stats);

int devide(candidates &input,candidates &result,modifier &options,statisticals &stats);

int output(std::vector<std::vector<unsigned> > &input,candidates &ref, char *out);

void findCluster(candidates &input, std::vector<std::vector<unsigned> > &dest, candidate::form indel, double coverage);

#endif  // PEMER_LITE_HEADER_H_
