/*=======================================================================
   Name:         Flexbar.cpp
   Authors:      Matthias Dodt and Johannes Roehr

   Description:  Flexbar - flexible barcode and adapter removal
   Version:      2.34
   Copyright:    GPL version 3

   SeqAn lib:    post release 1.4, revision 14066 on May 23, 2013
   TBB   lib:    version 4.0 update 5, stable release June 13, 2012
========================================================================*/

#include "Flexbar.h"
#include "Options.h"
#include "Enums.h"


int main(int argc, const char* argv[]){
	
	using namespace std;
	using namespace flexbar;
	
	using seqan::ArgumentParser;
	
	// performTest();
	
	const string version = "2.35 beta";
	const string date    = "May 25, 2013";
	
	ArgumentParser parser("flexbar");
	
	defineOptionsAndHelp(parser, version, date);
	parseCommandLine(parser, version, argc, argv);
	
	Options o;
	
	loadProgramOptions(o, parser);
	loadBarcodesAndAdapters(o);
	
	startComputation(o);
	printCompletedMessage(o);
	
	return 0;
}

