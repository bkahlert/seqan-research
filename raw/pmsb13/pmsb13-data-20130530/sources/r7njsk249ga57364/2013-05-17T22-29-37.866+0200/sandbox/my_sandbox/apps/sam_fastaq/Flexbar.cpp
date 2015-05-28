/*=======================================================================
   Name:         Flexbar.cpp
   Authors:      Matthias Dodt and Johannes Roehr

   Description:  Flexbar - flexible barcode and adapter removal
   Version:      2.33
   Copyright:    GPL version 3

   SeqAn lib:    post release 1.3.1, revision 13058 on October 9, 2012
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
	
	const string version = "2.33";
	const string date    = "March 18, 2013";
	
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

