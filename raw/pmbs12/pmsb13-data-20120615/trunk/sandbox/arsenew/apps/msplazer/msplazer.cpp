// ==========================================================================
//                                  msplazer
// ==========================================================================
// Copyright (c) 2011, Kathrin Trappe, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
//
// ==========================================================================
// Author: Kathrin Trappe <ktrappe@inf.fu-berlin.de>
// ==========================================================================

#include <seqan/basic.h>
#include <seqan/sequence.h>
//#include <seqan/graph_algorithms.h>

#include <seqan/misc/misc_cmdparser.h>

#include "../../../../core/apps/stellar/stellar.h"
//#include "../../../../core/apps/stellar/stellar_output.h"

//#include "stellar_routines.h"
//#include "create_stellarmatches_from_file.h"
#include "msplazer.h"
//#include "msplazer_algorithms.h"
//#include "msplazer_out.h"
#include "msplazer_main.h"

using namespace seqan;

// Program entry point
int main(int argc, char const ** argv)
{

	double start = sysTime();
	// command line parsing
	CommandLineParser parser("msplazer");

	_setParser(parser);
	if (!parse(parser, argc, argv)) {
		if (isSetShort(parser, 'h') || isSetShort(parser, 'v')) return 0; 
		shortHelp(parser, std::cerr);
		return 1;
	}

	/////////////////////////////////////////////////////////////////////
	//Initialisation
	MSplazerOptions msplazerOptions = MSplazerOptions();

	//Usage via fragmentStore
	/*Stellar can now handle FragmentStore input which is based on an Index using a StringSet based on Owner<ConcatDirect<>>. 
	  However, StellarMatch positions are wrong when processing multiple reads (i.e. after read 1).
	  Maybe bc. of consecutive positions due to the ConcatDirect format?
	*/
	/*
	   FragmentStore<void> fragments;
	   bool success = loadReads(fragments, queryFilename);
	   if(!success)
	   cout << "Unable to open file from " << queryFilename << endl;
	   for(unsigned i = 0; i < length(fragments.readSeqStore); ++i)
	   cout << "Read " << i << " : " << fragments.readSeqStore[i] << endl;
	   */	

	//stellar parameters
	//StellarOptions
	StellarOptions stellarOptions = StellarOptions();

	if (!_parseOptions(parser, stellarOptions, msplazerOptions)) {
		return 1;
	}

	msplazerOptions.databaseFile = stellarOptions.databaseFile;
	msplazerOptions.queryFile = stellarOptions.queryFile;
	
	_writeFileNames(stellarOptions);
	_writeParams(msplazerOptions);

    //msplazer wrapper function
	msplazer(stellarOptions, msplazerOptions);
	//msplazer(msplazerOptions) // If Stellar run is not necessary due to external match input


	std::cout << "TIME all " <<	(sysTime() - start) << "s" << std::endl;

	
	return 0;
}
