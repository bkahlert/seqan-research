// ==========================================================================
//                      create_stellarmatches_from_file
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
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Kathrin Trappe <kathrin.trappe@fu-berlin.de>
// ==========================================================================

#ifndef SANDBOX_KTRAPPE_APPS_MSPLAZER_MSPLAZER_OUT_H_
#define SANDBOX_KTRAPPE_APPS_MSPLAZER_MSPLAZER_OUT_H_

#include <iostream>
#include <fstream>

#include <seqan/file.h>      // For printing SeqAn Strings.

#include "../../../../core/apps/stellar/stellar.h"
#include "msplazer.h"

using namespace seqan;

//DotWriting call for read graphs
template < typename TSequence, typename TId, typename TMSplazerChain >
void _writeDotfiles(StringSet <QueryMatches< StellarMatch <TSequence, TId > > > & stellarMatches,
		StringSet<TSequence> const & queries,
		String< TMSplazerChain > & queryChains,
		MSplazerOptions const & msplazerOptions){

	for(unsigned i = 0; i < length(queryChains); ++i){
		if(!queryChains[i].isEmpty){
		//if(length(stellarMatches[i].matches) > 45){
			//String< char > fn = toString(msplazerOptions.outDir) + "read" + toString(i+1) + '_' + toString(msplazerOptions.queryFile) + '_' + toString(msplazerOptions.databaseFile) + ".dot";
			String< char > fn = toString(msplazerOptions.outDir) + "read" + toString(i+1) + '_' + toString(msplazerOptions.jobName) + ".dot";
			//std::cerr  << fn << std::endl;
			FILE* strmWrite = fopen(toCString(fn), "w");
			//write(strmWrite, queryChains[i].graph, DotDrawing());
			write(strmWrite, queryChains[i], stellarMatches[i].matches, length(queries[i]), DotDrawingMSplazer());

			fclose(strmWrite);
			//std::cerr << " completed dot write in: " << i << std::endl;
		}
	}
}

//Breakpoint writing call
template < typename TBreakpoint >
void _writeGlobalBreakpoints(String<TBreakpoint> const & globalBreakpoints,
		MSplazerOptions const & msplazerOptions,
		unsigned const & support,
		String<char> const & fileName){

	String< char > fn = toString(msplazerOptions.outDir) + toString(msplazerOptions.jobName) + toString(fileName);
	//std::cerr  << fn << std::endl;
	FILE* strmWrite = fopen(toCString(fn), "w");
	_streamWrite(strmWrite, "Global breakpoints found on best MSplazerChains sorted according to genome Ids\n");
	_streamWrite(strmWrite, "Database file: ");
	_streamWrite(strmWrite, msplazerOptions.databaseFile);
	_streamPut(strmWrite, '\n');
	_streamPut(strmWrite, BreakpointFileHeader());
	//print bps
	for(unsigned i = 0; i < length(globalBreakpoints); ++i){
		//if(globalBreakpoints[i].svtype != "none" && globalBreakpoints[i].support >= msplazerOptions.support)
		if(globalBreakpoints[i].svtype != "none" && globalBreakpoints[i].support >= support)
			_streamWrite(strmWrite, globalBreakpoints[i], i);
		//_streamPut(strmWrite, '\n');
	
	}
	fclose(strmWrite);
	std::cout << " completed writing " << fileName << std::endl;
}

///////////////////////////////////////////////////////////////////////////////
// Writes parameters from options object to std::cout
template<typename TOptions>
void
_writeParams(TOptions & options) {
//IOREV _notio_
	// Output user specified parameters
	if(options.outDir != "")
		std::cout << "Output directory        : " << options.outDir << std::endl;
	if(options.jobName != "")
		std::cout << "Job name        : " << options.jobName << std::endl;
	std::cout << "Thresholds:" << std::endl;
	std::cout << "  overlap threshold (oth)          : " << options.simThresh << std::endl;
	std::cout << "  gap threshold (gth)              : " << options.gapThresh << std::endl;
	std::cout << "  inital gap threshold (ith)       : " << options.initGapThresh << std::endl;
	std::cout << "Penalties:" << std::endl;
	std::cout << "  translocation penalty (tp)       : " << options.diffDBPen << std::endl;
	std::cout << "  inversion penalty (ip)           : " << options.diffStrandPen << std::endl;
	std::cout << "  order penalty (op)               : " << options.diffOrderPen << std::endl;

	std::cout << "  required read support (st)       : " << options.support << std::endl;

}

#endif  // #ifndef SANDBOX_MY_SANDBOX_APPS_MSPLAZER_MSPLAZER_OUT_H_
