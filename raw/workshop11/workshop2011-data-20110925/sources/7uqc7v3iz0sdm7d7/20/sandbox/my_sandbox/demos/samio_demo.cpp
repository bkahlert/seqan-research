// ==========================================================================
//                                 samio_demo
// ==========================================================================
// Copyright (c) 2006-2010, Knut Reinert, FU Berlin
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
// Author: Your Name <your.email@example.net>
// ==========================================================================

#include <iostream>
#include <fstream>
#include <seqan/bam_io.h>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/file.h>      // For printing SeqAn Strings.

using namespace seqan;

int main(int argc, char const ** argv)
{
    std::ofstream file2;
	CharString outfile = "demo.sam";
    file2.open(toCString(outfile), ::std::ios_base::out | ::std::ios_base::app);
    if (!file2.is_open()) {
 		std::cerr << "Could not open output file " << outfile << std::endl;
 	}else{
		typedef StringSet<CharString>      TNameStore;
		typedef NameStoreCache<TNameStore> TNameStoreCache;

		TNameStore refNameStore;
		appendValue(refNameStore, "Chr1");
		TNameStoreCache refNameStoreCache(refNameStore);
		BamIOContext<TNameStore> context(refNameStore, refNameStoreCache);
		int res;
		BamAlignmentRecord record;
		record.qName = "QueryName";
		record.flag = 0;
		record.rId = 0;
		record.pos = 10;
		record.mapQ = 255;
		record.bin = 0;
		resize(record.cigar, 5);
		record.cigar[0].count = 10;
		record.cigar[0].operation = 'M';
		record.cigar[1].count = 1;
		record.cigar[1].operation = 'D';
		record.cigar[2].count = 7;
		record.cigar[2].operation = 'M';
		record.cigar[3].count = 1;
		record.cigar[3].operation = 'I';
		record.cigar[4].count = 3;
		record.cigar[4].operation = 'M';
	
		record.rNextId = 0;
		record.pNext = 0;
		record.tLen = 0;
		record.seq = "ACTG";
		record.qual = "";
		record.tags = "";
		
		res = write2(file2, record, context, Sam());
		// res == 0 on success
		
		std::cout << "Writing report: " << res << std::endl;
	}

	
	return 0;
}