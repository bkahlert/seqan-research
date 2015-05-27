// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Iterate over BAM alignment, dump lines in SAM format, followed by the read
// alignment.
//
// USAGE: bam_print_alignments REF.fasta ALIGN.bam
// ==========================================================================

#include <iostream>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/file.h>      // For printing SeqAn Strings.
#include <seqan/align.h>
#include <seqan/bam_io.h>

#if SEQAN_HAS_ZLIB

int main(int argc, char const ** argv)
{
    using namespace seqan;

    // Check command line arguments.
    if (argc != 3)
    {
        std::cerr << "USAGE: bam_print_alignments REF.fasta FILE.bam" << std::endl;
        return 1;
    }

    // Read FASTA file.
    std::cerr << "Reading FASTA " << argv[1] << std::endl;
    StringSet<CharString> refNameStore;
    StringSet<Dna5String> seqs;
    std::fstream inSeq(argv[1], std::ios::binary | std::ios::in);
    if (!inSeq.good())
    {
        std::cerr << "Could not open FASTA file " << argv[1] << std::endl;
        return 1;
    }
    RecordReader<std::fstream, SinglePass<> > reader(inSeq);
    if (read2(refNameStore, seqs, reader, Fasta()) != 0)
        return 1;
    
    // Open BGZF stream.
    std::cerr << "Opening BAM " << argv[2] << std::endl;
    Stream<Bgzf> stream;
    if (!open(stream, argv[2], "r"))
    {
        std::cerr << "[ERROR] Could not open BAM file" << argv[1] << std::endl;
        return 1;
    }

    // Read Header.
    std::cerr << "Reading Header" << std::endl;
    NameStoreCache<StringSet<CharString> > refNameStoreCache(refNameStore);
    BamIOContext<StringSet<CharString> > context(refNameStore, refNameStoreCache);
    BamHeader header;
    if (readRecord(header, context, stream, Bam()) != 0)
    {
        std::cerr << "[ERROR] Could not read header from BAM file!" << std::endl;
        return 1;
    }

    // Stream through file, getting alignment and dumping it.
    std::cerr << "Reading Alignments..." << std::endl;
    Align<Dna5String> align;
    BamAlignmentRecord record;
    while (!atEnd(stream))
    {
        clear(record);
        if (readRecord(record, context, stream, Bam()) != 0)
        {
            std::cerr << "[ERROR] Error reading alignment!" << std::endl;
            return 1;
        }

        if (record.rId == BamAlignmentRecord::INVALID_REFID)
            continue;  // Skip * reference.

        // Convert BAM record to alignment.
        bamRecordToAlignment(align, seqs[record.rId], record);
        // Dump record as SAM and the alignment.
        write2(std::cout, record, context, Sam());
        std::cout << align << std::endl;
    }

    return 0;
}

#else  // #if SEQAN_HAS_ZLIB

int main()
{
    std::cerr << "zlib is required for bam_print_alignment demo." << std::endl;
    
    return 1;
}

#endif  // #if SEQAN_HAS_ZLIB
