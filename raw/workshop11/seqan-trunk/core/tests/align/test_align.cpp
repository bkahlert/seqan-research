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
// Author: Hannes Hauswedell <hauswedell@mi.fu-berlin.de>
// Author: Andreas Gogol-Doering <andreas.doering@mdc-berlin.de> 
// ==========================================================================

#include <seqan/basic.h>  
#include <seqan/file.h>   

#include "test_align_align.h"
#include "test_align_gaps.h"
#include "test_align_myers.h"
#include "test_align_matrix.h"
#include "test_local_align.h"

SEQAN_BEGIN_TESTSUITE("test_align") {
    // Call tests from test_align_align.cpp.
    SEQAN_CALL_TEST(test_align_align_char_string_array_gaps);
    SEQAN_CALL_TEST(test_align_align_dna_array_gaps);
    SEQAN_CALL_TEST(test_align_align_char_string_sumlist_gaps);
    SEQAN_CALL_TEST(test_align_align_dna_sumlist_gaps);

    SEQAN_CALL_TEST(test_align_gaps_base_char_string_array_gaps);
    SEQAN_CALL_TEST(test_align_gaps_base_char_string_sumlist_gaps);
    SEQAN_CALL_TEST(test_align_gaps_test_gaps_iterator);
    SEQAN_CALL_TEST(test_align_gaps_test_gap_manipulation_char_string_array_gaps);
    SEQAN_CALL_TEST(test_align_gaps_test_gap_manipulation_char_string_sumlist_gaps);
    SEQAN_CALL_TEST(test_align_gaps_test_trailing_gaps_char_string_array_gaps);
    SEQAN_CALL_TEST(test_align_gaps_test_count_characters_char_string_array_gaps);

    SEQAN_CALL_TEST(testLocalAlign);
    SEQAN_CALL_TEST(testLocalAlign2);
    SEQAN_CALL_TEST(testBandedLocalAlign);

    SEQAN_CALL_TEST(test_align_myers_test_short);
    SEQAN_CALL_TEST(test_align_myers_test_long);
    SEQAN_CALL_TEST(test_align_hirschberger);

    SEQAN_CALL_TEST(test_align_matrix);

    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/align/align_algorithms.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/align/align_base.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/align/align_cols_base.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/align/align_dynprog.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/align/align_hirschberg.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/align/align_iterator_base.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/align/align_local_dynprog.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/align/align_myers.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/align/align_trace.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/align/gaps_array.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/align/gaps_base.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/align/gaps_iterator_base.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/align/gaps_sumlist.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/align/hirschberg_set.h");
    SEQAN_VERIFY_CHECKPOINTS("core/include/seqan/align/matrix_base.h");
}
SEQAN_END_TESTSUITE
