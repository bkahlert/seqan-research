// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2012, Knut Reinert, FU Berlin
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
// Author: Andreas Gogol-Doering <andreas.doering@mdc-berlin.de>
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Storage-related alphabet interface part.  This means both
// construction type (simple, non-simple) and storage size.
// ==========================================================================

#ifndef SEQAN_CORE_INCLUDE_SEQAN_BASIC_ALPHABET_STORAGE_H_
#define SEQAN_CORE_INCLUDE_SEQAN_BASIC_ALPHABET_STORAGE_H_

#include <float.h>

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

template <typename TValue, typename TSpec>
class SimpleType;

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// TODO(holtgrew): Remove Ascii and Unicode alias. Also see #849.
typedef char Ascii;
//typedef unsigned char Byte;  // TODO(holtgrew): Disabling, remove together with Ascii and Unicode with #849
typedef wchar_t Unicode;

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction BitsPerValue
// ----------------------------------------------------------------------------

template <typename TValue>
struct BitsPerValue
{
    static const unsigned VALUE = sizeof(TValue) * 8;
    typedef unsigned Type;
};

template <typename TValue>
struct BitsPerValue<TValue const> : public BitsPerValue<TValue>
{};

// template <typename TValue>
// const unsigned BitsPerValue<TValue>::VALUE = ;

// ----------------------------------------------------------------------------
// Metafunction ValueSize
// ----------------------------------------------------------------------------

// TODO(holtgrew): Enable default implementation only for integers? Move to adapt builtins?

template <typename T>
struct ValueSize
{
    typedef __uint64  Type;
    static const Type VALUE = (__uint64(1) << BitsPerValue<T>::VALUE);
};

template <typename T>
struct ValueSize<T const> : ValueSize<T>
{};

// TODO(holtgrew): Use static assertion to make sure that ValueSize is never called on floating point numbers? Include assertion for __int64 and __uint64?

template <>
struct ValueSize<__int64>
{
    typedef __uint64  Type;
    static const Type VALUE = 0;
};

template <>
struct ValueSize<__uint64>
{
    typedef __uint64  Type;
    static const Type VALUE = 0;
};

template <>
struct ValueSize<double>
{
    typedef __uint64  Type;
    static const Type VALUE = 0;
};

template <>
struct ValueSize<float>
{
    typedef __uint64  Type;
    static const Type VALUE = 0;
};


// The internal value size is used for alphabets with piggyback qualities,
// for example Dna5Q.  Here, the public value size is 5 but the internal
// value size is 256.  

template <typename TValue> 
struct InternalValueSize_
        : public ValueSize<TValue>
{};

// ----------------------------------------------------------------------------
// Metafunction BytesPerValue
// ----------------------------------------------------------------------------

/**
.Metafunction.BytesPerValue:
..cat:Basic
..summary:Number of bytes needed to store a value.
..signature:BytesPerValue<T>::VALUE
..param.T:A class.
..returns.param.VALUE:Number of bytes needed to store $T$.
...default:$BitsPerValue / 8$, rounded up. For built-in types, this is the same as $sizeof(T)$.
..see:Metafunction.ValueSize
..see:Metafunction.BitsPerValue
..include:seqan/basic.h
*/

template <typename TValue>
struct BytesPerValue
{
    enum { VALUE = (BitsPerValue<TValue>::VALUE + 7) / 8 };
};

// ----------------------------------------------------------------------------
// Metafunction IntegralForValue
// ----------------------------------------------------------------------------

/**
.Metafunction.IntegralForValue:
..cat:Basic
..summary:Returns an itegral type that provides sufficient space to store a value.
..signature:IntegralForValue<T>::Type
..param.T:A class.
..returns.param.Type:An integral type that can store $T$ values.
..remarks:The type is the smallest unsigned integral type that has a size of at least @Metafunction.BytesPerValue@ bytes.
...tableheader:bytes|integral type
...table:1|$unsigned char$
...table:2|$unsigned short$
...table:3|$unsigned int$
...table:4|$unsigned int$
...table:5 and above|$__int64$
..remarks:Note that the returned integral type cannot store $T$ values, if $T$ takes more than 8 bytes, 
    since there exists no integral type that provides sufficient space to store types of this size.
..see:Metafunction.ValueSize
..see:Metafunction.BitsPerValue
..see:Metafunction.BytesPerValue
..include:seqan/basic.h
*/

template <int SIZE>
struct IntegralForValueImpl_
{
    typedef __int64 Type;
};

// TODO(holtgrew): Switch to __uint8, __uint16, __uint32?
template <>
struct IntegralForValueImpl_<1>
{
    typedef unsigned char Type;
};

template <>
struct IntegralForValueImpl_<2>
{
    typedef unsigned short Type;
};

template <>
struct IntegralForValueImpl_<3>
{
    typedef unsigned int Type;
};

template <>
struct IntegralForValueImpl_<4>
{
    typedef unsigned int Type;
};

template <typename TValue>
struct IntegralForValue : IntegralForValueImpl_<BytesPerValue<TValue>::VALUE>
{};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function ordValue()
// ----------------------------------------------------------------------------

// TODO(holtgrew): Enable only for integers, move to adapt builtins?

template <typename TValue>
inline typename ValueSize<TValue>::Type
ordValue(TValue const & c)
{
	return convert<unsigned>(static_cast<typename MakeUnsigned_<TValue>::Type const &>(c));
}

// The internal ord value is used for alphabets with piggyback qualities.

template <typename TValue>
inline typename ValueSize<TValue>::Type
_internalOrdValue(TValue const & c)
{
	return ordValue(c);
}

// ----------------------------------------------------------------------------
// Function valueSize<T>()
// ----------------------------------------------------------------------------

template <typename T>
inline typename ValueSize<T>::Type
valueSize()
{
    return +ValueSize<T>::VALUE;
}

}  // namespace seqan

#endif  // #ifndef SEQAN_CORE_INCLUDE_SEQAN_BASIC_ALPHABET_STORAGE_H_
