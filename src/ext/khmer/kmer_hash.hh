/*
This file is part of khmer, https://github.com/dib-lab/khmer/, and is
Copyright (C) 2010-2015, Michigan State University.
Copyright (C) 2015, The Regents of the University of California.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    * Redistributions in binary form must reproduce the above
      copyright notice, this list of conditions and the following
      disclaimer in the documentation and/or other materials provided
      with the distribution.

    * Neither the name of the Michigan State University nor the names
      of its contributors may be used to endorse or promote products
      derived from this software without specific prior written
      permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
LICENSE (END)

Contact: khmer-project@idyll.org
*/
#ifndef KMER_HASH_HH
#define KMER_HASH_HH

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>

#include "khmer.hh"

// test validity
#ifdef KHMER_EXTRA_SANITY_CHECKS
#   define is_valid_dna(ch) ((toupper(ch)) == 'A' || (toupper(ch)) == 'C' || \
			    (toupper(ch)) == 'G' || (toupper(ch)) == 'T')
#else
// NOTE: Assumes data is already sanitized as it should be by parsers.
//	     This assumption eliminates 4 function calls.
#   define is_valid_dna(ch) ((ch) == 'A' || (ch) == 'C' || \
			     (ch) == 'G' || (ch) == 'T')
#endif

// bit representation of A/T/C/G.
#ifdef KHMER_EXTRA_SANITY_CHECKS
#   define twobit_repr(ch) ((toupper(ch)) == 'A' ? 0LL : \
			    (toupper(ch)) == 'T' ? 1LL : \
			    (toupper(ch)) == 'C' ? 2LL : 3LL)
#else
// NOTE: Assumes data is already sanitized as it should be by parsers.
//	     This assumption eliminates 4 function calls.
#   define twobit_repr(ch) ((ch) == 'A' ? 0LL : \
			    (ch) == 'T' ? 1LL : \
			    (ch) == 'C' ? 2LL : 3LL)
#endif

#define revtwobit_repr(n) ((n) == 0 ? 'A' : \
                           (n) == 1 ? 'T' : \
                           (n) == 2 ? 'C' : 'G')

#ifdef KHMER_EXTRA_SANITY_CHECKS
#   define twobit_comp(ch) ((toupper(ch)) == 'A' ? 1LL : \
			    (toupper(ch)) == 'T' ? 0LL : \
			    (toupper(ch)) == 'C' ? 3LL : 2LL)
#else
// NOTE: Assumes data is already sanitized as it should be by parsers.
//	     This assumption eliminates 4 function calls.
#   define twobit_comp(ch) ((ch) == 'A' ? 1LL : \
			    (ch) == 'T' ? 0LL : \
			    (ch) == 'C' ? 3LL : 2LL)
#endif

// choose wisely between forward and rev comp.
#ifndef NO_UNIQUE_RC
#define uniqify_rc(f, r) ((f) < (r) ? (f) : (r))
#else
#define uniqify_rc(f,r)(f)
#endif


namespace khmer
{
// two-way hash functions.
HashIntoType _hash(const char * kmer, const WordLength k);
HashIntoType _hash(const char * kmer, const WordLength k,
                   HashIntoType& h, HashIntoType& r);
HashIntoType _hash_forward(const char * kmer, WordLength k);

std::string _revhash(HashIntoType hash, WordLength k);

// two-way hash functions, MurmurHash3.
HashIntoType _hash_murmur(const std::string& kmer);
HashIntoType _hash_murmur(const std::string& kmer,
                          HashIntoType& h, HashIntoType& r);
HashIntoType _hash_murmur_forward(const std::string& kmer);

/**
 * \class Kmer
 *
 * \brief Hold the hash values corresponding to a single k-mer.
 *
 * This class stores the forward, reverse complement, and
 * uniqified hash values for a given k-mer. It also defines
 * some basic operators and a utility function for getting
 * the string representation of the sequence. This is meant
 * to replace the original inelegant macros used for hashing.
 *
 * \author Camille Scott
 *
 * Contact: camille.scott.w@gmail.com
 *
 */
class Kmer
{

public:

    /// The forward hash
    HashIntoType kmer_f;
    /// The reverse (complement) hash
    HashIntoType kmer_r;
    /// The uniqified hash
    HashIntoType kmer_u;

    /** @param[in]   f forward hash.
     *  @param[in]   r reverse (complement) hash.
     *  @param[in]   u uniqified hash.
     */
    Kmer(HashIntoType f, HashIntoType r, HashIntoType u)
    {
        kmer_f = f;
        kmer_r = r;
        kmer_u = u;
    }

    /// @warning The default constructor builds an invalid k-mer.
    Kmer()
    {
        kmer_f = kmer_r = kmer_u = 0;
    }

    /// Allows complete backwards compatibility
    operator HashIntoType() const
    {
        return kmer_u;
    }

    bool operator< (const Kmer &other) const
    {
        return kmer_u < other.kmer_u;
    }

    std::string get_string_rep(WordLength K)
    {
        return _revhash(kmer_u, K);
    }

};

/**
 * \class KmerFactory
 *
 * \brief Build complete Kmer objects.
 *
 * The KmerFactory is a simple construct to emit complete
 * Kmer objects. The design decision leading to this class
 * stems from the issue of overloading the Kmer constructor
 * while also giving it a K size: you get ambiguous signatures
 * between the (kmer_u, K) and (kmer_f, kmer_r) cases. This
 * implementation also allows a logical architecture wherein
 * KmerIterator and Hashtable inherit directly from KmerFactory,
 * extending the latter's purpose of "create k-mers" to
 * "emit k-mers from a sequence" and "create and store k-mers".
 *
 * \author Camille Scott
 *
 * Contact: camille.scott.w@gmail.com
 *
 */
class KmerFactory
{
protected:
    WordLength _ksize;

public:

    KmerFactory(WordLength K): _ksize(K) {}

    /** @param[in]  kmer_u Uniqified hash value.
     *  @return A complete Kmer object.
     */
    Kmer build_kmer(HashIntoType kmer_u)
    {
        HashIntoType kmer_f, kmer_r;
        std:: string kmer_s = _revhash(kmer_u, _ksize);
        _hash(kmer_s.c_str(), _ksize, kmer_f, kmer_r);
        return Kmer(kmer_f, kmer_r, kmer_u);
    }

    /** Call the uniqify function and build a complete Kmer.
     *
     *  @param[in]  kmer_f Forward hash value.
     *  @param[in]  kmer_r Reverse complement hash value.
     *  @return A complete Kmer object.
     */
    Kmer build_kmer(HashIntoType kmer_f, HashIntoType kmer_r)
    {
        HashIntoType kmer_u = uniqify_rc(kmer_f, kmer_r);
        return Kmer(kmer_f, kmer_r, kmer_u);
    }

    /** Hash the given sequence and call the uniqify function
     *  on its results to build a complete Kmer.
     *
     *  @param[in]  kmer_s String representation of a k-mer.
     *  @return A complete Kmer object hashed from the given string.
     */
    Kmer build_kmer(std::string kmer_s)
    {
        HashIntoType kmer_f, kmer_r, kmer_u;
        kmer_u = _hash(kmer_s.c_str(), _ksize, kmer_f, kmer_r);
        return Kmer(kmer_f, kmer_r, kmer_u);
    }

    /** Hash the given sequence and call the uniqify function
     *  on its results to build a complete Kmer.
     *
     *  @param[in]  kmer_c The character array representation of a k-mer.
     *  @return A complete Kmer object hashed from the given char array.
     */
    Kmer build_kmer(const char * kmer_c)
    {
        HashIntoType kmer_f, kmer_r, kmer_u;
        kmer_u = _hash(kmer_c, _ksize, kmer_f, kmer_r);
        return Kmer(kmer_f, kmer_r, kmer_u);
    }
};

/**
 * \class KmerIterator
 *
 * \brief Emit Kmer objects generated from the given sequence.
 *
 * Given a string \f$S\f$ and a length \f$K > 0\f$, we define
 * the k-mers of \f$S\f$ as the set \f$S_{i..i+K} \forall i \in \{0..|S|-K+1\}\f$,
 * where \f$|S|\f$ is the length and \f$S_{j..k}\f$ is the half-open
 * substring starting at \f$j\f$ and terminating at \f$k\f$.
 *
 * KmerIterator mimics a python-style generator function which
 * emits the k-mers of the given sequence, in order, as Kmer objects.
 *
 * @warning This is not actually a valid C++ iterator, though it is close.
 *
 * \author Camille Scott
 *
 * Contact: camille.scott.w@gmail.com
 *
 */
class KmerIterator: public KmerFactory
{
protected:
    const char * _seq;

    HashIntoType _kmer_f, _kmer_r;
    HashIntoType bitmask;
    unsigned int _nbits_sub_1;
    unsigned int index;
    size_t length;
    bool initialized;
public:
    KmerIterator(const char * seq, unsigned char k);

    /** @param[in]  f The forward hash value.
     *  @param[in]  r The reverse complement hash value.
     *  @return The first Kmer of the sequence.
     */
    Kmer first(HashIntoType& f, HashIntoType& r);

    /** @param[in]  f The current forward hash value
     *  @param[in]  r The current reverse complement hash value
     *  @return The next Kmer in the sequence
     */
    Kmer next(HashIntoType& f, HashIntoType& r);

    Kmer first()
    {
        return first(_kmer_f, _kmer_r);
    }

    Kmer next()
    {
        return next(_kmer_f, _kmer_r);
    }

    /// @return Whether or not the iterator has completed.
    bool done()
    {
        return index >= length;
    }

    unsigned int get_start_pos() const
    {
        return index - _ksize;
    }

    unsigned int get_end_pos() const
    {
        return index;
    }
}; // class KmerIterator

}

#endif // KMER_HASH_HH
