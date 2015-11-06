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
#include <errno.h>
#include <sstream> // IWYU pragma: keep

#include "hashbits.hh"
#include "hashtable.hh"
#include "khmer_exception.hh"
#include "read_parsers.hh"

using namespace std;
using namespace khmer;
using namespace khmer:: read_parsers;

void Hashbits::save(std::string outfilename)
{
    if (!_counts[0]) {
        throw khmer_exception();
    }

    unsigned int save_ksize = _ksize;
    unsigned char save_n_tables = _n_tables;
    unsigned long long save_tablesize;
    unsigned long long save_occupied_bins = _occupied_bins;

    ofstream outfile(outfilename.c_str(), ios::binary);

    outfile.write(SAVED_SIGNATURE, 4);
    unsigned char version = SAVED_FORMAT_VERSION;
    outfile.write((const char *) &version, 1);

    unsigned char ht_type = SAVED_HASHBITS;
    outfile.write((const char *) &ht_type, 1);

    outfile.write((const char *) &save_ksize, sizeof(save_ksize));
    outfile.write((const char *) &save_n_tables, sizeof(save_n_tables));
    outfile.write((const char *) &save_occupied_bins,
                  sizeof(save_occupied_bins));

    for (unsigned int i = 0; i < _n_tables; i++) {
        save_tablesize = _tablesizes[i];
        unsigned long long tablebytes = save_tablesize / 8 + 1;

        outfile.write((const char *) &save_tablesize, sizeof(save_tablesize));

        outfile.write((const char *) _counts[i], tablebytes);
    }
    if (outfile.fail()) {
        throw khmer_file_exception(strerror(errno));
    }
    outfile.close();
}

/**
 * Loads @param infilename into Hashbits, with error checking on
 * file type and file version.  Populates _counts internally.
 */
void Hashbits::load(std::string infilename)
{
    ifstream infile;

    // configure ifstream to raise exceptions for everything.
    infile.exceptions(std::ifstream::failbit | std::ifstream::badbit |
                      std::ifstream::eofbit);

    try {
        infile.open(infilename.c_str(), ios::binary);
    } catch (std::ifstream::failure &e) {
        std::string err;
        if (!infile.is_open()) {
            err = "Cannot open k-mer graph file: " + infilename;
        } else {
            err = "Unknown error in opening file: " + infilename;
        }
        throw khmer_file_exception(err);
    }

    if (_counts) {
        for (unsigned int i = 0; i < _n_tables; i++) {
            delete[] _counts[i];
            _counts[i] = NULL;
        }
        delete[] _counts;
        _counts = NULL;
    }
    _tablesizes.clear();

    try {
        unsigned int save_ksize = 0;
        unsigned char save_n_tables = 0;
        unsigned long long save_tablesize = 0;
        unsigned long long save_occupied_bins = 0;
        char signature[4];
        unsigned char version, ht_type;

        infile.read(signature, 4);
        infile.read((char *) &version, 1);
        infile.read((char *) &ht_type, 1);
        if (!(std::string(signature, 4) == SAVED_SIGNATURE)) {
            std::ostringstream err;
            err << "Does not start with signature for a khmer file: 0x";
            for(size_t i=0; i < 4; ++i) {
                err << std::hex << (int) signature[i];
            }
            err << " Should be: " << SAVED_SIGNATURE;
            throw khmer_file_exception(err.str());
        } else if (!(version == SAVED_FORMAT_VERSION)) {
            std::ostringstream err;
            err << "Incorrect file format version " << (int) version
                << " while reading k-mer graph from " << infilename
                << "; should be " << (int) SAVED_FORMAT_VERSION;
            throw khmer_file_exception(err.str());
        } else if (!(ht_type == SAVED_HASHBITS)) {
            std::ostringstream err;
            err << "Incorrect file format type " << (int) ht_type
                << " while reading k-mer graph from " << infilename;
            throw khmer_file_exception(err.str());
        }

        infile.read((char *) &save_ksize, sizeof(save_ksize));
        infile.read((char *) &save_n_tables, sizeof(save_n_tables));
        infile.read((char *) &save_occupied_bins, sizeof(save_occupied_bins));

        _ksize = (WordLength) save_ksize;
        _n_tables = (unsigned int) save_n_tables;
        _occupied_bins = save_occupied_bins;
        _init_bitstuff();

        _counts = new Byte*[_n_tables];
        for (unsigned int i = 0; i < _n_tables; i++) {
            HashIntoType tablesize;
            unsigned long long tablebytes;

            infile.read((char *) &save_tablesize, sizeof(save_tablesize));

            tablesize = (HashIntoType) save_tablesize;
            _tablesizes.push_back(tablesize);

            tablebytes = tablesize / 8 + 1;
            _counts[i] = new Byte[tablebytes];

            unsigned long long loaded = 0;
            while (loaded != tablebytes) {
                infile.read((char *) _counts[i], tablebytes - loaded);
                loaded += infile.gcount();
            }
        }
        infile.close();
    } catch (std::ifstream::failure &e) {
        std::string err;
        if (infile.eof()) {
            err = "Unexpected end of k-mer graph file: " + infilename;
        } else {
            err = "Error reading from k-mer graph file: " + infilename;
        }
        throw khmer_file_exception(err);
    }
}

void Hashbits::update_from(const Hashbits &other)
{
    if (_ksize != other._ksize) {
        throw khmer_exception("both nodegraphs must have same k size");
    }
    if (_tablesizes != other._tablesizes) {
        throw khmer_exception("both nodegraphs must have same table sizes");
    }
    for (unsigned int table_num = 0; table_num < _n_tables; table_num++) {
        Byte * me = _counts[table_num];
        Byte * ot = other._counts[table_num];
        HashIntoType tablesize = _tablesizes[table_num];
        HashIntoType tablebytes = tablesize / 8 + 1;

        for (HashIntoType index = 0; index < tablebytes; index++) {
            me[index] |= ot[index];     // bitwise or
        }
    }
}

// vim: set sts=2 sw=2:
