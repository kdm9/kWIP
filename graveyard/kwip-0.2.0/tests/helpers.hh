/*
 * ============================================================================
 *
 *       Filename:  helpers.hh
 *    Description:  Helpers, looted from libqc++
 *        License:  GPL-3+
 *         Author:  Kevin Murray, spam@kdmurray.id.au
 *
 * ============================================================================
 */

#ifndef HELPERS_HH
#define HELPERS_HH


#include <cstdlib>
#include <cstdio>
#include <sstream>
#include <iostream>
#include <string>
#include <fstream>
#include <vector>


class TestConfig
{
public:
    TestConfig() :
        n_writable_files(0)
    {

    }

    ~TestConfig()
    {
        for (std::string &file: files_to_delete) {
            std::remove(file.c_str());
        }
    }

    std::string
    get_data_file(std::string file)
    {
        std::ostringstream fs;

        if (root.size() == 0) {
            _get_root();
        }

        fs << root << "/data/" << file;
        return fs.str();
    }

    std::string
    get_writable_file(std::string ext, bool keep=true)
    {
        std::ostringstream fs;

        if (root.size() == 0) {
            _get_root();
        }

        fs << root << "/data/" << ++n_writable_files << "." << ext;
        if (!keep) {
            files_to_delete.push_back(fs.str());
        }
        return fs.str();
    }

    static TestConfig *
    get_config()
    {
        return &instance;
    }

protected:
    void
    _get_root()
    {
        const char *envroot = std::getenv("LIBQCPP_DATA_ROOT");

        if (envroot != NULL) {
            root = envroot;
            return;
        }

        // Use full path if we can
        envroot = std::getenv("PWD");
        if (envroot != NULL) {
            root = envroot;
            return;
        }

        // Fall back on './'
        root = ".";
    }

    std::string root;
    std::vector<std::string> files_to_delete;
    size_t n_writable_files;

    static TestConfig instance;
};

extern TestConfig _config;

// returns true on identity, false otherwise
static inline bool
filecmp(const std::string &filename1,
        const std::string &filename2)
{
    char chr1, chr2;
    size_t idx = 0;
    std::ifstream fp1(filename1), fp2(filename2);

    if (!fp1 || !fp2) {
        std::cerr << "Couldn't open files\n";
        return false;
    }

    while (fp1.good() && fp2.good()) {
        chr1 = fp1.get();
        chr2 = fp2.get();
        if (chr1 != chr2) {
            std::cerr << "Chars differ at byte " << idx << " '" << (int)chr1 << "' '" << (int)chr2 << "'\n";
            return false;
        }
        idx++;
    }
    fp1.get();
    fp2.get();
    if (fp1.eof() && fp2.eof()) {
        return true;
    }
    std::cerr << "Not at eof " << fp1.eof() << " " << fp2.eof() << "\n";
    return false;
}

// returns true on identity, false otherwise
static inline bool
filestrcmp(const std::string &filename,
           const std::string &contents)
{
    std::ifstream fp(filename, std::ios_base::binary);
    size_t idx = 0;

    if (!fp) {
        return false;
    }

    while (fp.good() && idx < contents.size()) {
        char chr = fp.get();
        if (chr != contents[idx]) {
            return false;
        }
        idx++;
    }
    if (idx == contents.size() && fp.get() == EOF) {
        return true;
    }
    return false;
}

#endif /* HELPERS_HH */
