/**
 * Routines for loading/reading SDT files.
 */

#include <fstream>
#include <sstream>
#include <string>
#include <softlib/SFileException.h>
#include <softlib/SFile_SDT.h>

using namespace std;

/**
 * Load the entire SDT file and parse its contents.
 */
void SFile_SDT::Load() {
    string line, type;

    while (getline(sdtfile, line)) {
        if (line.empty())
            continue;

        stringstream ss(line);
        ss >> type;

        if (type == "@matrix") {
            string name;
            sfilesize_t rows, cols;

            if (!(ss >> name))
                throw SFileException("Expected name of matrix after '@matrix'.");

            if (!(ss >> rows) || !(ss >> cols))
                throw SFileException("Expected matrix dimensions (R C) after matrix name.");

            LoadMatrix(name, rows, cols);
        } else if (type == "@string") {
            string name;
            sfilesize_t length;

            if (!(ss >> name))
                throw SFileException("Expected name of string after '@string'.");

            if (!(ss >> length))
                throw SFileException("Expected string length after string name.");

            LoadString(name, length);
        } else
            throw SFileException("Expected data type specifier, but found token: '%s'.", type.c_str());
    }
}

/**
 * Load the matrix with the given name from this point
 * in the SDT file.
 *
 * name:  Name to give to matrix to load.
 * nrows: Number of rows of matrix.
 * ncols: Number of columns of matrix.
 */
void SFile_SDT::LoadMatrix(const string &name, const sfilesize_t nrows, const sfilesize_t ncols) {
    slibreal_t **vals = new slibreal_t*[nrows];
    vals[0] = new slibreal_t[nrows*ncols];

    for (sfilesize_t i = 0; i < nrows; i++) {
        string line;
        getline(sdtfile, line);
        stringstream ss(line);

        if (i > 0)
            vals[i] = vals[i-1] + ncols;

        for (sfilesize_t j = 0; j < ncols; j++) {
            if (!(ss >> vals[i][j]))
                throw SFileException("Invalid format of matrix '%s' at index (%u, %u).", name.c_str(), i, j);
        }
    }

    // Add variable
    SFile_SDT::matrix mat;
    mat.name = name;
    mat.value = vals;
    mat.nrows = nrows;
    mat.ncols = ncols;

    matrices.emplace(name, mat);
}

/**
 * Load the string with the given name from this point
 * in the SDT file.
 *
 * name:   Name to give to string to load.
 * length: Length of string to load (this parameter is not
 *         used but is present to prevent compiler warnings).
 */
void SFile_SDT::LoadString(const string &name, const sfilesize_t) {
    string line, value;
    int c;
    getline(sdtfile, line);

    stringstream ss(line);
    while (ss >> c)
        value += (char)c;

    // Add variable
    SFile_SDT::str str;
    str.name = name;
    str.value = value;

    strings.emplace(name, str);
}

