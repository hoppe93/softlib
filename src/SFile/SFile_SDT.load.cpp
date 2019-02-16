/**
 * Routines for loading/reading SDT files.
 */

#include <fstream>
#include <sstream>
#include <string>
#include <softlib/SFile_SDT.h>

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
            unsigned int rows, cols;
            char c;

            if (!(sdtfile >> name))
                throw SFileException("Expected name of matrix after '@matrix'.");

            if (!(sdtfile >> rows) || !stdfile.get(c) || c != 'x' || !(sdtfile >> cols))
                throw SFileException("Expected matrix dimensions (RxC) after matrix name.");

            LoadMatrix(name, rows, cols);
        } else if (type == "@string") {
            string name;
            unsigned int length;

            if (!(sdtfile >> name))
                throw SFileException("Expected name of string after '@string'.");

            if (!(sdtfile >> length))
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
void SFile_SDT::LoadMatrix(const string &name, const unsigned int nrows, const unsigned int ncols) {
    slibreal_t **vals = new slibreal_t*[nrows];
    vals[0] = new slibreal_t[nrows*ncols];

    for (unsigned int i = 0; i < nrows; i++) {
        string line;
        getline(sdtfile, line);
        stringstream ss(line);

        if (i > 0)
            vals[i] = vals[i-1];

        for (unsigned int j = 0; j < ncols; j++) {
            if (!(ss >> vals[i][j]))
                throw SFileException("Invalid format of matrix '%s' at index (%u ,%u).", name.c_str(), i, j);
        }
    }

    // TODO add variable
}

/**
 * Load the string with the given name from this point
 * in the SDT file.
 *
 * name:   Name to give to string to load.
 * length: Length of string to load (this parameter is not
 *         used but is present to prevent compiler warnings).
 */
void SFile_SDT::LoadString(const string &name, const unsigned int) {
    string str;
    getline(sdtfile, str);

    // TODO add variable
}

