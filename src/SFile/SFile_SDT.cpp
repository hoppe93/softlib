/**
 * Implementation of the 'Self-Descriptive Text' format,
 * designed specifically for SOFT.
 *
 * The original motivation for this format was to allow
 * SOFT to read magnetic field data exported from systems
 * without any capabilities of generating HDF5 or MAT
 * files (such as those at ASDEX-U).
 *
 * SPECIFICATION
 * An SDT file consists of blocks of variables, separated by
 * empty lines. All data is given as ASCII text.
 *
 * A variable definition consists of two or more lines. The
 * first line defines the size and name of the variable, while
 * the following lines contain the data of the variable. The
 * example below illustrates how a 3-by-2 matrix named 'abc'
 * should be specified:
 *
 * @matrix abc 3x2
 * 1 2
 * 3 4
 * 5 6
 *
 * The available types are
 *   
 *   @matrix   - Matrix (2-dimensional list of double-precision 
 *               floating point values)
 *   @string   - String
 */

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <softlib/SFile.h>
#include <softlib/SFile_SDT.h>

using namespace std;

/**
 * Open SDT file.
 *
 * filename: Name of file to open.
 * openmode: Mode to open the file in (read/write etc.)
 */
void SFile_SDT::Open(const string &filename, enum sfile_mode openmode) {
    mode = openmode;

    switch (mode) {
        case SFILE_MODE_READ:
            sdtfile.open(filename, ios::in);
            Load();
            Close();
            break;
        case SFILE_MODE_UPDATE:
            sdtfile.open(filename, ios::in | ios::out);
            break;
        case SFILE_MODE_WRITE:
            sdtfile.open(filename, ios::out | ios::trunc);
            break;
        default:
            throw SOFTLibException("Unrecognized option for opening the SDT file: %d.", mode);
    }
}

/**
 * Closes the SDT file.
 */
void SFile_SDT::Close() {
    sdtfile.close();
}

