#ifndef _SFILE_MAT_H
#define _SFILE_MAT_H

#ifdef OFFICIAL_MATLAB
#   include "SFile_MAT.Matlab.h"
#else
#   include "SFile_MAT.HDF5.h"
#endif

#endif/*_SFILE_MAT_H*/
