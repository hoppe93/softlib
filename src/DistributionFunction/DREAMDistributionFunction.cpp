/**
 * Interface for loading DREAM distribution functions.
 *
 * The 'Disruption Runaway Electron Analysis Model' (DREAM) is a 3D transport
 * solver developed for the study of runaway electrons during tokamak
 * disruptions. This module allows the user to load distribution functions from
 * DREAM output files.
 */

#include <string>
#include "softlib/DistributionFunction/DREAMDistributionFunction.h"
#include "softlib/SFile.h"


using namespace std;

/**
 * Constructor.
 *
 * filename:    Name of file to load distribution function from.
 * mf:          Magnetic field in which the distribution function lives.
 * distname:    Name of distribution function unknown to load (default: 'f_re' if existing, else 'f_hot').
 * logarithmic: If 'true', interpolates in the logarithm of f instead of f directly.
 * interptype:  Type of interpolation method to use for the momentum distributions.
 */
DREAMDistributionFunction::DREAMDistributionFunction(
    const string &filename, MagneticField2D *mf, const string &distname,
    bool logarithmic, int interptype
) {
    Load(filename, mf, distname, logarithmic, interptype);
}

/**
 * Destructor.
 */
DREAMDistributionFunction::~DREAMDistributionFunction() {
    delete [] this->rawdata->f;
    delete [] this->rawdata->r;
    delete [] this->rawdata->p;
    delete [] this->rawdata->xi;
    delete this->rawdata;
}

/**
 * Load a distribution function from the specified DREAM output file.
 *
 * filename:    Name of file to load distribution function from.
 * mf:          Magnetic field in which the distribution function lives.
 * distname:    Name of distribution function unknown to load (default: 'f_re' if existing, else 'f_hot').
 * logarithmic: If 'true', interpolates in the logarithm of f instead of f directly.
 * interptype:  Type of interpolation method to use for the momentum distributions.
 */
void DREAMDistributionFunction::Load(
    const string &filename, MagneticField2D *mf, const string &distname,
    bool logarithmic, int interptype
) {
    this->rawdata = this->__Load(filename, mf, distname);
    unsigned int
        nr = this->rawdata->nr,
        N = this->rawdata->np*this->rawdata->nxi;

    if (logarithmic) {
        for (unsigned int i = 0; i < nr; i++) {
            this->InsertMomentumSpaceDistributionLog(
                this->rawdata->np, this->rawdata->nxi, this->rawdata->r[i],
                this->rawdata->p, this->rawdata->xi, this->rawdata->f + i*N,
                interptype
            );
        }
    } else {
        for (unsigned int i = 0; i < nr; i++) {
            this->InsertMomentumSpaceDistribution(
                this->rawdata->np, this->rawdata->nxi, this->rawdata->r[i],
                this->rawdata->p, this->rawdata->xi, this->rawdata->f + i*N,
                interptype
            );
        }
    }
}

/**
 * Internal routine for loading a distribution function from a
 * DREAM output file.
 */
struct DREAMDistributionFunction::dreamdf_data *DREAMDistributionFunction::__Load(
    const string &filename, MagneticField2D *mf, const string &distname
) {
    SFile *sf = SFile::Create(filename, SFILE_MODE_READ);

    // First, determine the name of the distribution function.
    string dname, momgridname;
    if (distname == "") {   // Automatically detect distribution name
        if (sf->HasVariable("/eqsys/f_re"))
            dname = "/eqsys/f_re";
        else if (sf->HasVariable("/eqsys/f_hot"))
            dname = "/eqsys/f_hot";
        else
            throw SOFTLibException("DREAMDistributionFunction: No distribution function found in file '%s'. Try to explicitly give the name of the distribution function.", filename.c_str());
    } else if (distname[0] != '/') {    // Distribution name specified (but not full path)
        dname = "/eqsys/" + distname;
    } else
        dname = distname;

    if (!sf->HasVariable(dname))
        throw SOFTLibException("DREAMDistributionFunction: No distribution function found in file '%s'. Try to explicitly give the name of the distribution function.", filename.c_str());

    size_t l = dname.length();
    if (dname.substr(l-4) == "f_re")
        momgridname = "/grid/runaway";
    else if (dname.substr(l-5) == "f_hot")
        momgridname = "/grid/hottail";
    else
        throw SOFTLibException("DREAMDistributionFunction: Unrecognized type of distribution function specified: '%s'. Unable to determine which momentum grid to use.", dname.c_str());

    struct dreamdf_data *dat = new struct dreamdf_data;
    dat->distname = dname;

    /////////////////////////////
    // LOAD DATA
    /////////////////////////////
    sfilesize_t fsize[4], ndims, nr, np, nxi, nt;
    // f
    double *tf = sf->GetMultiArray_linear(dname, 4, ndims, fsize);
    double *tt = sf->GetDoubles1D("/grid/t", &nt);
    double *tr = sf->GetDoubles1D("/grid/r", &nr);
    double *tp = sf->GetDoubles1D(momgridname + "/p1", &np);
    double *tx = sf->GetDoubles1D(momgridname + "/p2", &nxi);
    
    if (fsize[1] != nr || fsize[2] != nxi || fsize[3] != np)
        throw SOFTLibException(
            "DREAMDistributionFunction: Invalid dimensions of distribution functions. The grid has size %llux%llux%llu, but f has size %llux%llux%llu.",
            nr, nxi, np, fsize[1], fsize[2], fsize[3]
        );

    sfilesize_t TIMEINDEX = (nt-1)*nr*np*nxi;
    dat->nr  = nr;
    dat->np  = np;
    dat->nxi = nxi;

    // Convert to slibreal_t (if necessary)
    if (std::is_same<slibreal_t, double>::value) {
        dat->f  = tf + TIMEINDEX;
        dat->r  = tr;
        dat->p  = tp;
        dat->xi = tx;
    } else {
        dat->f  = new slibreal_t[nr*nxi*np];
        dat->r  = new slibreal_t[nr];
        dat->p  = new slibreal_t[np];
        dat->xi = new slibreal_t[nxi];

        for (sfilesize_t i = 0; i < nr*nxi*np; i++)
            dat->f[i]  = (slibreal_t)(tf+TIMEINDEX)[i];
        for (sfilesize_t i = 0; i < nr; i++)
            dat->r[i]  = (slibreal_t)tr[i];
        for (sfilesize_t i = 0; i < np; i++)
            dat->p[i]  = (slibreal_t)tp[i];
        for (sfilesize_t i = 0; i < nxi; i++)
            dat->xi[i] = (slibreal_t)tx[i];

        delete [] tx;
        delete [] tp;
        delete [] tr;
        delete [] tf;
    }

    delete [] tt;

    // Shift radial grid (convert from normalized minor radius -> major radius)
    for (sfilesize_t i = 0; i < nr; i++) {
        dat->r[i] += mf->GetMagneticAxisR();

        if (i > 0 && dat->r[i] <= dat->r[i-1])
            throw SOFTLibException("DREAMDistributionFunction: The radial grid in the DREAM distribution function is not strictly increasing.");
    }

    return dat;
}

