/**
 * Implementation of module for generating 
 * Portable Pixmap forMat (PPM) images.
 */

#include <cstdio>
#include <fstream>
#include <iostream>
#include <softlib/ImageGenerator/ImageColormap.h>
#include <softlib/ImageGenerator/ImageGeneratorPPM.h>

using namespace std;


void ImageGeneratorPPM::SaveImage(
    const std::string& name, const slibreal_t *img,
    const unsigned int nrows, const unsigned int ncols,
    ImageColormap *cmap
) {
    ofstream file;
    file.open(name);

    /********************
     * Write PPM header *
     ********************/
    // Signature
    file << "P3" << endl;
    // Image size
    file << nrows << " " << ncols << endl;
    // Maximum color value (0-max)
    file << "255" << endl;

    /***************
     * Write image *
     ***************/
    for (unsigned int i = 0; i < nrows; i++) {
        for (unsigned int j = 0; j < ncols; j++) {
            // NOTE: The image should be output in transposed form!
            ImageColormap::scolor_t clr = cmap->Eval(img[j*nrows + i]);

            file << clr.red << " " << clr.green << " " << clr.blue << " ";
        }
        
        file << endl;
    }

    file.close();
}
