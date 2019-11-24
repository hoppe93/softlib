#ifndef _IMAGE_COLORMAP_H
#define _IMAGE_COLORMAP_H

#include <cmath>
#include <softlib/config.h>

class ImageColormap {
public:
    typedef struct {
        unsigned int red, green, blue;
    } scolor_t;

private:
    scolor_t *colors;
    unsigned int ncolors;
    slibreal_t vmin, vmax;

public:
    ImageColormap(
        const scolor_t *colors, unsigned int ncolors,
        slibreal_t vmin, slibreal_t vmax
    )
        : ncolors(ncolors), vmin(vmin), vmax(vmax) {
        
        this->colors = new scolor_t[ncolors];

        for (unsigned int i = 0; i < ncolors; i++) {
            this->colors[i].red   = colors[i].red;
            this->colors[i].green = colors[i].green;
            this->colors[i].blue  = colors[i].blue;
        }
    }

    scolor_t Eval(const slibreal_t v) {
        int l;
        slibreal_t i, f;

        i = ((v-vmin)/vmax * (ncolors-1));
        l = floor(i);

        scolor_t clr;
        if (l >= ((int)ncolors)-1) {
            clr.red   = colors[ncolors-1].red;
            clr.green = colors[ncolors-1].green;
            clr.blue  = colors[ncolors-1].blue;
        } else if (l < 0) {
            clr.red   = colors[0].red;
            clr.green = colors[0].green;
            clr.blue  = colors[0].blue;
        } else {
            f = i - (slibreal_t)l;

            // Linearly interpolate in colors
            clr.red   = (unsigned int)((1-f)*colors[l].red + f*colors[l+1].red);
            clr.green = (unsigned int)((1-f)*colors[l].green + f*colors[l+1].green);
            clr.blue  = (unsigned int)((1-f)*colors[l].blue + f*colors[l+1].blue);
        }

        return clr;
    }
};

#endif/*_IMAGE_COLORMAP_H*/
