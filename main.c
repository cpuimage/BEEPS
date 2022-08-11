
#include "browse.h"

#define USE_SHELL_OPEN

#define STB_IMAGE_STATIC
#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION

#include "stb_image.h"
#include "stb_image_write.h"

#include <math.h>
#include <stdio.h>
#include "timing.h"
#include <stdint.h>

#ifndef MIN
#define MIN(a, b)    ( (a) > (b) ? (b) : (a) )
#endif
#ifndef MAX
#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#endif


unsigned char *loadImage(const char *filename, int *Width, int *Height, int *Channels) {
    return (stbi_load(filename, Width, Height, Channels, 0));
}


void saveImage(const char *filename, int Width, int Height, int Channels, unsigned char *Output) {

    if (!stbi_write_jpg(filename, Width, Height, Channels, Output, 100)) {
        fprintf(stderr, "save JPEG fail.\n");
        return;
    }
#ifdef USE_SHELL_OPEN
    browse(filename);
#endif
}


void splitpath(const char *path, char *drv, char *dir, char *name, char *ext) {
    const char *end;
    const char *p;
    const char *s;
    if (path[0] && path[1] == ':') {
        if (drv) {
            *drv++ = *path++;
            *drv++ = *path++;
            *drv = '\0';
        }
    } else if (drv)
        *drv = '\0';
    for (end = path; *end && *end != ':';)
        end++;
    for (p = end; p > path && *--p != '\\' && *p != '/';)
        if (*p == '.') {
            end = p;
            break;
        }
    if (ext)
        for (s = end; (*ext = *s++);)
            ext++;
    for (p = end; p > path;)
        if (*--p == '\\' || *p == '/') {
            p++;
            break;
        }
    if (name) {
        for (s = p; s < end;)
            *name++ = *s++;
        *name = '\0';
    }
    if (dir) {
        for (s = path; s < p;)
            *dir++ = *s++;
        *dir = '\0';
    }
}

#define CLAMP255(x) (((x) > (255)) ? (255) : (((x) < (0)) ? (0) : (x)))

inline float calcWeight(const float weight, const float spatialContraDecay, const float diff) {
  return spatialContraDecay * expf(weight * (diff * diff)) * diff;
}

void BEEPS_Filter(const unsigned char *input, unsigned char *output, size_t width, size_t height, int channels,
                  size_t stride, float PhotometricStandardDeviation, float SpatialDecay, int RangeFilter) {

    float spatialContraDecay = 1.0f - SpatialDecay;
    float lambda = expf(-SpatialDecay);
    float rho =
            lambda + SpatialDecay / (2.0f - (spatialContraDecay + SpatialDecay));
    float inv_rho = 1.0f / rho;
    float weight = 0;
    int bpp = channels;
    if (bpp == 4)
        bpp--;
    // Gaussian
    if (RangeFilter == 0) {
        weight = -0.5f / (PhotometricStandardDeviation * PhotometricStandardDeviation);
    }
    // Hyperbolic Secant
    if (RangeFilter == 1) {
        float PI = 3.14159265358f;
        weight = -PI / (2.0f * PhotometricStandardDeviation);
    }
    // Euler
    if (RangeFilter == 2) {
        float euler = 2.718281828459f;
        weight = -powf((0.75f * euler) / (PhotometricStandardDeviation *
                                          PhotometricStandardDeviation),
                       1.0f / 3.0f);
        weight *= (PhotometricStandardDeviation < 0.0f)
                  ? (-1.0f)
                  : ((0.0 == PhotometricStandardDeviation) ? (0.0) : (1.0f));
    } 
    float *cache = (float *) calloc(stride * height, sizeof(float));
    if (cache == NULL)
        return;

    for (int y = 0; y < height; y++) {
        float *lineCache = cache + y * stride;
        const unsigned char *lineIn = input + y * stride;
        for (int x = 0; x < width; x++) {
            for (int c = 0; c < channels; ++c) {
                lineCache[c] = lineIn[c];
            }
            lineCache += channels;
            lineIn += channels;
        }
    }
    float *horizontal = cache;
    for (int y = 0; y < height; y++) {
        // forwards
        float *fNextHori = horizontal + y * stride + channels;
        float *fPrevHori = horizontal + y * stride;
        for (int x = 1; x < width; x++) {
            for (int c = 0; c < bpp; ++c) {
                fNextHori[c] -= calcWeight(weight, spatialContraDecay, fNextHori[c] - rho * fPrevHori[c]);
                fNextHori[c] *= inv_rho;
            }
            fNextHori += channels;
            fPrevHori += channels;
        }
        // backwards
        float *bPrevHori = horizontal + y * stride + (width - 2) * channels;
        float *bNextHori = horizontal + y * stride + (width - 1) * channels;
        for (int x = 1; x < width; x++) {
            for (int c = 0; c < bpp; ++c) {
                bPrevHori[c] -= calcWeight(weight, spatialContraDecay, bPrevHori[c] - rho * bNextHori[c]);
                bPrevHori[c] *= inv_rho;
            }
            bPrevHori -= channels;
            bNextHori -= channels;
        }
    }
    float *vertical = cache;
    for (int x = 0; x < width; x++) {
        // forwards
        float *fNextVert = vertical + x * channels + stride;
        float *fPrevVert = vertical + x * channels;
        for (int y = 1; y < height; y++) {
            for (int c = 0; c < bpp; ++c) {
                fNextVert[c] -= calcWeight(weight, spatialContraDecay, fNextVert[c] - rho * fPrevVert[c]);
                fNextVert[c] *= inv_rho;
            }
            fNextVert += stride;
            fPrevVert += stride;
        }
        // backwards
        float *bPrevVert = vertical + x * channels + (height - 2) * stride;
        float *bNextVert = vertical + x * channels + (height - 1) * stride;
        for (int y = 1; y < height; y++) {
            for (int c = 0; c < bpp; ++c) {
                bPrevVert[c] -= calcWeight(weight, spatialContraDecay, bPrevVert[c] - rho * bNextVert[c]);
                bPrevVert[c] *= inv_rho;
            }
            bPrevVert -= stride;
            bNextVert -= stride;
        }
    }
    for (int y = 0; y < height; y++) {
        const float *lineCache = cache + y * stride;
        unsigned char *lineOut = output + y * stride;
        for (int x = 0; x < width; x++) {
            for (int c = 0; c < channels; ++c) {
                lineOut[c] = (unsigned char) CLAMP255(lineCache[c]);
            }
            lineOut += channels;
            lineCache += channels;
        }
    }
    free(cache);
}

int main(int argc, char **argv) {
    printf("Bi-Exponential Edge-Preserving Blur Filter Implementation In C\n ");
    printf("blog:http://cpuimage.cnblogs.com/ \n ");
    if (argc < 2) {
        printf("usage: %s   image \n ", argv[0]);
        printf("eg: %s   d:\\image.jpg \n ", argv[0]);

        return (0);
    }
    char *in_file = argv[1];
    char drive[3];
    char dir[256];
    char fname[256];
    char ext[256];
    char out_file[1024];
    splitpath(in_file, drive, dir, fname, ext);
    sprintf(out_file, "%s%s%s_out.jpg", drive, dir, fname);

    int Width = 0;
    int Height = 0;
    int Channels = 0;
    unsigned char *inputImage = NULL;
    double startTime = now();
    inputImage = loadImage(in_file, &Width, &Height, &Channels);
    double nLoadTime = calcElapsed(startTime, now());
    printf("load time: %d ms.\n ", (int) (nLoadTime * 1000));
    if ((Channels != 0) && (Width != 0) && (Height != 0)) {
        unsigned char *outputImg = (unsigned char *) stbi__malloc(Width * Channels * Height * sizeof(unsigned char));
        if (!inputImage || !outputImg) {
            printf("load: %s fail!\n ", in_file);
            return -1;
        }
        float photometricStandardDeviation = 255.0f;
        float spatialDecay = 0.1;  // 0.1 [0.01, .250]
        int rangeFilter = 1; // def: 0  [0,2] [Gaussian|Hyperbolic Secant|Euler Constant]
        startTime = now();
        BEEPS_Filter(inputImage, outputImg, Width, Height, Channels, Width * Channels, photometricStandardDeviation,
                     spatialDecay, rangeFilter);
        double nProcessTime = calcElapsed(startTime, now());
        printf("process time: %d ms.\n ", (int) (nProcessTime * 1000));

        startTime = now();
        saveImage(out_file, Width, Height, Channels, outputImg);
        double nSaveTime = calcElapsed(startTime, now());
        printf("save time: %d ms.\n ", (int) (nSaveTime * 1000));

        stbi_image_free(outputImg);
        stbi_image_free(inputImage);

    } else {
        printf("load: %s fail!\n", in_file);
    }
    printf("press any key to exit. \n");
    getchar();
    return (EXIT_SUCCESS);
}
