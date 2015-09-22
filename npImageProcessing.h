#ifndef NPIMAGEPROCESSING_H
#define NPIMAGEPROCESSING_H

#include <visp/vpDebug.h>
#include <visp/vpException.h>
#include <visp/vpMath.h>
#include <visp/vpImage.h>
#include <visp/vpImageIo.h>
#include <visp/vpImageTools.h>
#include <visp/vpImageConvert.h>
#include <visp/vpImageFilter.h>
#include <visp/vpNoise.h>


class npImageProcessing
{
public:
    static void getGaussianBlur(vpImage<unsigned char> &Iblur, const vpImage<unsigned char> &Iin, const double sigma, const int size = 7, const bool normalize = true);
    static void getGaussianKernel(double *filter, unsigned int size, double sigma, bool normalize = true);
    static void getBlurImageCV(vpImage<unsigned char> &Iblur, const vpImage<unsigned char> &Iin, const double sigma);
    static void getNoisedImage(vpImage<unsigned char> &Inoised, const vpImage<unsigned char> &Iin,const double noise_mean, const double noise_sdv);
    static void getSobelImage(vpImage<unsigned char> &IsX, vpImage<unsigned char> &IsY, const vpImage<unsigned char> &Iin);
    static void cropImage(vpImage<unsigned char> &Icrop, const vpImage<unsigned char> &Iin,const vpImagePoint p1,const vpImagePoint p2);
    static void cropImage(vpImage<unsigned char> & Icrop, const vpImage<unsigned char> &Iin, const vpRect cropROI);

};

#endif
