/****************************************************************************
 *
 * Feature Luminance by perspective and parallel projection.
 *****************************************************************************/

#include <visp/vpMatrix.h>
#include <visp/vpHomogeneousMatrix.h>
#include <visp/vpDisplay.h>
#include <visp/vpPixelMeterConversion.h>
#include <visp/vpImageConvert.h>
#include <visp/vpImageFilter.h>
#include <visp/vpException.h>

#include <visp/vpDisplayGTK.h>
#include <visp/vpDisplayGDI.h>
#include <visp/vpDisplayOpenCV.h>
#include <visp/vpDisplayD3D.h>
#include <visp/vpDisplayX.h>

#include "npFeatureLuminance.h"
#include <math.h>
#include <iomanip>  // setprecision

using namespace std;

#define PI 3.1415926;

/*!
  \file vpFeatureLuminance.cpp
  \brief Class that defines the image luminance visual feature

  for more details see
  C. Collewet, E. Marchand, F. Chaumette. Visual
  servoing set free from image processing. In IEEE Int. Conf. on
  Robotics and Automation, ICRA'08, Pages 81-86, Pasadena, Californie,
  Mai 2008.
*/



/*!
  Initialize the memory space requested for vpFeatureLuminance visual feature.
*/
void
npFeatureLuminance::init()
{
    if (flags == NULL)
        flags = new bool[nbParameters];
    for (unsigned int i = 0; i < nbParameters; i++) flags[i] = false;

    //default value Z (1 meters)
    Z_fl = 1;

    firstTimeIn =0 ;

}


void
npFeatureLuminance::init(unsigned int _nbr, unsigned int _nbc, double _Z, projectionModel projModel, controlLawZ vfeature_Z)
{
    init() ;

    nbr = _nbr ;
    nbc = _nbc ;

    if((nbr < 2*bord) || (nbc < 2*bord)){
        throw vpException(vpException::dimensionError, "border is too important compared to number of row or column.");
    }

    // number of feature = nb column x nb lines in the images
    dim_s = (nbr-2*bord)*(nbc-2*bord) ;

    imIx.resize(nbr,nbc) ;
    imIy.resize(nbr,nbc) ;
    imG.resize(nbr,nbc) ;
    imGxy.resize(nbr,nbc) ;//for display gradient
    imNoise.resize(nbr,nbc) ;//for display noise (electron charging)

    s.resize(dim_s) ;
    sg.resize(dim_s);
    snv.resize(dim_s);

    if (pixInfo != NULL)
        delete [] pixInfo;

    pixInfo = new npLuminance[dim_s] ;

    Z_fl = _Z ;
    pjModel = projModel;

    vf_Z = vfeature_Z;

    S_c = 1;
    S_p = 0;

    denoise = false;
    nbrNoise = 0;
}

/*! 
  Default constructor that build a visual feature.
*/
npFeatureLuminance::npFeatureLuminance() : vpBasicFeature()
{
    nbParameters = 1;
    dim_s = 0 ;
    bord = 10 ; // default : 10
    flags = NULL;
    pixInfo = NULL;

    init() ;
}

/*! 
  Default destructor.
*/
npFeatureLuminance::~npFeatureLuminance()
{
    if (pixInfo != NULL) delete [] pixInfo ;
    if (flags != NULL) delete [] flags;
}



/*!
  Set the value of \f$ Z \f$ which represents the depth in the 3D camera frame.

  \param Z : \f$ Z \f$ value to set.
*/
void
npFeatureLuminance::setZ(const double Z)
{
    this->Z_fl = Z ;
    flags[0] = true;
}

void
npFeatureLuminance::setSigma(const double sigma)
{
    this->sigma = sigma ;
}


/*!
  Get the value of \f$ Z \f$ which represents the depth in the 3D camera frame.

  \return The value of \f$ Z \f$.
*/
double
npFeatureLuminance::getZ() const
{
    return Z_fl ;
}

vpMatrix
npFeatureLuminance::getLg()
{
    return Lg;
}

vpMatrix
npFeatureLuminance::getLdg()
{
    return Ldg;
}

vpColVector
npFeatureLuminance::getSg()
{
    return sg;
}

double
npFeatureLuminance::getGlobleSg()
{
    double S_g = 0;
    for (int m=0; m<sg.size();m++)
        S_g+= fabs(sg[m]);

   // cout << "S_g =" <<S_g <<endl;
    return S_g/(sg.size()-nbrNoise);
}

double
npFeatureLuminance::getGlobleSnv()
{
    double S_nv = 0;
    for (int m=0; m<snv.size();m++)
        S_nv+= fabs(snv[m]);

    //cout << " S_nv = " << S_nv << "/" << snv.size()-nbrNoise<< endl;

    return S_nv/double(snv.size()-nbrNoise);
}

double
npFeatureLuminance::getVariance()
{
    return variance;
}

double
npFeatureLuminance::getSpectralRatio()
{
    return spect_rat_mean;
}

double
npFeatureLuminance::getPhaseFactor()
{
    return phaseFactor;
}

vpMatrix
npFeatureLuminance::getVarianceMap()
{
    return varMap;
}

vpMatrix
npFeatureLuminance::getSigmaMap()
{
    return sigmaMap;
}

vpMatrix
npFeatureLuminance::getSpectralRatioMap()
{
    return spectMap;
}

void
npFeatureLuminance::setCameraParameters(vpCameraParameters &_cam)
{
    cam = _cam ;
}

/*!

  Build a luminance feature directly from the image
*/

void
npFeatureLuminance::matrixConvert(cv::Mat& src, vpMatrix& dest)
{
    cv::Size ksize;
    ksize = src.size();
    dest.resize(ksize.height,ksize.width);
    for(int i=0; i< ksize.height; i++)
        for(int j=0; j< ksize.width; j++)
            dest[i][j]=src.at<float>(i,j);
}

void
npFeatureLuminance::matrixConvert(cv::Mat& src, vpImage<double> & dest)
{
    cv::Size ksize;
    ksize = src.size();
    dest.init(ksize.height,ksize.width);
    for(int i=0; i< ksize.height; i++)
        for(int j=0; j< ksize.width; j++)
            dest[i][j]=src.at<float>(i,j);
}


void
npFeatureLuminance::matrixConvert(vpMatrix& src, cv::Mat& dest)
{
    cv::Size ksize;
    ksize.height = src.getRows();
    ksize.width = src.getCols();
    dest.create(ksize,CV_32F);
    for(int i=0; i< ksize.height; i++)
        for(int j=0; j< ksize.width; j++)
            dest.at<float>(i,j)=src[i][j];
}

void
npFeatureLuminance::matrixConvert(vpImage<double> &src, cv::Mat& dest)
{
    cv::Size ksize;
    ksize.height = src.getRows();
    ksize.width = src.getCols();
    dest.create(ksize,CV_32F);
    for(int i=0; i< ksize.height; i++)
        for(int j=0; j< ksize.width; j++)
            dest.at<float>(i,j)=src[i][j];
}



void
npFeatureLuminance::buildFrom(vpImage<unsigned char> &I)
{

    unsigned int l = 0;
    double Ix,Iy;

    double px = cam.get_px() ;
    double py = cam.get_py() ;

    Ic = I; // save image to Ic for computing sigma

    noiseLevel = I.getMaxValue()-1;
    //cout << "noise_level=" << noise_level << endl;

    if (firstTimeIn==0)
    {
        firstTimeIn=1 ;
        l =0 ;
        for (unsigned int i=bord; i < nbr-bord ; i++)
        {
            //   cout << i << endl ;
            for (unsigned int j = bord ; j < nbc-bord; j++)
            {

                double x=0,y=0;
                vpPixelMeterConversion::convertPointWithoutDistortion(cam,
                                                                      i, j,
                                                                      y, x)  ;

                pixInfo[l].x = x;
                pixInfo[l].y = y;
                pixInfo[l].Z_l = Z_fl ;

                l++;
            }
        }
    }


    for (int i=3; i < nbr-3 ; i++)
    {
        //   cout << i << endl ;
        for (int j = 3 ; j < nbc-3; j++)
        {
            // cout << dim_s <<" " <<l <<"  " <<i << "  " << j <<endl ;
            imIx[i][j] =  vpImageFilter::derivativeFilterX(I,i,j) ;//px *
            imIy[i][j] =  vpImageFilter::derivativeFilterY(I,i,j) ;//py *
            //cout << imIx[i][j] << endl;
        }
    }

    // Noramlized Variance
  if(vf_Z==NormalizedVariance || vf_Z==Gradient_Sigma_Variance)
  {
      //For Normalized variance
      double sum=0;
      double mean;
      for (unsigned int i=bord; i < nbr-bord ; i++)
          for (unsigned int j=bord; j < nbc-bord ; j++)
      {
            sum +=I[i][j];
      }

      mean = sum/(nbr*nbc);

      //cout << "mean = " << mean << endl;

      double mean_iv = 1/mean;
      double mean_iv_per = mean_iv/(nbr*nbc);

      l=0;

      for (unsigned int i=bord; i < nbr-bord ; i++)
          for (unsigned int j=bord; j < nbc-bord ; j++)
      {
            pixInfo[l].norVar = (I[i][j]-mean)*(I[i][j]-mean)*mean_iv;

            pixInfo[l].Iz = 2*mean_iv*mean_iv*(I[i][j]-mean); //? should not be good !
            //sg[l] = fabs((I[i+1][j]- I[i][j]));//F-1 absolute gradiant
            //sg[l] = (I[i+1][j]- I[i][j])*(I[i+1][j]- I[i][j]);//F-2 Squared gradiant
            //sg[l] = (I[i+2][j]- I[i][j])*(I[i+2][j]- I[i][j]);//F-3 Brenner gradiant
            //sg[l] = (I[i][j]-mean)*(I[i][j]-mean);//F-10 variance
            snv[l] = (I[i][j]-mean)*(I[i][j]-mean)*mean_iv;//F-11 Normalized variance
            //sg[l] = I[i][j]*I[i+1][j]- I[i][j]*I[i+2][j];//F-12 AutoCorrelation
            //sg[l] = (I[i][j]*I[i+1][j])-mean*mean;//F-13 Standard-deviation-based correlation
            //sg[l] = I[i][j]*I[i][j];// F-18 image power

            //cout << "snv[l]=" << snv[l] << endl;
            l++;
      }

      /*
      for (int u=0; u < filter_size ; u++)
          for (int v = 0 ; v < filter_size; v++)
          {
              int u2=u*u;
              int v2=v*v;
              double exp_uv = exp(-(u2+v2)*sigma_2_inv);
              Iuv[u][v] = -m*Z_2_inv*((u2+v2)*sigma_2_inv-1)*pi_inv*sigma_3_inv*exp_uv;
          }*/
  }

    //For luminance
    l= 0 ;
    for (unsigned int i=bord; i < nbr-bord ; i++)
    {
        //   cout << i << endl ;
        for (unsigned int j = bord ; j < nbc-bord; j++)
        {

            // cout << dim_s <<" " <<l <<"  " <<i << "  " << j <<endl ;
            Ix =  px * vpImageFilter::derivativeFilterX(I,i,j) ;
            Iy =  py * vpImageFilter::derivativeFilterY(I,i,j) ;


            pixInfo[l].I  =  I[i][j] ;
            //cout << "pixInfo[" << l <<"].I = " <<pixInfo[l].I  << endl;
            s[l]  =  I[i][j] ;
            pixInfo[l].Ix  = Ix;
            pixInfo[l].Iy  = Iy;

            l++;
        }
    }

    // for image gradient
    l= 0 ;
    for (unsigned int i=bord; i < nbr-bord ; i++)
    {
        //   cout << i << endl ;
        for (unsigned int j = bord ; j < nbc-bord; j++)
        {
            sg[l] = ( vpMath::sqr(vpImageFilter::derivativeFilterX(I,i,j)) + vpMath::sqr(vpImageFilter::derivativeFilterY(I,i,j)) ) ;//sqrt

            // detection of noise (pixels with maximum gris level)
            if ((int)Ic[i][j]>noiseLevel && denoise)
            {
                imNoise[i][j]=255;
                noisePosition.push_back(l);
                sg[l]=0;
            }
            else
                imNoise[i][j]=0;

            //cout << "sg[" << l << "] = " << sg[l] <<endl;
            l++;
        }
    }
}

double npFeatureLuminance::computeImageGradient(const vpImage<unsigned char> &I)
{
    double g=0;
    for (unsigned int i=bord; i < I.getRows()-bord ; i++)
    {
        for (unsigned int j = bord ; j < I.getCols()-bord; j++)
        {
            g += fabs( vpMath::sqr(vpImageFilter::derivativeFilterX(I,i,j)) + vpMath::sqr(vpImageFilter::derivativeFilterY(I,i,j)) ) ;//sqrt
        }
    }
    //cout  << "g = " << g << endl;
    return g / ((I.getRows()-2*bord)*(I.getCols()-2*bord));
}

void
npFeatureLuminance::computeGradientJacobian(const vpImage <unsigned char> &Id)
{

    //    cout << "Id size: " << Id.getHeight() << "x" << Id.getWidth() <<endl;
    double Ix,Iy;
    double px = cam.get_px() ;
    double py = cam.get_py() ;

    vpMatrix Iuv; // derivative of gaussian kernel
    //vpMatrix I2uv; // second derivaive of gaussian kernel
    int filter_size = 15; // should be odd!
    double m = -10;//a coefficient which depends on the sensor(diamter of aperture, distance between image plane and lens, focus distance...);

    Iuv.resize(filter_size,filter_size);

    double sigma_2_inv = 0.5/(sigma*sigma);
    double sigma_3_inv = 1/(sigma*sigma*sigma);

    const double pi_inv = 1/PI;
    double Z_2_inv = 1/(Z_fl*Z_fl);

    int u,v;

    for (int i=0; i < filter_size ; i++)
        for (int j = 0 ; j < filter_size; j++)
        {
            u = i - (filter_size-1)/2;
            v = j - (filter_size-1)/2;
            int u2=u*u;
            int v2=v*v;
            double exp_uv = exp(-(u2+v2)*sigma_2_inv);
            Iuv[i][j] = -m*Z_2_inv*((u2+v2)*sigma_2_inv-1)*pi_inv*sigma_3_inv*exp_uv;
        }

 /*   cout <<  "Iuv :" << endl;
    Iuv.print(std::cout,6);
*/

     vpImage<double> ImS;
    ImS.init(filter_size,filter_size);

/*    vpMatrix ImS;
    ImS.resize(filter_size,filter_size);
*/
    unsigned int l= 0 ;

 /* ----------without optimization ---------------*/

    for (int i=bord; i < nbr-bord ; i++)
    {
        for (int j = bord ; j < nbc-bord; j++)
        {
           double Ixs = 0;
           double Iys = 0;

          //----Convolution----
            for (int ii= 0; ii < filter_size ; ii++)
              for (int jj= 0; jj < filter_size ; jj++)
              {
                  u = ii - (filter_size-1)/2;
                  v = jj - (filter_size-1)/2;
                if( (i+ii >= 3) && (i+ii < nbr-3) && (j+jj >= 3) && (j+jj < nbc-3))
                {
                  ImS[ii][jj] = Id[i+ii][j+jj]*Iuv[ii][jj];// i+/-ii isn't the same thing
                }
              }

          //cout <<  "ImS :" << endl;
           // ImS.print(std::cout,6);

            for (int u = 3; u < filter_size-3; u++)
              for (int v = 3; v < filter_size-3; v++)
              {
                  Ixs += vpImageFilter::derivativeFilterX(ImS,u,v) ;//ImS[u+1][v]-ImS[u][v]; //
                  Iys += vpImageFilter::derivativeFilterY(ImS,u,v) ;//ImS[u][v+1]-ImS[u][v];
              }

            pixInfo[l].Ixs = Ixs;
            pixInfo[l].Iys = Iys;

           // cout << "I[" << i << "][" << j <<"] = " << pixInfo[l].I <<endl;
         //   cout << "Ixs, Iys [" << i << "][" << j <<"] = " << Ixs << "\t" << Iys <<endl;

            Ix =  imIx[i][j] ;
            Iy =  imIy[i][j] ;
            imGxy[i][j] = vpMath::sqr(Ix) + vpMath::sqr(Iy) ;

            // Calcul de Z
            pixInfo[l].Ix_g  = Ix * px;
            pixInfo[l].Iy_g  = Iy * py;

            l++;
        }
    }

}

void
npFeatureLuminance::computeVarianceFourier(const vpImage <unsigned char> &Id)
{
    cv::Mat cId, cI, cFId,cFI, cHanning; // desired/current image in space and frequence domain
    cId.create(nbr,nbc,CV_32F);
    cI.create(nbr,nbc,CV_32F);
    vpImageConvert::convert(Id,cId);
    vpImageConvert::convert(Ic,cI);

    //nbr = Id.getRows();
    //nbc = Id.getCols();

    cv::Mat cI_float, cId_float;
    cI.convertTo(cI_float,CV_32F);
    cId.convertTo(cId_float,CV_32F);

    // applying hanning window
 /*   cv::Size sHanning(9,9);
    cv::createHanningWindow(cHanning, sHanning,CV_32F);
    cv::filter2D(cI_float,cI_float,-1,cHanning);
    cv::filter2D(cId_float,cId_float,-1,cHanning);
*/
    cv::vector<cv::Mat> planes(2);

    cv::dft(cId_float,cFId, cv::DFT_SCALE|cv::DFT_COMPLEX_OUTPUT);

    cv::split(cFId, planes);
    cv::Mat cFId_re=planes[0];
    cv::Mat cFId_im=planes[1];

    planes.clear();
    cv::dft(cI_float,cFI,cv::DFT_SCALE|cv::DFT_COMPLEX_OUTPUT);
    cv::split(cFI, planes);
    cv::Mat cFI_re=planes[0];
    cv::Mat cFI_im = planes[1];

    double nbr_inv = 1.0/nbr;
    double nbc_inv = 1.0/nbc;
    double u,v;
    double var_sum = 0;
    double spec_rat_sum = 0;
    vpMatrix var; // variance
    var.resize(nbr,nbc);

    vpMatrix spect_rat; // spectral ratio
    spect_rat.resize(nbr,nbc);

    vpMatrix sigmaM;
    sigmaM.resize(nbr,nbc);
    double pi = PI;
    int nbrErr = 0;
    int nbrErr_spec = 0;
    double pi_inv = 1.0 / (2.0 * pi* pi);
    int radius = 30; // range for computing sigma
    int boader = 2;

    phaseFactor = 0.;

    for(int i = boader; i < nbr-boader; i++)
        for(int j = boader;  j < nbc-boader; j++)
        {
            // ------------------------------variance---------------------------------
            if (i < nbr/2)
                u=(i+1)*nbr_inv;
            else
                u=(i+1)*nbr_inv-1;
            if (j < nbc/2)
                v=(j+1)*nbc_inv;
            else
                v=(j+1)*nbc_inv-1;
            double u2v2_inv = 1/(u*u+v*v);
            //cout << "u2v2_inv = " << u2v2_inv << endl;
            //cout << "log(fabs(cFId_real.at<float>(i,j))) = " << log(fabs(cFId_real.at<float>(i,j))) << endl;
            //var[i][j] = 4 * (log(fabs(cFId_re.at<float>(i,j)))-log(fabs(cFI_re.at<float>(i,j))))* u2v2_inv;
            //var[i][j] = pi_inv * (log(fabs(cFId_re.at<float>(i,j)))-log(fabs(cFI_re.at<float>(i,j))))* u2v2_inv;
            var[i][j] = pi_inv*(log(cFId_re.at<float>(i,j)*cFId_re.at<float>(i,j) + cFId_im.at<float>(i,j)* cFId_im.at<float>(i,j))\
                                      -log(cFI_re.at<float>(i,j)*cFI_re.at<float>(i,j) + cFI_im.at<float>(i,j)* cFI_im.at<float>(i,j)))* u2v2_inv;          

            if (isinf(var[i][j]) || isnan(var[i][j]) || var[i][j]<0 || var[i][j]>50)//
            {
                //cout << "inf !!! : var[" << i << "][" << j <<"] = " << var[i][j] << endl;
                var[i][j] = 0; // remove Inf and NaN, (and negative value)
                //if ((i<radius && j < radius) || (i>nbr-radius && j < radius) || (i<radius && j > nbc -radius) ||(i>nbr-radius && j > nbc -radius))
                if ((i*i+j*j< radius*radius) || ((nbr-i)*(nbr-i)+j*j< radius*radius) || (i*i+(nbc-j)*(nbc-j)< radius*radius) || ((nbr-i)*(nbr-i)+(nbc-j)*(nbc-j)< radius*radius))
                nbrErr ++;
            }

            // ------------------------------spectral ratio---------------------------------
            double spec_den = cFId_re.at<float>(i,j)*cFId_re.at<float>(i,j) + cFId_im.at<float>(i,j)*cFId_im.at<float>(i,j);
            double spec_num_re = cFI_re.at<float>(i,j)*cFId_re.at<float>(i,j) + cFI_im.at<float>(i,j)*cFId_im.at<float>(i,j);
            double spec_num_im = cFId_re.at<float>(i,j)*cFI_im.at<float>(i,j) - cFI_re.at<float>(i,j)*cFId_im.at<float>(i,j);

            spect_rat[i][j] = (spec_num_re*spec_num_re + spec_num_im*spec_num_im)/(spec_den*spec_den); //
           // cout << "spec_rat[" << i << "][" << j <<"] = " << spect_rat[i][j] << endl;

            if (isinf(spect_rat[i][j]) || isnan(spect_rat[i][j]) || spect_rat[i][j]<0 )//
            {
                cout << "inf !!! : spec_rat[" << i << "][" << j <<"] = " << spect_rat[i][j] << endl;
                spect_rat[i][j] = 0; // remove Inf and NaN, (and negative value)
                //if ((i<radius && j < radius) || (i>nbr-radius && j < radius) || (i<radius && j > nbc -radius) ||(i>nbr-radius && j > nbc -radius))
                if ((i*i+j*j< radius*radius) || ((nbr-i)*(nbr-i)+j*j< radius*radius) || (i*i+(nbc-j)*(nbc-j)< radius*radius) || ((nbr-i)*(nbr-i)+(nbc-j)*(nbc-j)< radius*radius))
                nbrErr_spec ++;
            }

            //if ((i<radius && j < radius) || (i>nbr-radius && j < radius) || (i<radius && j > nbc -radius) ||(i>nbr-radius && j > nbc -radius))
            if ((i*i+j*j< radius*radius) || ((nbr-i)*(nbr-i)+j*j< radius*radius) || (i*i+(nbc-j)*(nbc-j)< radius*radius) || ((nbr-i)*(nbr-i)+(nbc-j)*(nbc-j)< radius*radius))
            {
                spec_rat_sum += (spect_rat[i][j]);//fabs

                var_sum += (var[i][j]);//fabs

                // test phase as a cost function
                if(cFI_re.at<float>(i,j) != 0 || cFId_re.at<float>(i,j) != 0)
                    phaseFactor += fabs(fabs(cFI_im.at<float>(i,j)/cFI_re.at<float>(i,j))-fabs(cFId_im.at<float>(i,j)/cFId_re.at<float>(i,j) ));
            }

            sigmaM[i][j] = (var[i][j]>0);//sqrt(fabs(var[i][j]));

        }

    varMap = var;
    sigmaMap = sigmaM;
    spectMap = spect_rat;

    //variance = var_sum/(nbr*nbc-nbrErr);

    variance = var_sum/(4*(radius-boader)*(radius-boader)-nbrErr);
    spect_rat_mean = spec_rat_sum/(4*(radius-boader)*(radius-boader)-nbrErr_spec);
}




/*!

  Compute and return the interaction matrix \f$ L_I \f$. The computation is made
  thanks to the values of the luminance features \f$ I \f$
*/
void
npFeatureLuminance::interaction(vpMatrix &L)
{
    double x,y,Ix,Iy,z,Zinv;

    double  Ix_g, Iy_g, Ixs, Iys;

    if(pjModel==perspective)
    {
        L.resize(dim_s,6) ;
        for(unsigned int m = 0; m< L.getRows(); m++)
        {
            Ix = pixInfo[m].Ix;
            Iy = pixInfo[m].Iy;

            x = pixInfo[m].x ;
            y = pixInfo[m].y ;
            Zinv =  1 / pixInfo[m].Z_l;

            {
                L[m][0] = Ix * Zinv;
                L[m][1] = Iy * Zinv;
                L[m][2] = -(x*Ix+y*Iy)*Zinv;
                L[m][3] = -Ix*x*y-(1+y*y)*Iy;
                L[m][4] = (1+x*x)*Ix + Iy*x*y;
                L[m][5]  = Iy*x-Ix*y;
            }
        }
    }
    else if(pjModel==parallel)
    {
        L.resize(dim_s,6) ;
        for(unsigned int m = 0; m< L.getRows(); m++)
        {
            Ix = pixInfo[m].Ix;
            Iy = pixInfo[m].Iy;

            x = pixInfo[m].x ;
            y = pixInfo[m].y ;
            z = pixInfo[m].Z_l ;

            {
                L[m][0] = Ix;
                L[m][1] = Iy ;
                /*L[m][2] = -Iy * z;
            L[m][3] = Ix * z;
            L[m][4]  = Iy*x-Ix*y;*/

                L[m][2] =  0 ;
                L[m][3] = -Iy * z;
                L[m][4] = Ix * z;
                L[m][5]  = Iy*x-Ix*y;
            }
        }
    }
    else if(pjModel==parallelZ)
    {
        L.resize(dim_s,6) ;
        Lg.resize(dim_s,6);
        Ldg.resize(dim_s,6);
        Jg.resize(dim_s);

        std::list<int>::iterator it;
        it=noisePosition.begin();

        //CostFunctionSg=0;

        for(unsigned int m = 0; m< L.getRows(); m++)
        {
            Zinv =  1 / pixInfo[m].Z_l;
            Ix = pixInfo[m].Ix;
            Iy = pixInfo[m].Iy;

            x = pixInfo[m].x ;
            y = pixInfo[m].y ;
            z = pixInfo[m].Z_l;

         /*   Ixx = pixInfo[m].Ixx;
            Ixy = pixInfo[m].Ixy;
            Iyx = pixInfo[m].Iyx;
            Iyy = pixInfo[m].Iyy;*/

            Ix_g = pixInfo[m].Ix_g;
            Iy_g = pixInfo[m].Iy_g;

            Ixs = pixInfo[m].Ixs;
            Iys = pixInfo[m].Iys;

            //Ixss = pixInfo[m].Ixss;
            //Iyss = pixInfo[m].Iyss;

            //cout << "Ix and Ix_g: "<<Ix << "  " << Ix_g << endl;

         /*   double A = ( Ixx*Ix_g+Iyx*Iy_g );
            double B = ( Ixy*Ix_g+Iyy*Iy_g );*/
            double D = (Ix*Ixs +Iy*Iys); // maximize image gradient as cost function => mimimize derivative of image gradient
            //    double E = (Ixs*Ixs+Ix*Ixss +Iys*Iys+Iy*Iyss); // derivative of D;

            int noiseId = *it; //here is the electron charging, remove it!
            //cout << "noiseId= " << noiseId << endl;
            if (m==noiseId && denoise)
            {
                D=0;
                it++;
            }

            /*
          if(m%20000==0)
          {
              cout << "D[" << m << "]=" << D <<endl;
              cout << "E[" << m << "]=" << E <<endl;
          }
*/
            //cout << "Ix,Ixs, Iy,Iys= " << Ix << "\t" << pixInfo[m].Ixs << "\t"<< Iy << "\t"<< Iys << "\t" << endl;

           // double C= -( Ix+Iy )*pixInfo[m].norVar;
            double C = pixInfo[m].Iz;

            {
                L[m][0] = Ix;
                L[m][1] = Iy ;
                L[m][2] = 0;//(A+B);
                L[m][3] = -Iy * z;
                L[m][4] = Ix * z;
                L[m][5]  = Iy*x-Ix*y;
            }

            if(vf_Z== ImageGradient || vf_Z== Gradient_Sigma_Variance)
            {
                Lg[m][0] = 0;
                Lg[m][1] = 0;
                Lg[m][2] = 2*D;//(A*x+B*y)*Zinv;//2*(A+B)//2*D;
                Lg[m][3] = 0;
                Lg[m][4] = 0;
                Lg[m][5] = 0;

                Jg[m] = 2*D;


                /*
              Ldg[m][0] = 0;
              Ldg[m][1] = 0;
              Ldg[m][2] = 2*E;
              Ldg[m][3] = 0;
              Ldg[m][4] = 0;
              Ldg[m][5] = 0;
*/
                //Lz[m] = Lg[m][2];
                //CostFunctionSg += (A+B);
                //cout << "CostFunctionLz="<<CostFunctionSg << endl;
            }
            else if(vf_Z==NormalizedVariance || vf_Z== Gradient_Sigma_Variance)
            {
                Lg[m][0] = 0;
                Lg[m][1] = 0;
                Lg[m][2] = C;
                Lg[m][3] = 0;
                Lg[m][4] = 0;
                Lg[m][5] = 0;

                Jg[m] = 2*C;
            }
        }
    }
    else
    {
        L.resize(dim_s,6) ;
        for(unsigned int m = 0; m< L.getRows(); m++)
        {
            Ix = pixInfo[m].Ix;
            Iy = pixInfo[m].Iy;

            x = pixInfo[m].x ;
            y = pixInfo[m].y ;
            Zinv =  1 / pixInfo[m].Z_l;

            {
                L[m][0] = Ix * Zinv;
                L[m][1] = Iy * Zinv;
                L[m][2] = -(x*Ix+y*Iy)*Zinv;
                L[m][3] = -Ix*x*y-(1+y*y)*Iy;
                L[m][4] = (1+x*x)*Ix + Iy*x*y;
                L[m][5]  = Iy*x-Ix*y;
            }
        }
    }
    //cout << "L=" << L << endl;
}

/*!
  Compute and return the interaction matrix \f$ L_I \f$. The computation is made
  thanks to the values of the luminance features \f$ I \f$
*/

vpMatrix  npFeatureLuminance::interaction(const unsigned int /* select */)
{
    /* static */ vpMatrix L  ; // warning C4640: 'L' : construction of local static object is not thread-safe
    interaction(L) ;
    return L ;
}


/*!
  Compute the error \f$ (I-I^*)\f$ between the current and the desired

  \param s_star : Desired visual feature.
  \param e : Error between the current and the desired features.

*/
void
npFeatureLuminance::error(const vpBasicFeature &s_star,
                          vpColVector &e)
{
    e.resize(dim_s) ;

    for (unsigned int i =0 ; i < dim_s ; i++)
    {
        e[i] = s[i] - s_star[i] ;
    }
}

/*!
  Compute the error \f$ (I-I^*)\f$ between the current and the desired image gradient

  \param s_star : Desired visual feature.
  \param e : Error between the current and the desired features.
*/
void
npFeatureLuminance::sgError(const vpColVector &sg_star,
                            vpColVector &eg)
{
    eg.resize(dim_s) ;

    for (unsigned int i =0 ; i < dim_s ; i++)
    {
        eg[i] = sg[i] - sg_star[i] ;
        // eg[i] = sg[i] ; // use maximazing cost function C=delta(x)^2+delta(y)^2
    }
}



/*!
  Compute the error \f$ (I-I^*)\f$ between the current and the desired

  \param s_star : Desired visual feature.
  \param select : Not used.

*/
vpColVector
npFeatureLuminance::error(const vpBasicFeature &s_star,
                          const unsigned int /* select */)
{
    /* static */ vpColVector e ; // warning C4640: 'e' : construction of local static object is not thread-safe

    error(s_star, e) ;

    return e ;

}


/*!

  Not implemented.

 */
void
npFeatureLuminance::print(const unsigned int /* select */) const
{
    static int firsttime =0 ;

    if (firsttime==0)
    {
        firsttime=1 ;
        vpERROR_TRACE("not implemented") ;
        // Do not throw and error since it is not subject
        // to produce a failure
    }
}


/*!

  Not implemented.

 */
void
npFeatureLuminance::display(const vpCameraParameters & /* cam */,
                            const vpImage<unsigned char> & /* I */,
                            const vpColor &/* color */,
                            unsigned int /* thickness */) const
{
    static int firsttime =0 ;

    if (firsttime==0)
    {
        firsttime=1 ;
        vpERROR_TRACE("not implemented") ;
        // Do not throw and error since it is not subject
        // to produce a failure
    }
}

/*!

  Not implemented.

 */
void
npFeatureLuminance::display(const vpCameraParameters & /* cam */,
                            const vpImage<vpRGBa> & /* I */,
                            const vpColor &/* color */,
                            unsigned int /* thickness */) const
{
    static int firsttime =0 ;

    if (firsttime==0)
    {
        firsttime=1 ;
        vpERROR_TRACE("not implemented") ;
        // Do not throw and error since it is not subject
        // to produce a failure
    }
}


/*!
  Create an object with the same type.

  \code
  vpBasicFeature *s_star;
  vpFeatureLuminance s;
  s_star = s.duplicate(); // s_star is now a vpFeatureLuminance
  \endcode

*/
npFeatureLuminance *npFeatureLuminance::duplicate() const
{
    npFeatureLuminance *feature = new npFeatureLuminance ;
    return feature ;
}


/*
 * Local variables:
 * c-basic-offset: 2
 * End:
 */
