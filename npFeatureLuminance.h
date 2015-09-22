/****************************************************************************
 *
 * $Id: vpFeatureLuminance.h 3653 2012-03-28 12:43:23Z fspindle $
 *
 * This file is part of the ViSP software.
 * Copyright (C) 2005 - 2012 by INRIA. All rights reserved.
 * 
 * This software is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * ("GPL") version 2 as published by the Free Software Foundation.
 * See the file LICENSE.txt at the root directory of this source
 * distribution for additional information about the GNU GPL.
 *
 * For using ViSP with software that can not be combined with the GNU
 * GPL, please contact INRIA about acquiring a ViSP Professional 
 * Edition License.
 *
 * See http://www.irisa.fr/lagadic/visp/visp.html for more information.
 * 
 * This software was developed at:
 * INRIA Rennes - Bretagne Atlantique
 * Campus Universitaire de Beaulieu
 * 35042 Rennes Cedex
 * France
 * http://www.irisa.fr/lagadic
 *
 * If you have questions regarding the use of this file, please contact
 * INRIA at visp@inria.fr
 * 
 * This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
 * WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
 *
 *
 * Description:
 *   Luninance based feature .
 *
 * Authors:
 * Eric Marchand
 *
 *****************************************************************************/

#ifndef npFeatureLuminance_h
#define npFeatureLuminance_h

#include <visp/vpMatrix.h>
#include <visp/vpBasicFeature.h>
#include <visp/vpImage.h>
#include <list>



/*!
  \file vpFeatureLuminance.h
  \brief Class that defines the image luminance visual feature

  for more details see
  C. Collewet, E. Marchand, F. Chaumette. Visual
  servoing set free from image processing. In IEEE Int. Conf. on
  Robotics and Automation, ICRA'08, Pages 81-86, Pasadena, Californie,
  Mai 2008.
*/


/*!
  \class vpLuminance
  \brief Class that defines the luminance and gradient of a point

  \sa vpFeatureLuminance
*/


class VISP_EXPORT npLuminance
{
 public:
  double x, y;   // point coordinates (in meter)
  double I ; // pixel intensity
  double Ix,Iy ; // pixel gradient
  double Z_l; // pixel depth
  double Ixx,Iyy ; // pixel gradient 2
  double Ixy,Iyx ; // pixel gradient 2
  double Ix_g, Iy_g; // pixel gradient to control z
  double Ixs, Iys; // pixel gradient dsdx/dsdy by x/y and sigma (using image gradient to control z)
  double Ixss, Iyss; // pixel gradient ds2dx/ds2dy by x/y and sigma (using image gradient to control z)

  double norVar; //normalized variance
  double Isig; //pixel gradient  dI/dsig, using normalized variance to control z

  double Iz; // ds/dz for controling z

};


/*!
  \class vpFeatureLuminance
  \brief Class that defines the image luminance visual feature

  For more details see
  C. Collewet, E. Marchand, F. Chaumette. Visual
  servoing set free from image processing. In IEEE Int. Conf. on
  Robotics and Automation, ICRA'08, Pages 81-86, Pasadena, Californie,
  Mai 2008.
*/

class VISP_EXPORT npFeatureLuminance : public vpBasicFeature
{
 protected:
  //! FeaturePoint depth (required to compute the interaction matrix)
  //! default Z = 1m
  double Z_fl ;
  double sigma; //sdv of gaussian filter
  double variance; // variance in gaussian filter, ie sigma^2, used for sigma based model
  double A; // coefficent sigma = A* (Z-Z*)

  double phaseFactor; // use phase infor as a cost function : !!! THIS DOES NOT WORK

  //! Number of rows.
  unsigned int nbr ;
  //! Number of column.
  unsigned int nbc ;
  //! Border size.
  unsigned int bord ;

  //! image gradient as additional visual feature
  vpColVector sg;

  //! Normalized variance as additional visual feature
  vpColVector snv;

  //! (old) interaction matrix for image gradient
  vpMatrix Lg;

  //! Jacobian (interaction matrix) for image gradient
  vpColVector Jg;

  //! interaction matrix for normalized variance
  vpMatrix Lnv;

  //! Derivative of interaction matrix for image gradient
  vpMatrix Ldg;

  //! variance map (matrix of variance computed from fourier transform)
  vpMatrix varMap;

  //! sigma map (matrix of sigma (standard deviation) computed from fourier transform)
  vpMatrix sigmaMap;

  //! spectral ratio map computed from fourier transform)
  vpMatrix spectMap;

  //! test mean value of spectral ratio
  double spect_rat_mean;


  //! Store the image (as a vector with intensity and gradient I, Ix, Iy) 
  npLuminance *pixInfo ;
  int  firstTimeIn  ;

  bool denoise; //flag for remove the charging of electroin (white part, 255)
  int noiseLevel; // tolerance of noise level, 254 maximum
  std::list<int> noisePosition; // list to save the position of the noise

  vpImage<unsigned char> Ic; // current image
  //vpImage<unsigned char> Id; // desired image

 public:

  vpCameraParameters cam ;
  vpImage<double> imG, imIx, imIy, imGxy, imNoise;

  double S_c; //current feature
  double S_p; //previous feature

  int nbrNoise; // number of pixels that is noise

public: 
  typedef enum {
       perspective,
       parallel,
       parallelZ
   }projectionModel;

  projectionModel pjModel;

public:
  typedef enum {
       ImageGradient,
       NormalizedVariance,
       SigmaGaussian,
       Gradient_Sigma_Variance,//compute the three cost function, but not interaction matrix, only for test cost function
       Regression
   }controlLawZ;

  controlLawZ vf_Z;


  void init() ;
  void init(unsigned int _nbr, unsigned int _nbc, double _Z, projectionModel projModel, controlLawZ vfeature_Z = ImageGradient) ;

  npFeatureLuminance() ;
 
  //! Destructor.
  virtual ~npFeatureLuminance()  ;

 public:

  void buildFrom(vpImage<unsigned char> &I) ;

  void setCameraParameters(vpCameraParameters &_cam)  ;

  void setZ(const double Z_fl) ;
  void setSigma(const double sigma);
  double getZ() const  ;

  vpMatrix getLg();
  vpColVector getJg() {return Jg;}
  vpColVector getSg();
  vpMatrix getLdg();
  double getGlobleSg();
  double getGlobleSnv();
  double getVariance();
  double getSpectralRatio();
  double getPhaseFactor();
  vpMatrix getVarianceMap();
  vpMatrix getSigmaMap();
  vpMatrix getSpectralRatioMap();

  double computeImageGradient(const vpImage<unsigned char> &I);
  void computeVarianceFourier(const vpImage<unsigned char> &Id);
  void computeGradientJacobian(const vpImage<unsigned char> &Id);

  vpMatrix  interaction(const unsigned int select = FEATURE_ALL);
  void      interaction(vpMatrix &L );

  vpColVector error(const vpBasicFeature &s_star, const unsigned int select = FEATURE_ALL)  ;
  void error(const vpBasicFeature &s_star, vpColVector &e)  ;
  //! Compute the error between a visual features and zero
  vpColVector error(const unsigned int select = FEATURE_ALL)  ;

  void sgError(const vpColVector &sg_star,
             vpColVector &eg)  ;

  void print(const unsigned int select = FEATURE_ALL ) const ;

  npFeatureLuminance *duplicate() const ;

  void display(const vpCameraParameters &cam,
               const vpImage<unsigned char> &I,
               const vpColor &color=vpColor::green, unsigned int thickness=1) const ;
  void display(const vpCameraParameters &cam,
               const vpImage<vpRGBa> &I,
               const vpColor &color=vpColor::green, unsigned int thickness=1) const ;

  void  matrixConvert(vpMatrix &src, cv::Mat &dest);
  void  matrixConvert(vpImage<double> &src, cv::Mat &dest);
  void  matrixConvert(cv::Mat& src, vpMatrix& dest);
  void  matrixConvert(cv::Mat& src, vpImage<double> & dest);


} ;

#endif
