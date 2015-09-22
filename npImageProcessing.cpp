#include "npImageProcessing.h"

void npImageProcessing::getGaussianBlur(vpImage<unsigned char> &Iblur, const vpImage<unsigned char> &I, const double sigma, const int size, bool normalize)
{

    double *fg=new double[(size+1)/2] ;
    getGaussianKernel(fg, size, sigma, normalize) ;
   // std::cout << "fg[0] = " << fg[0] <<std::endl;
    if(fg[0]<0)
    {
        vpImage<unsigned char> Ib  = I;
        Iblur = Ib;
    }
    else
    {
    vpImage<double> GIx;
    vpImage<double> GI ;
    vpImageFilter::filterX(I, GIx,fg,size);
    vpImageFilter::filterY(GIx, GI,fg,size);

    Iblur.init(I.getHeight(),I.getWidth());

    for(int i=0;i<Iblur.getNumberOfPixel();i++)
            Iblur.bitmap[i] = floor(GI.bitmap[i]+0.5);

    GIx.destroy();
    delete[] fg;
    }
}

void npImageProcessing::getGaussianKernel(double *filter, unsigned int size, double sigma, bool normalize)
{
  if (size%2 != 1)
    throw (vpImageException(vpImageException::incorrectInitializationError,
          "Bad Gaussian filter size"));

  int middle = (int)(size-1)/2;
  if (sigma<= 0)
  {
      filter[0] = -1;
      for( int i=1; i<= middle; i++)
      {
        filter[i] = 0;
      }
  }
  else
  {

      double sigma2 = vpMath::sqr(sigma);
      for( int i=0; i<= middle; i++)
      {
        filter[i] = (1./(sigma*sqrt(2.*M_PI)))*exp(-(i*i)/(2.*sigma2));
      }
      if (normalize) {
        //renormalization
        double sum=0;
        for(int i=1; i<=middle; i++)
        {
          sum += 2*filter[i] ;
        }
        sum += filter[0];

        for(int i=0; i<=middle; i++)
        {
          filter[i] = filter[i]/sum;
        }
      }
  }
}


void npImageProcessing::getBlurImageCV(vpImage<unsigned char> &Iblur, const vpImage<unsigned char> &Iin, const double sigma)
{
    cv::Mat I,Ib;
    vpImageConvert::convert(Iin,I);
    cv::Size ksize;
    ksize.width = 15;
    ksize.height = 15;

    if(sigma>0)
    {
        cv::GaussianBlur(I, Ib, ksize,sigma);
        vpImageConvert::convert(Ib,Iblur);
    }
    else
        Iblur = Iin;
}

void npImageProcessing::getNoisedImage(vpImage<unsigned char> &Inoised, const vpImage<unsigned char> &Iin, const double noise_mean,const double noise_sdv)
{

    Inoised = Iin;
    //  Inoised.init(Itexture.getRows(),Itexture.getCols(),0);

    vpGaussRand noise(noise_sdv, noise_mean);
    for(int i=0;i<Iin.getNumberOfPixel();i++)
        {
            double gauss = noise();
            double noised =(double) Iin.bitmap[i] + gauss;
            if (noised < 0)
                Inoised.bitmap[i] = 0;
            else if (noised > 255)
                Inoised.bitmap[i] = 255;
            else
                Inoised.bitmap[i] = noised;
        }
}

void npImageProcessing::getSobelImage(vpImage<unsigned char> &IsX, vpImage<unsigned char> &IsY, const vpImage<unsigned char> &Iin)
{
    cv::Mat I;

/*  vpImage<unsigned char> Ib;
    getBlurImage(Ib,Iin,1);
*/
    vpImageConvert::convert(Iin,I);

    /// Generate grad_x and grad_y
      cv::Mat grad_x, grad_y;

      /// Gradient X
      //Scharr( src_gray, grad_x, ddepth, 1, 0, scale, delta, BORDER_DEFAULT );
      cv::Sobel( I, grad_x, -1, 1, 0);

      /// Gradient Y
      //Scharr( src_gray, grad_y, ddepth, 0, 1, scale, delta, BORDER_DEFAULT );
      cv::Sobel( I, grad_y, -1, 0, 1);

    vpImageConvert::convert(grad_x,IsX);
    vpImageConvert::convert(grad_y,IsY);

}

void npImageProcessing::cropImage(vpImage<unsigned char> &Icrop,const vpImage<unsigned char> &Iin,const vpImagePoint p1,const vpImagePoint p2)
{
    cv::Mat I;
    vpImageConvert::convert(Iin,I);

    double y1,x1,H,W;

    x1 = std::min(p1.get_j(),p2.get_j());
    y1 = std::min(p1.get_i(),p2.get_i());

    W = fabs(p1.get_j()-p2.get_j());
    H = fabs(p1.get_i()-p2.get_i());

    if((W * H == 0))
        std::cerr << "check crop image size" << std::endl;

    cv::Rect ROI(x1, y1, W, H);

    cv::Mat Ic = I(ROI);
    vpImageConvert::convert(Ic,Icrop);

}

void npImageProcessing::cropImage(vpImage<unsigned char> & Icrop, const vpImage<unsigned char> &Iin,const vpRect cropROI)
{
    cv::Mat I;
    vpImageConvert::convert(Iin,I);

//    std::cout << "cropROI=" << cropROI << std::endl;
    cv::Rect ROI(cropROI.getLeft(), cropROI.getTop(), cropROI.getWidth(), cropROI.getHeight());

    cv::Mat Ic2 = I(ROI);
    cv::Mat Ic;
    Ic2.copyTo(Ic);

//    //test
// /*    std::cout << Ic.size() << std::endl;

//     cv::imshow("Image", I);
//     cv::waitKey(1000);
//     cv::imshow("crop", Ic);
//     cv::waitKey(1000);
//*/
    vpImageConvert::convert(Ic,Icrop);


}
