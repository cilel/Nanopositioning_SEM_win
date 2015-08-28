/// include std and socket

#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>

/// include visp
#include <visp/vpImage.h>
#include <visp/vpImageIo.h>
#include <visp/vpImageConvert.h>

#include <visp/vpConfig.h>
#include <visp/vpDebug.h>
#include <visp/vpDisplayX.h>
#include <visp/vpDisplayD3D.h>
#include <visp/vpDisplayGTK.h>
#include <visp/vpDisplayGDI.h>
#include <visp/vpDisplayOpenCV.h>
/// include opencv
#include <cv.hpp>

#define BUFFER_IMG 49153 //: buffer ayant la taille de l'image en nb de pixels:786433 Buffer pour l'image

using namespace std;
using namespace cv;

int
main(int argc, const char ** argv)
{
// SOCKET sock;    // déclaration du socket pour la communication par TCP/IP
int sock; // here I use int instead of Socket since I can't find a class SOCKET
char bufferImg[BUFFER_IMG];    //déclaration du buffer pour stocker l'image

// !! I think there maybe some socket configuration to do, eg, port and ip adress. These should be added.

// to display the image from visp
#if defined VISP_HAVE_X11
vpDisplayX d;
cout << "x11" << endl;
#elif defined VISP_HAVE_GDI
vpDisplayGDI d;
cout << "GDI" << endl;
#elif defined VISP_HAVE_GTK
vpDisplayGTK d;
cout << "GTK" << endl;
#else
cout << "nothing to display" << endl;
#endif

vpImage<unsigned char> I;
I.resize(360,360);

#if defined(VISP_HAVE_X11) || defined(VISP_HAVE_GDI) || defined(VISP_HAVE_GTK)
    {
            d.init(I, 1000, 100, "input image") ;
            vpDisplay::display(I) ;
            vpDisplay::flush(I) ;
        }
#endif

// here we use a loop to show continues image
long iter=0;
while(1)
{

     // get image from c#
    send(sock,"start",5,0);         // demander une acquisition
    recv(sock, bufferImg, BUFFER_IMG, 0);    // réception de l'image
    Mat image(360, 360,CV_8UC1, bufferImg);    // passage de l'image à OpenCV

    // convert the Mat to VpImage
    vpImageConvert::convert(image,I);

    // display the image
    vpDisplay::display(I) ;
    vpDisplay::flush(I) ;

iter++;

}
    vpDisplay::getClick(I);
}

