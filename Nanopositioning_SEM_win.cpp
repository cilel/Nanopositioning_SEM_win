/****************************************************************************/

// The main fonction of Nanopositioning

/****************************************************************************/
//#include <math.h>

#include <visp/vpDebug.h>

#include <visp/vpImage.h>
#include <visp/vpImageIo.h>
#include <visp/vpImageTools.h>

#include <visp/vpCameraParameters.h>
#include <visp/vpTime.h>
#include <visp/vpVelocityTwistMatrix.h>

#include <visp/vpMath.h>
#include <visp/vpHomogeneousMatrix.h>
#include <visp/vpDisplayGTK.h>
#include <visp/vpDisplayGDI.h>
#include <visp/vpDisplayOpenCV.h>
#include <visp/vpDisplayD3D.h>
#include <visp/vpDisplayX.h>

#include "npFeatureLuminance.h"
#include "npRegression.h"
#include "npImageProcessing.h"

#include <visp/vpParseArgv.h>

//#include <stdlib.h>

#include <visp/vpParseArgv.h>
#include <visp/vpIoTools.h>
#include <visp/vpPlot.h>

#include <visp/vpNoise.h>
#include <visp/vpImageFilter.h>
#include <visp/vpImageConvert.h>
#include <visp/vpExponentialMap.h>

#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>


#define PORT 1085
#define SRV_IP "127.0.0.1"
#define UDP_VIDEO_IP "127.0.0.1"

#define BUFFER_IMG 49153 //: buffer ayant la taille de l'image en nb de pixels:786433 Buffer pour l'image

#define MU 0.04

// List of allowed command line options
#define GETOPTARGS	"cdi:n:hp:bu:af"

using namespace std;

typedef enum {
    perspective,
    parallel, //donot control translation along Z
    parallelZ,// control translation along Z
}projectionModel;

projectionModel pjModel;

typedef enum {
    Gauss_statistic,// gauss noise
    Gauss_dynamic // gauss noise in real time
}noiseModel;

noiseModel nsModel;


/*!

  Set the program options.

  \param argc : Command line number of parameters.
  \param argv : Array of command line parameters.
  \param ipath : Input image path.
  \param click_allowed : Mouse click activation.
  \param display : Display activation.
  \param niter : Number of iterations.

  \return false if the program has to be stopped, true otherwise.

*/
bool getOptions(int argc, const char **argv, std::string &ipath,
                bool &click_allowed, bool &display, int &niter, bool &add_noise,bool &blur, bool &denoise, bool &savevideo, double &scale)
{
    const char *optarg;
    int c;
    while ((c = vpParseArgv::parse(argc, argv, GETOPTARGS, &optarg)) > 1) {

        switch (c) {
        case 'c': click_allowed = false; break;
        case 'd': display = false; break;
        case 'i': ipath = optarg; break;
        case 'n': niter = atoi(optarg); break;
        case 'a': add_noise = true;  break;
        case 'b': blur = true;  break;//for parallelz
        case 'f': denoise = true;  break;
        case 'v': savevideo = true;  break;
        case 'p':
            if(!strcmp( optarg, "PERS" ))
                pjModel = perspective;
            else if (!strcmp( optarg,"PARA" ))
                pjModel = parallel;
            else if (!strcmp( optarg,"PARZ" ))
                pjModel = parallelZ;
            else
                pjModel = perspective;
            break;
        case 'u':
            if(!strcmp( optarg, "dm" ))
                scale = 10;
            else if (!strcmp( optarg,"cm" ))
                scale = 100;
            else if (!strcmp( optarg,"mm" ))
                scale = 1000;
            else if (!strcmp( optarg,"um" ))
                scale = 1000000;
            else if (!strcmp( optarg,"nm" ))
                scale = 1000000000;
            else // m
                scale = 1;
            break;
        default:
            return false; break;
        }
    }

    if ((c == 1) || (c == -1)) {
        // standalone param or error
        std::cerr << "ERROR: " << std::endl;
        std::cerr << "  Bad argument " << optarg << std::endl << std::endl;
        return false;
    }

    return true;
}

void diep(const char *s)
{
    perror(s);
    exit(1);
}

int send_wMe(vpHomogeneousMatrix wMe, double scale)
{

    vpTranslationVector T;
    wMe.extract(T);
    vpThetaUVector R;
    wMe.extract(R);

    //cout << "Translation=\n"<<T << "\nRotation=\n" << R << endl;

    /* diep(), #includes and #defines like in the server */

    double Pose_send[6];
    Pose_send[0]=T[0]*1000/scale;//convert meter to milimeter
    Pose_send[1]=T[1]*1000/scale;
    Pose_send[2]=T[2]*1000/scale;
    Pose_send[3]=R[0];
    Pose_send[4]=R[1];
    Pose_send[5]=R[2];

    cout << "Pose_send(wMe)=\n" << Pose_send[0]<< " " << Pose_send[1]<< " " << Pose_send[2]<< " " << Pose_send[3]<< " " << Pose_send[4]<< " " << Pose_send[5] << endl;

    struct sockaddr_in si_other;
    int s, i, slen=sizeof(si_other);

    if ((s=socket(AF_INET, SOCK_DGRAM, IPPROTO_UDP))==-1)
        diep("socket");

    memset((char *) &si_other, 0, sizeof(si_other));
    si_other.sin_family = AF_INET;
    si_other.sin_port = htons(PORT);
    if (inet_aton(SRV_IP, &si_other.sin_addr)==0) {
        fprintf(stderr, "inet_aton() failed\n");
        exit(1);
    }

    if (sendto(s, Pose_send, sizeof(Pose_send), 0, (struct sockaddr *) &si_other, slen)==-1)
        diep("sendto()");
    close(s);

    vpTime::wait(600);

    return 0;

}

void getFilteredImage(vpImage<unsigned char> &Iout, vpImage<unsigned char> Iin)
{
    cv::Mat Ii,Io;
    vpImageConvert::convert(Iin,Ii);
    cv::medianBlur(Ii, Io, 3);

    /*  cv::Size ksize;
    ksize.width = 9;
    ksize.height = 9;
    cv::GaussianBlur(Io, Io, ksize, 0.5);*/
    vpImageConvert::convert(Io,Iout);
}

void getDividedImage(vpImage<unsigned char> &Iout, vpImage<unsigned char> Iin)
{
    vpImage<unsigned char> Io;
    Iin.halfSizeImage(Io);
   // Iin.quarterSizeImage(Io);
    Iout=Io;
}

int
main(int argc, const char ** argv)
{
    std::string opt_ipath;
    std::string ipath;
    std::string filename;
    std::string readFileFlag;
    bool opt_click_allowed = true;
    bool opt_display = true;
    bool denoise = false;
    bool savevideo = true;
    bool divide_image = false;
    pjModel = parallelZ;  // here I choose parallel projection model on considering Z
    nsModel = Gauss_dynamic;

    // SOCKET sock;    // déclaration du socket pour la communication par TCP/IP
    int sock; // here I use int instead of Socket since I can't find a class SOCKET
    char bufferImg[BUFFER_IMG];    //déclaration du buffer pour stocker l'image

    npFeatureLuminance::controlLawZ clZ = npFeatureLuminance::ImageGradient;

    double scale = 1; // specific unit, 1 means in meter, 1e3 means in mm...

    // a folder "Result" TO BE created outside of built folder
    ofstream filecMo,fileVelociy,fileResidu,fileodMo,fileSg;
    filecMo.open ("../Result/Trajectory.txt");
    fileVelociy.open("../Result/Velocity.txt");
    fileResidu.open("../Result/Residual.txt");
    fileodMo.open("../Result/odMo.txt");
    fileSg.open("../Result/Sg.txt");

    // Read the command line options
/*    if (getOptions(argc, argv, opt_ipath, opt_click_allowed,
                   opt_display, opt_niter,add_noise,blur,denoise,savevideo,scale) == false) {
        return (-1);
    }
*/

    double Z = 0.01;
    double sigma =1;

    // Get the option values
    if (!opt_ipath.empty())
        ipath = opt_ipath;


    // Test if an input path is set
    if (opt_ipath.empty() && ipath.empty()){
        std::cerr << std::endl
                  << "ERROR:" << std::endl;
        std::cerr << "  Use -i <visp image path> option or set VISP_INPUT_IMAGE_PATH "
                  << std::endl
                  << "  environment variable to specify the location of the " << std::endl
                  << "  image path where test images are located." << std::endl << std::endl;
        exit(-1);
    }

    vpHomogeneousMatrix cMo,cMod,wMe,eMo,cMw,wMcR,wMc,wMo,Tr,cMe;

    // -------------------------Relative sample position, TO BE modified  -------------------------
    wMcR.buildFrom(0*0.001*scale,0*0.001*scale,2*0.01*scale,vpMath::rad(0),vpMath::rad(0),vpMath::rad(0)); // camera position
    wMe.buildFrom(0*0.001*scale,0*0.001*scale,0*0.01*scale,0,0,0); // center of robot, here the world frame is defined as the init pose of robot
    wMo.buildFrom(0*0.001*scale,0*0.001*scale,0*0.01*scale,0,0,0); // center of sample

    cout << "scale=" << scale << endl;
    cout << "wMe_desired=\n" << wMe << endl;

    Tr.buildFrom(0,0,0,vpMath::rad(180),0,0);
    wMc = wMcR*Tr;

    cout << "camera pose wMc=\n" << wMc << endl;

    cMw = wMc.inverse();
    eMo = wMe.inverse() * wMo;
    //eMo.setIdentity();

    cMod = cMw * wMo;
    cMe = cMw * wMe;

    send_wMe(wMe,scale);// to ensure the init pose of plateform

    vpHomogeneousMatrix cModR;

    cModR = wMcR.inverse()*wMo;

/*    //test
    vpThetaUVector R_cMod;
    cModR.extract(R_cMod);
    cout<< "Rot_cMod in deg=" << vpMath::deg(R_cMod[0])<< "\t"<< vpMath::deg(R_cMod[1])<< "\t"<< vpMath::deg(R_cMod[2])<<"\n";
*/

    cout << "cMod=\n"<< cMod << endl;

    // ----------------move to desired pose -----------------
    vpColVector vm;//velocity to be sent for the init pose, in meter or specified unit
    vm.resize(6);
    // vm[0]= 10e-6*scale; // Tx, 10e-6: 10 um
    // vm[1]= 5e-6*scale; // Ty
    //vm[2]= 5e-6*scale; // Tz
    // vm[3]=vpMath::rad(-0.1); // Rx
    // vm[4]=vpMath::rad(0.1); // Ry
   //  vm[5]=vpMath::rad(2); // Rz

    wMe = wMe * vpExponentialMap::direct(vm,1);
    cMo = cMw * wMe * eMo ;

    vpHomogeneousMatrix wMe_ogn = wMe; // save original wMe

    cout << "vpExponentialMap::direct(vm,1)=\n" << vpExponentialMap::direct(vm,1) << endl;

    cout << "wMe_desired=\n" << wMe << endl;
    cout << "cMe_desired=\n" << cMe << endl;
    cout << "cMo_desired=\n" << cMo << endl;

    send_wMe(wMe,scale);

    vpImage<unsigned char> I,Id;
    filename = ipath;
    const char * filename_c = filename.c_str();

    vpTime::wait(50);

    // -----------------------------get image from c#--------------------------
   send(sock,"start",5,0);         // demander une acquisition
   recv(sock, bufferImg, BUFFER_IMG, 0);    // réception de l'image
   cv::Mat image(360, 360,CV_8UC1, bufferImg);    // passage de l'image à OpenCV

   // convert the Mat to VpImage
   vpImageConvert::convert(image,Id);

    if (divide_image)
        getDividedImage(Id,Id);
    if(denoise)
        getFilteredImage(Id,Id);

    int Iw, Ih;//size of image
    Iw = Id.getWidth();
    Ih = Id.getHeight();

    //-----------------------------Calibration parameters, TO BE modified----------------------
    vpCameraParameters cam;
    if((pjModel == parallel) || ( pjModel == parallelZ))
        cam.initPersProjWithoutDistortion(100000/scale, 100000/scale, (int)Iw/2,(int)Ih/2);
    else
        cam.initPersProjWithoutDistortion(4000, 4000, (int)Iw/2, (int)Ih/2);


    // desired visual feature built from the image
    npFeatureLuminance sId ;
    if(pjModel == parallel)
        sId.init( Id.getHeight(), Id.getWidth(), cMod[2][3], npFeatureLuminance::parallel) ;
    else if(pjModel == parallelZ)
        sId.init( Id.getHeight(), Id.getWidth(), cMod[2][3], npFeatureLuminance::parallelZ,clZ) ;
    else
        sId.init( Id.getHeight(), Id.getWidth(), cMod[2][3], npFeatureLuminance::perspective) ;
    sId.setCameraParameters(cam) ;
    sId.setSigma(1);//just a value, not useful
    sId.buildFrom(Id);
    double Sgd_sum=sId.getGlobleSg();
   // double s_d = sId.computeImageGradient(Id);
    double Snvd = sId.getGlobleSnv();

    // display the image
#if defined VISP_HAVE_X11
    vpDisplayX d;
#elif defined VISP_HAVE_GDI
    vpDisplayGDI d;
#elif defined VISP_HAVE_GTK
    vpDisplayGTK d;
#elif defined VISP_HAVE_OPENCV
    vpDisplayOpenCV d;
#endif

#if defined(VISP_HAVE_X11) || defined(VISP_HAVE_GDI) || defined(VISP_HAVE_GTK) || defined(VISP_HAVE_OPENCV)
    if (opt_display) {
        d.init(Id, 1620, 10, "Photometric VS desired feature : s*") ;
        vpDisplay::display(Id);
        vpDisplay::flush(Id);
    }
    if (opt_display && opt_click_allowed) {
        std::cout << "Click in the image to continue..." << std::endl;
        vpDisplay::getClick(Id) ;
    }
#endif

    /*----------------------Here move the plateform for moving back --------------*/


    vm.resize(6);
    // vm[0]= 10e-6*scale; // Tx
     vm[1]= 5e-6*scale; // Ty
    //vm[2]= 5e-6*scale; // Tz
    // vm[3]=vpMath::rad(-0.1); // Rx
    // vm[4]=vpMath::rad(0.1); // Ry
   //  vm[5]=vpMath::rad(2); // Rz

    // vpHomogeneousMatrix wMe_tmp =  wMe * vpExponentialMap::direct(vm,1);

    //  cout << "wMe_tmp=\n" <<  wMe_tmp << endl;

    wMe = wMe_ogn * vpExponentialMap::direct(vm,1);
    cMo = cMw * wMe * eMo ;

    cout << "vpExponentialMap::direct(vm,1)=\n" << vpExponentialMap::direct(vm,1) << endl;

    cout << "wMe_first=\n" << wMe << endl;
    cout << "cMe_first=\n" << cMe << endl;
    cout << "cMo_first=\n" << cMo << endl;

    send_wMe(wMe,scale);

    I.resize(0,0);

     // -----------------------------get image from c#--------------------------
    send(sock,"start",5,0);         // demander une acquisition
    recv(sock, bufferImg, BUFFER_IMG, 0);    // réception de l'image
    cv::Mat img_temp(360, 360,CV_8UC1, bufferImg);    // passage de l'image à OpenCV
    image = img_temp;
    // convert the Mat to VpImage
    vpImageConvert::convert(img_temp,I);

    if (divide_image)
        getDividedImage(I,I);
    if(denoise)
        getFilteredImage(I,I);


#if defined(VISP_HAVE_X11) || defined(VISP_HAVE_GDI) || defined(VISP_HAVE_GTK) || defined(VISP_HAVE_OPENCV)
    if (opt_display) {
        d.init(I, 1620, 10, "s") ;
        vpDisplay::display(I) ;
        vpDisplay::flush(I) ;
    }
    if (opt_display && opt_click_allowed) {
        std::cout << "Click in the image to continue..." << std::endl;
        vpDisplay::getClick(I) ;
    }
#endif  

    vpImage<unsigned char> Idiff ;
    Idiff = I ;

    vpImageTools::imageDifference(I,Id,Idiff) ;

    // Affiche de l'image de difference
#if defined VISP_HAVE_X11
    vpDisplayX d1;
#elif defined VISP_HAVE_GDI
    vpDisplayGDI d1;
#elif defined VISP_HAVE_GTK
    vpDisplayGTK d1;
#elif defined VISP_HAVE_OPENCV
    vpDisplayOpenCV d1;
#endif
#if defined(VISP_HAVE_X11) || defined(VISP_HAVE_GDI) || defined(VISP_HAVE_GTK) || defined(VISP_HAVE_OPENCV)
    if (opt_display) {
        d1.init(Idiff, 1640+(int)I.getWidth(), 10, "s-s* ") ;
        vpDisplay::display(Idiff) ;
        vpDisplay::flush(Idiff) ;
    }
#endif
    // create the robot (here a simulated free flying camera)
    // ------------------------------------------------------
    // Visual feature, interaction matrix, error
    // s, Ls, Lsd, Lt, Lp, etc
    // ------------------------------------------------------

    // current visual feature built from the image
    // (actually, this is the image...)

    npFeatureLuminance sI ;
    if(pjModel == parallel)
        sI.init( I.getHeight(), I.getWidth(), cMod[2][3], npFeatureLuminance::parallel) ;
    else if(pjModel == parallelZ)
        sI.init( I.getHeight(), I.getWidth(), cMod[2][3], npFeatureLuminance::parallelZ,clZ) ;
    else
        sI.init( I.getHeight(), I.getWidth(), cMod[2][3], npFeatureLuminance::perspective) ;
    sI.setCameraParameters(cam) ;
    sI.setSigma(sigma);
    sI.buildFrom(I) ;


    vpImage<unsigned char> Idip;

#if defined VISP_HAVE_X11
    vpDisplayX ds;
#elif defined VISP_HAVE_GDI
    vpDisplayGDI ds;
#elif defined VISP_HAVE_GTK
    vpDisplayGTK ds;
#elif defined VISP_HAVE_OPENCV
    vpDisplayOpenCV ds;
#endif

    if(pjModel == parallelZ)
    {
        vpImageConvert::convert(sI.imGxy,Idip);
#if defined(VISP_HAVE_X11) || defined(VISP_HAVE_GDI) || defined(VISP_HAVE_GTK)|| defined(VISP_HAVE_OPENCV)
        {
            ds.init(Idip, 2680, 10, "image gradient") ;
            vpDisplay::display(Idip);
            vpDisplay::flush(Idip);
        }
#endif
    }


    // Matrice d'interaction, Hessien, erreur,...
    vpMatrix Lsd;   // matrice d'interaction a la position desiree
    vpMatrix Hsd;  // hessien a la position desiree
    vpMatrix H ; // Hessien utilise pour le levenberg-Marquartd
    vpColVector error ; // Erreur I-I*, photometric information

    // Compute the interaction matrix
    // link the variation of image intensity to camera motion
    // here it is computed at the desired position

    vpMatrix Lgsd; // interaction matrix (using image gradient) of the desired position
    vpMatrix Hgsd; // hessien of the desired position (using image gradient)
    vpMatrix Hg;  // Hessien for  levenberg-Marquartd (using image gradient)
    vpColVector sg_error; // error sg-sg*, image gradientmm


    vpVelocityTwistMatrix cVw;
    cVw.buildFrom(cMw);
    cout << "cVw=\n" << cVw <<endl;

    vpMatrix Js;// visual feature Jacobian
    vpMatrix Jn;// robot Jacobian
    vpMatrix diagHsd;// diag(Hsd)

    //For parallelZ, image gradient
    vpMatrix Jgs;// visual feature Jacobian
    vpMatrix Jgn;// robot Jacobian
    vpMatrix diagHgsd;// diag(Hsd)
    vpMatrix Jgsd;

    double vard = 0; // square of desired sigma

    sId.interaction(Lsd) ;
    cout << "Size of Lsd:" << Lsd.getRows() << "x" << Lsd.getCols() <<endl;
    Lgsd = sId.getLg();
    cout << "Size of Lgsd:" << Lgsd.getRows() << "x" << Lgsd.getCols() <<endl;

    if(pjModel==parallel )
    {
        Jn.resize(6,5);
        Jn[0][0]=1;
        Jn[1][1]=1;
        Jn[3][2]=1;
        Jn[4][3]=1;
        Jn[5][4]=1;


        //sI.interaction(Lsd) ;
        Js=-Lsd*cVw*Jn;
        //cout << "Lsd*cVw=\n" << Lsd*cVw << endl;

        //cout << "Jn=\n" << Jn << endl;
        cout << "Size of Jn:" << Jn.getRows() << "x" << Jn.getCols() <<endl;

        cout << "Js[0][0]=" << Js[0][0] <<endl;
        cout << "Size of Js:" << Js.getRows() << "x" << Js.getCols() <<endl;

        // Compute the Hessian H = L^TL
        Hsd = Js.AtA() ;

        //cout << "Hsd=\n" << Hsd <<endl;
        // Compute the Hessian diagonal for the Levenberg-Marquartd
        // optimization process
        unsigned int n = 5 ;
        diagHsd.resize(n,n) ;
        diagHsd.eye(n);
        for(unsigned int i = 0 ; i < n ; i++) diagHsd[i][i] = Hsd[i][i];
    }
    else if (pjModel==parallelZ)
    {
        Jn.resize(6,5);
        Jn[0][0]=1;
        Jn[1][1]=1;
        Jn[3][2]=1;
        Jn[4][3]=1;
        Jn[5][4]=1;

        Js=-Lsd*cVw*Jn;
        cout << "Js[0][0]=" << Js[0][0] <<endl;
        cout << "Size of Jn:" << Jn.getRows() << "x" << Jn.getCols() <<endl;
        cout << "Size of Js:" << Js.getRows() << "x" << Js.getCols() <<endl;


        // Compute the Hessian H = L^TL
        Hsd = Js.AtA() ;

        // Compute the Hessian diagonal for the Levenberg-Marquartd
        // optimization process
        unsigned int n = 5 ;
        diagHsd.resize(n,n) ;
        diagHsd.eye(n);
        for(unsigned int i = 0 ; i < n ; i++) diagHsd[i][i] = Hsd[i][i];

        // use fourrier based method to compute variance (sigma^2) of gaussian kernel
        sI.computeVarianceFourier(I);
        vard = sI.getVariance();

        cout << "vard = " << vard <<endl;

    }
    else //perspective
    {
        Jn.resize(6,6);
        Jn.setIdentity();
        Js=-Lsd*cVw*Jn;

        cout << "Size of Js:" << Js.getRows() << "x" << Js.getCols() <<endl;
        //cout << "Js=\n" << Js <<endl;

        // Compute the Hessian H = L^TL
        Hsd = Js.AtA() ;

        cout << "Hsd=\n" << Hsd <<endl;

        // Compute the Hessian diagonal for the Levenberg-Marquartd
        // optimization process
        unsigned int n = 6 ;
        diagHsd.resize(n,n) ;
        diagHsd.eye(n);
        for(unsigned int i = 0 ; i < n ; i++) diagHsd[i][i] = Hsd[i][i];
    }


    vpPlot graphy(4, 600, 800, 2420, 1300, "Nanopositioning");
    vpPlot graphy2(4, 600, 800, 020, 1300, "Nanopositioning");

    graphy.initGraph(0,3);
    graphy.initGraph(1,3);
    graphy.initGraph(2,3);
    graphy.initGraph(3,3);

    graphy.setTitle(0,"odMo: m");
    graphy.setTitle(1,"odMo: deg");
    graphy.setTitle(2,"Velocity: m/s ");
    graphy.setTitle(3,"Velocity: deg/s");


    char unit[40];
    graphy2.initGraph(0,1);
    graphy2.initGraph(1,1);
    graphy2.initGraph(2,1);
    graphy2.setTitle(0,"Norm Error");
    if (pjModel==parallelZ)
    {
        graphy2.setTitle(1,"Gradient (cost function)");
        graphy2.setTitle(2,"regression");
    }
    else
    {
        graphy2.setTitle(1,"Trajectory of object");
        strncpy( unit, "x", 40 );
        graphy2.setUnitX(1,unit);
        strncpy( unit, "y", 40 );
        graphy2.setUnitY(1,unit);
        strncpy( unit, "z", 40 );
        graphy2.setUnitZ(1,unit);
        graphy2.setColor(0,0,vpColor::red);
    }


    char legend[40];
    strncpy( legend, "Norm Error", 40 );
    graphy2.setLegend(0,0,legend);
    strncpy( legend, "Tx", 40 );
    graphy.setLegend(0,0,legend);
    strncpy( legend, "Ty", 40 );
    graphy.setLegend(0,1,legend);
    strncpy( legend, "Tz", 40 );
    graphy.setLegend(0,2,legend);
    strncpy( legend, "Rx", 40 );
    graphy.setLegend(1,0,legend);
    strncpy( legend, "Ry", 40 );
    graphy.setLegend(1,1,legend);
    strncpy( legend, "Rz", 40 );
    graphy.setLegend(1,2,legend);
    strncpy( legend, "Tx", 40 );
    graphy.setLegend(2,0,legend);
    strncpy( legend, "Ty", 40 );
    graphy.setLegend(2,1,legend);
    strncpy( legend, "Tz", 40 );
    graphy.setLegend(2,2,legend);
    strncpy( legend, "Rx", 40 );
    graphy.setLegend(3,0,legend);
    strncpy( legend, "Ry", 40 );
    graphy.setLegend(3,1,legend);
    strncpy( legend, "Rz", 40 );
    graphy.setLegend(3,2,legend);

    // ------------------------------------------------------
    // Control law

    vpColVector e ;// velocity to be multiply by lamda
    vpColVector v ; // camera velocity send to the robot
    vpColVector eg; // velocity of z axis
    vpColVector vg ; // camera velocity of z axis
    double vgd;// camera velocity of z axis, double

    // ----------------------------------------------------------
    // Minimisation

    double mu,mu1,mu2,mu3 ;  // mu = 0 : Gauss Newton ; mu != 0  : LM
    double lambdaTZ  = 1;// gain for Translation on Z : lambdaTZ*lambda
    double lambda,lambda1, lambda2, lambda3 ; //gain

    mu1       =  0.0001;
    mu2    =  0;
    mu3    =  0;
    lambda1= 1;
    lambda2 = 0.4*lambda2;
    lambda3   = 0.1*lambda2;

    mu = mu1;
    lambda = lambda1;

    // ----------------------------------------------------------
    int iter   = 0;
    int iter1 = 150 ; // swicth to Gauss Newton after iterGN iterations
    int iter2 = 250 ;

    double normError = 1000; // norm error = |I-I*|
    double normError_p = 0; // previous norm error
    double threshold=0.02;// condition of convergence
    //double convergence_threshold = 0.001;

    cout<< "threshold=" << threshold << endl;

    // vpHomogeneousMatrix edMe,e/*dMw ;
    vpHomogeneousMatrix odMo ;

    fileResidu <<"#projection model:" << pjModel << "\t init pose:" << vm.t() << endl;

    double Sg_previsous; // recode image gradient at 1st iteration to compute the sign of veloctiy
    int sign=1; // sign for velocity


    /******************************* Driving direction ************************************/

    std::cout << "------------------- Driving direction -------------------------" << std::endl ;

    std::cout << "///////// plus /////////////" << std::endl ;
    vm.resize(6);
    vm[2]=20e-6;

    wMe = wMe * vpExponentialMap::direct(vm,1);
    cMo = cMw * wMe * eMo;
    send_wMe(wMe,scale);

    // -----------------------------get image from c#--------------------------
   send(sock,"start",5,0);         // demander une acquisition
   recv(sock, bufferImg, BUFFER_IMG, 0);    // réception de l'image
   cv::Mat imgTempP(360, 360,CV_8UC1, bufferImg);    // passage de l'image à OpenCV
   image = imgTempP;
   // convert the Mat to VpImage
   vpImageConvert::convert(imgTempP,I);

    cout  <<"sigma=" << sigma << endl;
    cout << "cMo[2][3]=" << cMo[2][3] << endl;

    if (divide_image)
        getDividedImage(I,I);
    if(denoise)
        getFilteredImage(I,I);

    double s_plus = sI.computeImageGradient(I);
    double z_plus = cMo[2][3];

    std::cout << "///////// minus /////////////" << std::endl ;
    vm.resize(6);
    vm[2]=-40e-6;

    wMe = wMe * vpExponentialMap::direct(vm,1);
    cMo = cMw * wMe * eMo;
    send_wMe(wMe,scale);

    // -----------------------------get image from c#--------------------------
   send(sock,"start",5,0);         // demander une acquisition
   recv(sock, bufferImg, BUFFER_IMG, 0);    // réception de l'image
   cv::Mat imgTempM(360, 360,CV_8UC1, bufferImg);    // passage de l'image à OpenCV
   image = imgTempM;
   // convert the Mat to VpImage
   vpImageConvert::convert(imgTempM,I);

    cout  <<"sigma=" << sigma << endl;
    cout << "cMo[2][3]=" << cMo[2][3] << endl;

    if (divide_image)
        getDividedImage(I,I);
    if(denoise)
        getFilteredImage(I,I);

    double s_minus = sI.computeImageGradient(I);
    double z_minus = cMo[2][3];


    std::cout << "///////// origin /////////////" << std::endl ;
    vm.resize(6);
    vm[2]=20e-6;

    wMe = wMe * vpExponentialMap::direct(vm,1);
    cMo = cMw * wMe * eMo;
    send_wMe(wMe,scale);

    // -----------------------------get image from c#--------------------------
   send(sock,"start",5,0);         // demander une acquisition
   recv(sock, bufferImg, BUFFER_IMG, 0);    // réception de l'image
   cv::Mat imgTempO(360, 360,CV_8UC1, bufferImg);    // passage de l'image à OpenCV
   image = imgTempO;
   // convert the Mat to VpImage
   vpImageConvert::convert(imgTempO,I);

    cout  <<"sigma=" << sigma << endl;
    cout << "cMo[2][3]=" << cMo[2][3] << endl;

    if (divide_image)
        getDividedImage(I,I);
    if(denoise)
        getFilteredImage(I,I);

    double s_0 = sI.computeImageGradient(I);
    double z_0 = cMo[2][3];

    cout << "s_plus and minus: " <<  s_plus << " | " << s_minus << endl;

    if (s_minus>s_plus)
        sign = -1;

    /*if (var_minus < var_plus)
        sign_var = -1;
*/
    cout << "sign=" << sign << endl;

    vpColVector G_reg(3);
    vpColVector Z_reg(3);

    G_reg[0] = s_minus;
    G_reg[1] = s_0;
    G_reg[2] = s_plus;
    Z_reg[0] = z_minus;
    Z_reg[1] = z_0;
    Z_reg[2] = z_plus;


    double JG;
    double var; // square of sigma

    //regression
    vpColVector p; // parameters of rational function
    npRegression reg(Z_reg,G_reg,npRegression::Cauchy) ;

    // using image gradient, TO BE modified to adapt the real case
    double Sg_threshold = 1e-2;
    double Sg_gain = 1e-5;

    /******************************* Here the loop begins ************************************/
    do
    {
        std::cout << "--------------------------------------------" << iter++ << std::endl ;

        // -----------------------------get image from c#--------------------------
       send(sock,"start",5,0);         // demander une acquisition
       recv(sock, bufferImg, BUFFER_IMG, 0);    // réception de l'image
       cv::Mat imgTemp(360, 360,CV_8UC1, bufferImg);    // passage de l'image à OpenCV
       image = imgTemp;
       // convert the Mat to VpImage
       vpImageConvert::convert(imgTemp,I);

        cout << "cMo[2][3]=" << cMo[2][3] << endl;
        if (divide_image)
            getDividedImage(I,I);
        if(denoise)
            getFilteredImage(I,I);


#if defined(VISP_HAVE_X11) || defined(VISP_HAVE_GDI) || defined(VISP_HAVE_GTK) 
        if (opt_display) {
            vpDisplay::display(I) ;
            vpDisplay::flush(I) ;
        }
#endif
        vpImageTools::imageDifference(I,Id,Idiff) ;
#if defined(VISP_HAVE_X11) || defined(VISP_HAVE_GDI) || defined(VISP_HAVE_GTK) 
        if (opt_display) {
            vpDisplay::display(Idiff) ;
            vpDisplay::flush(Idiff) ;
        }
#endif

        // Compute current visual feature
        // sI.set_Z(cMo[m][3]);

        //sI.setSigma(sigma);
        sI.buildFrom(I);

        //using computed sigma
/*        sI.computeVarianceFourier(Id);
        double variance = sI.getVariance();
        double computedSigma = sqrt(fabs(variance));
        sI.setSigma(computedSigma);
        sI.buildFrom(I);
        sI.computeGradientJacobian(Id);
*/
        // compute current error
        sI.error(sId,error);

        double Sg_sum=0; // sum of square of norm of image gradient, cost function
        double Snv = 0; // globle normalized variance, cost function

        if (pjModel==parallelZ)
        {

            Sg_sum = sI.getGlobleSg();
            Snv = sI.getGlobleSnv();

            fileSg << iter << "\t" << Sg_sum << "\t" << (cMod[2][3]-cMo[2][3])*1e6;

            sI.sgError(sId.getSg(),sg_error) ;

        }

        normError_p = normError;

        normError = sqrt(error.sumSquare()/error.getRows());

        graphy2.plot(0,0,iter,normError);

        cout << "|e| "<< normError <<endl ;
        fileResidu << iter << "\t" << normError << endl;

        // ---------- Levenberg Marquardt method --------------
        if(iter > iter2)
        {
            mu = mu3 ;
            lambda = lambda3;
        }
        else if(iter > iter1)
        {
            mu = mu2 ;
            lambda = lambda2;
        }

        if(iter == 1)
        {
            vpImageIo::writePNG(Idiff,"../Result/img_init.png");
            fileResidu << "#\t mu:"<< mu << "\t lambda:" << lambda  << endl;
            Sg_previsous = Sg_sum;
        }
        else if (iter == iter1 || iter == iter2)
            fileResidu << "#\t mu:"<< mu << "\t lambda:" << lambda  << endl;


        // Compute the levenberg Marquartd term
        H = ((mu * diagHsd) + Hsd).inverseByLU();
        //	compute the control law
        e = H * Js.t() *error ;

        v =  -lambda*e;

        if(pjModel==parallelZ)
        {
            Jn.resize(6,5);
            Jn[0][0]=1;
            Jn[1][1]=1;
            Jn[3][2]=1;
            Jn[4][3]=1;
            Jn[5][4]=1;

            vpMatrix Lg;//interaction matrix for image gradient
            //vpMatrix Ldgs;
            vpMatrix Ld_temp; // temporar matrix for compute Lgs

            //vpMatrix Jdgs;

            sI.interaction(Ld_temp);

            Lg = sI.getLg();
            //cout << "Ldg: " << Ldg << endl;

            Jgn.resize(6,1);
            Jgn[2][0]=1;

            Jgs=-Lg*cVw*Jgn; //  derivative of visual feature (jacobian)

            // cout << "Jgs[0][0]="<< Jgs[0][0] << endl;

            vpRowVector Jg;
            Jg.resize(1);

            int nbrJgs=Jgs.getRows()-sI.nbrNoise;

            for (int m=0; m<Jgs.getRows();m++)
            {
                Jg[0]+=(Jgs[m][0]);//fabs
                // Jdg[0]+=(Jdgs[m][0]);
            }
            Jg[0] = Jg[0]/nbrJgs;

            //cout << "nbrJgs=" << nbrJgs << endl;

            //graphy2.plot(0,0,iter,Jg[0]);

            cout << "Jg=" << Jg << endl;
            // cout << "Jdg=" << Jdg << endl;

            //-----------------------------------

            cout << "Sg_sum and Sgd_sum:" <<  Sg_sum << " " << Sgd_sum << endl;
            if(iter== 1)
                JG = Jg[0];

            cout << "JG = " << JG << endl;

            Hgsd = Jgs.AtA();
            diagHgsd.resize(1,1);
            diagHsd[0][0] = Hgsd[0][0];


            double alpha = 1e-2;

        if (clZ == npFeatureLuminance::Regression)// regression
        {
            if(iter ==1)
            {
                //reg.linearRegression();
                reg.ridgeRegression();
            }
            else
            {
                reg.addXY(cMo[2][3],sI.computeImageGradient(I));
                cout << "G = "<< sI.computeImageGradient(I) <<endl;
                reg.ridgeRegression();
            }
            p = reg.getParameters();
            cout << "rational parameters:" << p.t() << endl;

            double z = cMo[2][3];

            double J_reg_inv = reg.jacobianInv(z);
            cout << "J_reg_inv : " <<J_reg_inv << endl;

            if (iter < 10)
                vgd = -sign * 10e-6 * (Sg_sum-Sgd_sum);
            else
                vgd = -alpha * J_reg_inv * (Sg_sum-Sgd_sum);

            fileSg<< "\t" << Jg[0] << "\t" << Lg.sumSquare() << "\t" << var <<endl;

            cout << "vgd=" << vgd << endl;
        }
        else if (clZ == npFeatureLuminance::ImageGradient)
        {
            if (fabs(Sg_sum-Sgd_sum)>Sg_threshold)
                vgd = -sign * Sg_gain * (Sg_sum-Sgd_sum);
            else
                vgd = 0;
        }


        }

        if(pjModel==parallel)
        {
            vpColVector vc=v;
            v.resize(6);
            v[5]=vc[4];
            v[4]=vc[3];//vc[3]
            v[3]=vc[2];//vc[2]
            v[2]=0;
            v[1]=vc[1];//vc[1]
            v[0]=vc[0];//vc[0]
        }
        else if (pjModel==parallelZ)
        {
            if (iter <= 2000) // here we first try 2 DOF case, Tx and Ty
            {
                vpColVector vc=v;
                v.resize(6);
           //     v[5]=vc[4];//0;
           //    v[4]=vc[3];
           //     v[3]=vc[2];
                //v[2]=vgd;
                v[1]=vc[1];
                v[0]=vc[0];
            }
            else
            {
 /*               vpColVector vc=v;
                v.resize(6);
                v[5]=vc[4];
                v[4]=vc[3];
                v[3]=vc[2];
               v[2]=vgd;
                v[1]=vc[1];
                v[0]=vc[0];
*/
            }

        }

   for(int i=0;i<3;i++)
            graphy.plot(2,i,iter,v[i]/scale);
        for(int i=0;i<3;i++)
            graphy.plot(3,i,iter,vpMath::deg(v[i+3]));

        //convert m/s, rad/s to um/s, rad/s
        vpColVector vPlot=v;
        for(int i=0;i<3;i++)
            vPlot[i]=v[i]*1e6;

        fileVelociy << iter << "\t" << vPlot.t() << endl;
        //
        cout << "v=" << v.t() << endl;
        cout << "lambda = " << lambda << "  mu = " << mu ;
        cout << " |Tc| = " << sqrt(v.sumSquare()) << endl;

        wMe = wMe * vpExponentialMap::direct(v,0.04);

    //    cout << "wMe2_current=\n" << wMe << endl;

        send_wMe(wMe,scale);

        cMo = cMw * wMe * eMo;

        filecMo << iter;
        for(int m=0;m<3;m++)
            filecMo << "\t" << cMo[m][3]*1e6; //convert m to um
        filecMo << endl;

        cout << "cMo_new=\n" << cMo << endl;

        odMo = cMod.inverse() * cMo;

        vpTranslationVector TodMo;
        odMo.extract(TodMo);
        vpThetaUVector RodMo;
        odMo.extract(RodMo);
        vpColVector TodMoPlot;
        TodMoPlot.resize(3);
        vpColVector RodMoPlot;
        RodMoPlot.resize(3);

        /*--------Here is plot in realtime and save data for gnuplot-------*/
        //convert m, rad to um, deg
        for(int i=0;i<3;i++)
            TodMoPlot[i]= TodMo[i]*1e6;
        for(int i=0;i<3;i++)
            RodMoPlot[i]= vpMath::deg(RodMo[i]);

        fileodMo << iter << "\t" << TodMoPlot.t() << RodMoPlot.t() << endl;

        for(int i=0;i<3;i++)
            graphy.plot(0,i,iter,TodMo[i]/scale);
        for(int i=0;i<3;i++)
            graphy.plot(1,i,iter,vpMath::deg(RodMo[i]));

        if(pjModel==parallelZ)
        {
            graphy2.plot(1,0,iter,Sg_sum);
        }
        else
            graphy2.plot(1,0,cMo[0][3]/scale,cMo[1][3]/scale,cMo[2][3]/scale);

        // cout << "cMoTT:" << cMo[0][3]/scale<<cMo[1][3]/scale<<cMo[2][3]/scale << endl;


        //simulation of regression
        if(pjModel==parallelZ && clZ == npFeatureLuminance::Regression)
        {
            vpColVector G_r(100);
            vpColVector Z_r(100);

            double Z_min = cMod[2][3] - 200e-6;

            for (int i = 0; i < 100 ; i++)
            {    
               Z_r[i] = Z_min + i*4e-6;
               G_r[i] = (p[0]*Z_r[i] *Z_r[i]+p[1]*Z_r[i]+p[2])/(Z_r[i]*Z_r[i]+p[3]*Z_r[i]+p[4]);
            }

            Z_min = Z_r.getMinValue();
            graphy2.initGraph(2,1);
            //graphy2.initRange(2,Z_min,Z_max,G_min,G_max);

            for (int i = 0; i < 100 ; i++)
                graphy2.plot(2,0,Z_r[i],G_r[i]);
        }

        /*------------Here begin to save the image for video--------*/
        if(savevideo==true)
        {
            char img_filename_c[80];
            sprintf(img_filename_c, "../Result/video/img_current_%d.png", iter);
            string img_filename(img_filename_c);
            vpImageIo::writePNG(I,img_filename);
            sprintf(img_filename_c, "../Result/video/img_diff_%d.png", iter);
            string img_filename_d(img_filename_c);
            vpImageIo::writePNG(Idiff,img_filename_d);
/*
            if(pjModel == parallelZ)
            {
                sprintf(img_filename_c, "../Result/video/img_gradient_%d.png", iter);
                string img_filename_g(img_filename_c);
                vpImageIo::writePNG(Idip,img_filename_g);
            }*/
        }
        /*-----------Here end to save the image for video---------*/

        if(iter > 200)
            vpImageIo::writePNG(Idiff,"../Result/img_end.png");

    }
    // while(normError > threshold  && iter < opt_niter);// && !(vpMath::equal(normError,normError_p, convergence_threshold) && normError < 1.5*threshold));
    while(1) ;

    vpImageIo::writePNG(Idiff,"../Result/img_end.png");

    filecMo.close();
    cout << "===========================END==============================" << endl;

    while(1)
        graphy2.plot(1,0,cMo[0][3]/scale,cMo[1][3]/scale,cMo[2][3]/scale);


    //vpDisplay::getClick(graphy.I);
    vpDisplay::getClick(graphy2.I);


}

