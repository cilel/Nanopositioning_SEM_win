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
#include <visp/vpImageConvert.h>
#include <visp/vpHomogeneousMatrix.h>


#define PORT 1085
#define SRV_IP "127.0.0.1" // to be changed to the server ip

using namespace std;
using namespace cv;

void diep(const char *s)
{
    perror(s);
    exit(1);
}

int send_velocity(vpColVector &v)
{

    double velocity_send[6];
    velocity_send[0]=v[0];
    velocity_send[1]=v[1];
    velocity_send[2]=v[2];
    velocity_send[3]=v[0];
    velocity_send[4]=v[1];
    velocity_send[5]=v[2];

    cout << "velocity_sent=\n" << velocity_send[0]<< " " << velocity_send[1]<< " " << velocity_send[2]<< " " << velocity_send[3]<< " " << velocity_send[4]<< " " << velocity_send[5] << endl;

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

    if (sendto(s, velocity_send, sizeof(velocity_send), 0, (struct sockaddr *) &si_other, slen)==-1)
        diep("sendto()");
    close(s);

    vpTime::wait(600);

    return 0;

}


int
main(int argc, const char ** argv)
{
    vpColVector velocity(6); // here we should check the UNIT for the velocity, in meter or mm or um?
    velocity[0] = 1; // translation along x
    velocity[1] = 0; // translation along y
    velocity[2] = 0; // translation along z
    velocity[3] = 0; // rotation along x
    velocity[4] = 0; // rotation along y
    velocity[5] = 0; // rotation along z

    while(1) // loop
    {
       // here can put some changement of velocity

       int rst = send_velocity(velocity);
       if( rst!=0)
       {
            cerr << "error in sent_velocity" <<endl;
       }
    }
}
