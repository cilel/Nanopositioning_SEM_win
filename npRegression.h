#include <visp/vpDebug.h>
#include <visp/vpException.h>
#include <visp/vpMath.h>
#include <visp/vpHomogeneousMatrix.h>
#include <visp/vpHomography.h>

class npRegression
{
public:
  typedef enum {
       Gaussian,
       Cauchy, // particular case of Rational, 1/(p0x^2+p1x+p2)
       Polynomial, //Quadratic (p0x^2+p1x+p2)
       Rational // Quadratic (p0x^2+p1x+p2)/(x^2+p3x+p4)
   }regressionFunction;
   regressionFunction rgModel;

public:

    void init(const vpColVector &x0,const vpColVector &y0,const regressionFunction rgM);
    void init();
    npRegression();
    npRegression(const vpColVector &x0,const vpColVector &y0,const regressionFunction rgM);
    void setRegressionFunction(const regressionFunction rgM);

    void linearRegression();
    void ridgeRegression(const double mu = 1e-3);
    double jacobianInv(double x0);
    double jacobian(double x0);
    static double jacobian(double x0, vpColVector p0, regressionFunction rg = Rational);
    double hessianInv(double x0);
    double hessian(double x0);
    static double hessian(double x0,vpColVector p0, regressionFunction rg= Rational);
    double reponse(double x0);
    double getMaxPosition();

    void setXY(const vpColVector &x0, const vpColVector &y0);
    void addXY(const double x0,const  double y0);

    vpColVector getParameters();

    void nonlinearRegression(const vpColVector &G, const vpColVector &Z, double &b, double &c);

private:
    vpColVector x;
    vpColVector y;
    vpColVector p; // paramaters





};
