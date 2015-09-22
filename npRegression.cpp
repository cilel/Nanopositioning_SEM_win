#include "npRegression.h"

void  npRegression::init(const vpColVector &x0,const vpColVector &y0,const  regressionFunction rgM)
{
    x = x0 ;
    y = y0 ;
    rgModel = rgM;
}

void  npRegression::init()
{
    rgModel = Rational;
}

npRegression::npRegression()
{
    init() ;
}

npRegression::npRegression(const vpColVector &x0,const vpColVector &y0,const  regressionFunction rgM)
{
    x = x0 ;
    y = y0 ;
    rgModel = rgM;
}

void npRegression::setXY(const vpColVector &x0, const vpColVector &y0)
{
    if (x.getRows() != y.getRows())
        std::cerr << "input vector x and y shoud have the same size" << std::endl;
    else
    {
        x = x0 ;
        y = y0 ;
    }
}

void npRegression::addXY(const double x0, const double y0)
{
    unsigned nbr = y.getRows();
    y.resize(nbr+1,false);
    y[nbr] = y0;
    x.resize(nbr+1,false);
    x[nbr] = x0;
}

void npRegression::setRegressionFunction(const regressionFunction rgM)
{
    rgModel = rgM;
}


// compute regression,
void npRegression::linearRegression()
{
    unsigned num = y.getRows();

    vpMatrix A;
    vpMatrix b; // A*p =b;
    switch (rgModel)
    {
    case Rational:
        A.resize(num,5);
        b.resize(num,1);
        for(int i = 0;i<num;i++)
        {
            b[i][0] = y[i]*x[i]*x[i];

            A[i][0] = x[i]*x[i];
            A[i][1] = x[i];
            A[i][2] = 1;
            A[i][3] = -y[i]*x[i];
            A[i][4] = -y[i];
        }
        break;

    case Polynomial:
        A.resize(num,3);
        b.resize(num,1);
        for(int i = 0;i<num;i++)
        {
            b[i][0] = y[i];

            A[i][0] = x[i]*x[i];
            A[i][1] = x[i];
            A[i][2] = 1;
        }
        break;

    case Cauchy:
        A.resize(num,3);
        b.resize(num,1);
        for(int i = 0;i<num;i++)
        {
            b[i][0] = 1;

            A[i][0] = y[i]*x[i]*x[i];
            A[i][1] = y[i]*x[i];
            A[i][2] = y[i];
        }
        break;
    default:
        A.resize(num,5);
        b.resize(num,1);
        for(int i = 0;i<num;i++)
        {
            b[i][0] = y[i]*x[i]*x[i];

            A[i][0] = x[i]*x[i];
            A[i][1] = x[i];
            A[i][2] = 1;
            A[i][3] = -y[i]*x[i];
            A[i][4] = -y[i];
        }
        break;

    }

 //   cout << "A:" << A.t() << endl;
 //   cout << "B:\n" << B << endl;

    p= A.pseudoInverse()*b;
}

// compute regression,
void npRegression::ridgeRegression(const double mu)
{
    unsigned num = y.getRows();

    vpMatrix A;
    vpMatrix b; //  A*p = b;
    switch (rgModel)
    {
    case Rational:
        A.resize(num,5);
        b.resize(num,1);
        for(int i = 0;i<num;i++)
        {
            b[i][0] = y[i]*x[i]*x[i];

            A[i][0] = x[i]*x[i];
            A[i][1] = x[i];
            A[i][2] = 1;
            A[i][3] = -y[i]*x[i];
            A[i][4] = -y[i];
        }
        break;

    case Polynomial:
        A.resize(num,3);
        b.resize(num,1);
        for(int i = 0;i<num;i++)
        {
            b[i][0] = y[i];

            A[i][0] = x[i]*x[i];
            A[i][1] = x[i];
            A[i][2] = 1;
        }
        break;

    case Cauchy:
        A.resize(num,3);
        b.resize(num,1);
        for(int i = 0;i<num;i++)
        {
            b[i][0] = 1;

            A[i][0] = y[i]*x[i]*x[i];
            A[i][1] = y[i]*x[i];
            A[i][2] = y[i];
        }
        break;
    default:
        A.resize(num,5);
        b.resize(num,1);
        for(int i = 0;i<num;i++)
        {
            b[i][0] = y[i]*x[i]*x[i];

            A[i][0] = x[i]*x[i];
            A[i][1] = x[i];
            A[i][2] = 1;
            A[i][3] = -y[i]*x[i];
            A[i][4] = -y[i];
        }
        break;

    }

 //   cout << "A:" << A.t() << endl;
 //   cout << "B:\n" << B << endl;

    vpMatrix AtA = A.AtA() ;
    unsigned int n = A.getCols() ;
    vpMatrix diagAtA(n,n) ;
    diagAtA.eye(n);
    for(unsigned int i = 0 ; i < n ; i++) diagAtA[i][i] = AtA[i][i];

    vpMatrix H = AtA + mu*diagAtA;

    p=H.inverseByLU()*A.t()*b;
}

//compute inverse of derivativeb from regression
double npRegression::jacobianInv(double x0)
{
    double J_reg_inv;
    double J_den,J_num;
    switch (rgModel)
    {
    case Rational:
         J_den = (x0*x0+p[3]*x0+p[4])*(x0*x0+p[3]*x0+p[4]);
         J_num = p[1]*(p[4]-x0*x0)-p[2]*(p[3]+2*x0)+p[0]*x0*(p[3]*x0+2*p[4]);
        J_reg_inv = J_den/J_num;
        break;
    case Polynomial:
        J_reg_inv = 1/(2*p[0]*x0+p[1]);
        break;
    case Cauchy:
         J_den = (p[0]*x0*x0+p[1]*x0+p[2])*(p[0]*x0*x0+p[1]*x0+p[2]);
         J_num = -(2*p[0]*x0+p[1]);
        J_reg_inv = J_den/J_num;
        break;
    default:
         J_den = (x0*x0+p[3]*x0+p[4])*(x0*x0+p[3]*x0+p[4]);
         J_num = p[1]*(p[4]-x0*x0)-p[2]*(p[3]+2*x0)+p[0]*x0*(p[3]*x0+2*p[4]);
        J_reg_inv = J_den/J_num;
        break;
    }
    return J_reg_inv;
}

//compute derivativeb from regression
double npRegression::jacobian(double x0)
{
    double J_reg;
    double J_den,J_num;
    switch (rgModel)
    {
    case Rational:
         J_den = (x0*x0+p[3]*x0+p[4])*(x0*x0+p[3]*x0+p[4]);
         J_num = p[1]*(p[4]-x0*x0)-p[2]*(p[3]+2*x0)+p[0]*x0*(p[3]*x0+2*p[4]);
        J_reg = J_num/J_den;
        break;
    case Polynomial:
        J_reg = (2*p[0]*x0+p[1]);
        break;
    case Cauchy:
         J_den = (p[0]*x0*x0+p[1]*x0+p[2])*(p[0]*x0*x0+p[1]*x0+p[2]);
         J_num = -(2*p[0]*x0+p[1]);
        J_reg = J_num/J_den;
        break;
    default:
         J_den = (x0*x0+p[3]*x0+p[4])*(x0*x0+p[3]*x0+p[4]);
         J_num = p[1]*(p[4]-x0*x0)-p[2]*(p[3]+2*x0)+p[0]*x0*(p[3]*x0+2*p[4]);
        J_reg = J_num/J_den;
        break;
    }
    return J_reg;
}

//compute derivativeb from regression
double npRegression::jacobian(double x0,vpColVector p0,regressionFunction rg)
{
    double J_reg;
    double J_den,J_num;
    switch (rg)
    {
    case Rational:
         J_den = (x0*x0+p0[3]*x0+p0[4])*(x0*x0+p0[3]*x0+p0[4]);
         J_num = p0[1]*(p0[4]-x0*x0)-p0[2]*(p0[3]+2*x0)+p0[0]*x0*(p0[3]*x0+2*p0[4]);
        J_reg = J_num/J_den;
        break;
    case Polynomial:
        J_reg = (2*p0[0]*x0+p0[1]);
        break;
    case Cauchy:
         J_den = (p0[0]*x0*x0+p0[1]*x0+p0[2])*(p0[0]*x0*x0+p0[1]*x0+p0[2]);
         J_num = -(2*p0[0]*x0+p0[1]);
        J_reg = J_num/J_den;
        break;
    default:
         J_den = (x0*x0+p0[3]*x0+p0[4])*(x0*x0+p0[3]*x0+p0[4]);
         J_num = p0[1]*(p0[4]-x0*x0)-p0[2]*(p0[3]+2*x0)+p0[0]*x0*(p0[3]*x0+2*p0[4]);
        J_reg = J_num/J_den;
        break;
    }
    return J_reg;
}

//second derivative (inverse)
double npRegression::hessianInv(double x0)
{
    double H_reg_inv;
    double H_den,H_num;
    switch (rgModel)
    {
    case Rational:
         H_den = -(x0*x0+p[3]*x0+p[4])*(x0*x0+p[3]*x0+p[4])*(x0*x0+p[3]*x0+p[4]);
         H_num = 2*(p[0]*p[3]-p[1])*x0*x0*x0+6*(p[0]*p[4]-p[2])*x0*x0+6*(p[1]*p[4]-p[2]*p[3])+\
                 2*(p[2]*p[4]+p[1]*p[3]*p[4]-p[2]*p[3]*p[3]-p[0]*p[4]*p[4]);
        H_reg_inv = H_den/H_num;
        break;
    case Polynomial:
        H_reg_inv = 1/(2*p[0]);
        break;
    case Cauchy:
         H_den = (p[0]*x0*x0+p[1]*x0+p[2])*(p[0]*x0*x0+p[1]*x0+p[2])*(p[0]*x0*x0+p[1]*x0+p[2]);
         H_num = 6*p[0]*p[0]*x0*x0+6*p[0]*p[1]*x0-2*(p[1]*p[1]-p[0]*p[2]);
        H_reg_inv = H_den/H_num;
        break;
    default:
         H_den = (x0*x0+p[3]*x0+p[4])*(x0*x0+p[3]*x0+p[4]);
         H_num = p[1]*(p[4]-x0*x0)-p[2]*(p[3]+2*x0)+p[0]*x0*(p[3]*x0+2*p[4]);
        H_reg_inv = H_den/H_num;
        break;
    }
    return H_reg_inv;
}

//second derivative
double npRegression::hessian(double x0)
{
    double H_reg;
    double H_den,H_num;
    switch (rgModel)
    {
    case Rational:
         H_den = -(x0*x0+p[3]*x0+p[4])*(x0*x0+p[3]*x0+p[4])*(x0*x0+p[3]*x0+p[4]);
         H_num = 2*(p[0]*p[3]-p[1])*x0*x0*x0+6*(p[0]*p[4]-p[2])*x0*x0+6*(p[1]*p[4]-p[2]*p[3])+\
                 2*(p[2]*p[4]+p[1]*p[3]*p[4]-p[2]*p[3]*p[3]-p[0]*p[4]*p[4]);
        H_reg = H_num/H_den;
        break;
    case Polynomial:
        H_reg = 2*p[0];
        break;
    case Cauchy:
         H_den = (p[0]*x0*x0+p[1]*x0+p[2])*(p[0]*x0*x0+p[1]*x0+p[2])*(p[0]*x0*x0+p[1]*x0+p[2]);
         H_num = 6*p[0]*p[0]*x0*x0+6*p[0]*p[1]*x0-2*(p[1]*p[1]-p[0]*p[2]);
        H_reg = H_num/H_den;
        break;
    default:
         H_den = (x0*x0+p[3]*x0+p[4])*(x0*x0+p[3]*x0+p[4]);
         H_num = p[1]*(p[4]-x0*x0)-p[2]*(p[3]+2*x0)+p[0]*x0*(p[3]*x0+2*p[4]);
        H_reg = H_num/H_den;
        break;
    }
    return H_reg;
}

//second derivative
double npRegression::hessian(double x0, vpColVector p0,regressionFunction rg)
{
    double H_reg;
    double H_den,H_num;
    switch (rg)
    {
    case Rational:
         H_den = -(x0*x0+p0[3]*x0+p0[4])*(x0*x0+p0[3]*x0+p0[4])*(x0*x0+p0[3]*x0+p0[4]);
         H_num = 2*(p0[0]*p0[3]-p0[1])*x0*x0*x0+6*(p0[0]*p0[4]-p0[2])*x0*x0+6*(p0[1]*p0[4]-p0[2]*p0[3])+\
                 2*(p0[2]*p0[4]+p0[1]*p0[3]*p0[4]-p0[2]*p0[3]*p0[3]-p0[0]*p0[4]*p0[4]);
        H_reg = H_num/H_den;
        break;
    case Polynomial:
        H_reg = 2*p0[0];
        break;
    case Cauchy:
         H_den = (p0[0]*x0*x0+p0[1]*x0+p0[2])*(p0[0]*x0*x0+p0[1]*x0+p0[2])*(p0[0]*x0*x0+p0[1]*x0+p0[2]);
         H_num = 6*p0[0]*p0[0]*x0*x0+6*p0[0]*p0[1]*x0-2*(p0[1]*p0[1]-p0[0]*p0[2]);
        H_reg = H_num/H_den;
        break;
    default:
         H_den = (x0*x0+p0[3]*x0+p0[4])*(x0*x0+p0[3]*x0+p0[4]);
         H_num = p0[1]*(p0[4]-x0*x0)-p0[2]*(p0[3]+2*x0)+p0[0]*x0*(p0[3]*x0+2*p0[4]);
        H_reg = H_num/H_den;
        break;
    }
    return H_reg;
}

// get maximum value (ymax) position (x)
double npRegression::getMaxPosition()
{
    double x1,x2,h1,h2,j1m,j1p;
    switch (rgModel)
    {
    case Polynomial:
    case Cauchy:
        return -0.5*(p[1]/p[0]);
    case Rational:
        x1= (p[0]*p[4]-p[2]+sqrt((p[0]*p[4]-p[2])*(p[0]*p[4]-p[2])-\
             (p[0]*p[3]-p[1])*(p[1]*p[4]-p[2]*p[3])))/(p[1]-p[0]*p[3]);

        x2= (p[0]*p[4]-p[2]-sqrt((p[0]*p[4]-p[2])*(p[0]*p[4]-p[2])-\
                (p[0]*p[3]-p[1])*(p[1]*p[4]-p[2]*p[3])))/(p[1]-p[0]*p[3]);
        j1m = jacobian(x1-1.);
        j1p = jacobian(x1+1.);
        std::cout << "x1,x2 = " << x1 << " , " << x2 <<std::endl;
        //std::cout << "j1,j2 = " << j1m << " , " << j1p <<std::endl;
   /*     h1 = hessian(x1);
          h2 = hessian(x2);
        std::cout << "h1,h2 = " << h1 << " , " << h2 <<std::endl;
        */
        if (j1m>0 && j1p<0)
        {
            return x1;
        }
        else if (j1m<0 && j1p>0)
        {
            return x2;
        }
        else
            return 0;

    default:
        return 0;
    }
}


// compute a value of y for a given x from regression function
double npRegression::reponse(double x0)
{
    switch (rgModel)
    {
    case Rational:
        return (p[0]*x0*x0+p[1]*x0+p[2])/(x0*x0+p[3]*x0+p[4]);
    case Polynomial:
        return (p[0]*x0*x0+p[1]*x0+p[2]);
    case Cauchy:
        return 1/(p[0]*x0*x0+p[1]*x0+p[2]);
    default:
        return 0;
    }
}

vpColVector npRegression::getParameters()
{
    return p;
}

// compute regression of G(Z), in case of Cauchy (1/(1+c*(x-b)^2)), G* and G are normalized: dividing by G*
// do not work well...
void npRegression::nonlinearRegression(const vpColVector &G, const vpColVector &Z, double &b, double &c)
{
   // cout << "G= " << G << endl;
    //cout << "Z= " << Z << endl;

    double r = 10; // square of residual error
    unsigned iter = 0;
    unsigned iterMax = 4000;

    //Levenberg-Marquardt
    double mu = 1e-4;
    double lambda = 1e-2;

    vpColVector Ge(G.size());// estimated Gradient
    vpColVector err(G.size()); // residual error

    vpColVector p; // function parameters
    p.resize(2);
    p[0] = b;
    p[1] = c;

    vpColVector dp;

    while (r>1e-3 && iter < iterMax)
    {

        vpMatrix J;
        J.resize(G.size(),2);

        for(int i = 0; i< G.size(); i++)
        {
            double zminusb2 = (Z[i]-p[0])*(Z[i]-p[0]);
            double finv = 1+p[1]*zminusb2;
            double finv2 = finv*finv;

            // estimate Gradient
            Ge[i] = 1/finv;

            // compute Jacobian
            J[i][0] = 2*p[1]*(Z[i]-p[0]) / finv2;
            J[i][1] = -zminusb2 / finv2;
        }

         err = Ge- G;
         r = err.sumSquare();

        // compute the pseudo inverse of the interaction matrix
        vpMatrix Jp ;

        //Levenberg-Marquardt
        vpMatrix JtJ = J.AtA();
       // cout << "JtJ= " << JtJ << endl;

        vpMatrix diagJtJ(JtJ.getRows(), JtJ.getCols()) ; diagJtJ = 0 ;
        for (int i=0 ; i < JtJ.getCols() ; i++)
            diagJtJ[i][i] = JtJ[i][i] ;

       // cout << "diagJtJ= " << diagJtJ << endl;

        Jp = (JtJ + mu*diagJtJ).inverseByLU()*J.t(); //inverseByQR

//       Jp = J.pseudoInverse();
        //cout << "Jp= " << Jp << endl;

        dp = -lambda*Jp*err;

        p += dp;

       // cout << "p= " << p << endl;
       // cout << "r= " << r << endl;

        iter ++;
    }
    b = p[0];
    c = p[1];
}





