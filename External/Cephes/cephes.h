#ifndef CCEPHES
#define CCEPHES
#ifdef __cplusplus
extern "C" {
#endif
/*
Wrapper Header to access distribution functions from the Cephes Mathematical Library by Stephen L. Moshier
We use the Double precision arithmetic function library. http://www.moshier.net/double.zip
*/
namespace Cephes {
    /*airy.c: Airy Function Wrapper*/
    extern int airy(double,double*,double*,double*,double*);

    /*fac.c: Factorial Function Wrapper*/
    extern double fac(int);

    /*isnan.c wrapper*/
    extern int isnan(double);
    extern int signbit(double);
    extern int isfinite(double);

    /*ndtr.c: Normal Distribution Function Wrapper*/
    extern double ndtr(double);

    /*ndtri.c: Inverse of Normal Distribution Function Wrapper*/
    extern double ndtri(double);

    /*gamma.c wrapper*/
    extern double gamma(double);
    extern double lgam(double);

    /*gdtr.c: Gamma distribution Function Wrapper*/
    extern double gdtr(double,double,double);
    extern double gdtrc(double,double,double);

    /*igam.c: Incomplete Gamma Integral and Complemented Incomplete Gamma Integral Wrapper*/
    extern double igamc(double,double);
    extern double igam(double,double);

    /*igami.c: Inverse of Complemented Incomplete Gamma Integral*/
    extern double igami(double,double);

    /*rgamma.c: Reciprocal Gamma Function Wrapper*/
    extern double rgamma(double);

    /*incbet.c: Incomplete Beta Integral Wrapper*/
    extern double incbet(double,double,double);

    /*beta.c: Beta function Wrapper*/
    extern double beta(double,double);
    extern double lbeta(double,double);

    /*btdtr.c: Beta Distribution Wrapper*/
    extern double btdtr(int,int,double);

    /*bdtr.c: Binomial Distribution Wrapper*/
    extern double bdtrc(int,int,double);
    extern double bdtr(int,int,double);
    extern double bdtri(int,int,double);

    /*incbi.c: Inverse of Incomplete Beta Integral Wrapper*/
    extern double incbi(double,double,double);

    /*chbevl.c: Evaluate Chebyshev series*/
    extern double chbevl(double, double[], int);

    /*chdtr.c: Chi-squared distribution*/
    extern double chdtrc(double,double);
    extern double chdtr(double,double);
    extern double chdtri(double,double);

    /*dawsn.c: Dawson's Integral*/
    extern double dawsn(double);

    /*drand.c: Pseudorandom number generator*/
    extern int drand(double*);

    /*hyperg.c: Confluent hypergeometric function*/
    extern double hyperg(double,double,double);
    extern double hyp2f0(double,double,double,int,double*);

    /*hyp2f1.c: Gauss Hypergeometric Function F21*/
    extern double hyp2f1(double,double,double,double);

    /*psi.c: Psi (digamma) Function Wrapper*/
    extern double psi(double);

    /*stdtr.c: Student's T Distribution Wrapper*/
    extern double stdtr(int,double);
    extern double stdtri(int,double);

    /*pdtr.c: Poisson Distribution Wrapper*/
    extern double pdtrc(int,double);
    extern double pdtr(int,double);
    extern double pdtri(int,double);

    /*fdtr.c: F Distribution Wrapper*/
    extern double fdtrc(int,int,double);
    extern double fdtr(int,int,double);
    extern double fdtri(int,int,double);

    /*fresnl.c: Fresnel Integral*/
    extern int fresnl(double,double*,double*);

    /*nbdtr.c: Negative binomial distribution Wrapper*/
    extern double nbdtrc(int,int,double);
    extern double nbdtr(int,int,double);
    extern double nbdtri(int,int,double);

}; //Namespace

#ifdef __cplusplus
}
#endif
#endif
