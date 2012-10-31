/**
 * pkolmogorov2x and psmirnov2x
 * 
 * Taken from R version 2.15.1 (under GPL)
 **/

#include <R.h>
#include <Rmath.h>		/* constants */

#include "kolmogorov_smirnov.h"

static double K(int n, double d);
static void m_multiply(double *A, double *B, double *C, int m);
static void m_power(double *A, int eA, double *V, int *eV, int m, int n);



/* Two-sided two-sample */
void
psmirnov2x(double *x, Sint *m, Sint *n)
{
    double md, nd, q, *u, w;
    Sint i, j;

    if(*m > *n) {
	i = *n; *n = *m; *m = i;
    }
    md = (double) (*m);
    nd = (double) (*n);
    /*
       q has 0.5/mn added to ensure that rounding error doesn't
       turn an equality into an inequality, eg abs(1/2-4/5)>3/10 

    */
    q = (0.5 + floor(*x * md * nd - 1e-7)) / (md * nd);
    u = (double *) R_alloc(*n + 1, sizeof(double));

    for(j = 0; j <= *n; j++) {
	u[j] = ((j / nd) > q) ? 0 : 1;
    }
    for(i = 1; i <= *m; i++) {
	w = (double)(i) / ((double)(i + *n));
	if((i / md) > q)
	    u[0] = 0;
	else
	    u[0] = w * u[0];
	for(j = 1; j <= *n; j++) {
	    if(fabs(i / md - j / nd) > q) 
		u[j] = 0;
	    else
		u[j] = w * u[j] + u[j - 1];
	}
    }
    *x = u[*n];
}

/* The two-sided one-sample 'exact' distribution */
void
pkolmogorov2x(double *x, Sint *n)
{
    /* x is input and output. */

    *x = K(*n, *x);
}

static double
K(int n, double d)
{
    /* Compute Kolmogorov's distribution.
       Code published in
	 George Marsaglia and Wai Wan Tsang and Jingbo Wang (2003),
	 "Evaluating Kolmogorov's distribution".
	 Journal of Statistical Software, Volume 8, 2003, Issue 18.
	 URL: http://www.jstatsoft.org/v08/i18/.
    */

   int k, m, i, j, g, eH, eQ;
   double h, s, *H, *Q;

   /* 
      The faster right-tail approximation is omitted here.
      s = d*d*n; 
      if(s > 7.24 || (s > 3.76 && n > 99)) 
          return 1-2*exp(-(2.000071+.331/sqrt(n)+1.409/n)*s);
   */
   k = (int) (n * d) + 1;
   m = 2 * k - 1;
   h = k - n * d;
   H = (double*) Calloc(m * m, double);
   Q = (double*) Calloc(m * m, double);
   for(i = 0; i < m; i++)
       for(j = 0; j < m; j++)
	   if(i - j + 1 < 0)
	       H[i * m + j] = 0;
	   else
	       H[i * m + j] = 1;
   for(i = 0; i < m; i++) {
       H[i * m] -= pow(h, i + 1);
       H[(m - 1) * m + i] -= pow(h, (m - i));
   }
   H[(m - 1) * m] += ((2 * h - 1 > 0) ? pow(2 * h - 1, m) : 0);
   for(i = 0; i < m; i++)
       for(j=0; j < m; j++)
	   if(i - j + 1 > 0)
	       for(g = 1; g <= i - j + 1; g++)
		   H[i * m + j] /= g;
   eH = 0;
   m_power(H, eH, Q, &eQ, m, n);
   s = Q[(k - 1) * m + k - 1];
   for(i = 1; i <= n; i++) {
       s = s * i / n;
       if(s < 1e-140) {
	   s *= 1e140;
	   eQ -= 140;
       }
   }
   s *= pow(10., eQ);
   Free(H);
   Free(Q);
   return(s);
}

static void
m_multiply(double *A, double *B, double *C, int m)
{
    /* Auxiliary routine used by K().
       Matrix multiplication.
    */
    int i, j, k;
    double s;
    for(i = 0; i < m; i++)
	for(j = 0; j < m; j++) {
	    s = 0.;
	    for(k = 0; k < m; k++)
		s+= A[i * m + k] * B[k * m + j];
	    C[i * m + j] = s;
	}
}

static void
m_power(double *A, int eA, double *V, int *eV, int m, int n)
{
    /* Auxiliary routine used by K().
       Matrix power.
    */
    double *B;
    int eB , i;

    if(n == 1) {
	for(i = 0; i < m * m; i++)
	    V[i] = A[i];
	*eV = eA;
	return;
    }
    m_power(A, eA, V, eV, m, n / 2);
    B = (double*) Calloc(m * m, double);
    m_multiply(V, V, B, m);
    eB = 2 * (*eV);
    if((n % 2) == 0) {
	for(i = 0; i < m * m; i++)
	    V[i] = B[i];
	*eV = eB;
    }
    else {
	m_multiply(A, B, V, m);
	*eV = eA + eB;
    }
    if(V[(m / 2) * m + (m / 2)] > 1e140) {
	for(i = 0; i < m * m; i++)
	    V[i] = V[i] * 1e-140;
	*eV += 140;
    }
    Free(B);
}
