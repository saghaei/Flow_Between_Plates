/* This library uses Thomas method
   for solving tri-diagonal matrix.
   Tri-Diogonal Matrix Algorithm (TDMA) */

#include <iostream>
#include <iomanip>
using namespace std;
double* triDiagonal(unsigned, double, double, double);

double* triDiagonal(unsigned n, double mid, double side, double* tmp) {
double *a  = new double[n];
double *b  = new double[n];
double *c  = new double[n];
double *y  = new double[n];
double *gm = new double[n+1];
double *bt = new double[n+1];
static double *x  = new double[n+1];

a[0] = 0;
c[n-1] = 0;
gm[0] = 0;
bt[0] = 0;
bt[n] = 0;
gm[n] = 0;
x[n]  = 0;

for(int i=1 ; i<n ; i++)
a[i] = side;

for(int i=0 ; i<n ; i++)
b[i] = mid;

for(int i=0 ; i<n-1 ; i++)
c[i] = side;

for(int i=0 ; i<n ; i++)
y[i] = tmp[i];

for(int i=0 ; i<=n ; i++) {
gm[i+1] = -c[i] / (a[i]*gm[i]+b[i]);
bt[i+1] = (y[i]-a[i]*bt[i]) / (a[i]*gm[i]+b[i]);
}

for(int i=n ; i>=1 ; i--)
x[i-1] = gm[i]*x[i]+bt[i];

delete[] a, b, c, y, gm, bt;

return x;

delete[] x;
}
