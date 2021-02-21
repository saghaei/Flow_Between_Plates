/* This piece of a code simulates the flow between two
   parallel plates with infinite size while one plate
   is moving with constant velocity.                   */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include "triDiagonal.hpp"
using namespace std;

void FTCS(double**, const double, const size_t);
void DuFortFrankle(double**, const double, const size_t);
void Laasonen(double**, const double, const size_t);
void CrankNicolson(double**, const double, const size_t);

const size_t im=41;       // Number of vertical Nodes
const double nu=0.000217; // Viscousity [m^2/s]
const double u0=40;       // Lower surface speed [m/s]
const double dy=0.001;    // Spacing between nodes

int main() {
const size_t nm=541;      // Number of iterations (time elapse)
const double dt=0.002;    // Time step
unsigned short int mtd;   // Solving method

// Allocating memory for array
double **u = new double*[im];
for(size_t r=0 ; r<im ; r++)
u[r] = new double[nm];

// Applying I.C.
for(size_t r=1 ; r<im-1 ; r++)
u[r][0]=0;

// Applying B.C.
for(size_t c=0 ; c<nm ; c++)
u[0][c]=u0;
for(size_t c=0 ; c<nm ; c++)
u[im-1][c]=0;

cout << "Fluid bounded by two parallel plates." << endl;
cout << "\n1) Forward Time Central Space\n2) DuFort-Frankle\n";
cout << "3) Laasonen\n4) Crank-Nicolson\n";
do {
cout << "\nEnter solving method: ";
cin >> mtd;
} while(mtd<1 || mtd>4);

switch (mtd) {
case 1:
FTCS(u,dt,nm);
break;
case 2:
DuFortFrankle(u,dt,nm);
break;
case 3:
Laasonen(u,dt,nm);
break;
case 4:
CrankNicolson(u,dt,nm);
break;
default:
cout << "Program execution ended." << endl;
exit(1);
break;
}

// Writing to file
fstream output("data", ios::out);
if(!output) {
cerr << "File could not be created." << endl;
exit(1);
}
for(size_t r=0 ; r<im ; r++) {
for(size_t c=0 ; c<nm ; c++)
output << setw(4) << setprecision(1) << fixed << left << u[r][c] << '\t';
output << endl;
}

// Releasing allocated memory
for(size_t c=0 ; c<im ; c++)
delete[] u[c];
delete[] u;

return 0;
}

void FTCS(double** u, const double dt, const size_t nm) {
for(size_t n=0 ; n<nm-1 ; n++)
for(size_t i=1 ; i<im-1 ; i++)
u[i][n+1] = u[i][n] + ((nu*dt)/(dy*dy)) * (u[i+1][n] - 2*u[i][n] + u[i-1][n]);
}

void DuFortFrankle(double** u, const double dt, const size_t nm) {
for(size_t n=0 ; n<nm-1 ; n++)
for(size_t i=1 ; i<im-1 ; i++)
u[i][n+1] = ((dy*dy-2*nu*dt)/(dy*dy+2*nu*dt)) * u[i][n-1] + ((2*nu*dt)/(dy*dy+2*nu*dt)) * (u[i+1][n] + u[i-1][n]);
}

void Laasonen(double** u, const double dt, const size_t nm) {
const double k = (nu*dt)/(dy*dy);
static double* tmp = new double[im-2];
double *temp;
for(size_t i=0 ; i<im-2 ; i++)
tmp[i] = u[i+1][0];
tmp[0] += k*u[0][1];
tmp[im-3] += k*u[im-1][1];
for(size_t n=1 ; n<nm ; n++) {
temp = triDiagonal(im-2,(2*k+1),(-k),tmp);
for(size_t i=0 ; i<im-2 ; i++)
u[i+1][n] = temp[i];
for(size_t i=0 ; i<im-2 ; i++)
tmp[i] = u[i+1][n];
tmp[0] += k*u[0][n];
tmp[im-3] += k*u[im-1][n];
}
delete[] temp, tmp;
}

void CrankNicolson(double** u, const double dt, const size_t nm) {
const double k = (nu*dt)/(2*dy*dy);
double* u_half = new double[im-2];
double* temp;
for(size_t i=0 ; i<im-2 ; i++)
u_half[i] = (1-2*k)*u[i+1][0] + k*(u[i][0]+u[i+2][0]);
u_half[0] += k*u[0][1];
u_half[im-3] += k*u[im-2][1];
for(size_t n=1 ; n<nm ; n++) {
temp = triDiagonal(im-2,(2*k+1),(-k),u_half);
for(size_t i=0 ; i<im-2 ; i++)
u[i+1][n] = temp[i];
for(size_t r=0 ; r<im-2 ; r++)
u_half[r] = (1-2*k)*u[r+1][n] + k*(u[r][n]+u[r+2][n]);
u_half[0] += k*u[0][n];
u_half[im-3] += k*u[im-2][n];
}
delete[] temp, u_half;
}
