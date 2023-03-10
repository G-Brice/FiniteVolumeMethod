#include <iostream>
#include <cmath>
#include <Eigen/Core>
#include <fstream>

using namespace std;
using namespace Eigen;

/******************************************************************************/
/* Parameters of the problem                                                  */
/******************************************************************************/

const int N = 200;   // Number of points of the grid
double h    = 1./N;   // Space step
double a    = 1.;     // Velocity in the transport equation
double T    = 2.;     // Final time
double dt   = h/2./a; // Time step

/******************************************************************************/
/* Definition of the function g to compute the numerical flux                 */
/******************************************************************************/

// Function to define the numerical flux F_{i+1/2} as g(U_i, U_{i+1}):
double g(double u, double v) {
  return (a>0)*a*u + (a<0)*a*v; // Upwind
  //return a * (u+v) / 2 + (u-v); // LF
  //return a * (u+v) / 2 - 0.5* dt/h*a*a*(v-u); // LW
}

/******************************************************************************/
/* Definition of the initial condition u0                                     */
/******************************************************************************/

const double PI = 4*atan(1); // Constant pi

/* Initial condition */
double u0(double x) {
  //return sin(2*PI*x);
  return (x>0.5);
}

int main() {

  ofstream myfile;
  myfile.open ("data.txt");

  VectorXd U0(N+1);   // U^{n}
  VectorXd U1(N+1);   // U^{n+1}
  VectorXd Xmid(N+1); // Midpoints of the mesh
  VectorXd Xint(N+1); // Interface points of the mesh

/******************************************************************************/
/* Initializations                                                            */
/******************************************************************************/

  Xmid(0) = 0;
  for (int i=1; i<=N; i++) {
    Xmid(i) = Xmid(i-1)+h;
  }

  Xint(0) = h/2;
  for (int i=1; i<=N; i++) {
    Xint(i) = Xint(i-1)+h;
  }

  for (int i=0; i<=N; i++) {
    U0(i) = u0(Xmid(i));
  }

  myfile << N << endl;

  myfile << Xmid << endl;

  double t = 0; // Current time

  myfile << t << endl << U0 << endl;

/******************************************************************************/
/* Time step, update of the solution                                          */
/******************************************************************************/

  while (t <= T) {

    t += dt;
    for (int i=1; i<N; i++) {
      U1(i) = U0(i) - dt/h * (g(U0(i), U0(i+1)) - g(U0(i-1), U0(i)));
    }

    // Periodic boundary conditions:
    U1(0) = U0(0) - dt/h* (g(U0(0), U0(1)) - g(U0(N-1), U0(0)));
    U1(N) = U1(0);

    // Homogeneous Neumann boundary conditions:
    //U1(0) = U1(1);
    //U1(N) = U1(N-1);

    U0 = U1;

    myfile << t << endl << U0 << endl;
  }

  myfile.close();

  return 0;
}
