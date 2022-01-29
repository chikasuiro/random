#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <fstream>
using namespace std;

constexpr double zeta = 0.05;
constexpr double omega0 = 1.0;
constexpr double largeomega = 6.0;
constexpr double h = 20.0;

constexpr double dt = 1e-3;

double f_largex(double largev) {
  return largev;
}

double f_largev(double largex, double largev, double theta, double g) {
  return -zeta*largev - pow(omega0,2)*largex - g*pow(largex,3) + h*cos(theta);
}

double f_theta() {
  return largeomega;
}

int main(int argc, char *argv[]) {
  if (argc < 3) return -1;
  const double g = atof(argv[1]);
  ofstream outfile(argv[2]);
  outfile << "# t X V theta" << endl;

  double largex = 1;
  double largev = 1;
  double theta  = 1;

  for (int i = 0; i < 1e5; ++i) {
    outfile << scientific
            << i << " " << largex << " " << largev << " " << theta
            << endl;
    
    double k1[3], k2[3], k3[3], k4[3];

    k1[0] = f_largex(largev);
    k1[1] = f_largev(largex, largev, theta, g);
    k1[2] = f_theta();

    k2[0] = f_largex(largev+0.5*dt*k1[1]);
    k2[1] = f_largev(largex+0.5*dt*k1[0], largev+0.5*dt*k1[1], theta+0.5*dt*k1[2], g);
    k2[2] = f_theta();

    k3[0] = f_largex(largev+0.5*dt*k2[1]);
    k3[1] = f_largev(largex+0.5*dt*k2[0], largev+0.5*dt*k2[1], theta+0.5*dt*k2[2], g);
    k3[2] = f_theta();

    k4[0] = f_largex(largev+dt*k3[1]);
    k4[1] = f_largev(largex+dt*k3[0], largev+dt*k3[1], theta+dt*k3[2], g);
    k4[2] = f_theta();

    largex += dt/6.0 * (k1[0]+k4[0]+2.0*(k2[0]+k3[0]));
    largev += dt/6.0 * (k1[1]+k4[1]+2.0*(k2[1]+k3[1]));
    theta  += dt/6.0 * (k1[2]+k4[2]+2.0*(k2[2]+k3[2]));
  }

  outfile.close();

  return 0;
}
