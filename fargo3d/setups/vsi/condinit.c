#include "fargo3d.h"

real DensityDisk(int j, int k) {
  real r = Ymed(j) * sin(Zmed(k));
  real z = Ymed(j) * cos(Zmed(k));
  real p = -SIGMASLOPE;
  real q = 2*FLARINGINDEX - 1;
  real h = ASPECTRATIO*pow(r/R0,FLARINGINDEX);
  real densityDisk = 0.0;

  densityDisk = SIGMA0 / sqrt(2.0*M_PI) / (R0*ASPECTRATIO) * pow(r/R0, p) * \
    exp(1.0/(h*h) * (r / sqrt(r*r + z*z) - 1));
  return densityDisk;
}

real DensityCorona(int j) {
  real r = Ymed(j);
  real densityCorona = DENSITYCONTRAST * SIGMA0 / sqrt(2.0 * M_PI) / (R0 * ASPECTRATIO); 
  densityCorona *= pow(R0 / r, 1.0 / (GAMMA - 1.0));
  return densityCorona;
}

real EnergyCorona(int j) {
  real r = Ymed(j);
  real energyCorona = DENSITYCONTRAST * SIGMA0 / sqrt(2.0 * M_PI) / (R0 * ASPECTRATIO);
  energyCorona *= 1.0 / GAMMA * G * MSTAR / R0 * pow(R0 / r, GAMMA / (GAMMA - 1.0));
  return energyCorona;
}

void Init() {
  int i,j,k;
  real *v1;
  real *v2;
  real *v3;
  real *e;
  real *e0;
  real *rho;
  real h;
  
  rho = Density->field_cpu;
  e   = Energy->field_cpu;
  e0   = Energy_initial->field_cpu;
  v1  = Vx->field_cpu;
  v2  = Vy->field_cpu;
  v3  = Vz->field_cpu;

  srand(42);

  for (k=0; k<Nz+2*NGHZ; k++) {
    for (j=0; j<Ny+2*NGHY; j++) {

      real r = Ymed(j) * sin(Zmed(k));
      real z = Ymed(j) * cos(Zmed(k));
      real r3 = r*r*r;
      real omega = sqrt(G*MSTAR/r3);

      for (i=NGHX; i<Nx+NGHX; i++) {
        real p = -SIGMASLOPE;
        real q = 2*FLARINGINDEX - 1;
        h = ASPECTRATIO*pow(r/R0,FLARINGINDEX);
        real cs = h * r * omega;
        real del_v1 = 0.1 * cs * (rand() - RAND_MAX / 2) / RAND_MAX;
        real del_v2 = 0.1 * cs * (rand() - RAND_MAX / 2) / RAND_MAX;
        real del_v3 = 0.1 * cs * (rand() - RAND_MAX / 2) / RAND_MAX;


        // printf("%e\n", del_v1);
        v2[l] = del_v2;
        v3[l] = del_v3;

        real densDisk = DensityDisk(j, k);
        real densCorona = DensityCorona(j);
        #ifdef ADIABATIC
        real eDisk = densDisk*h*h*G*MSTAR/r/(GAMMA-1.0);
        #endif
        #ifdef ISOTHERMAL
        real eDisk = h*sqrt(G*MSTAR/r);
        #endif

        real eCorona = EnergyCorona(j);

        if (eDisk > eCorona) {
            rho[l] = densDisk;
            e[l] = eDisk;
            v1[l] = omega*r;
            v1[l] *= sqrt((p+q)*h*h + (1+q) - q*r / sqrt(r*r + z*z));
        } else {
            rho[l] = densCorona;
            e[l] = eCorona;
            v1[l] = 0.0;
        }

        e0[l] = e[l] / rho[l] * (GAMMA - 1.0);

        v1[l] -= OMEGAFRAME*r;
        v1[l] += del_v1;

      }
    }
  }
}

void CondInit() {
  
  const int id_gas = 0;
  const int feedback = NO;
  Fluids[id_gas] = CreateFluid("gas",GAS);
  SelectFluid(id_gas);
  Init();
}
