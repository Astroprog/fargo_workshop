#include "fargo3d.h"

void Init() {
  
  int i,j,k;
  
  real *sigma  = Density->field_cpu;
  real *vphi   = Vx->field_cpu;
  real *vr     = Vy->field_cpu;
  real *cs     = Energy->field_cpu;

  i = j = k = 0;
  
  for (j=0; j<Ny+2*NGHY; j++) {
    for (i=0; i<Nx+2*NGHX; i++) {
      
      real vkmed = sqrt(G*MSTAR/ymed(j));
      
      sigma[l] = SIGMA0*pow(ymed(j)/R0,-SIGMASLOPE);
      cs[l]    = vkmed*ASPECTRATIO;
      vphi[l]  = vkmed*sqrt(1.0 - ASPECTRATIO*ASPECTRATIO*(1.0 + SIGMASLOPE));

      real nu  = ALPHA*cs[l]*ASPECTRATIO*ymin(j);
      vr[l]    = -1.5*nu/ymin(j);
    }
  } 
}

void CondInit() {
   Fluids[0] = CreateFluid("gas",GAS);
   SelectFluid(0);
   Init();
}
