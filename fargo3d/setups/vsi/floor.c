//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

void Floor_cpu() {

//<USER_DEFINED>
  INPUT(Density);
  INPUT(Energy);
  INPUT(Vx);
  INPUT(Vy);
  INPUT(Vz);
  OUTPUT(Density);
  OUTPUT(Energy);
  OUTPUT(Vx);
  OUTPUT(Vy);
  OUTPUT(Vz);
//<\USER_DEFINED>


//<EXTERNAL>
  real* dens = Density->field_cpu;
  real* e    = Energy->field_cpu;
  real* vx   = Vx->field_cpu;
  real* vy   = Vy->field_cpu;
  real* vz   = Vz->field_cpu;
  int pitch  = Pitch_cpu;
  int stride = Stride_cpu;
  int size_x = Nx+2*NGHX;
  int size_y = Ny+2*NGHY;
  int size_z = Nz+2*NGHZ;
  int fluidtype = Fluidtype;
  real mstar_cgs =  MSTAR_CGS;
  real g_cgs = G_CGS;
  real r0 = R0;
  real r0_cgs = R0_CGS;
  real density_corona0 = 1e-2 * DENSITYCONTRAST * SIGMA0 / sqrt(2.0 * M_PI) / (R0 * ASPECTRATIO);
  real temp_conv = (GAMMA - 1) * MSTAR_CGS * G_CGS / R0_CGS / 1.380649e-16 * 2.4 * 1.6726231e-24;
  real eps = EPSILON;
//<\EXTERNAL>

//<INTERNAL>
  int i;
  int j;
  int k;
  int ll;
  real density_corona;
  real temperature;
//<\INTERNAL>

//<CONSTANT>
// real xmin(Nx+1);
// real ymin(Ny+2*NGHY+1);
// real zmin(Nz+2*NGHZ+1);
// real DENSITYCONTRAST(1);
// real TEMPERATUREFLOOR(1);
// real SIGMA0(1);
// real ASPECTRATIO(1);
// real GAMMA(1);
//<\CONSTANT>
  
//<MAIN_LOOP>

  i = j = k = 0;

#ifdef Z
  for (k=0; k<size_z; k++) {
#endif
#ifdef Y
    for (j=0; j<size_y; j++) {
#endif
#ifdef X
      for (i=0; i<size_x; i++ ) {
#endif
//<#>
	ll = l;

  if (fluidtype == GAS) {
    density_corona = density_corona0 * pow(r0 / ymed(j), 1.0 / (GAMMA - 1.0));
    if (dens[ll] < density_corona) {

      vx[ll] *= dens[ll] / density_corona;
      vy[ll] *= dens[ll] / density_corona;
      vz[ll] *= dens[ll] / density_corona;
      dens[ll] = density_corona;
    }

    temperature = temp_conv * e[ll] / dens[ll];
    if (temperature < TEMPERATUREFLOOR) {
      e[ll] = TEMPERATUREFLOOR * dens[ll] / temp_conv;
    }
  } else {

    // density_corona = eps * density_corona0 * pow(r0 / ymed(j), 1.0 / (GAMMA - 1.0));
    // if (dens[ll] < density_corona) {

    //   vx[ll] *= dens[ll] / density_corona;
    //   vy[ll] *= dens[ll] / density_corona;
    //   vz[ll] *= dens[ll] / density_corona;
    //   dens[ll] = density_corona;
    // }

  }

//<\#>
#ifdef X
      }
#endif
#ifdef Y
    }
#endif
#ifdef Z
  }
#endif
//<\MAIN_LOOP>
}
