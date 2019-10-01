#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(void)
{
//  const int Lmax = 500;
//  const int Lmax = 1000;
  const int Lmax = 1500;
//  const int Lmax = 2000;
//  const int Lmax = 2500;
//  const int Lmax = 10000;
  const double pi = M_PI;

  double invLmax = 1.0/Lmax;
  int i;
  int ix,iy;
  double kx,ky;
  double Qx,Qy;
  double parQ;
  double parS, parJ, parJp;
  double mag;
  double ene;
  double tmp;

//  Qx = pi;
//  Qy = pi;
//  parS = 0.5;
//  parS = 1.0;
  parS = 1.5;
//  parS = 100.0;
  parJ = 1.0;
  parJp = 0.0;
//  parJp = 0.49;

  double dispJQ;
  double **dispJ;
  dispJ = (double**)malloc(Lmax*sizeof(double*));
  dispJ[0] = (double*)malloc(Lmax*Lmax*sizeof(double));
  for(i=0; i<Lmax; i++){
    dispJ[i] = dispJ[0] + i*Lmax;
  }
  double **dispJp;
  dispJp = (double**)malloc(Lmax*sizeof(double*));
  dispJp[0] = (double*)malloc(Lmax*Lmax*sizeof(double));
  for(i=0; i<Lmax; i++){
    dispJp[i] = dispJp[0] + i*Lmax;
  }
  double **dispJm;
  dispJm = (double**)malloc(Lmax*sizeof(double*));
  dispJm[0] = (double*)malloc(Lmax*Lmax*sizeof(double));
  for(i=0; i<Lmax; i++){
    dispJm[i] = dispJm[0] + i*Lmax;
  }
  double **dispA;
  dispA = (double**)malloc(Lmax*sizeof(double*));
  dispA[0] = (double*)malloc(Lmax*Lmax*sizeof(double));
  for(i=0; i<Lmax; i++){
    dispA[i] = dispA[0] + i*Lmax;
  }
  double **dispW;
  dispW = (double**)malloc(Lmax*sizeof(double*));
  dispW[0] = (double*)malloc(Lmax*Lmax*sizeof(double));
  for(i=0; i<Lmax; i++){
    dispW[i] = dispW[0] + i*Lmax;
  }

  printf("# L J Jp S mag E_cla E_qua Q/pi\n");

//for(parJp=0.0; parJp<=0.5+1.0e-10; parJp+=0.025){
for(parS=0.5; parS<=2.5+1.0e-10; parS+=0.5){

  parQ = pi;
//  parQ = acos(-parJ/parJp*0.5);
  Qx = parQ;
  Qy = parQ;

  for(ix=0; ix<Lmax; ix++){
    for(iy=0; iy<Lmax; iy++){
      // J_k
      kx = 2.0*pi*ix*invLmax;
      ky = 2.0*pi*iy*invLmax;
      dispJ[ix][iy] =
        + parJ*( cos(kx)+cos(ky) )
        + parJp*( cos(kx+ky) );
      // J_Q+k
      kx = Qx + 2.0*pi*ix*invLmax;
      ky = Qy + 2.0*pi*iy*invLmax;
      dispJp[ix][iy] =
        + parJ*( cos(kx)+cos(ky) )
        + parJp*( cos(kx+ky) );
      // J_Q-k
      kx = Qx - 2.0*pi*ix*invLmax;
      ky = Qy - 2.0*pi*iy*invLmax;
      dispJm[ix][iy] =
        + parJ*( cos(kx)+cos(ky) )
        + parJp*( cos(kx+ky) );
    }
  }
  // J_Q
  kx = Qx;
  ky = Qy;
  dispJQ =
    + parJ*( cos(kx)+cos(ky) )
    + parJp*( cos(kx+ky) );

  for(ix=0; ix<Lmax; ix++){
    for(iy=0; iy<Lmax; iy++){
      dispA[ix][iy] =
        + 0.25*(dispJp[ix][iy] + dispJm[ix][iy])
        + 0.5*dispJ[ix][iy]
        - dispJQ;
//      dispW[ix][iy] =
//        2.0*parS*sqrt( (dispJ[ix][iy] - dispJQ)*(0.5*(dispJp[ix][iy] + dispJm[ix][iy]) - dispJQ) );
      tmp = (dispJ[ix][iy] - dispJQ)*(0.5*(dispJp[ix][iy] + dispJm[ix][iy]) - dispJQ);
      if(tmp < 0){
        printf("# error: for parJp=%f, dispW(%d,%d)^2=%f < 0 !\n",parJp,ix,iy,tmp);
        if(fabs(tmp) < 1.0e-10) tmp = 1.0e-32;
      }
      dispW[ix][iy] = 2.0*parS*sqrt(tmp);
    }
  }

  mag = parS + 0.5;
  for(ix=0; ix<Lmax; ix++){
    for(iy=0; iy<Lmax; iy++){
      if(dispW[ix][iy] > 1.0e-10){
        mag += -parS/Lmax/Lmax * dispA[ix][iy]/dispW[ix][iy];
//        mag += -0.5/Lmax/Lmax * dispA[ix][iy]/dispW[ix][iy];
      }
    }
  }

  ene = dispJQ*parS*parS;
  for(ix=0; ix<Lmax; ix++){
    for(iy=0; iy<Lmax; iy++){
      ene += 1.0/Lmax/Lmax * (
//        + parS*(dispW[ix][iy] - dispA[ix][iy])
        + 0.5*dispW[ix][iy] - parS*dispA[ix][iy]
//        + 0.5*(dispW[ix][iy] - dispA[ix][iy])
//        + dispW[ix][iy] // alpha_k^dag alpha_k = 0 for GS
        );
    }
  }

/*
  for(ix=0; ix<Lmax; ix++){
    for(iy=0; iy<Lmax; iy++){
      kx = 2.0*pi*ix*invLmax;
      ky = 2.0*pi*iy*invLmax;
      printf("%d %d %f %f %f %f %f\n",ix,iy,kx,ky,dispJ[ix][iy],dispA[ix][iy],dispW[ix][iy]);
    }
    printf("\n");
  }
  printf("\n");
  printf("# %d %f\n",Lmax,mag);
*/

//  printf("%d %f %f %f %f\n",Lmax,parJp,mag,dispJQ*parS*parS,ene);
  printf("%d %f %f %f %f %f %f %f\n",Lmax,parJ,parJp,parS,mag,dispJQ*parS*parS,ene,parQ/pi);

}

  free(dispJ[0]);
  free(dispJ);
  free(dispJp[0]);
  free(dispJp);
  free(dispJm[0]);
  free(dispJm);
  free(dispA[0]);
  free(dispA);
  free(dispW[0]);
  free(dispW);

  return 0;
}
