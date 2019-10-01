#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main
(
  int argc,
  char *argv[]
)
{
//  const int Lmax = 500;
//  const int Lmax = 1000;
  const int Lmax = 1500;
//  const int Lmax = 2000;
//  const int Lmax = 2500;
//  const int Lmax = 10000;
  const int Lz = 2;
  const double pi = M_PI;

  double invLmax = 1.0/Lmax;
  double invLz = 1.0/Lz;
  int i,j;
  int ix,iy,iz;
  double kx,ky,kz;
  double Qx,Qy,Qz;
  double parQ;
  double parS, parJ, parJp, parJd;
  double mag;
  double ene;
  double tmp;

//  Qx = pi;
//  Qy = pi;
//  Qz = pi;

  parS = 0.5;
//  parS = 1.0;
//  parS = 1.5;
//  parS = 100.0;

  parJ = 1.0;

//  parJp = 0.0;
//  parJp = 0.49;
  parJp = 0.50;
//  parJp = 0.51;
//  parJp = 1.0;

  parJd = 0.0;
//  parJd = 0.5;
//  parJd = 1.0;
//  parJd = 2.0;
//  parJd = 14.0;
//  parJd = 15.0;
//  parJd = 20.0;
//  parJd = 100.0;

  double dispJQ;
  double ***dispJ;
  dispJ = (double***)malloc(Lmax*sizeof(double**));
  dispJ[0] = (double**)malloc(Lmax*Lmax*sizeof(double*));
  dispJ[0][0] = (double*)malloc(Lmax*Lmax*Lz*sizeof(double));
  for(i=0; i<Lmax; i++){
    dispJ[i] = dispJ[0] + i*Lmax;
    for(j=0; j<Lmax; j++){
      dispJ[i][j] = dispJ[0][0] + (i*Lmax+j)*Lz;
    }
  }
  double ***dispJp;
  dispJp = (double***)malloc(Lmax*sizeof(double**));
  dispJp[0] = (double**)malloc(Lmax*Lmax*sizeof(double*));
  dispJp[0][0] = (double*)malloc(Lmax*Lmax*Lz*sizeof(double));
  for(i=0; i<Lmax; i++){
    dispJp[i] = dispJp[0] + i*Lmax;
    for(j=0; j<Lmax; j++){
      dispJp[i][j] = dispJp[0][0] + (i*Lmax+j)*Lz;
    }
  }
  double ***dispJm;
  dispJm = (double***)malloc(Lmax*sizeof(double**));
  dispJm[0] = (double**)malloc(Lmax*Lmax*sizeof(double*));
  dispJm[0][0] = (double*)malloc(Lmax*Lmax*Lz*sizeof(double));
  for(i=0; i<Lmax; i++){
    dispJm[i] = dispJm[0] + i*Lmax;
    for(j=0; j<Lmax; j++){
      dispJm[i][j] = dispJm[0][0] + (i*Lmax+j)*Lz;
    }
  }
  double ***dispA;
  dispA = (double***)malloc(Lmax*sizeof(double**));
  dispA[0] = (double**)malloc(Lmax*Lmax*sizeof(double*));
  dispA[0][0] = (double*)malloc(Lmax*Lmax*Lz*sizeof(double));
  for(i=0; i<Lmax; i++){
    dispA[i] = dispA[0] + i*Lmax;
    for(j=0; j<Lmax; j++){
      dispA[i][j] = dispA[0][0] + (i*Lmax+j)*Lz;
    }
  }
  double ***dispW;
  dispW = (double***)malloc(Lmax*sizeof(double**));
  dispW[0] = (double**)malloc(Lmax*Lmax*sizeof(double*));
  dispW[0][0] = (double*)malloc(Lmax*Lmax*Lz*sizeof(double));
  for(i=0; i<Lmax; i++){
    dispW[i] = dispW[0] + i*Lmax;
    for(j=0; j<Lmax; j++){
      dispW[i][j] = dispW[0][0] + (i*Lmax+j)*Lz;
    }
  }

  char opt;
  /* prameter settings */
  for(i = 0; i < argc; ++i){
    if(*argv[i] == '-'){
      opt = *(argv[i]+1);
      if(opt == 's'){
        parS = atof(argv[i+1]);
      }else if(opt == 'p'){
        parJp = atof(argv[i+1]);
      }else if(opt == 'd'){
        parJd = atof(argv[i+1]);
      }
    }
  }

  printf("# L J Jp Jd S mag E_cla E_qua Q/pi\n");

//for(parJp=0.0; parJp<=0.5+1.0e-10; parJp+=0.025){
//for(parS=0.5; parS<=2.5+1.0e-10; parS+=0.5){
//for(parJd=0.0; parJd<=20.0+1.0e-10; parJd+=0.1){
//for(parJd=0.0; parJd<=80.0+1.0e-10; parJd+=5.0){

  parQ = pi;
//  parQ = acos(-parJ/parJp*0.5);
  Qx = parQ;
  Qy = parQ;
  Qz = pi;

  for(ix=0; ix<Lmax; ix++){
    for(iy=0; iy<Lmax; iy++){
      for(iz=0; iz<Lz; iz++){
        // J_k
        kx = 2.0*pi*ix*invLmax;
        ky = 2.0*pi*iy*invLmax;
        kz = 2.0*pi*iz*invLz;
        dispJ[ix][iy][iz] =
          + 2.0*parJ*( cos(kx)+cos(ky) )
          + 2.0*parJp*( cos(kx+ky) )
          + parJd*cos(kz);
        // J_Q+k
        kx = Qx + 2.0*pi*ix*invLmax;
        ky = Qy + 2.0*pi*iy*invLmax;
        kz = Qz + 2.0*pi*iz*invLz;
        dispJp[ix][iy][iz] =
          + 2.0*parJ*( cos(kx)+cos(ky) )
          + 2.0*parJp*( cos(kx+ky) )
          + parJd*cos(kz);
      // J_Q-k
        kx = Qx - 2.0*pi*ix*invLmax;
        ky = Qy - 2.0*pi*iy*invLmax;
        kz = Qz - 2.0*pi*iz*invLz;
        dispJm[ix][iy][iz] =
          + 2.0*parJ*( cos(kx)+cos(ky) )
          + 2.0*parJp*( cos(kx+ky) )
          + parJd*cos(kz);
      }
    }
  }
  // J_Q
  kx = Qx;
  ky = Qy;
  kz = Qz;
  dispJQ =
    + 2.0*parJ*( cos(kx)+cos(ky) )
    + 2.0*parJp*( cos(kx+ky) )
    + parJd*cos(kz);

  for(ix=0; ix<Lmax; ix++){
    for(iy=0; iy<Lmax; iy++){
      for(iz=0; iz<Lz; iz++){
        dispA[ix][iy][iz] =
          + 0.25*(dispJp[ix][iy][iz] + dispJm[ix][iy][iz])
          + 0.5*dispJ[ix][iy][iz]
          - dispJQ;
//        dispW[ix][iy][iz] =
//          2.0*parS*sqrt( (dispJ[ix][iy][iz] - dispJQ)*(0.5*(dispJp[ix][iy][iz] + dispJm[ix][iy][iz]) - dispJQ) );
        tmp = (dispJ[ix][iy][iz] - dispJQ)*(0.5*(dispJp[ix][iy][iz] + dispJm[ix][iy][iz]) - dispJQ);
        if(tmp < 0){
          printf("# error: for parJp=%f and parJd=%f, dispW(%d,%d,%d)^2=%f < 0 !\n",parJp,parJd,ix,iy,iz,tmp);
          if(fabs(tmp) < 1.0e-10) tmp = 1.0e-32;
        }
        dispW[ix][iy][iz] = 2.0*parS*sqrt(tmp);
      }
    }
  }

  mag = parS + 0.5;
  for(ix=0; ix<Lmax; ix++){
    for(iy=0; iy<Lmax; iy++){
      for(iz=0; iz<Lz; iz++){
        if(dispW[ix][iy][iz] > 1.0e-10){
          mag += -parS/Lmax/Lmax/Lz * dispA[ix][iy][iz]/dispW[ix][iy][iz];
//          mag += -0.5/Lmax/Lmax_Lz * dispA[ix][iy][iz]/dispW[ix][iy][iz];
        }
      }
    }
  }

  ene = dispJQ*parS*parS;
  for(ix=0; ix<Lmax; ix++){
    for(iy=0; iy<Lmax; iy++){
      for(iz=0; iz<Lz; iz++){
        ene += 1.0/Lmax/Lmax/Lz * (
//          + parS*(dispW[ix][iy][iz] - dispA[ix][iy][iz])
          + 0.5*dispW[ix][iy][iz] - parS*dispA[ix][iy][iz]
//          + 0.5*(dispW[ix][iy][iz] - dispA[ix][iy][iz])
//          + dispW[ix][iy][iz] // alpha_k^dag alpha_k = 0 for GS
          );
      }
    }
  }

/*
  for(ix=0; ix<Lmax; ix++){
    for(iy=0; iy<Lmax; iy++){
      for(iz=0; iz<Lz; iz++){
        kx = 2.0*pi*ix*invLmax;
        ky = 2.0*pi*iy*invLmax;
        kz = 2.0*pi*iz*invLz;
//        printf("%d %d %f %f %f %f %f\n",ix,iy,kx,ky,dispJ[ix][iy],dispA[ix][iy],dispW[ix][iy]);
        printf("%d %d %d %f %f %f %f %f %f\n",ix,iy,iz,kx,ky,kz,dispJ[ix][iy][iz],dispA[ix][iy][iz],dispW[ix][iy][iz]);
      }
      printf("\n");
    }
    printf("\n");
  }
  printf("\n");
  printf("# %d %f\n",Lmax,mag);
*/

//  printf("%d %f %f %f %f\n",Lmax,parJp,mag,dispJQ*parS*parS,ene);
//  printf("%d %f %f %f %f %f %f %f\n",Lmax,parJ,parJp,parS,mag,dispJQ*parS*parS,ene,parQ/pi);
  printf("%d %f %f %f %f %f %f %f %f\n",Lmax,parJ,parJp,parJd,parS,mag,dispJQ*parS*parS,ene,parQ/pi);

//}

  free(dispJ[0][0]);
  free(dispJ[0]);
  free(dispJ);
  free(dispJp[0][0]);
  free(dispJp[0]);
  free(dispJp);
  free(dispJm[0][0]);
  free(dispJm[0]);
  free(dispJm);
  free(dispA[0][0]);
  free(dispA[0]);
  free(dispA);
  free(dispW[0][0]);
  free(dispW[0]);
  free(dispW);

  return 0;
}
