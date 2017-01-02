#include "lbm.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


int main(void){

  lbm_data lbm;

  lbm_init(&lbm, 1.);

  lbm_display(&lbm);

  return 0;

}

void lbm_init(lbm_data *lbm, double tmax){

  lbm->omega[0] = 4./9;
  lbm->omega[1] = 1./9;
  lbm->omega[2] = 1./9;
  lbm->omega[3] = 1./9;
  lbm->omega[4] = 1./9;
  lbm->omega[5] = 1./36;
  lbm->omega[6] = 1./36;
  lbm->omega[7] = 1./36;
  lbm->omega[8] = 1./36;

  lbm->vit_dir[0][0] = 0;
  lbm->vit_dir[0][1] = 0;

  lbm->vit_dir[1][0] = 1;
  lbm->vit_dir[1][1] = 0;

  lbm->vit_dir[2][0] = 0;
  lbm->vit_dir[2][1] = 1;

  lbm->vit_dir[3][0] = -1;
  lbm->vit_dir[3][1] = 0;

  lbm->vit_dir[4][0] = 0;
  lbm->vit_dir[4][1] = -1;

  lbm->vit_dir[5][0] = 1;
  lbm->vit_dir[5][1] = 1;

  lbm->vit_dir[6][0] = -1;
  lbm->vit_dir[6][1] = 1;

  lbm->vit_dir[7][0] = -1;
  lbm->vit_dir[7][1] = -1;

  lbm->vit_dir[8][0] = 1;
  lbm->vit_dir[8][1] = -1;

  lbm->tmax = tmax;
  lbm->tnow = 0;

  lbm->dt = DX / sqrt(3.) / C;

  double w[M];
  double f[NV];
  for(int i = 0; i < N1; i++){
    for(int j = 0; j < N2; j++){
      double x = i * DX;
      double y = j * DX;
      double rho, u1, u2;
      cond_init(x, y, lbm->tnow, &rho, &u1, &u2);
      w[0] = rho;
      w[1] = rho * u1;
      w[2] = rho * u2;
      macro2micro(lbm, w, f);
      for(int k = 0; k < NV; k++){
	lbm->fnow[k][i][j] = f[k];
      }
    }
  }
  
  for(int k = 0; k < NV; k++){
     for(int i = 0; i < N1; i++){
        for(int j = 0; j < N2; j++){
	  lbm->fnext[k][i][j] = 0;
	}
     }
  }
}

void lbm_display(lbm_data *lbm){

  printf("tmax=%f tnow=%f dt=%f \n",lbm->tmax, lbm->tnow, lbm->dt);
  
  printf("poids=");
  for(int k = 0; k < NV; k++){
    printf("%f ", lbm->omega[k]);
  }
  printf("\n");
  
  for(int k = 0; k < NV; k++){
    printf("vitesse %d = %d %d \n",k,
	   lbm->vit_dir[k][0],
	   lbm->vit_dir[k][1]);
  }

   for(int k = 0; k < NV; k++){
     for(int i = 0; i < N1; i++){
        for(int j = 0; j < N2; j++){
	  printf("k=%d i=%d j=%d fnow=%f fnext=%f \n",
		 k, i, j,
		 lbm->fnow[k][i][j],
		 lbm->fnext[k][i][j]);
	}
     }
   }

   // gnuplot 
   // plot 'lbm_rho.dat' matrix with image
   // set pm3d map
   // splot 'lbm_rho.dat' matrix
   // set pm3d interpolate 2,2
   // splot 'lbm_rho.dat' matrix
   
   // display in a gnuplot file
   FILE * gnufile;
   gnufile = fopen("lbm_rho.dat", "w" );
   for(int j = 0; j < N2; j++){
     for(int i = 0; i < N1; i++){
	  double f[NV];
	  double w[M];
	  for(int k = 0; k < NV; k++){
	    f[k] = lbm->fnow[k][i][j];
	  }
	  micro2macro(lbm, f, w);
	  fprintf(gnufile, "%f ", w[0]);
	}
	fprintf(gnufile, "\n");
     }
     fclose(gnufile);
     system("gnuplot gnuplot_com");
}

void micro2macro(lbm_data *lbm, double *f, double* w){


  w[0]=0;  // rho
  w[1]=0;  // rho u1
  w[2]=0;  // rho u2
  
  for(int k = 0; k < NV; k++){
    w[0] += lbm->omega[k] * f[k];
    printf("w=%f f=%f k=%d \n",w[0],f[k],k);
    w[1] += lbm->omega[k] * lbm->vit_dir[k][0]  * sqrt(3.) * C * f[k];
    w[2] += lbm->omega[k] * lbm->vit_dir[k][1]  * sqrt(3.) * C * f[k];
  }

}

void macro2micro(lbm_data *lbm, double *w, double *feq){
  
  double r = w[0];
  double u1 = w[1] / r;
  double u2 = w[2] / r;
  double uu = u1 * u1 + u2 * u2;
  
  for(int k = 0; k < NV; k++){
    double vu = lbm->vit_dir[k][0] * sqrt(3.) * C * u1 +
      lbm->vit_dir[k][1] * sqrt(3.) * C * u2;
    double vuc = vu / C / C;
    feq[k] = r * (1 + vuc + vuc * vuc - uu / 2 / C / C);
    printf("k=%d feq=%f uu=%f vuc=%f \n",k,feq[k],uu,vuc); 
  }
}

// condition initiale
void cond_init(double x, double y, double t, double *rho, double *u1, double *u2){

  *rho = atan(x-1) * sin(2 * 3.14 * y) + 10;
  *u1 = 0;
  *u2 = 0;

}


