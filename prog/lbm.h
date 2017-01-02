#ifndef _LBM_H
#define _LBM_H


// nombre de variables conservées
#define M 3

// nombre de vitesses
#define NV 9

// nombre de dimensions
#define DIM 2

// vitesse du son
#define C 0.6

// taille de la grille
#define N1 40
#define N2 20

#define L2 1.
#define L1 (L2 * N1 / N2)

// pas en espace
#define DX (L2 / N2)


typedef struct lbm_data{

  // fonction de distribution now
  double fnow[NV][N1][N2];
  // fonction de distribution next
  double fnext[NV][N1][N2];

  // poids
  double omega[NV];

  // directions des vitesses
  int vit_dir[NV][DIM];
  
  // pas en temps
  double dt;

  // temps courant
  double tnow;

  // temps final
  double tmax;
  
} lbm_data;

void lbm_init(lbm_data *lbm, double tmax);

// condition initiale
void cond_init(double x, double y, double t, double *rho, double *u1, double *u2);

void lbm_display(lbm_data *lbm);

// calcule une fonction de distribution équilibre
// en fonction de rho, rho u1, rho u2 (stocké dans w)
void macro2micro(lbm_data *lbm, double *w, double* feq);

// calcule rho, rho u1, rho u2 (stocké dans w) à partir
// de f[k]
void micro2macro(lbm_data *lbm, double *f, double *w);

#endif
