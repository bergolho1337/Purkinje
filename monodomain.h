#ifndef MONODOMAIN_H
#define MONODOMAIN_H

#include <iostream>
#include <iomanip>
#include <cmath>
#include "purkinje.h"

// Define um ponteiro para funcao que recebe f(t,y) como parametro e retorna um 'double'
typedef
	double (*Func) (double t, double *y);

void setInitialConditions (Graph *g, double *y0, int n_eq);
void solveODE (Graph *g, Func *func, double t, double Dt, int n_eq);
void solvePDE (Graph *g, Func *func, double t, double Cm, double beta, double Dt, double Dx, int n_eq);
void solvePDE2 (Graph *g, Func *func, double t, double Cm, double beta, double Dt, double Dx, int n_eq);
void checkCFL (double r);
void nextTimestep (Graph *g, int n_eq);

void writeGraphic (Graph *g, double t, int plot, int n_eq, int num_method);
void plotGraphic (int plot, int num_method);

// IMPLICITO
void makeMatrix_A (Graph *g, double beta, double Cm, double Dt, double Dx);
void makeMatrix_A2 (Graph *g, double beta, double Cm, double Dt, double Dx);
void makeVector_b (Graph *g, double beta, double Cm, double Dt, double Dx);
void solvePDE_Imp (Graph *g);
void solveODE_Imp (Graph *g, Func *func, double t, double Cm, double Dt, int n_eq);
void solveODE_Imp2 (Graph *g, Func *func, double t, double Cm, double Dt, int n_eq);
double* LU (int n);

#endif
