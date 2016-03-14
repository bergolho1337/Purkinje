/* ===================================================================================
	MITCHELL & SHAEFFER MODEL
    A Two-Current Model for the Dynamics of Cardiac Membrane (MITCHELL and SCHAEFFER)
   ===================================================================================
*/

#include <cmath>

// Define o tipo Func, igual esta no 'main.cpp'
typedef
	double (*Func) (double t, double *y);

// PARAMETROS DO MODELO

const char model_name[40] = "Mitchell & Shaeffer";

// Duas equacoes para resolver
const int num_equations = 2;

// Constantes relacionadas a equacao da voltagem (dV/dt)
const double t_close = 150.0;
const double t_open = 120.0;
const double t_out = 6.0;
const double t_in = 0.3;

// Estimulo inicial: Paper --> 0.1 mV
const double v_stim = 1.0;

// Mudanca de voltagem
const double v_gate = 0.13;

// Controle do tipo de estimulo
// true = (S1)  --> Single estimulus
// false = (S2) --> Double estimulus
bool single = true;

// Condicoes iniciais: variaveis gate
double h0 = 1.0;

// Parametros da equacao do monodominio
const double beta = 0.55556;	// (mm^-1)
const double Cm = 0.1;		// (mF/mm^2)

// Definicoes das funcoes do modelo
double dvdt (double t, double *y)
{
	return ( ((y[1]*pow(y[0],2)*(1-y[0]))/t_in) - (y[0]/t_out) );
}

double dhdt (double t, double *y)
{
	if (y[0] < v_gate)
		return ( (1-y[1])/t_open );
	else
		return ( -y[1]/t_close );
}

// **************************************************************************************************
//	[!] Essas funcoes devem estar implementadas em qualquer cabecalho de modelo celular !!!
// **************************************************************************************************
double* getInitialConditions ()
{
	double *y = new double[num_equations];
	y[0] = v_stim;
	y[1] = h0;
	return (y);
}

Func* getFunctions ()
{
	Func *f = new Func[num_equations];
	f[0] = dvdt;
	f[1] = dhdt;
	return (f);
}





