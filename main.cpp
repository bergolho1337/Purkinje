/* ----------------------------------------------------------------------------------------
	MODELO MACROSCOPICO DAS FIBRAS DE PURKINJE
   ----------------------------------------------------------------------------------------
	1) Geracao automatica das fibras eh feita usando o L-System modificado;
	2) Modelo celular utilizado eh o de Mitchell Shaeffer (Two Current Model);
	3) Resolucao da equacao do monodominio nas celulas eh feita por volumes finitos;
   ---------------------------------------------------------------------------------------- 
	Try example: ./macroscopic 500 1.8 0.05
*/

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include "purkinje.h"
#include "monodomain.h"
#include "mitchell.h"

using namespace std;

// Define um ponteiro para funcao que recebe f(t,y) como parametro e retorna um 'double'
typedef
	double (*Func) (double t, double *y);

int main (int argc, char *argv[])
{
	clock_t begin, end;
	double Dt, Dx;
	double t_max;
	double elapsed;
	Func *f;
	double *y_old;
	int plot = 19;
	int num_method;

	if (argc-1 < 4)
	{
		cout << "================= MACROSCOPIC PURKINJE MODEL ====================" << endl;
		cout << "Usage:> ./macroscopic <t_max> <Dx> <Dt> <number_method_PDE> [in.txt]" << endl;
		cout << "t_max = Maximum time of the simulation. (ms)" << endl;
		cout << "Dx = Size of the discretization on space. (mm)" << endl;
		cout << "Dt = Size of the discretization on time. (ms)" << endl;
		cout << "number_method_PDE = Number of the method that solves the PDE" << endl;
		cout << "\t1 = FTCS (Forward in Time Centered in Space)" << endl;
		cout << "\t2 = BTCS (Backward in Time Centered in Space)" << endl;
		cout << "\t3 = Roda os dois metodos para a comparacao" << endl;
		cout << "[in.txt] = An optional file that can be pass to construct the Purkinje structure" << endl;
	}
	else
	{
		int N, M;
		double t;
		Graph *g;

		// Dispara o cronometro
		begin = clock();
		// Cria diretorio para armazenar os .VTK da animacao
		system("mkdir VTK");
		
		// Inicializar o modelo celular
		f = getFunctions();
		y_old = getInitialConditions();

		// Receber os parametros de entrada
		t_max = atof(argv[1]);
		Dx = atof(argv[2]);
		Dt = atof(argv[3]);
		num_method = atof(argv[4]);

		if (argc-1 == 5)
		{
			FILE *file = fopen(argv[5],"r");
			if (!file)
			{
				cout << "[-] ERROR! File \"" << argv[5] << "\" not found!" << endl;
				exit(1);
			}
			// Criando grafo a partir da entrada
			g = new Graph(file,y_old,num_equations,Dx);
			setInitialConditions(g,y_old,num_equations);
			// Calcula o numero de subintervalos de discretizacao no espaco e no tempo
			M = nearbyint(t_max/Dt);
			N = nearbyint(g->getXmax()/Dx);
			cout << "[!] Solving monodomain equation ..." << endl;
			cout << "\tSubintervals in space: N = " << N << endl;
			cout << "\tSubintervals in time: M = " << M << endl;
			switch (num_method)
			{
				case 1: cout << "\tUsing FTCS ..." << endl;
					break;
					
				case 2: cout << "\tUsing BTCS ..." << endl;
						// Se o método for implícito monta a matriz A e gera a decomposição LU
						if (num_method == 2)
							makeMatrix_A(g,beta,Cm,Dt,Dx);
							break;
			}
			// Loop no tempo
			for (int j = 0; j < M; j++)
			{
				t = j*Dt;
				writeGraphic(g,t,plot,num_equations,num_method);
				if (M > 1000)
				{
					if (j % 10 == 0)
						g->writeVTK(j);
				}
				else
					g->writeVTK(j);
				// Desvia para qual metodo foi selecionado
				switch (num_method)
				{
					// FTCS
					case 1: {
							// Para cada celula resolver a EDO primeiro, loop no espaco esta dentro desse procedimento
							solveODE(g,f,t,Dt,num_equations);
							// Agora cada célula possui um V* e calculamos o V_n+1 de cada célula pela EDP
							solvePDE2(g,f,t,Cm,beta,Dt,Dx,num_equations);
							break;
						}
					// BTCS
					case 2: {
							solveODE_Imp(g,f,t,Cm,Dt,num_equations);
							makeVector_b(g,beta,Cm,Dt,Dx);
							solvePDE_Imp(g);
							break;
						}
				}
				// Avanca no tempo
				nextTimestep(g,num_equations);
			}
		}
		// Nao foi passado arquivo de entrada, logo gerar automaticamente as fibras
		else
		{
			g = new Graph(y_old,num_equations,Dx);
			setInitialConditions(g,y_old,num_equations);

			// Calcula o numero de subintervalos de discretizacao no espaco e no tempo
			M = nearbyint(t_max/Dt);
			N = nearbyint(g->getXmax()/Dx);
			cout << "[!] Solving monodomain equation ..." << endl;
			cout << "\tSubintervals in space: N = " << N << endl;
			cout << "\tSubintervals in time: M = " << M << endl;
			switch (num_method)
			{
				case 1: cout << "\tUsing FTCS ..." << endl;
					break;
				case 2: cout << "\tUsing BTCS ..." << endl;
					break;
			}
			// Se o método for implícito monta a matriz A e gera a decomposição LU
			if (num_method == 2)
				makeMatrix_A(g,beta,Cm,Dt,Dx);

			// Loop no tempo
			for (int j = 0; j < M; j++)
			{
				//cout << "\nIteration " << j << endl;
				t = j*Dt;
				writeGraphic(g,t,plot,num_equations,num_method);
				if (M > 1000)
				{
					if (j % 10 == 0)
						g->writeVTK(j);
				}
				else
					g->writeVTK(j);
				switch (num_method)
				{
					// FTCS
					case 1: {
							// Para cada celula resolver a EDO primeiro, loop no espaco esta dentro desse procedimento
							solveODE(g,f,t,Dt,num_equations);
							// Agora cada célula possui um V* e calculamos o V_n+1 de cada célula pela EDP
							solvePDE2(g,f,t,Cm,beta,Dt,Dx,num_equations);
							break;
						}
					// BTCS
					case 2: {
							// Resolve EDO primeiro
							solveODE_Imp(g,f,t,Cm,Dt,num_equations);
							makeVector_b(g,beta,Cm,Dt,Dx);
							solvePDE_Imp(g);
							break;
						}
				}
				// Avanca no tempo
				nextTimestep(g,num_equations);
			}
			//g->writeVTK(0);
		}		
		plotGraphic(plot,num_method);
		g->printGraph();		

		end = clock();
		elapsed = (double)(end-begin) / (double)CLOCKS_PER_SEC;
		cout << "Time elapsed: " << elapsed << " s" << endl;

		delete [] f;
		delete [] y_old;

		cout << "[+] Simulation completed with sucess!" << endl;
	}
}	
