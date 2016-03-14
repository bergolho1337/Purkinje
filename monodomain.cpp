#include "monodomain.h"

using namespace std;

double **A;					// Matriz do metodo implicito
double *b;					// Vetor do metodo implicito
double *pivot;				// Vetor dos pivos

// Atribuir a condicao inicial do problema
void setInitialConditions (Graph *g, double *y0, int n_eq)
{
	int i;
	cout << "[!] Setting the initial conditions." << endl;
	Node *ptr = g->getListNodes();
	// Na primeira celula setamos as condicoes iniciais de estimulo
	for (i = 0; i < n_eq; i++)
		ptr->getCell()->y[i] = y0[i];
	ptr = ptr->getNext();
	// Nas proximas celulas setamos o estimulo no estado de relaxamento.
	while (ptr != NULL)
	{
		for (i = 0; i < n_eq; i++)
		{
			if (i == 0)
				ptr->getCell()->y[i] = 0.0;
			else
				ptr->getCell()->y[i] = y0[i];
		}
		//cout << "Cell " << ptr->getId() << " v_stim = " << ptr->getCell()->y[0] << " h0 = " << ptr->getCell()->y[1] << endl;
		ptr = ptr->getNext();
	}
}

// Resolve o sistema não-linear de EDOs utilizando o método Euler explícito (Reação)
void solveODE (Graph *g, Func *func, double t, double Dt, int n_eq)
{
	Node *ptr;
	Cell *c;
	int i;
	double f[n_eq];
	
	ptr = g->getListNodes();
	// Loop no espaco
	while (ptr != NULL)
	{
		c = ptr->getCell();
		for (i = 0; i < n_eq; i++)
		{
			// Calcular o potencial transmembranico intermediário -> V_{i/2} = V*
			if (i == 0)
			{
				f[i] = func[i](t,c->y)*Dt;
				c->V_star = c->y[i] + f[i];
			}
			// Calcula já o valor das variaveis gate para o próximo passo de tempo -> h_{i+1}
			else
			{
				f[i] = func[i](t,c->y)*Dt;
				c->y_new[i] = c->y[i] + f[i];
			}
		}
		ptr = ptr->getNext();
	}
}

// Resolve a EDP relacionada por meio do método FTCS. ** Considerando que a condutividade é igual para todos os volumes ** (Difusao)
void solvePDE (Graph *g, Func *func, double t, double Cm, double beta, double Dt, double Dx, int n_eq)
{
	Node *ptr;
	Cell *c;
	double f, r;

	ptr = g->getListNodes();
	while (ptr != NULL)
	{
		c = ptr->getCell();
		// Verifica quantas arestas estao ligados ao volume de controle atual
		switch (ptr->getNumEdges())
		{
			case 1: {
					// Se for a primeira celula, so resolver a EDO
					if (ptr->getId() == 1)
					{
						f = func[0](t,c->y)*Dt; 
						c->y_new[0] = c->y[0] + f;
					}
					// Senao a celula eh folha, logo aplica-se condicao de Neumann dv/dt=0
					else
					{
						r = (c->sigma*Dt)/(Dx*Dx*Cm*beta);
						//cout << r << endl;
						if (r >= 0.5)
						{
							cout << "[-] ERROR! CFL condition was violated! r = " << r << endl;
							exit(1);
						}
						Cell *c2;
						c2 = ptr->getEdges()->getDest()->getCell();
						c->y_new[0] = r*(c2->V_star - c->V_star) + c->V_star; 
					}
					break;
				}
			case 2: {
					// No interior da malha, deve ser resolvida usando FTCP (Forward Time, Centered Space)
					r = (c->sigma*Dt)/(Dx*Dx*Cm*beta);
					if (r >= 0.5)
					{
						cout << "[-] ERROR! CFL condition was violated! r = " << r << endl;
						exit(1);
					}
					Cell *c2, *c3;
					c2 = ptr->getEdges()->getDest()->getCell();
					c3 = ptr->getEdges()->getNext()->getDest()->getCell();
					c->y_new[0] = r*(c2->V_star - 2*c->V_star + c3->V_star) + c->V_star;
					break;
				}
			case 3: {
					// Senao celula de bifurcacao
					r = (c->sigma*Dt)/(Dx*Dx*Cm*beta);
					if (r >= 0.5)
					{
						cout << "[-] ERROR! CFL condition was violated! r = " << r << endl;
						exit(1);
					}
					Cell *c2, *c3, *c4;
					c2 = ptr->getEdges()->getDest()->getCell();
					c3 = ptr->getEdges()->getNext()->getDest()->getCell();
					c4 = ptr->getEdges()->getNext()->getNext()->getDest()->getCell();
					c->y_new[0] = r*(c2->V_star + c3->V_star + c4->V_star - 3*c->V_star) + c->V_star;
					break;
				 }
			default: {
					// Celula de bifurcacao caso geral
					r = (c->sigma*Dt)/(Dx*Dx*Cm*beta);
					if (r >= 0.5)
					{
						cout << "[-] ERROR! CFL condition was violated! r = " << r << endl;
						exit(1);
					}
					double V_star = 0.0;
					Edge *ptrl = ptr->getEdges();
					while (ptrl != NULL)
					{
						V_star += ptrl->getDest()->getCell()->V_star;
						ptrl = ptrl->getNext();
					}
					c->y_new[0] = r*(V_star - (ptr->getNumEdges()+1)*c->V_star) + c->V_star;
				 }
		}
		ptr = ptr->getNext();
	}
}

// Resolve a EDP relacionada por meio do método FTCS. ** Caso geral ** (Difusao)
void solvePDE2 (Graph *g, Func *func, double t, double Cm, double beta, double Dt, double Dx, int n_eq)
{
	Node *ptr;
	Cell *c;
	double f, r;

	ptr = g->getListNodes();
	while (ptr != NULL)
	{
		c = ptr->getCell();
		// Verifica quantas arestas estao ligados ao volume de controle atual
		switch (ptr->getNumEdges())
		{
			case 1: {
					// Se for a primeira celula, so resolver a EDO
					if (ptr->getId() == 1)
					{
						f = func[0](t,c->y)*Dt; 
						c->y_new[0] = c->y[0] + f;
					}
					// Senao a celula eh folha, logo aplica-se condicao de Neumann dv/dt=0
					else
					{
						r = (c->sigma*Dt)/(Dx*Dx*Cm*beta);
						//cout << r << endl;
						if (r >= 0.5)
						{
							cout << "[-] ERROR! CFL condition was violated! r = " << r << endl;
							exit(1);
						}
						Cell *c2;
						c2 = ptr->getEdges()->getDest()->getCell();
						c->y_new[0] = r*(c2->V_star - c->V_star) + c->V_star; 
					}
					break;
				}
			default: {
					// Celula de bifurcacao caso geral
					Cell *c1, *c2;
					double V_star = 0.0;
					double sigma_med;
					Edge *ptrl = ptr->getEdges();
					c1 = ptr->getCell();
					while (ptrl != NULL)
					{
						c2 = ptrl->getDest()->getCell();
						sigma_med = (c1->sigma + c2->sigma)/2.0;
						r = (Dt*sigma_med)/(Dx*Dx*Cm*beta);
						checkCFL(r);
						// Corrente entrando no volume
						if (ptrl->getId() < ptr->getId())
							V_star += r*(c2->V_star - c1->V_star);
						// Corrente saindo do volume
						else
						{
							// Checa se uma corrente de reentrada esta prestes a acontecer, nesse caso liberar a corrente alterando a condutividade
							if (c->y[0] < 0.1 && c2->y[0] > 0.3)
							{
								r = (Dt*sigma_c)/(Dx*Dx*Cm*beta);
								checkCFL(r);
							}
							V_star -= r*(c1->V_star - c2->V_star);
						}
						ptrl = ptrl->getNext();
					}
					c->y_new[0] = V_star + c1->V_star;
				 }
		}
		ptr = ptr->getNext();
	}
}

// Constrói a matriz do método implícito e já gera a decomposição LU.
void makeMatrix_A (Graph *g, double beta, double Cm, double Dt, double Dx)
{
	Node *ptr;
	Edge *ptrl;
	int i;	
	int total_nodes;
	double alfa, sum;
	total_nodes = g->getTotalNodes();
	ptr = g->getListNodes();
	alfa = beta*Cm*Dx*Dx/Dt;	
	// Aloca memoria
	A = new double*[total_nodes]();
	for (i = 0; i < total_nodes; i++)
		A[i] = new double[total_nodes]();
	pivot = new double[total_nodes];
	
	while (ptr != NULL)
	{
		//cout << "Node " << ptr->getId() << endl;
		ptrl = ptr->getEdges();
		sum = 0.0;
		while (ptrl != NULL)
		{
			//cout << "Edge " << ptrl->getId() << endl;
			sum += ptrl->getDest()->getCell()->sigma;
			A[ptr->getId()-1][ptrl->getId()-1] = -ptrl->getDest()->getCell()->sigma;
			ptrl = ptrl->getNext();
		}
		sum += alfa;
		// Condição de Neumann nas folhas
		if (ptr->getNumEdges() == 1 && ptr->getId() != 1)
			A[ptr->getId()-1][ptr->getId()-1] = ptr->getCell()->sigma;
		else
			A[ptr->getId()-1][ptr->getId()-1] = sum;
		ptr = ptr->getNext(); 
	}
	// Decomposição LU
	int j, k, p;
	double Amax, t, m, r, Mult;
	// 1 PASSO: Transformar a matriz A do problema em duas matrizes triangulares L e U. 	
	for (i = 0; i < total_nodes; i++)
		pivot[i] = i;
	for (j = 0; j < total_nodes-1; j++)
	{
		// Escolher pivot
		p = j;
		Amax = abs(A[j][j]);
		// Verifica na coluna a ser eliminada qual elemento possui o maior valor absoluto, este elemento será o pivô.		
		for (k = j+1; k < total_nodes; k++)
		{
			if (abs(A[k][j]) > Amax)
			{
				Amax = abs(A[k][j]);
				p = k;
			}
		}
		// Se (p != j) então deve-se trocar de linhas 
		if (p != j)
		{
			for (k = 0; k < total_nodes; k++)
			{
				t = A[j][k];
				A[j][k] = A[p][k];
				A[p][k] = t;	
			}
			m = pivot[j];
			pivot[j] = pivot[p];
			pivot[p] = m;
		}
		if (abs(A[j][j]) != 0)
		{
			// Eliminação de Gauss
			r = 1 / A[j][j];
			for (i = j+1; i < total_nodes; i++)
			{
				Mult = A[i][j]*r;
				A[i][j] = Mult;
				for (k = j+1; k < total_nodes; k++)
					A[i][k] = A[i][k] - Mult*A[j][k];
			}
		}
	}
}

void makeMatrix_A2 (Graph *g, double beta, double Cm, double Dt, double Dx)
{
	Node *ptr;
	Edge *ptrl;
	int i;	
	int total_nodes;
	double alfa, sum;
	total_nodes = g->getTotalNodes();
	ptr = g->getListNodes();
	alfa = beta*Cm*Dx*Dx/Dt;	
	// Aloca memoria
	A = new double*[total_nodes]();
	b = new double[total_nodes]();
	for (i = 0; i < total_nodes; i++)
		A[i] = new double[total_nodes]();
	
	while (ptr != NULL)
	{
		ptrl = ptr->getEdges();
		sum = 0.0;
		while (ptrl != NULL)
		{
			sum += ptrl->getDest()->getCell()->sigma;
			A[ptr->getId()-1][ptrl->getId()-1] = -ptrl->getDest()->getCell()->sigma;	
			ptrl = ptrl->getNext();
		}
		sum += alfa;
		// Condição de Neumann nas folhas
		if (ptr->getNumEdges() == 1 && ptr->getId() != 1)
			A[ptr->getId()-1][ptr->getId()-1] = ptr->getCell()->sigma;
		else
			A[ptr->getId()-1][ptr->getId()-1] = sum;
		ptr = ptr->getNext(); 
	}
}

void makeVector_b (Graph *g, double beta, double Cm, double Dt, double Dx)
{
	Node *ptr;
	double alfa;
	b = new double[g->getTotalNodes()];
	alfa = beta*Cm*Dx*Dx/Dt; 
	ptr = g->getListNodes();
	//cout << "Vector b" << endl;
	while (ptr != NULL)
	{
		// Condição de Neumann
		if (ptr->getNumEdges() == 1 && ptr->getId() != 1)
			b[ptr->getId()-1] = 0.0;
		else
			b[ptr->getId()-1] = alfa*ptr->getCell()->V_star;
		//cout << b[ptr->getId()-1] << endl;
		ptr = ptr->getNext();
	}
	//cout << endl;
	//exit(1);
}

void solvePDE_Imp (Graph *g)
{
	Node *ptr;
	// Resolve o sistema do metodo implicito
	double *Vnext = LU(g->getTotalNodes());
	// Seta os valores de V_n+1 nos volumes de controle
	ptr = g->getListNodes();
	//cout << "Vnext" << endl;
	//for (int i = 0; i < g->getTotalNodes(); i++)
	//	cout << Vnext[i] << endl;
	//cout << "============" << endl;
	while (ptr != NULL)
	{
		ptr->getCell()->y_new[0] = Vnext[ptr->getId()-1];
		//cout << "Cell " << ptr->getId() << " " << ptr->getCell()->y_new[0] << endl;
		ptr = ptr->getNext();
	}
	//cout << "========================" << endl;
	delete [] Vnext;
	delete [] b;
	//exit(1);
}

void solveODE_Imp (Graph *g, Func *func, double t, double Cm, double Dt, int n_eq)
{
	Node *ptr;
	Cell *c;
	int i;
	double f[n_eq];
	
	// Na primeira célula mantemos o estímulo sendo gerado: V^n --> V^{n+1}
	ptr = g->getListNodes();
	//cout << "V*\tVn" << endl;
	// Nas outras células considera-se reação-difusão
	while (ptr != NULL)
	{
		c = ptr->getCell();
		for (i = 0; i < n_eq; i++)
		{
			// Calcular o potencial transmembranico intermediario: V^{n} --> V^{*}
			if (i == 0)
			{
				f[i] = func[i](t,c->y)*Dt;
				c->V_star = c->y[i] + f[i];
				//cout << "Cell " << ptr->getId() << " " << c->V_star << "\t" << c->y[0] << endl;
			}
			// Calcula já o valor das variaveis gate para o próximo passo de tempo: h_{i} --> h_{i+1}
			else
			{
				f[i] = func[i](t,c->y)*Dt;
				c->y_new[i] = c->y[i] + f[i];
			}
		}
		ptr = ptr->getNext();
	}
	//exit(1);
}

void solveODE_Imp2 (Graph *g, Func *func, double t, double Cm, double Dt, int n_eq)
{
	Node *ptr;
	Cell *c;
	int i;
	double f[n_eq];
	
	// Na primeira célula mantemos o estímulo sendo gerado: V^n --> V^{n+1}
	ptr = g->getListNodes();
	//cout << "V*\tVn" << endl;
	// Nas outras células considera-se reação-difusão
	while (ptr != NULL)
	{
		c = ptr->getCell();

		// Checa se a celula esta estimulada
		if (c->y[0] > 0.2)
			ptr->setStimulated(true);
		else
			ptr->setStimulated(false);
		
		for (i = 0; i < n_eq; i++)
		{
			// Calcular o potencial transmembranico intermediario: V^{n} --> V^{*}
			if (i == 0)
			{
				f[i] = func[i](t,c->y)*(Dt/2.0);
				c->V_star = c->y[i] + f[i];
				//cout << "Cell " << ptr->getId() << " " << c->V_star << "\t" << c->y[0] << endl;
			}
			// Calcula já o valor das variaveis gate para o próximo passo de tempo: h_{i} --> h_{i+1}
			else
			{
				f[i] = func[i](t,c->y)*Dt;
				c->y_new[i] = c->y[i] + f[i];
			}
		}
		ptr = ptr->getNext();
	}

}

void writeGraphic (Graph *g, double t, int plot, int n_eq, int num_method)
{
	int i;
	char filename[20];
	sprintf(filename,"data%d.dat",num_method);
	FILE *file = fopen(filename,"a");
	Node *ptr = g->getListNodes();
	while (ptr != NULL)
	{
		if (ptr->getId() == plot)
			break;
		ptr = ptr->getNext();
	}
	fprintf(file,"%e ",t);
	for (i = 0; i < n_eq; i++)
		fprintf(file," %e",ptr->getCell()->y[i]);
	fprintf(file,"\n");
	fclose(file);
}

void plotGraphic (int plot, int num_method)
{
	char filename[20];
	sprintf(filename,"data%d.dat",num_method);
	FILE *arq = fopen("graph.plt","w+");
	fprintf(arq,"set grid\n");
	fprintf(arq,"set terminal png\n");
	fprintf(arq,"set output \"monodomain.png\"\n");
	fprintf(arq,"set title \"x = %d\"\n",plot);
	fprintf(arq,"plot \"%s\" using 1:2 title \"Vm\" w l\n",filename);
	fclose(arq);
	if (!system("gnuplot graph.plt"))
		cout << "[+] Grafico plotado com sucesso! Metodo " << num_method << endl;
	else
		cout << "[-] ERRO! Plotando grafico!" << endl;
}

// Resolve um sistema linear em que a matriz A já foi decomposta em 
double* LU (int n)
{
	int i, j, k;
	double *y = new double[n];
	double soma;
	k = pivot[0];
	y[0] = b[k];
	// Realizar substituições sucessivas para resolver o sistema triangular inferior: Ly = b
	for (i = 1; i < n; i++)
	{
		soma = 0;
		for (j = 0; j <= i-1; j++)
			soma += A[i][j]*y[j]; 
		k = pivot[i];
		y[i] = b[k] - soma;	
	}
	//cout << "On LU" << endl;
	// Realizar substituições retroativas para resolver o sistema triangular superior: Ux = y 			
	double *x = new double[n];
	x[n-1] = y[n-1] / A[n-1][n-1];
	for (i = n-2; i >= 0; i--)
	{
		soma = 0;
		for (j = i+1; j < n; j++)
			soma += A[i][j]*x[j];
		x[i] = (y[i] - soma) / A[i][i];
	}	

	delete [] y;

	return (x);
}

void nextTimestep (Graph *g, int n_eq)
{
	int i;
	Node *ptr = g->getListNodes();
	while (ptr != NULL)
	{
		for (i = 0; i < n_eq; i++)
		{
			ptr->getCell()->y[i] = ptr->getCell()->y_new[i];
		}
		//cout << endl;
		//cout << ptr->getCell()->y[0] << "\t" << ptr->getCell()->y_new[0] << endl;
		ptr = ptr->getNext();
	}
}

void checkCFL (double r)
{
	if (r >= 0.5)
	{
		cout << "[-] ERROR! CFL condition was violated! r = " << r << endl;
		exit(1);
	}
}