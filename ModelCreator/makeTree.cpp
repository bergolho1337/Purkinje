/*
	PROGRAMA USADO PARA CRIAR OS ARQUIVOS DE ENTRADA USADOS NAS SIMULACOES
*/

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>

using namespace std;

#define PI 3.14159265

const double Dx = 1.8;					// Tamanho de um segmento
const double sigma_c_N = 0.2;			// Condutividade normal
const double sigma_c_P = 0.0;			// Condutividade com problema
const int num_segments = 20;			// Numero de segementos por ramo

// Modelo 1: atividade normal das fibras Vs atividade paranormal 
void makeModel_1 ()
{
	int i;
	double x, y, dx, dy;
	FILE *file = fopen("in1.txt","w+");
	// Numero de nos e numero de arestas
	fprintf(file,"%d %d\n",num_segments*3,num_segments*3-2);
	// Segmento 1
	for (i = 0; i < num_segments; i++)
	{
		x = i*Dx;
		fprintf(file,"%e %e %e %e\n",x,0.0,0.0,sigma_c_N); 
	}
	dx = cos(PI/4.0)*Dx;
	dy = sin(PI/4.0)*Dx;
	// Segmento 2
	for (i = 0; i < num_segments; i++)
	{
		x += dx;	
		y += dy;
		fprintf(file,"%e %e %e %e\n",x,y,0.0,sigma_c_N);
	}
	x = i*Dx;
	y = 0.0;
	dx = cos(-PI/4.0)*Dx;
	dy = sin(-PI/4.0)*Dx;
	// Segmento 3
	for (i = 0; i < num_segments; i++)
	{
		x += dx;	
		y += dy;
		if (i < num_segments/2 || i > num_segments/2)
			fprintf(file,"%e %e %e %e\n",x,y,0.0,sigma_c_N);
		else
			fprintf(file,"%e %e %e %e\n",x,y,0.0,sigma_c_P);
	}
	// Arestas
	for (i = 0; i < num_segments-1; i++)
		fprintf(file,"%d %d\n",i+1,i+2);
	fprintf(file,"%d %d\n",i+1,i+2);
	fprintf(file,"%d %d\n",i,i+num_segments+2);
	for (i = num_segments; i < 2*num_segments-1; i++)
		fprintf(file,"%d %d\n",i+1,i+2);
	for (i = 2*num_segments; i < 3*num_segments-1; i++)
		fprintf(file,"%d %d\n",i+1,i+2);
	fclose(file);
}

// Modelo de teste para corrente de reentrada
void makeModel_2 ()
{
	int i;
	double x, y, dx, dy;
	double X1, Y1;
	double X2, Y2;
	FILE *file = fopen("in2.txt","w+");
	// Numero de nos e numero de arestas
	fprintf(file,"%d %d\n",num_segments*5,num_segments*5-1);
	// Segmento 1
	for (i = 0; i < num_segments; i++)
	{
		x = i*Dx;
		fprintf(file,"%e %e %e %e\n",x,0.0,0.0,sigma_c_N); 
	}
	dx = cos(PI/4.0)*Dx;
	dy = sin(PI/4.0)*Dx;
	// Segmento 2
	for (i = 0; i < num_segments; i++)
	{
		x += dx;	
		y += dy;
		fprintf(file,"%e %e %e %e\n",x,y,0.0,sigma_c_N);
	}
	x = i*Dx;
	y = 0.0;
	dx = cos(-PI/4.0)*Dx;
	dy = sin(-PI/4.0)*Dx;
	// Segmento 3
	for (i = 0; i < num_segments; i++)
	{
		x += dx;	
		y += dy;
		if (i < num_segments/3 || i > num_segments/3)
			fprintf(file,"%e %e %e %e\n",x,y,0.0,sigma_c_N);
		else
			fprintf(file,"%e %e %e %e\n",x,y,0.0,sigma_c_P);
	}
	// Segmento 4
	x -= 3*dx;
	y -= 3.5*dy;
	for (i = 0; i < 2*num_segments; i++)
	{
		y -= dy;
		fprintf(file,"%e %e %e %e\n",x,y,0.0,sigma_c_N/60.0);
	}
	// Arestas
	for (i = 0; i < num_segments-1; i++)
		fprintf(file,"%d %d\n",i+1,i+2);
	fprintf(file,"%d %d\n",i+1,i+2);
	fprintf(file,"%d %d\n",i,i+num_segments+2);
	for (i = num_segments; i < 2*num_segments-1; i++)
		fprintf(file,"%d %d\n",i+1,i+2);
	for (i = 2*num_segments; i < 3*num_segments-1; i++)
		fprintf(file,"%d %d\n",i+1,i+2);

	for (i = 3*num_segments; i < 5*num_segments-8; i++)
		fprintf(file,"%d %d\n",i+1,i+2);
	fprintf(file,"%d %d\n",i+1,2*num_segments-3);
	fprintf(file,"%d %d\n",3*num_segments-4,3*num_segments+1);
	fclose(file);
}


int main (int argc, char *argv[])
{
	int numModel;
	if (argc-1 < 1)
	{
		cout << "Usage:> ./makeTree <number_of_model>" << endl;
		exit(1);
	}
	else
	{
		numModel = atoi(argv[1]);
		switch (numModel)
		{
			case 1: makeModel_1();
				break;
			case 2: makeModel_2();
				break;
		}
		cout << "[+] Model file created with sucess!" << endl;
	}
}


