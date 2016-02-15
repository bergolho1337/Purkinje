#ifndef PURKINJE_H
#define PURKINJE_H

#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <ctime>
#include <random>

#define Iter 6				// Numero de iteracoes de crescimento
#define l_bra 1.8			// Tamanho do ramo (mm)
#define w_1 0.75			// Peso especificado pelo usuário, ligado ao quanto a fibra deve desviar
#define t_cube 0.5			// Dimensão do mini-cubo de Sobel
#define L 200.0				// Lado do cubo que envolve a estrutura da árvore (domínio)
#define sigma_c 0.2			// Condutividade padrao citoplasmatica (mS/mm)

using namespace std;

class Edge;
class Node;

// Estrutura de uma célula na fibra. É utilizada para armazenar dados da solução da equação do monodomínio.
struct Cell
{
	int type;			// Tipo de célula.  
	double *y;			// Resolução da equação no instante n. 	
	double *y_new;			// Resolução da equação no instante n+1
	double sigma;			// Condutividade    		
	double V_star;			// V*
}typedef Cell;

// Estrutura de uma aresta no grafo
class Edge
{
private:
	int id;				// Identificador do no destino
	double size;			// Tamanho da aresta, distancia euclidiana
	Edge *next;			// Ponteiro para a proxima aresta
	Node *dest;			// Ponteiro para o no destino
public:
	Edge ();
	Edge (int id, double size);
// Inline
	void setId (int id) { this->id = id; }
	void setSize (double size) { this->size = size; }
	void setNext (Edge *next) { this->next = next; }
	void setDest (Node *dest) { this->dest = dest; }
	int getId () { return (id); }
	double getSize () { return (size); }
	Edge* getNext () { return (next); }
	Node* getDest () { return (dest); }
};

// Estrutura de nó do grafo
class Node
{
private:
	int id;				// Identificador do nó
	double x, y, z;			// Coordenadas (x,y,z)
	int num_edges;			// Contador do número de arestas
	bool grow;			// Controlador se o no esta ou nao em crescimento 
	Cell *cell;			// Estrutura células contém dados sobre a solução da equação do monodomínio
	Node *next;			// Ponteiro para o próximo nó na lista de nós
	Edge *edges;			// Ponteiro para a lista de arestas
	
public:
// Construtores
	Node ();
	Node (int id, double x, double y, double z, bool grow);
	Node (int id, double x, double y, double z, double sigma);
// Inline
	void setId (int id) { this->id = id; }
	void setX (double x) { this->x = x; }
	void setY (double y) { this->y = y; }
	void setZ (double z) { this->z = z; }
	void setNumEdges (int num_edges) { this->num_edges = num_edges; }
	void setGrow (bool grow) { this->grow = grow; }
	void setCell (Cell *cell) { this->cell = cell; }
	void setNext (Node *next) { this->next = next; }
	void setEdges (Edge *edges) { this->edges = edges; }
	int getId () { return (id); }
	double getX () { return (x); }
	double getY () { return (y); }
	double getZ () { return (z); }
	bool getGrow () { return (grow); }
	int getNumEdges () { return (num_edges); }
	Cell* getCell () { return (cell); }
	Node* getNext () { return (next); }
	Edge* getEdges () { return (edges); }
	
};

// Estrutura do grafo
class Graph
{
private:
	Node *listNodes;			// Ponteiro para a lista de nós
	Node *lastNode;				// Ponteiro para último nó da lista de nós
	int total_nodes;			// Contador de nós
	int total_edges;			// Contador de arestas
	double Xmax;				// Valor máximo do domínio no eixo x

// Funções privadas
	double calcNorm (double x1, double y1, double z1, double x2, double y2, double z2)	// Calcula a distancia euclidiana
	{
		return (sqrt(pow((x1-x2),2) + pow((y1-y2),2) + pow((z1-z2),2)));
	}
	bool isInsideCube (Node *ptr, double x, double y, double z, double width, double lenght, double height)
	{
		return (ptr->getX() <= x + lenght && ptr->getX() >= x && ptr->getY() >= y - width && ptr->getY() <= y && ptr->getZ() >= z - height && ptr->getZ() <= z);
	}

public:
// Construtores
	Graph (double *y0, int n_eq, double Dx);
	Graph (FILE *file, double *y0, int n_eq, double Dx);
// Inline
	void setListNodes (Node *n) { listNodes = n; }
	void setLastNode (Node *n) { lastNode = n; }
	void setTotalNodes (int tn) { total_nodes = tn; }
	void setTotalEdges (int te) { total_edges = te; }
	void setXmax (int xm) { Xmax = xm; }
	Node* getListNodes () { return (listNodes); }
	Node* getLastNode () { return (lastNode); }
	int getTotalNodes () { return (total_nodes); }
	int getTotalEdges () { return (total_edges); }
	double getXmax () { return (Xmax); }
// Funções
	void makeRoot (double Dx, int n_eq);
	void makeRoot2 (double Dx, int n_eq);
	void initializeRandomGenerator ();
	Node* searchNode (int id);
	void insertNodeGraph (Node *p);
	void insertEdgeGraph (int id_1, int id_2);
	void growBranch (Node *p);
	void calculateDirection (double *d_gra, Node *p, double teta);
	void mallocSobel ();
	void Sobel (Node *p, int type);
	double calculateSizeBranch ();
	void makeBranch (Node *p);
	bool checkCollision (Node *p, double x, double y, double z);
	bool checkLimits (Node *p, double x, double y, double z);
	void printSobel ();
	void printGradient ();
	void printGraph ();
	void writeVTK (int k);
}; 

#endif
