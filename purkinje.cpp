#include "purkinje.h"
#include "queue.h"

Queue *queue;						// Fila que controla os nos em crescimento
int num_eq;						// Numero de equacoes do modelo celular
int cont_segments = 1;					// Contador de segmentos criados para cada ramo (Maximo = 5)
double Toler = 0.3*l_bra;				// Tolerãncia (distância) para juntar segmentos muito próximos. (é o "k" do artigo)
double rand_number[1000];				// Vetor de numeros aleatorios, seguindo uma distribuicao normal
double v[3];						// Vetor da direcao de crescimento do ramo
double d_gra[3];					// Vetor do gradiente de distancia
double sobel[3][3][3];					// Filtro de Sobel
// Eixo y: i++, j++, k++
// Eixo z: k--, i++, j++
// Eixo x: j++, i++, k++ 


Node::Node () { }
Edge::Edge () { }

// Construtor de um no
Node::Node (int id, double x, double y, double z, bool grow)
{
	this->id = id;
	this->x = x;
	this->y = y;
	this->z = z;
	this->num_edges = 0;
	this->grow = grow;
	this->cell = new Cell();
	this->cell->type = 0;
	this->cell->y = new double[num_eq];
	this->cell->y_new = new double[num_eq];
	this->cell->sigma = sigma_c;
	this->cell->V_star = 0.0;
	this->edges = NULL;
	this->next = NULL;
}

Node::Node (int id, double x, double y, double z, double sigma)
{
	this->id = id;
	this->x = x;
	this->y = y;
	this->z = z;
	this->num_edges = 0;
	this->grow = false;
	this->cell = new Cell();
	this->cell->type = 0;
	this->cell->y = new double[num_eq]();
	this->cell->y_new = new double[num_eq]();
	this->cell->sigma = sigma;
	this->cell->V_star = 0.0;
	this->edges = NULL;
	this->next = NULL;
}

// Construtor de uma aresta
Edge::Edge (int id, double size)
{
	this->id = id;
	this->size = size;
	this->next = NULL;
	this->dest = NULL;
}

// Calcula o tamanho do segmento do ramo usando o vetor de *rand_number
double Graph::calculateSizeBranch ()
{
	int i;
	double number;
	// Selecionar um elemento do vetor de aleatorios
	i = rand() % 1000;
	number = rand_number[i];
	return (l_bra + number);
	//return (l_bra);
} 

// Checa se o ponto esta fora dos limites do dominio, segue a formula:
// d' = d - (d.n)n
bool Graph::checkLimits (Node *p, double x, double y, double z)
{
	double l = L/2.0;
	bool tick = true;
	// Fora dos limites no eixo x, fazer a projeção usando a normal da superfície: nx = (1,0,0)
	if (x < -l || x > l)
	{
		v[0] = v[0] - (v[0]*1);
		v[1] = v[1] - (v[1]*0);
		v[2] = v[2] - (v[2]*0);
		tick = false;
	}
	// Fora dos limites no eixo y, fazer a projeção usando a normal da superfície: ny = (0,1,0)
	else if (y < -l || y > l)
	{
		v[0] = v[0] - (v[0]*0);
		v[1] = v[1] - (v[1]*1);
		v[2] = v[2] - (v[2]*0);
		tick = false;
	}
	// Fora dos limites no eixo z, fazer a projeção usando a normal da superfície: nz = (0,0,1)
	else if (z < -l || z > l)
	{
		v[0] = v[0] - (v[0]*0);
		v[1] = v[1] - (v[1]*0);
		v[2] = v[2] - (v[2]*1);
		tick = false;
	}
	if (!tick)
	{
		cout << "[!] Point has been set out of the domain! Making a projection to avoid errors ..." << endl;
		// Calcula a posição do novo ponto que atende o domínio
		x = p->getX() + calculateSizeBranch()*v[0];
		y = p->getY() + calculateSizeBranch()*v[1];
		z = p->getZ() + calculateSizeBranch()*v[2];

		Node *ptr_new = new Node(++total_nodes,x,y,z,false);
		insertNodeGraph(ptr_new);
		insertEdgeGraph(p->getId(),ptr_new->getId());
		insertEdgeGraph(ptr_new->getId(),p->getId());
	}
	return (tick);
}

// Checa criterio de colisao com outro no na arvore. (true = Liberado) e (false = Recusado)
bool Graph::checkCollision (Node *p, double x, double y, double z)
{
	double dist_prox, dist;
	Node *ptr, *p_prox;
	ptr = getListNodes();		// Ponteiro para percorrer a lista de nos
	p_prox = NULL;			// Ponteiro para o no mais proximo ate entao
	dist_prox = 0.0;		// Distancia euclidiana do no mais proximo
	while (ptr != NULL)
	{
		dist = calcNorm(x,y,z,ptr->getX(),ptr->getY(),ptr->getZ());
		// Verifica se existe algum nó na árvore muito próximo do que irá nascer. A proximidade é avaliada por tolerância (Toler)
		if ( dist <= Toler )
		{
			if (p_prox == NULL)
			{
				p_prox = ptr;
				dist_prox = dist;
			}
			else
			{
				// Atualiza o ponto para ser o mais proximo possivel
				if (dist < dist_prox)
				{
					p_prox = ptr;
					dist_prox = dist;
				}
			}
		}
		ptr = ptr->getNext();
	}
	// Nao foi encontrado nenhum no proximo do ponto, logo liberar o crescimento
	if (p_prox == NULL)
		return (true);
	// Foi encontrado pelo menos um no proximo, sendo assim liga-se os dois por meio de uma aresta e paramos o crescimento do ramo
	else
	{
		insertEdgeGraph(p->getId(),p_prox->getId());
		insertEdgeGraph(p_prox->getId(),p->getId());
		return (false);
	}
}

// Gerar os ramos da arvore. Sendo que a direcao de crescimento ja foi calculada e esta no vetor *v
void Graph::makeBranch (Node *p)
{
	double x, y, z;
	// Ponto de crescimento
	x = p->getX() + calculateSizeBranch()*v[0];
	y = p->getY() + calculateSizeBranch()*v[1];
	z = p->getZ() + calculateSizeBranch()*v[2];
	
	//cout << "d" << endl;
	//cout << "x = " << x << endl;
	//cout << "y = " << y << endl;
	//cout << "z = " << z << endl;
	// Estamos liberados para gerar o segmento ? O ponto de crescimento esta muito proximo de alguem ou ultrapassa o dominio ?
	if (checkCollision(p,x,y,z) && checkLimits(p,x,y,z))
	{
		Node *ptr_new = new Node(++total_nodes,x,y,z,true);
		insertNodeGraph(ptr_new);
		insertEdgeGraph(p->getId(),ptr_new->getId());
		insertEdgeGraph(ptr_new->getId(),p->getId());
		cont_segments++;
		
		// Deve-se repetir 4 vezes para completar o crescimento do ramo.
		if (cont_segments <= 4)
			Sobel(ptr_new,1);
		else
			queue->Enqueue(ptr_new);
	}
}

// Calcula a direcao de crescimento pela formula: d = d_ori + w1*d_gra
void Graph::calculateDirection (double *d_gra, Node *p, double teta)
{
	double norm;

	// Captura o no anterior ao de crescimento para calcular a direcao original de crescimento: d_ori
	Node *ptr = searchNode(p->getEdges()->getId());
	v[0] = p->getX() - ptr->getX();
	v[1] = p->getY() - ptr->getY();
	v[2] = p->getZ() - ptr->getZ();
	norm = sqrt(pow(v[0],2) + pow(v[1],2) + pow(v[2],2));
	// Vetor unitario: d_ori
	v[0] = v[0]/norm;
	v[1] = v[1]/norm;
	v[2] = v[2]/norm;
	
	// Substituir na formula de direcao
	if (teta >= 0)
	{
		v[0] = v[0] + w_1*d_gra[0];
		v[1] = v[1] + w_1*d_gra[1];
		v[2] = v[2] + w_1*d_gra[2];
	}
	else
	{
		v[0] = v[0] - w_1*d_gra[0];
		v[1] = v[1] - w_1*d_gra[1];
		v[2] = v[2] - w_1*d_gra[2];
	}
	// Vetor unitário da direção crescimento
	norm = sqrt(pow(v[0],2) + pow(v[1],2) + pow(v[2],2));

	v[0] = v[0]/norm;
	v[1] = v[1]/norm;
	v[2] = v[2]/norm;

	// Criar o ramo na direção calculada, repete-se a geração do ramo 5 vezes dentro desse procedimento
	makeBranch(p);
	// Resetar o contador de ramos
	cont_segments = 1;
}

void Graph::Sobel (Node *p, int type)
{
	bool tick;
	int i, j, k;
	double x, y, z;
	double width, lenght, height;
	double norm = 0.0;
	Node *ptr = getListNodes();
	// Dimensoes do mini-cubo
	width = lenght = height = t_cube;
	// Seta as posicoes x, y, z do primeiro mini-cubo a ser testado
	x = p->getX() - ((width*3)/2);
	y = p->getY() - ((lenght*3)/2);
	z = p->getZ() - ((height*3)/2);
	// Aumenta-se o tamanho dos mini-cubo caso o gradiente resulte em nulo
	while (norm == 0.0)
	{
		width++; lenght++; height++;
		x = p->getX() - ((width*3)/2);
		y = p->getY() - ((lenght*3)/2);
		z = p->getZ() - ((height*3)/2);
		// Itera-se em relacao a todos os mini-cubos
		for (i = 0; i < 3; i++)
		{
			z = p->getZ() + ((height*3)/2);
			for (j = 0; j < 3; j++)
			{
				y = p->getY() + ((lenght*3)/2);
				for (k = 0; k < 3; k++)
				{
					ptr = getListNodes();
					tick = false;
					// Nao checar o cubo central
					if (i == 1 && j == 1 && k == 1)
						sobel[i][j][k] = 0;
					// Para os outros ...
					else
					{
						// Testar para todos os nos do grafo se algum deles esta dentro do mini-cubo atual
						while (ptr != NULL)
						{
							// Se acharmos pelo menos um no, marcamos o mini-cubo
							if (isInsideCube(ptr,x,y,z,width,lenght,height))
							{
								tick = true;
								break;
							}
							ptr = ptr->getNext();
						}
						if (tick)
							sobel[i][j][k] = 1;
						else
							sobel[i][j][k] = 0;
					}
					y -= lenght;
				}
				z -= height;
			}
			x += width;
		}
		// Convolucao
		d_gra[0] = (-sobel[0][0][0] + sobel[0][0][2] - 2*sobel[0][1][0] + 2*sobel[0][1][2] - sobel[0][2][0] + sobel[0][2][2]) +
		   (-2*sobel[1][0][0] + 2*sobel[1][0][2] - 4*sobel[1][1][0] + 4*sobel[1][1][2] - 2*sobel[1][2][0] + 2*sobel[1][2][2]) +
                   (-sobel[2][0][0] + sobel[2][0][2] - 2*sobel[2][1][0] + 2*sobel[2][1][2] - sobel[2][2][0] + sobel[2][2][2]);

		d_gra[1] = (sobel[0][0][2] + 2*sobel[1][0][2] + sobel[2][0][2] - sobel[0][2][2] - 2*sobel[1][2][2] - sobel[2][2][2]) +
		   (2*sobel[0][0][1] + 4*sobel[1][0][1] + 2*sobel[2][0][1] - 2*sobel[0][2][1] - 4*sobel[1][2][1] - 2*sobel[2][2][1]) +
                   (sobel[0][0][0] + 2*sobel[1][0][0] + sobel[2][0][0] - sobel[0][2][0] - 2*sobel[1][2][0] - sobel[2][2][0]);

		d_gra[2] = (-sobel[0][2][2] - 2*sobel[0][2][1] - sobel[0][2][0] - 2*sobel[1][2][2] - 4*sobel[1][2][1] - 2*sobel[1][2][0] - sobel[2][2][2] - 2*sobel[2][2][1] - sobel[2][2][0]) + (sobel[0][0][2] + 2*sobel[0][0][1] + sobel[0][0][0] + 2*sobel[1][0][2] + 4*sobel[1][0][1] + 2*sobel[1][0][0] + sobel[2][0][2] + 2*sobel[2][0][1] + sobel[2][0][0]);

		//cout << d_gra[0] << " " << d_gra[1] << " " << d_gra[2] << endl;
		norm = sqrt(pow(d_gra[0],2)+pow(d_gra[1],2)+pow(d_gra[2],2));
	}
	// Invertendo o valor do gradiente de distância e calculando sua direção pelo vetor unitário
	d_gra[0] = -d_gra[0]/norm;
	d_gra[1] = -d_gra[1]/norm;
	d_gra[2] = -d_gra[2]/norm;

	switch (type)
	{
		// Inicio do ramo
		case 0:
		{
			calculateDirection(d_gra,p,1);
			Sobel(p,2);
			calculateDirection(d_gra,p,-1);
			break;
		}
		// Meio do ramo
		case 1:
		{
			calculateDirection(d_gra,p,0);
			break;
		}
		// Apos crescer um lado do ramo
		default:
			break;
	}
}

void Graph::growBranch (Node *p)
{
	Sobel(p,0);
}

// Busca um no pelo seu identificador 
Node* Graph::searchNode (int id)
{
	Node *ptr = getListNodes();
	while (ptr != NULL)
	{
		if (ptr->getId() == id)
			return (ptr);
		ptr = ptr->getNext();
	}
	cout << "[-] ERROR! Node \"" << id << "\" not found!" << endl;
	return (NULL);
}

// Insere um noh no grafo
void Graph::insertNodeGraph (Node *p)
{
	if (p != NULL)
	{
		if (lastNode == NULL)
		{
			setListNodes(p);
			setLastNode(p);
		}
		else
		{
			lastNode->setNext(p);
			setLastNode(p);
		} 
	}
	else
		cout << "[-] ERROR! Pointer *p is not valid!" << endl;
}

// Insere uma aresta no grafo
void Graph::insertEdgeGraph (int id_1, int id_2)
{
	Node *ptr1, *ptr2;
	ptr1 = searchNode(id_1);
	ptr2 = searchNode(id_2);

	Edge *edge = new Edge(id_2,calcNorm(ptr1->getX(),ptr1->getY(),ptr1->getZ(),ptr2->getX(),ptr2->getY(),ptr2->getZ()));
	if (ptr1->getEdges() == NULL)
	{
		edge->setDest(ptr2);
		ptr1->setEdges(edge);
	}
	else
	{
		Edge *ptrl = ptr1->getEdges();
		while (ptrl->getNext() != NULL)
			ptrl = ptrl->getNext();
		edge->setDest(ptr2);
		ptrl->setNext(edge);
	}
	// Incrementa o contador de arestas do no origem
	ptr1->setNumEdges(ptr1->getNumEdges()+1);
	// Incrementa o contador global de arestas
	total_edges++;
}

void Graph::makeRoot (double Dx, int n_eq)
{
	Node *ptr1, *ptr2;

	cout << "\t[!] Making root" << endl;
	total_nodes = total_edges = 0;
	num_eq = n_eq;
	listNodes = NULL;
	lastNode = NULL;

	ptr1 = new Node(++total_nodes,0.0,0.0,0.0,false);
	ptr2 = new Node(++total_nodes,0.0,0.0,l_bra,true);
	queue = new Queue();
	insertNodeGraph(ptr1);
	insertNodeGraph(ptr2);
	insertEdgeGraph(1,2);
	insertEdgeGraph(2,1);
	// Enfileira no em crescimento
	queue->Enqueue(ptr2);
	//queue->printQueue();
}


// Inicializa o vetor de numeros aleatorios usado no calculo do tamanho dos ramos
void Graph::initializeRandomGenerator ()
{
	int i;
	srand(time(NULL));
	default_random_engine generator;
	normal_distribution<double> distribution(0.0,l_bra*0.4);
	for (i = 0; i < 1000; i++)
		rand_number[i] = distribution(generator);
	
}

// Construtor do grafo. Constroi uma arvore de Purkinje automaticamente pelo L-System 
Graph::Graph (double *y0, int n_eq, double Dx)
{
	int k, cont;
	Node *ptr;
	cout << "[!] Generating automatic Purkinje structure ..." << endl;
	initializeRandomGenerator();
	makeRoot(Dx,n_eq);
	for (k = 0; k < Iter; k++)
	{
		cout << "\t[!] Iteration " << k+1 << "! Growing branches ..." << endl;
		cont = queue->getInTheQueue();
		while (cont > 0)
		{
			ptr = queue->Dequeue();
			growBranch(ptr);
			cont--;
		}
	}
	setXmax(total_nodes*Dx);
	cout << "[+] Purkinje structure was build with sucess!" << endl;
	cout << "[+] There are " << total_nodes << " nodes" << endl;
	cout << "[+] There are " << total_edges << " edges" << endl;
	//printGraph();
	
}

// Constroi uma arvore de Purkinje a partir de um arquivo de entrada
Graph::Graph (FILE *file, double *y0, int n_eq, double Dx)
{
	cout << "[!] Generating Purkinje structure from file ..." << endl;
	Node *ptr;
	int V, E, i, u, v;
	double x, y, z, sigma;
	// Inicializa parametros
	total_nodes = total_edges = 0;
	num_eq = n_eq;
	listNodes = NULL;
	lastNode = NULL;

	fscanf(file,"%d %d",&V,&E);
	for (i = 0; i < V; i++)
	{
		fscanf(file,"%lf %lf %lf %lf",&x,&y,&z,&sigma);
		ptr = new Node(++total_nodes,x,y,z,sigma);
		insertNodeGraph(ptr);
	}
	for (i = 0; i < E; i++)
	{
		fscanf(file,"%d %d",&u,&v);
		insertEdgeGraph(u,v);
		insertEdgeGraph(v,u);
	}
	setXmax(total_nodes*Dx);
	cout << "[+] Purkinje structure was build with sucess!" << endl;
	cout << "[+] There are " << total_nodes << " nodes" << endl;
	cout << "[+] There are " << total_edges << " edges" << endl;
	//printGraph();
}

// Imprime o grafo na tela
void Graph::printGraph ()
{
	Node *ptr;
	Edge *ptrl;
	ptr = getListNodes();
	cout << "======================= PRINTING GRAPH ================================" << endl;
	while (ptr != NULL)
	{
		cout << "|| " << ptr->getId() << " " << ptr->getX() << " " << ptr->getY() << " " << ptr->getZ() << " " << ptr->getNumEdges() << " " << ptr->getGrow() << " ||";
		ptrl = ptr->getEdges();
		while (ptrl != NULL)
		{
			cout << " --> || " << ptrl->getId() << " " << ptrl->getSize() << " " << ptrl->getDest()->getX() << " " << ptrl->getDest()->getY() << " " << ptrl->getDest()->getZ() << " ||";
			ptrl = ptrl->getNext();
		}
		cout << endl;
		ptr = ptr->getNext();
	}
	cout << "=======================================================================" << endl;
}

// Escreve o arquivo VTK com base no grafo atual "k"
void Graph::writeVTK (int k)
{
	char filename[50];
	FILE *fileVTK;
	Node *ptr;
	Edge *ptrl;

	ptr = getListNodes();
	sprintf(filename,"./VTK/monodomain%d.vtk",k);
	fileVTK = fopen(filename,"w+");
	fprintf(fileVTK,"# vtk DataFile Version 3.0\n");
	fprintf(fileVTK,"Monodomain\n");
	fprintf(fileVTK,"ASCII\n");
	fprintf(fileVTK,"DATASET POLYDATA\n");
	// Nos
	fprintf(fileVTK,"POINTS %d float\n",total_nodes);
	while (ptr != NULL)
	{
		fprintf(fileVTK,"%e %e %e\n",ptr->getX(),ptr->getY(),ptr->getZ());
		ptr = ptr->getNext();
	}
	// Arestas
	fprintf(fileVTK,"LINES %d %d\n",total_edges,total_edges*3);
	ptr = getListNodes();
	while (ptr != NULL)
	{
		ptrl = ptr->getEdges();
		while (ptrl != NULL)
		{
			fprintf(fileVTK,"2 %d %d\n",ptr->getId()-1,ptrl->getId()-1);
			ptrl = ptrl->getNext();
		}
		ptr = ptr->getNext();
	}
	// Solucao da equacao do monodominio
	fprintf(fileVTK,"POINT_DATA %d\n",getTotalNodes());
	fprintf(fileVTK,"SCALARS vm float 1\n");
	fprintf(fileVTK,"LOOKUP_TABLE default\n");
	ptr = getListNodes();
	while (ptr != NULL)
	{
		fprintf(fileVTK,"%e\n",ptr->getCell()->y[0]);
		ptr = ptr->getNext();
	}			
	fclose(fileVTK);
}
