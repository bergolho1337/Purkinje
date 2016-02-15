#include "queue.h"

QNode::QNode (QNode *next, QNode *prev, Node *p) 
{ 
	this->next = next;
	this->prev = prev;
	this->node = p;
}

// Enfileira um no em crescimento 
void Queue::Enqueue (Node *p)
{
	if (head == NULL)
	{
		QNode *ptr = new QNode(NULL,NULL,p);
		setHead(ptr);
		setLast(ptr);
	}
	else
	{
		QNode *ptr, *aux;
		aux = getLast();
		ptr = new QNode(NULL,aux,p);
		aux->setNext(ptr);
		setLast(ptr);
	}
	in_the_queue++;
}

// Verificar se a fila esta vazia
bool Queue::isEmpty ()
{
	return (in_the_queue = 0 ? true : false);
}

// Desinfileirar um no em crescimento. *** Respeita a ordenacao da fila
Node* Queue::Dequeue ()
{
	Node *ptr;
	if (head == last)
	{
		ptr = head->getNode();
		setLast(NULL);
		setHead(NULL);
	}
	else
	{
		ptr = head->getNode();
		setHead(head->getNext());
		head->setPrev(NULL);
	}
	in_the_queue--;
	return (ptr);
}

// Desinfileirar um no em crescimento. *** Escolhe um elemento aleatorio
Node* Queue::Dequeue2 (int cont)
{
	Node *ptr;
	QNode *p1, *p2;
	if (head == last)
	{
		ptr = head->getNode();
		setLast(NULL);
		setHead(NULL);
	}
	else
	{
		// Sorteia um elemento da fila
		int choosen = rand() % in_the_queue;
		p1 = getHead();
		while (choosen > 0)
		{
			p1 = p1->getNext();
			choosen--;
		}
		// Captura o no escolhido
		ptr = p1->getNode();
		// Acerta a fila
		if (p1 == head)
		{
			p2 = p1->getNext();
			p2->setPrev(NULL);
			setHead(p2);
		}
		else if (p1 == last)
		{
			p2 = p1->getPrev();
			p2->setNext(NULL);
			setLast(p2);
		}
		else
		{
			p1->getPrev()->setNext(p1->getNext());
			p1->getNext()->setPrev(p1->getPrev());
		}
	}
	in_the_queue--;
	return (ptr);
}

// Imprime os elementos na fila
void Queue::printQueue ()
{
	QNode *ptr = getHead();
	cout << "============== GROWING QUEUE =========================================" << endl;	
	cout << "head --> last" << endl;	
	cout << "In the queue = " << getInTheQueue() << endl;	
	while (ptr != NULL)
	{
		cout << ptr->getNode()->getId() << " ";
		ptr = ptr->getNext();
	}
	cout << "\n======================================================================" << endl;
	cout << endl;
}


