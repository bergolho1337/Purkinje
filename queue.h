#ifndef QUEUE_H
#define QUEUE_H

#include <iostream>
#include "purkinje.h"

using namespace std;

class QNode
{
private:
	QNode *next;				// Ponteiro para o proximo elemento da fila
	QNode *prev;				// Ponteiro para o elemento anterior da fila
	Node *node;				// Ponteiro para o no em crescimento da arvore
public:
	QNode (QNode *next, QNode *prev, Node *p);
// Inline
	void setNext (QNode *p) { next = p; }
	void setPrev (QNode *p) { prev = p; }
	void setNode (Node *p) { node = p; }
	QNode* getNext () { return (next); }
	QNode* getPrev () { return (prev); }
	Node* getNode () { return (node); }
};

class Queue
{
private:
	QNode *head;				// Ponteiro para a cabeca da fila
	QNode *last;				// Ponteiro para o ultimo da fila
	int in_the_queue;			// Contador de nos em crescimento na fila
public:
	Queue () { this->head = NULL; this->last = NULL; in_the_queue = 0; }
// Inline
	void setHead (QNode *p) { head = p; }
	void setLast (QNode *p) { last = p; }
	QNode* getHead () { return (head); }
	QNode* getLast () { return (last); }
	int getInTheQueue () { return (in_the_queue); }
	
	void Enqueue (Node *p);
	Node* Dequeue ();
	Node* Dequeue2 (int cont);
	bool isEmpty ();
	void printQueue ();
};

#endif
