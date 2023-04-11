#define MAX_NO_NODES 20
#define MAX_NO_BEES 30
#define EPS 0.001
#include <iostream>
#include <cstdio>
#include <unordered_set>
#include <set>
#include <vector>

using namespace std;

// input parameters: spectrum, size of a graph
double C[MAX_NO_NODES];
int n;

/* input parameters: number of bees and constructive moves per forward pass */
//int bees, NC;

// number of edges, will be calculated after input
int m;

FILE* fin;
FILE* fout;
FILE* matr_spectrum; /* input and output files */

// binary array(composed of real symmetric matrix), representing solution for each bee
double As[MAX_NO_BEES][MAX_NO_NODES * MAX_NO_NODES];
double it_bestAs[MAX_NO_NODES * MAX_NO_NODES];

// array containing spectral distances from c vector for each bee
double y[MAX_NO_BEES];
double it_bestY;

// global minimum spectral distance
double ygmin;

// graph with global minimum spectral distance
double Agmin[MAX_NO_NODES * MAX_NO_NODES];