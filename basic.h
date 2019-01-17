#ifndef BASIC_KMEANS_H
#define BASIC_KMEANS_H

#include <random> // For uniform random -> laplace
#include <cmath> // for logarithm
#include <string.h>

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))
#define ABS(a) (((a)>0)?(a):-(a))
#define SIGN(a) (((a)>0)?1:-1)

#define MAXDEPTH 5

class points;
class Node;

// For noise insertion
double laplace(double mu, double b, int seed = 0);
double dist(double* A, double* B, int dim);
void printPoints(points *A, int k = 10);

class points{
public:
    int size;
    int dim;
    double *weights;
    double *d;
    double **p;
    void fillZero();
    double distPoint(int i, int j);
    void copy(const points *A);
    void add(const points *A, int idx, int target);
    void assign(double* coord, int idx); // assign coord into data->p[idx]
    int includeIn(int idx, int M, double *domain) const;
    points(){}
    points(int length, int dimension){
        this->size = length;
        this->dim = dimension;
        this->d = new double[this->size * this->dim];
        this->p = new double*[this->size];
        this->weights = new double[this->size];
        for(int i = 0; i < this->size; ++i){
            this->p[i] = &(this->d[i * this->dim]);
            this->weights[i] = 1.0; 
        }
    }
    ~points(){
        delete[] this->p;
        delete[] this->d;
        delete[] this->weights;
    }
};


#endif