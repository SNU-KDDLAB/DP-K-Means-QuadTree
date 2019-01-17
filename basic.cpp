
#include <iostream>
#include <iomanip>
#include "basic.h"

using std::cout;
using std::endl;
using std::setw;
double laplace(double mu, double b, int seed){
    static std::mt19937 gen(seed);
    std::uniform_real_distribution<double> uniDist(-0.5,0.5);
    double uni = uniDist(gen);
    double laplace = mu  - b * SIGN(uni) * log(1- 2*ABS(uni));
    return laplace;
};

double dist(double* A, double* B, int dim){
    double temp = 0.0;
    for(int i = 0; i < dim; ++i)
        temp += pow(A[i]-B[i], 2);
    return temp;
}

void printPoints(points *A, int k){
    int size = MIN(k, A->size);
    for (int i = 0; i< A->size; ++i){
        for(int j = 0; j< A->dim; ++j)
            cout<<setw(8)<<A->p[i][j]<<" ";
        cout<<endl;
    }
}

// class points
double points::distPoint(int i, int j ){
    return pow(dist(this->p[i], this->p[j], this->dim),0.5);
}

// return the index of bucket that this->p[idx] is in (domain) with M buckets
int points::includeIn(int idx, int M, double *domain) const{
    int loc = 0;
    for(int i = 0; i<this->dim; ++i){
        double a = M * double((this->p[idx][i] - domain[2*i])/ (domain[2*i+1]-domain[2*i]));
        if (int(a) == M)    a=M-1;
        loc += int(a) * pow(M, i);
    }
    return loc;
}

void points::fillZero(){
    for(int i = 0; i< size*dim; ++i){
        this->d[i]=0.0; 
    }
}
void points::copy(const points *A){
    this->size = A->size;
    this->dim = A->dim;
    memcpy(this->d, A->d, sizeof(double) * this->size * this->dim);
    //for(int i = 0; i < this->size; ++i)
    //   this->p[i] = &(this->d[i * this->dim]); 
    
}
void points::add(const points *A, int idx, int target){ // add A->p[idx] to this->p[target]
    for(int i = 0; i< dim; ++i)
        this->p[target][i] += (A->weights[idx] * A->p[idx][i]);
}






