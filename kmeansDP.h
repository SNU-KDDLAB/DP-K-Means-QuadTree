
#ifndef KMEANSDP_H
#define KMEANSDP_H

#include "basic.h"

#define DBL_MAX std::numeric_limits<double>::max()
#define DBL_MIN std::numeric_limits<double>::min()

#define MAX_ITER 10

// For measuring Error
double calcSSE(const points *data, const points *center);

// Initialize Center
void initCenter(points *center, double *domain, int dim, int seed);

void kmeansIter(const points *data, points *center, double *domain = NULL, double epsilon = -1.);

// 0. Original K-means Clustering
double kmeans(const points *data, const points *initial, points *outCenter, int maxIter = MAX_ITER);

// 1. Naive Approach
double kmeansNaiveDP(const points *data, const points *initial, points *outCenter, double *domain, double epsilon, int maxIter = MAX_ITER);

// 2. EUG K-means  (CODASPY 2016, Differentially Private K-Means Clustering)
double kmeansEUGDP(const points *data, const points *initial, points *outCenter, double *domain, double epsilon);

// Proposed Alg.
// 3. Differential Privacy K-means Alg using Quad Tree 
// Quad Tree building
double kmeansDP(const points *data, const points *initial, points *outCenter,double *domain, double epsilon, int seed, int maxH, double T, int noiseType, double ratio);




#endif