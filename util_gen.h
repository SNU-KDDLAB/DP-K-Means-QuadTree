#ifndef UTIL_GEN_H
#define UTIL_GEN_H

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))
#define ABS(a) (((a)>0)?(a):-(a))
#define SIGN(a) (((a)>0)?1:-1)

enum class Noise{ uni = 1, geo, rel };

//double geoNoise(double eps, int depth, int tot_depth, int seed = 0);
double uniformNoise(double eps, int tot_depth, int seed = 0);
double geoBudget(double eps, int depth, int tot_depth);
double uniformBudget(double eps, int tot_depth);
double relativeBudget(double eps, int depth, int tot_depth, int dims);

double findBudget(Noise noiseType, double eps, int depth, int tot_depth, int dims);

#endif