////////////////////////////////
//   Differential K-means     //
//  Written by HJ. Koo @ SNU  //
// Email: hjkoo@kdd.snu.ac.kr //
////////////////////////////////

#include <iostream> //cout, endl
#include <fstream> // file read
#include <limits> // double max
#include <stdlib.h> // Random
#include <assert.h> // Assert
#include <iomanip>

#include "basic.h"
#include "kmeansDP.h"

//#define DEBUG 1
using std::setw;
using std::cout; using std::endl;
using std::ifstream;

// Prototypes
void readData(const char* path, points *&data, double *&domain);
void performExp(points *data, double* domain, int k, double epsilon, int seed , double* result, int maxH, double T, double ratio);

int main(int argc, char* args[]){
	// Load data
	points *data;
	double* domain;
	if(argc<2)	readData("test", data, domain);
	else	readData(args[1], data, domain);
	
	// Check the domain of the data
	for(int i = 0; i< data->dim; ++i)
		cout<<i<<"-th dim:["<<domain[2*i]<<", "<<domain[2*i+1]<<"]"<<endl;
	cout<<endl;

	int sizeExp = 5+3;
	double *result = new double[sizeExp];
	
	// Experiment parameter
	int seed = 1;
	int n_cluster_min = 5;
	int n_cluster_max = 6;
	int n_cluster_step = 1;
	double epsilon_min = 0.05;
	double epsilon_max = 0.06;

	// Tree
	double T_min = data->size/1000.;
	double T_max = T_min+1;

	// ratio for tree & noise
	double ratio_min = 0.05;
	double ratio_max = 1.0;
	//int minH = 4*log(data->size)/(log(t)*data->dim);
	int minH = log(data->size)/(data->dim)-1;
	//int minH = 3;
	int maxH = minH+1;


	int iterExp = 30;	
	// Perform Experiment	
	cout<<"kmeans vs NaiveDP vs EUGDP vs Proposed (uni, rel)"<<endl;
	for(double ratio = ratio_min; ratio<ratio_max; ratio+=0.05){
		for(int t =T_min; t < T_max; t*= 2){
			for (int h = minH; h<maxH; ++h){
				for (double epsilon = epsilon_min; epsilon < epsilon_max; epsilon *= 10){
					for(int k= n_cluster_min; k < n_cluster_max; k+= n_cluster_step){

						double per1, per2;
						double avg1 = 0.0, avg2 = 0.0;
						for(int i = 0; i< iterExp;++i){
							// Initialize the result
							for(int temp = 0;temp<sizeExp;++temp)	result[temp]=0.;
							
							performExp(data, domain, k, epsilon, seed+i, result, h, t, ratio);
							per1 = ((result[2] -result[3])/result[0] * 100.); //uni
							per2 = ((result[2] -result[4])/result[0] * 100.); //rel
							avg1 += per1;
							avg2 += per2;
							cout<<per1<<", "<<per2<<" ";
						}
						avg1 /= iterExp;
						avg2 /= iterExp;

						// Output the result
						cout<<"k = "<<k<<" ratio = "<<ratio<<" eps = "<<epsilon<<" maxH = "<< h <<" T = "<<setw(3)<<t<<"  ";
						/*
						for(int i = 0;i<sizeExp-3;++i)
							cout<< setw(12) << result[i]/iterExp <<" ";
						cout<<" EUG "<<setw(5)<<result[5]/iterExp <<" K_ "<<setw(5)<< result[6]/iterExp;
						double per1 = ((result[2] -result[3])/result[0] * 100.); //uni
						int per2 = int((result[2] -result[4])/result[0] * 100.); //rel
						*/
						cout<<"   uni_total: "<<setw(3)<<avg1<<"  "<<"rel_total: "<<setw(3)<<avg2<<"  "<<endl;

					}
				}
			}
		}
	}

	// de-allocation
	delete[] result;
	delete[] domain;
	delete data;

	cout<<"    Done!   "<<endl;
	return 0;
}

void performExp(points *data, double* domain, int k, double epsilon, int seed, double* result, int maxH, double T, double ratio){

	points *initial = new points(k, data->dim);
	initCenter(initial, domain, data->dim, seed);
#ifdef  DEBUG
	cout<<"Initial Center"<<endl;
	printPoints(initial);
#endif
	// For k-means result
	points *center1 = new points(k, data->dim);

	// Perform k-means clusterings
	double SSE1 = kmeans(data, initial, center1);
#ifdef  DEBUG
	cout<<"Kmeans 1"<<endl;
	printPoints(center1);	
#endif

	double SSE2 = kmeansNaiveDP(data, initial, center1, domain, epsilon);
#ifdef  DEBUG
	cout<<"Kmeans 2"<<endl;
	printPoints(center1);
#endif

	int n_bucket = kmeansEUGDP(data, initial, center1, domain, epsilon);
	double SSE3 = calcSSE(data, center1);
#ifdef  DEBUG
	cout<<"Kmeans 3"<<endl;
	printPoints(center1);	
#endif

	int n_bucket1 = kmeansDP(data, initial, center1, domain, epsilon, seed, maxH, T, 1, ratio);
	double SSE4 = calcSSE(data, center1);
	#ifdef  DEBUG
	cout<<"Kmeans 4"<<endl;
	printPoints(center1);
#endif
	int n_bucket2 = kmeansDP(data, initial, center1, domain, epsilon, seed, maxH, T, 2, ratio);
	double SSE5 = calcSSE(data, center1);
#ifdef  DEBUG
	cout<<"Kmeans 5"<<endl;
	printPoints(center1);
#endif

	// Evaluation
	//cout<<"SSE @ k = "<<k<<" "<<SSE1<<", "<<SSE2<<", "<<SSE3<<", "<<SSE4<<endl;		
	result[0] += SSE1/data->size;
	result[1] += SSE2/data->size;
	result[2] += SSE3/data->size;
	result[3] += SSE4/data->size;
	result[4] += SSE5/data->size;
	result[5] += n_bucket;
	result[6] += n_bucket1;
	result[7] += n_bucket2;
	delete initial;
	delete center1;
}

void readData(const char* path, points *&data, double *&domain){
	cout<<"Reading the data! :  "<<path<<endl;
	ifstream inputFile(path);
	if(inputFile.is_open()){
		int size, dim, count = 0;
		inputFile>>size;
		inputFile>>dim;
		data = new points(size, dim);
		domain = new double[dim*2];
		for(int i = 0; i<dim;++i){
			domain[2*i]=DBL_MAX;
			domain[2*i+1]=DBL_MIN;
		}
		// Data load
		while(inputFile>>data->d[count]){
			domain[2*(count%dim)] = MIN(data->d[count], domain[2*(count%dim)]);
			domain[2*(count%dim)+1] = MAX(data->d[count], domain[2*(count%dim)+1]);
			++count;
		}
		double area = 1.0;
		for(int i = 0; i<dim;++i)
			area *= (domain[2*i+1]-domain[2*i]);
		cout<<"Size, Dim = "<<size<<", "<<dim<<endl;
		cout<<"Density = "<<size/area<<endl;
		cout<<endl;
		inputFile.close();
	}
	else
		cout<<"Unable to open file: "<<path<<endl;
	return;
}