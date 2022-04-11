#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

#include "vtk.h"
#include "data.h"

struct Constants constants = {
	.c = 299792458,
	.mu = 4.0 * M_PI * 1.0e-7,
	.eps = 1.0 / (constants.c * constants.c * constants.mu),
	.cfl = 0.6363961031
};

struct Specifics specifics;
struct Arrays arrays;
struct Cuda cuda_consts;

double *** host_E;
double *** host_B;

// Time to run for / or number of steps
double T = 0.0001;
int steps = 0;

/**
 * @brief Allocate a 2D CUDA array that is addressable using square brackets
 * 
 * @param m The first dimension of the array
 * @param n The second dimension of the array
 * @return double** A 2D array
 */
void alloc_2d_cuda_array(int m, int n, double **array, size_t *pitch) {
  	cudaMallocPitch((void **)array, pitch, n*sizeof(double), m);
	*pitch = (*pitch) / sizeof(double);
	// double *tmp = (double *) calloc((*pitch)*m, sizeof(double));
	// cudaMemcpy(*array, tmp, (*pitch)*m*sizeof(double), cudaMemcpyHostToDevice);
}

/**
 * @brief Free a 2D CUDA array
 * 
 * @param array The 2D array to free
 */
void free_2d_cuda_array(double *array) {
	cudaFree(array);
}

/**
 * @brief Allocate a 3D CUDA array that is addressable using square brackets
 * 
 * @param m The first dimension of the array
 * @param n The second dimension of the array
 * @param o The third dimension of the array
 * @return double*** A 3D array
 */
void alloc_3d_cuda_array(int m, int n, int o, double **array, size_t *pitch) {
	cudaMallocPitch((void **)array, pitch, n*o*sizeof(double), m);
	*pitch = (*pitch) / sizeof(double);
	printf("Pitch is %ld\n", *pitch);
	double *tmp = (double *) calloc((*pitch)*m, sizeof(double));
	cudaMemcpy(*array, tmp, (*pitch)*m*sizeof(double), cudaMemcpyHostToDevice);
}

/**
 * @brief Free a 3D CUDA array
 * 
 * @param array The 3D array to free
 */
void free_3d_cuda_array(double* array) {
	cudaFree(array);
}

/**
 * @brief Allocate a 3D C array that is addressable using square brackets
 * 
 * @param m The first dimension of the array
 * @param n The second dimension of the array
 * @param o The third dimension of the array
 * @return double*** A 3D array
 */
double ***alloc_3d_array(int m, int n, int o) {
	double ***x;
	x = (double***) malloc(m*sizeof(double **));
	x[0] = (double **) malloc(m*n*sizeof(double *));
	x[0][0] = (double *) calloc(m*n*o,sizeof(double));
	for (int i = 1; i < m; i++) {
		x[i] = &x[0][i*n];
	}
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			if (i == 0 && j == 0) continue;
			x[i][j] = &x[0][0][i*n*o + j*o];
		}
	}
	return x;
}

/**
 * @brief Free a 3D C array
 * 
 * @param array The 3D array to free
 */
void free_3d_array(double*** array) {
	free(array[0][0]);
	free(array[0]);
	free(array);
}
