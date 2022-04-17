#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <string.h>

#include "vtk.h"
#include "data.h"
#include "setup.h"
#include "args.h"

/**
 * @brief Set up some default values before arguments have been loaded
 * 
 */
void set_defaults() {
	lengthX = 1.0;
	lengthY = 1.0;

	X = 4000;
	Y = 4000;
	
	T = 1.6e-9;

	set_default_base();
}

void load_originals(char *filename) {
	char *tail = strrchr(filename, '/');
	char newfile[100] = "original";
	strcat(newfile, tail);
	FILE * fp;
    char * line = NULL;
    size_t len = 0;
    ssize_t read;
	fp = fopen(newfile, "r");
	int reading = 0;
	int iter = 0;
	const char s[2] = " ";
	char *token;
	if (fp == NULL)
		exit(EXIT_FAILURE);
	while ((read = getline(&line, &len, fp)) != -1) {
		int j = iter / 101;
		int i = iter % 101;
		if (reading == 1 && strcmp(line, "VECTORS B_field float\n") != 0) {
			/* get the first token */
			token = strtok(line, s);
			O_E[i][j][0] = atof(token);
			/* walk through other tokens */
			while( token != NULL ) {
				token = strtok(NULL, s);
				O_E[i][j][1] = atof(token);
				break;	// Breaking because in O_E we only read first 2 dimensions
			}
			iter++;
		}
		if (reading == 2) {
			/* get the last char string */
			char *last = strrchr(line, ' ');
			O_B[i][j][2] = atof(last);
			iter++;
		}
		if (strcmp(line, "VECTORS E_field float\n") == 0) {
			reading = 1;
			iter = 0;
		}
		if (strcmp(line, "VECTORS B_field float\n") == 0) {
			reading = 2;
			iter = 0;
		}
	}
	fclose(fp);
}

/**
 * @brief Set up some of the values required for computation after arguments have been loaded
 * 
 */
void setup() {
	dx = lengthX / X;
	dy = lengthY / Y;

	dt = cfl * (dx > dy ? dx : dy) / c;
	
	if (steps == 0) // only set this if steps hasn't been specified
		steps = (int) (T / dt);
}

/**
 * @brief Allocate all of the arrays used for computation
 * 
 */
void allocate_arrays() {
	Ex_size_x = X; Ex_size_y = Y+1;
	Ex = alloc_2d_array(X, Y+1);
	Ey_size_x = X+1; Ey_size_y = Y;
	Ey = alloc_2d_array(X+1, Y);
	
	Bz_size_x = X; Bz_size_y = Y;
	Bz = alloc_2d_array(X, Y);
	
	E_size_x = X+1; E_size_y = Y+1; E_size_z = 3;
	E = alloc_3d_array(E_size_x, E_size_y, E_size_z);

	B_size_x = X+1; B_size_y = Y+1; B_size_z = 3;
	B = alloc_3d_array(B_size_x, B_size_y, B_size_z);
}

/**
 * @brief Free all of the arrays used for the computation
 * 
 */
void free_arrays() {
	free_2d_array(Ex);
	free_2d_array(Ey);
	free_2d_array(Bz);
	free_3d_array(E);
	free_3d_array(B);
}

/**
 * @brief Set up a guassian to curve around the centre
 * 
 */
void problem_set_up() {
	#pragma omp parallel
	{
		#pragma omp for collapse(2)
		for (int i = 0; i < Ex_size_x; i++ ) {
			for (int j = 0; j < Ex_size_y; j++) {
				double xcen = lengthX / 2.0;
				double ycen = lengthY / 2.0;
				double xcoord = (i - xcen) * dx;
				double ycoord = j * dy;
				double rx = xcen - xcoord;
				double ry = ycen - ycoord;
				double rlen = sqrt(rx*rx + ry*ry);
				double tx = (rlen == 0) ? 0 : ry / rlen;
				double mag = exp(-400.0 * (rlen - (lengthX / 4.0)) * (rlen - (lengthX / 4.0)));
				Ex[i][j] = mag * tx;
			}
		}
		#pragma omp for collapse(2)
		for (int i = 0; i < Ey_size_x; i++ ) {
			for (int j = 0; j < Ey_size_y; j++) {
				double xcen = lengthX / 2.0;
				double ycen = lengthY / 2.0;
				double xcoord = i * dx;
				double ycoord = (j - ycen) * dy;
				double rx = xcen - xcoord;
				double ry = ycen - ycoord;
				double rlen = sqrt(rx*rx + ry*ry);
				double ty = (rlen == 0) ? 0 : -rx / rlen;
				double mag = exp(-400.0 * (rlen - (lengthY / 4.0)) * (rlen - (lengthY / 4.0)));
				Ey[i][j] = mag * ty;
			}
		}
	}
}

void compare_line(int len, char **buf, double mags[]) {
	char *token;
	char *ptr;
	if (len != -1) {
		token = strtok(*buf, " ");
		int cnt = 0;
		while (token)
		{
			double value = strtod(token, &ptr);
			// value = round(value * round_by) / round_by; // Rounding to account only first 15 decimal points
			printf("%.15f v %.15f\n", value, mags[cnt]);
			double diff = fabs(value - mags[cnt]);
			total_error += diff;
			token = strtok(NULL, " ");
			cnt++;
		}
		
	}
}
