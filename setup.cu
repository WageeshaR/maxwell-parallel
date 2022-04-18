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
	specifics.lengthX = 1.0;
	specifics.lengthY = 1.0;

	specifics.X = 4000;
	specifics.Y = 4000;
	
	T = 1.6e-9;

	set_default_base();
}

/**
 * @brief Set up some of the values required for computation after arguments have been loaded
 * 
 */
void setup() {
	specifics.dx = specifics.lengthX / specifics.X;
	specifics.dy = specifics.lengthY / specifics.Y;

	specifics.dt = constants.cfl * (specifics.dx > specifics.dy ? specifics.dx : specifics.dy) / constants.c;
	
	if (steps == 0) // only set this if steps hasn't been specified
		steps = (int) (T / specifics.dt);
}

/**
 * @brief Allocate all of the arrays used for computation
 * 
 */
void allocate_arrays() {
	arrays.Ex_size_x = specifics.X; arrays.Ex_size_y = specifics.Y+1;
	alloc_2d_cuda_array(specifics.X, specifics.Y+1, &arrays.Ex, &arrays.ex_pitch);
	arrays.Ey_size_x = specifics.X+1; arrays.Ey_size_y = specifics.Y;
	alloc_2d_cuda_array(specifics.X+1, specifics.Y, &arrays.Ey, &arrays.ey_pitch);
	
	arrays.Bz_size_x = specifics.X; arrays.Bz_size_y = specifics.Y;
	alloc_2d_cuda_array(specifics.X, specifics.Y, &arrays.Bz, &arrays.bz_pitch);
	
	arrays.E_size_x = specifics.X+1; arrays.E_size_y = specifics.Y+1; arrays.E_size_z = 3;
	alloc_3d_cuda_array(arrays.E_size_x, arrays.E_size_y, arrays.E_size_z, &arrays.E, &arrays.e_pitch);
	host_E = alloc_3d_array(arrays.E_size_x, arrays.E_size_y, arrays.E_size_z);

	arrays.B_size_x = specifics.X+1; arrays.B_size_y = specifics.Y+1; arrays.B_size_z = 3;
	alloc_3d_cuda_array(arrays.B_size_x, arrays.B_size_y, arrays.B_size_z, &arrays.B, &arrays.b_pitch);
	host_B = alloc_3d_array(arrays.B_size_x, arrays.B_size_y, arrays.B_size_z);
}

/**
 * @brief Free all of the arrays used for the computation
 * 
 */
void free_arrays() {
	free_2d_cuda_array(arrays.Ex);
	free_2d_cuda_array(arrays.Ey);
	free_2d_cuda_array(arrays.Bz);
	free_3d_cuda_array(arrays.E);
	free_3d_cuda_array(arrays.B);
	free_3d_array(host_E);
	free_3d_array(host_B);
}

/**
 * @brief Set up a guassian to curve around the centre
 * 
 */
__global__ void problem_set_up(Arrays arrays, Specifics specifics) {
	double xcen = specifics.lengthX / 2.0;
	double ycen = specifics.lengthY / 2.0;


    for (int i = 0; i < arrays.Ex_size_x; i++ ) {
        for (int j = 0; j < arrays.Ex_size_y; j++) {
			double xcoord = (i - xcen) * specifics.dx;
			double ycoord = j * specifics.dy;
			double rx = xcen - xcoord;
			double ry = ycen - ycoord;
			double rlen = sqrt(rx*rx + ry*ry);
			double tx = (rlen == 0) ? 0 : ry / rlen;
			double mag = exp(-400.0 * (rlen - (specifics.lengthX / 4.0)) * (rlen - (specifics.lengthY / 4.0)));
			arrays.Ex[i * arrays.ex_pitch + j] = mag * tx;
		}
	}
    for (int i = 0; i < arrays.Ey_size_x; i++ ) {
        for (int j = 0; j < arrays.Ey_size_y; j++) {
            double xcoord = i * specifics.dx;
            double ycoord = (j - ycen) * specifics.dy;
            double rx = xcen - xcoord;
            double ry = ycen - ycoord;
            double rlen = sqrt(rx*rx + ry*ry);
            double ty = (rlen == 0) ? 0 : -rx / rlen;
			double mag = exp(-400.0 * (rlen - (specifics.lengthY / 4.0)) * (rlen - (specifics.lengthY / 4.0)));
            arrays.Ey[i*arrays.ey_pitch + j] = mag * ty;
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
			double e_value = strtod(token, &ptr);
			total_error += fabs(e_value - mags[cnt]);
			token = strtok(NULL, " ");
			cnt++;
		}
		
	}
}