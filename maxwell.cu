#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

#include "args.h"
#include "vtk.h"
#include "data.h"
#include "setup.h"

/**
 * @brief Update the magnetic and electric fields. The magnetic fields are updated for a half-time-step. The electric fields are updated for a full time-step.
 * 
 */
__global__ void update_fields(Constants constants, Specifics specifics, Arrays arrays) {
	// for (int i = 0; i < arrays.Bz_size_x; i++) {
	// 	for (int j = 0; j < arrays.Bz_size_y; j++) {
	// 		arrays.Bz[i][j] = arrays.Bz[i][j] - (specifics.dt / specifics.dx) * (arrays.Ey[i+1][j] - arrays.Ey[i][j])
	// 			                + (specifics.dt / specifics.dy) * (arrays.Ex[i][j+1] - arrays.Ex[i][j]);
	// 	}
	// }

	// for (int i = 0; i < arrays.Ex_size_x; i++) {
	// 	for (int j = 1; j < arrays.Ex_size_y-1; j++) {
	// 		arrays.Ex[i][j] = arrays.Ex[i][j] + (specifics.dt / (specifics.dy * constants.eps * constants.mu)) * (arrays.Bz[i][j] - arrays.Bz[i][j-1]);
	// 	}
	// }

	// for (int i = 1; i < arrays.Ey_size_x-1; i++) {
	// 	for (int j = 0; j < arrays.Ey_size_y; j++) {
	// 		arrays.Ey[i][j] = arrays.Ey[i][j] - (specifics.dt / (specifics.dx * constants.eps * constants.mu)) * (arrays.Bz[i][j] - arrays.Bz[i-1][j]);
	// 	}
	// }
}

/**
 * @brief Apply boundary conditions
 * 
 */
__global__ void apply_boundary(Arrays arrays) {
	// for (int i = 0; i < arrays.Ex_size_x; i++) {
	// 	arrays.Ex[i][0] = -arrays.Ex[i][1];
	// 	arrays.Ex[i][arrays.Ex_size_y-1] = -arrays.Ex[i][arrays.Ex_size_y-2];
	// }

	// for (int j = 0; j < arrays.Ey_size_y; j++) {
	// 	arrays.Ey[0][j] = -arrays.Ey[1][j];
	// 	arrays.Ey[arrays.Ey_size_x-1][j] = -arrays.Ey[arrays.Ey_size_x-2][j];
	// }
}

/**
 * @brief Resolve the Ex, Ey and Bz fields to grid points and sum the magnitudes for output
 * 
 * @param E_mag The returned total magnitude of the Electric field (E)
 * @param B_mag The returned total magnitude of the Magnetic field (B) 
 */
__global__ void resolve_to_grid(double *E_mag, double *B_mag, Arrays arrays) {
	// *E_mag = 0.0;
	// *B_mag = 0.0;

	// for (int i = 1; i < arrays.E_size_x-1; i++) {
	// 	for (int j = 1; j < arrays.E_size_y-1; j++) {
	// 		arrays.E[i][j][0] = (arrays.Ex[i-1][j] + arrays.Ex[i][j]) / 2.0;
	// 		arrays.E[i][j][1] = (arrays.Ey[i][j-1] + arrays.Ey[i][j]) / 2.0;
	// 		//E[i][j][2] = 0.0; // in 2D we don't care about this dimension

	// 		*E_mag += sqrt((arrays.E[i][j][0] * arrays.E[i][j][0]) + (arrays.E[i][j][1] * arrays.E[i][j][1]));
	// 	}
	// }
	
	// for (int i = 1; i < arrays.B_size_x-1; i++) {
	// 	for (int j = 1; j < arrays.B_size_y-1; j++) {
	// 		//B[i][j][0] = 0.0; // in 2D we don't care about these dimensions
	// 		//B[i][j][1] = 0.0;
	// 		arrays.B[i][j][2] = (arrays.Bz[i-1][j] + arrays.Bz[i][j] + arrays.Bz[i][j-1] + arrays.Bz[i-1][j-1]) / 4.0;

	// 		*B_mag += sqrt(arrays.B[i][j][2] * arrays.B[i][j][2]);
	// 	}
	// }
}

/**
 * @brief The main routine that sets up the problem and executes the timestepping routines
 * 
 * @param argc The number of arguments passed to the program
 * @param argv An array of the arguments passed to the program
 * @return int The return value of the application
 */
int main(int argc, char *argv[]) {
	set_defaults();
	parse_args(argc, argv);
	setup();

	printf("Running problem size %f x %f on a %d x %d grid.\n", specifics.lengthX, specifics.lengthY, specifics.X, specifics.Y);
	
	if (verbose) print_opts();
	
	allocate_arrays();

	problem_set_up<<<1,1>>>(arrays, specifics);
	cudaDeviceSynchronize();

	// start at time 0
	double t = 0.0;
	int i = 0;
	while (i < steps) {
		// apply_boundary(arrays);
		// update_fields(constants, specifics, arrays);

		t += specifics.dt;

		// if (i % output_freq == 0) {
		// 	double E_mag, B_mag;
		// 	resolve_to_grid(&E_mag, &B_mag, arrays);
		// 	printf("Step %8d, Time: %14.8e (dt: %14.8e), E magnitude: %14.8e, B magnitude: %14.8e\n", i, t, specifics.dt, E_mag, B_mag);

		// 	if ((!no_output) && (enable_checkpoints))
		// 		write_checkpoint(i);
		// }

		i++;
	}

	// double E_mag, B_mag;
	// resolve_to_grid(&E_mag, &B_mag, arrays);

	// printf("Step %8d, Time: %14.8e (dt: %14.8e), E magnitude: %14.8e, B magnitude: %14.8e\n", i, t, specifics.dt, E_mag, B_mag);
	printf("Simulation complete.\n");

	// if (!no_output) 
	// 	write_result();

	// free_arrays();

	exit(0);
}


