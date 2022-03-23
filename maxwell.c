#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <omp.h>

#include "args.h"
#include "vtk.h"
#include "data.h"
#include "setup.h"

/**
 * @brief Update the magnetic and electric fields. The magnetic fields are updated for a half-time-step. The electric fields are updated for a full time-step.
 * 
 */
void update_fields() {
	#pragma omp parallel for schedule(static) collapse(2)
	for (int i = 0; i < Bz_size_x; i++) {
		for (int j = 0; j < Bz_size_y; j++) {
			Bz[i][j] = Bz[i][j] - (dt / dx) * (Ey[i+1][j] - Ey[i][j])
				                + (dt / dy) * (Ex[i][j+1] - Ex[i][j]);
		}
	}

	#pragma omp parallel for schedule(static) collapse(2)
	for (int i = 0; i < Ex_size_x; i++) {
		for (int j = 1; j < Ex_size_y-1; j++) {
			Ex[i][j] = Ex[i][j] + (dt / (dy * eps * mu)) * (Bz[i][j] - Bz[i][j-1]);
		}
	}

	#pragma omp parallel for schedule(static) collapse(2)
	for (int i = 1; i < Ey_size_x-1; i++) {
		for (int j = 0; j < Ey_size_y; j++) {
			Ey[i][j] = Ey[i][j] - (dt / (dx * eps * mu)) * (Bz[i][j] - Bz[i-1][j]);
		}
	}
}

/**
 * @brief Apply boundary conditions
 * 
 */
void apply_boundary() {
	#pragma omp parallel for schedule(static)
	for (int i = 0; i < Ex_size_x; i++) {
		Ex[i][0] = -Ex[i][1];
		Ex[i][Ex_size_y-1] = -Ex[i][Ex_size_y-2];
	}

	#pragma omp parallel for schedule(static)
	for (int j = 0; j < Ey_size_y; j++) {
		Ey[0][j] = -Ey[1][j];
		Ey[Ey_size_x-1][j] = -Ey[Ey_size_x-2][j];
	}
}

/**
 * @brief Resolve the Ex, Ey and Bz fields to grid points and sum the magnitudes for output
 * 
 * @param E_mag The returned total magnitude of the Electric field (E)
 * @param B_mag The returned total magnitude of the Magnetic field (B) 
 */
void resolve_to_grid(double *E_mag, double *B_mag) {
	double local_E_mag = 0.0;
	double local_B_mag = 0.0;

	#pragma omp parallel for schedule(static) collapse(2) reduction(+:local_E_mag)
	for (int i = 1; i < E_size_x-1; i++) {
		for (int j = 1; j < E_size_y-1; j++) {
			E[i][j][0] = (Ex[i-1][j] + Ex[i][j]) / 2.0;
			E[i][j][1] = (Ey[i][j-1] + Ey[i][j]) / 2.0;
			//E[i][j][2] = 0.0; // in 2D we don't care about this dimension

			local_E_mag += sqrt((E[i][j][0] * E[i][j][0]) + (E[i][j][1] * E[i][j][1]));
		}
	}
	
	#pragma omp parallel for schedule(static) collapse(2) reduction(+:local_B_mag)
	for (int i = 1; i < B_size_x-1; i++) {
		for (int j = 1; j < B_size_y-1; j++) {
			//B[i][j][0] = 0.0; // in 2D we don't care about these dimensions
			//B[i][j][1] = 0.0;
			B[i][j][2] = (Bz[i-1][j] + Bz[i][j] + Bz[i][j-1] + Bz[i-1][j-1]) / 4.0;

			local_B_mag += sqrt(B[i][j][2] * B[i][j][2]);
		}
	}

	*E_mag = local_E_mag;
	*B_mag = local_B_mag;
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
	if (enable_comparison == 1) {
		if (X != 100 || Y != 100) {
			printf("To enable comparison X and Y must be 100\n");
			printf("Exiting with status 1\n");
			exit(1);
		}
	}
	setup();
	omp_set_num_threads(omp_num_threads);
	total_error = 0;

	printf("Running problem size %f x %f on a %d x %d grid.\n", lengthX, lengthY, X, Y);
	
	if (verbose) print_opts();
	
	allocate_arrays();

	problem_set_up();

	// start at time 0
	double t = 0.0;
	int i = 0;
	double start, end;
	start = omp_get_wtime();
	while (i < steps) {
		apply_boundary();
		update_fields();

		t += dt;

		if (i % output_freq == 0) {
			double E_mag, B_mag;
			resolve_to_grid(&E_mag, &B_mag);
			printf("Step %8d, Time: %14.8e (dt: %14.8e), E magnitude: %14.8e, B magnitude: %14.8e\n", i, t, dt, E_mag, B_mag);

			if ((!no_output) && (enable_checkpoints))
				write_checkpoint(i);
		}

		i++;
	}
	double E_mag, B_mag;
	resolve_to_grid(&E_mag, &B_mag);

	printf("Step %8d, Time: %14.8e (dt: %14.8e), E magnitude: %14.8e, B magnitude: %14.8e\n", i, t, dt, E_mag, B_mag);
	printf("Simulation complete.\n");
	end = omp_get_wtime();
	printf("Elapsed wall clock time is %fs\n", end - start);
	if (enable_comparison == 1)
		printf("Total error is %.12e\n", total_error);

	if (!no_output) 
		write_result();

	free_arrays();

	exit(0);
}


