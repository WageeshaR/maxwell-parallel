#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <omp.h>
#include <string.h>

#include "args.h"
#include "vtk.h"
#include "data.h"
#include "setup.h"

/**
 * @brief Update the magnetic and electric fields. The magnetic fields are updated for a half-time-step. The electric fields are updated for a full time-step.
 * 
 */
void update_fields() {
	#pragma omp parallel
	{
		#pragma omp for collapse(2)
		for (int i = 0; i < Bz_size_x; i++) {
			for (int j = 0; j < Bz_size_y; j++) {
				Bz[i][j] = Bz[i][j] - (dt / dx) * (Ey[i+1][j] - Ey[i][j])
									+ (dt / dy) * (Ex[i][j+1] - Ex[i][j]);
			}
		}

		#pragma omp for collapse(2)
		for (int i = 0; i < Ex_size_x; i++) {
			for (int j = 1; j < Ex_size_y-1; j++) {
				Ex[i][j] = Ex[i][j] + (dt / (dy * eps * mu)) * (Bz[i][j] - Bz[i][j-1]);
			}
		}

		#pragma omp for collapse(2)
		for (int i = 1; i < Ey_size_x-1; i++) {
			for (int j = 0; j < Ey_size_y; j++) {
				Ey[i][j] = Ey[i][j] - (dt / (dx * eps * mu)) * (Bz[i][j] - Bz[i-1][j]);
			}
		}
	}
}

/**
 * @brief Apply boundary conditions
 * 
 */
void apply_boundary() {
	#pragma omp parallel 
	{
		#pragma omp for
		for (int i = 0; i < Ex_size_x; i++) {
			Ex[i][0] = -Ex[i][1];
			Ex[i][Ex_size_y-1] = -Ex[i][Ex_size_y-2];
		}

		#pragma omp for
		for (int j = 0; j < Ey_size_y; j++) {
			Ey[0][j] = -Ey[1][j];
			Ey[Ey_size_x-1][j] = -Ey[Ey_size_x-2][j];
		}
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

	#pragma omp parallel
	{
		#pragma omp for collapse(2) reduction(+:local_E_mag)
		for (int i = 1; i < E_size_x-1; i++) {
			for (int j = 1; j < E_size_y-1; j++) {
				E[i][j][0] = (Ex[i-1][j] + Ex[i][j]) / 2.0;
				E[i][j][1] = (Ey[i][j-1] + Ey[i][j]) / 2.0;
				//E[i][j][2] = 0.0; // in 2D we don't care about this dimension

				local_E_mag += sqrt((E[i][j][0] * E[i][j][0]) + (E[i][j][1] * E[i][j][1]));
			}
		}
		
		#pragma omp for collapse(2) reduction(+:local_B_mag)
		for (int i = 1; i < B_size_x-1; i++) {
			for (int j = 1; j < B_size_y-1; j++) {
				//B[i][j][0] = 0.0; // in 2D we don't care about these dimensions
				//B[i][j][1] = 0.0;
				B[i][j][2] = (Bz[i-1][j] + Bz[i][j] + Bz[i][j-1] + Bz[i-1][j-1]) / 4.0;

				local_B_mag += sqrt(B[i][j][2] * B[i][j][2]);
			}
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
	setup();
	omp_set_num_threads(omp_num_threads);
	total_error = 0;
	
	FILE *comp_file;
	char comp_file_name[1024];
	char *buffer = (char *) malloc(32 * sizeof(char));
	size_t bufsize = 32;
	if (comp_mode != 0) {
		sprintf(comp_file_name_base, "comp/comp_%d_%d%%s", X, Y);
		if (comp_mode == 1) {
			sprintf(comp_file_name, comp_file_name_base, "_mag.cmp");
			comp_file = fopen(comp_file_name, "r");
			if (comp_file == NULL) {
				printf("Unable to find comparison file %s\n", comp_file_name);
				exit(1);
			}
		}
	}

	printf("Running problem size %f x %f on a %d x %d grid.\n", lengthX, lengthY, X, Y);
	
	if (verbose) print_opts();
	
	double start, end;
	start = omp_get_wtime();

	allocate_arrays();
	problem_set_up();

	// start at time 0
	double t = 0.0;
	int i = 0;
	int comp_line_len = 0;
	double round_by = pow(10, 15);
	while (i < steps) {
		apply_boundary();
		update_fields();

		t += dt;
		if (comp_mode == 1)
			comp_line_len = getline(&buffer, &bufsize, comp_file);

		if (i % output_freq == 0) {
			double E_mag, B_mag;
			resolve_to_grid(&E_mag, &B_mag);
			if (comp_mode == 1) {
				char *token;
				char *ptr;
				double mags[2] = { E_mag, B_mag };
				if (comp_line_len != -1) {
					token = strtok(buffer, " ");
					int cnt = 0;
					while (token)
					{
						double value = strtod(token, &ptr);
						value = round(value * round_by) / round_by; // Rounding to account only first 15 deciman points
						double diff = abs(value - mags[cnt]);
						total_error += diff;
						cnt++;
					}
					
				}
			}
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

	if (!no_output) 
		write_result();

	free_arrays();

	exit(0);
}


