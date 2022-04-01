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
	int tidx = blockIdx.x*blockDim.x + threadIdx.x;
	int tidy = blockIdx.y*blockDim.y + threadIdx.y;
	int cpx_ex = arrays.Ex_size_x / (gridDim.x * blockDim.x); // Cells per cuda thread in X direction (keeping this same for all fields)
	int cpy_ey = arrays.Ey_size_y / (gridDim.y * blockDim.y); // Cells per cuda thread in Y direction

	for (int i = cpx_ex * tidx; i < cpx_ex * (tidx + 1); i++) {
		for (int j = cpy_ey * tidy; j < cpy_ey * (tidy + 1); j++) {
			arrays.Bz[i * arrays.Bz_size_y + j] = arrays.Bz[i * arrays.Bz_size_y + j] - (specifics.dt / specifics.dx) * (arrays.Ey[(i+1) * arrays.Ey_size_y + j] - arrays.Ey[i * arrays.Ey_size_y + j])
				                					+ (specifics.dt / specifics.dy) * (arrays.Ex[i * arrays.Ex_size_y + j + 1] - arrays.Ex[i * arrays.Ex_size_y + j]);
			// printf("%14.8e\n", arrays.Bz[i * arrays.Bz_size_y + j]);
		}
	}

	for (int i = cpx_ex * tidx; i < cpx_ex * (tidx + 1); i++) {
		for (int j = cpy_ey * tidy; j < cpy_ey * (tidy + 1); j++) {
			if (tidy == 0 && j == 0)
				continue;
			arrays.Ex[i * arrays.Ex_size_y + j] = arrays.Ex[i * arrays.Ex_size_y + j]
													+ (specifics.dt / (specifics.dy * constants.eps * constants.mu)) * (arrays.Bz[i * arrays.Bz_size_y + j] - arrays.Bz[i * arrays.Bz_size_y + j - 1]);
		}
	}

	for (int i = cpx_ex * tidx; i < cpx_ex * (tidx + 1); i++) {
		for (int j = cpy_ey * tidy; j < cpy_ey * (tidy + 1); j++) {
			if (tidx == 0 && i == 0)
				continue;
			arrays.Ey[i * arrays.Ey_size_y + j] = arrays.Ey[i * arrays.Ey_size_y + j]
													- (specifics.dt / (specifics.dx * constants.eps * constants.mu)) * (arrays.Bz[i * arrays.Bz_size_y + j] - arrays.Bz[(i - 1) * arrays.Bz_size_y + j]);
		}
	}
}

/**
 * @brief Apply boundary conditions
 * 
 */
__global__ void apply_boundary(Arrays arrays) {
	int tidx = blockIdx.x*blockDim.x + threadIdx.x;
	int tidy = blockIdx.y*blockDim.y + threadIdx.y;
	int cpx_ex = arrays.Ex_size_x / (gridDim.x * blockDim.x); // Cells per cuda thread in X direction
	int cpy_ey = arrays.Ey_size_y / (gridDim.y * blockDim.y); // Cells per cuda thread in Y direction
	
	for (int i = tidx * cpx_ex; i < cpx_ex * (tidx + 1); i++) {
		if (tidy == 0) {
			arrays.Ex[i * arrays.Ex_size_y] = -arrays.Ex[i * arrays.Ex_size_y + 1];
			// printf("%14.8e\n", -arrays.Ex[i * arrays.Ex_size_y + 1]);
		}
		if (tidy == gridDim.y * blockDim.y - 1)
			arrays.Ex[i * arrays.Ex_size_y + arrays.Ex_size_y - 1] = -arrays.Ex[i * arrays.Ex_size_y + arrays.Ex_size_y - 2];
	}

	for (int j = tidy * cpy_ey; j < cpy_ey * (tidy + 1); j++) {
		if (tidx == 0)
			arrays.Ey[j] = -arrays.Ey[arrays.Ey_size_y + j];
		if (tidx == gridDim.x * blockDim.x - 1)
			arrays.Ey[arrays.Ey_size_y * (arrays.Ey_size_x - 1) + j] = -arrays.Ey[arrays.Ey_size_y * (arrays.Ey_size_x - 2) + j];
	}
}

/**
 * @brief Resolve the Ex, Ey and Bz fields to grid points and sum the magnitudes for output
 * 
 * @param E_mag The returned total magnitude of the Electric field (E)
 * @param B_mag The returned total magnitude of the Magnetic field (B) 
 */
__global__ void resolve_to_grid(double *E_mag, double *B_mag, Arrays arrays) {
	int tidx = blockIdx.x*blockDim.x + threadIdx.x;
	int tidy = blockIdx.y*blockDim.y + threadIdx.y;
	int cpx_x = (arrays.E_size_x - 1) / (gridDim.x * blockDim.x); // Cells per cuda thread in X direction
	int cpy_y = (arrays.E_size_y - 1) / (gridDim.y * blockDim.y); // Cells per cuda thread in Y direction
	E_mag[tidy * gridDim.x * blockDim.x + tidx] = 0.0;
	B_mag[tidy * gridDim.x * blockDim.x + tidx] = 0.0;

	for (int i = cpx_x * tidx; i < cpx_x * (tidx + 1); i++) {
		for (int j = cpy_y * tidy; j < cpy_y * (tidy + 1); j++) {
			if ((tidx == 0 && i == 0) || (tidy == 0 && j == 0))
				continue;
			arrays.E[(i * arrays.E_size_y + j) * arrays.E_size_z] = (arrays.Ex[(i-1) * arrays.Ex_size_y + j] + arrays.Ex[i * arrays.Ex_size_y + j]) / 2.0;
			arrays.E[(i * arrays.E_size_y + j) * arrays.E_size_z + 1] = (arrays.Ey[i * arrays.Ey_size_y + j - 1] + arrays.Ey[i * arrays.Ey_size_y + j]) / 2.0;

			E_mag[tidy * gridDim.x * blockDim.x + tidx] += sqrt((arrays.E[(i * arrays.E_size_y + j) * arrays.E_size_z] * arrays.E[(i * arrays.E_size_y + j) * arrays.E_size_z])
																+ (arrays.E[(i * arrays.E_size_y + j) * arrays.E_size_z + 1] * arrays.E[(i * arrays.E_size_y + j) * arrays.E_size_z + 1]));
		}
	}
	
	for (int i = cpx_x * tidx; i < cpx_x * (tidx + 1); i++) {
		for (int j = cpy_y * tidy; j < cpy_y * (tidy + 1); j++) {
			if ((tidx == 0 && i == 0) || (tidy == 0 && j == 0))
				continue;
			arrays.B[(i * arrays.B_size_y + j) * arrays.B_size_z + 2] = (arrays.Bz[(i-1) * arrays.Bz_size_y + j] + arrays.Bz[i * arrays.Bz_size_y + j]
																		+ arrays.Bz[i * arrays.Bz_size_y + j - 1] + arrays.Bz[(i-1) * arrays.Bz_size_y + j - 1]) / 4.0;

			B_mag[tidy * gridDim.x * blockDim.x + tidx] += sqrt(arrays.B[(i * arrays.B_size_y + j) * arrays.B_size_z + 2] * arrays.B[(i * arrays.B_size_y + j) * arrays.B_size_z + 2]);
		}
	}
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
	dim3 blockShape = dim3(10,10);
	dim3 gridShape = dim3(5,5);
	int total_threads = gridShape.x * blockShape.x * gridShape.y * blockShape.y;
	problem_set_up<<<1,1>>>(arrays, specifics);
	double *E_mag_vec = (double *) calloc(total_threads, sizeof(double));
	double *B_mag_vec = (double *) calloc(total_threads, sizeof(double));
	double *d_E_mag_vec, *d_B_mag_vec;
	cudaMalloc(&d_E_mag_vec, total_threads * sizeof(double));
	cudaMalloc(&d_B_mag_vec, total_threads * sizeof(double));

	// start at time 0
	double t = 0.0;
	int i = 0;
	while (i < steps) {
		apply_boundary<<<gridShape, blockShape>>>(arrays);
		update_fields<<<gridShape, blockShape>>>(constants, specifics, arrays);
		t += specifics.dt;

		if (i % output_freq == 0) {
			double E_mag = 0;
			double B_mag = 0;
			resolve_to_grid<<<gridShape, blockShape>>>(d_E_mag_vec, d_B_mag_vec, arrays);
			cudaDeviceSynchronize();
			cudaMemcpy(E_mag_vec, d_E_mag_vec, total_threads * sizeof(double), cudaMemcpyDeviceToHost);
			cudaMemcpy(B_mag_vec, d_B_mag_vec, total_threads * sizeof(double), cudaMemcpyDeviceToHost);
			for (int u = 0; u < total_threads; u++) {
				E_mag += E_mag_vec[u];
				B_mag += B_mag_vec[u];
			}
			printf("Step %8d, Time: %14.8e (dt: %14.8e), E magnitude: %14.8e, B magnitude: %14.8e\n", i, t, specifics.dt, E_mag, B_mag);

			if ((!no_output) && (enable_checkpoints)) {
				cudaMemcpy(&host_E[0][0][0], arrays.E, arrays.E_size_x * arrays.E_size_y * arrays.E_size_z * sizeof(double), cudaMemcpyDeviceToHost);
				write_checkpoint(i);
			}
		}

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


