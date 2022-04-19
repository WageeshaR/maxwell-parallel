#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <sys/time.h>

#include "args.h"
#include "vtk.h"
#include "data.h"
#include "setup.h"

/**
 * @brief Update the magnetic and electric fields. The magnetic fields are updated for a half-time-step. The electric fields are updated for a full time-step.
 * 
 */
__global__ void update_fields(Constants constants, Specifics specifics, Arrays arrays) {
	/** 
	 * @param tidx x id of current thread
	 * @param tidy y id of current thread
	 * @param cpx_x number of cells in X direction for this thread
	 * @param cpx_y number of cells in Y direction for this thread
	 * each memory access are adjusted with the pitched alignment
	 */
	int tidx = blockIdx.x*blockDim.x + threadIdx.x;
	int tidy = blockIdx.y*blockDim.y + threadIdx.y;
	int cpx_ex = arrays.Ex_size_x / (gridDim.x * blockDim.x);
	int cpy_ey = arrays.Ey_size_y / (gridDim.y * blockDim.y);

	for (int i = cpx_ex * tidx; i < cpx_ex * (tidx + 1); i++) {
		for (int j = cpy_ey * tidy; j < cpy_ey * (tidy + 1); j++) {
			arrays.Bz[i * arrays.bz_pitch + j] = arrays.Bz[i * arrays.bz_pitch + j] - (specifics.dt / specifics.dx) * (arrays.Ey[(i+1) * arrays.ey_pitch + j] - arrays.Ey[i * arrays.ey_pitch + j])
				                					+ (specifics.dt / specifics.dy) * (arrays.Ex[i * arrays.ex_pitch + j + 1] - arrays.Ex[i * arrays.ex_pitch + j]);
		}
	}

	for (int i = cpx_ex * tidx; i < cpx_ex * (tidx + 1); i++) {
		for (int j = cpy_ey * tidy; j < cpy_ey * (tidy + 1); j++) {
			if (tidy == 0 && j == 0)
				continue;
			arrays.Ex[i * arrays.ex_pitch + j] = arrays.Ex[i * arrays.ex_pitch + j]
													+ (specifics.dt / (specifics.dy * constants.eps * constants.mu)) * (arrays.Bz[i * arrays.bz_pitch + j] - arrays.Bz[i * arrays.bz_pitch + j - 1]);
		}
	}

	for (int i = cpx_ex * tidx; i < cpx_ex * (tidx + 1); i++) {
		for (int j = cpy_ey * tidy; j < cpy_ey * (tidy + 1); j++) {
			if (tidx == 0 && i == 0)
				continue;
			arrays.Ey[i * arrays.ey_pitch + j] = arrays.Ey[i * arrays.ey_pitch + j]
													- (specifics.dt / (specifics.dx * constants.eps * constants.mu)) * (arrays.Bz[i * arrays.bz_pitch + j] - arrays.Bz[(i - 1) * arrays.bz_pitch + j]);
		}
	}
}

/**
 * @brief Apply boundary conditions
 * 
 */
__global__ void apply_boundary(Arrays arrays) {
	/** 
	 * @param tidx x id of current thread
	 * @param tidy y id of current thread
	 * @param cpx_x number of cells in X direction for this thread
	 * @param cpx_y number of cells in Y direction for this thread
	 * each memory access are adjusted with the pitched alignment
	 */
	int tidx = blockIdx.x*blockDim.x + threadIdx.x;
	int tidy = blockIdx.y*blockDim.y + threadIdx.y;
	int cpx_ex = arrays.Ex_size_x / (gridDim.x * blockDim.x);
	int cpy_ey = arrays.Ey_size_y / (gridDim.y * blockDim.y);
	
	for (int i = tidx * cpx_ex; i < cpx_ex * (tidx + 1); i++) {
		if (tidy == 0) {
			arrays.Ex[i * arrays.ex_pitch] = -arrays.Ex[i * arrays.ex_pitch + 1];
		}
		if (tidy == gridDim.y * blockDim.y - 1)
			arrays.Ex[i * arrays.ex_pitch + arrays.Ex_size_y - 1] = -arrays.Ex[i * arrays.ex_pitch + arrays.Ex_size_y - 2];
	}

	for (int j = tidy * cpy_ey; j < cpy_ey * (tidy + 1); j++) {
		if (tidx == 0)
			arrays.Ey[j] = -arrays.Ey[arrays.ey_pitch + j];
		if (tidx == gridDim.x * blockDim.x - 1)
			arrays.Ey[arrays.ey_pitch * (arrays.Ey_size_x - 1) + j] = -arrays.Ey[arrays.ey_pitch * (arrays.Ey_size_x - 2) + j];
	}
}

/**
 * @brief Resolve the Ex, Ey and Bz fields to grid points and sum the magnitudes for output
 * 
 * @param E_mag The returned total magnitude of the Electric field (E)
 * @param B_mag The returned total magnitude of the Magnetic field (B) 
 */
__global__ void resolve_to_grid(double *E_mag, double *B_mag, Arrays arrays) {
	/** 
	 * @param tidx x id of current thread
	 * @param tidy y id of current thread
	 * @param cpx_x number of cells in X direction for this thread
	 * @param cpx_y number of cells in Y direction for this thread
	 * each memory access are adjusted with the pitched alignment
	 */
	int tidx = blockIdx.x*blockDim.x + threadIdx.x;
	int tidy = blockIdx.y*blockDim.y + threadIdx.y;
	int cpx_x = (arrays.E_size_x - 1) / (gridDim.x * blockDim.x);
	int cpy_y = (arrays.E_size_y - 1) / (gridDim.y * blockDim.y);
	E_mag[tidy * gridDim.x * blockDim.x + tidx] = 0.0;
	B_mag[tidy * gridDim.x * blockDim.x + tidx] = 0.0;

	for (int i = cpx_x * tidx; i < cpx_x * (tidx + 1); i++) {
		for (int j = cpy_y * tidy; j < cpy_y * (tidy + 1); j++) {
			if ((tidx == 0 && i == 0) || (tidy == 0 && j == 0))
				continue;
			arrays.E[i * arrays.e_pitch + j * arrays.E_size_z] = (arrays.Ex[(i-1) * arrays.ex_pitch + j] + arrays.Ex[i * arrays.ex_pitch + j]) / 2.0;
			arrays.E[i * arrays.e_pitch + j * arrays.E_size_z + 1] = (arrays.Ey[i * arrays.ey_pitch + j - 1] + arrays.Ey[i * arrays.ey_pitch + j]) / 2.0;

			E_mag[tidy * gridDim.x * blockDim.x + tidx] += sqrt((arrays.E[i * arrays.e_pitch + j * arrays.E_size_z] * arrays.E[i * arrays.e_pitch + j * arrays.E_size_z])
																+ (arrays.E[i * arrays.e_pitch + j * arrays.E_size_z + 1] * arrays.E[i * arrays.e_pitch + j * arrays.E_size_z + 1]));
		}
	}
	
	for (int i = cpx_x * tidx; i < cpx_x * (tidx + 1); i++) {
		for (int j = cpy_y * tidy; j < cpy_y * (tidy + 1); j++) {
			if ((tidx == 0 && i == 0) || (tidy == 0 && j == 0))
				continue;
			arrays.B[i * arrays.b_pitch + j * arrays.B_size_z + 2] = (arrays.Bz[(i-1) * arrays.bz_pitch + j] + arrays.Bz[i * arrays.bz_pitch + j]
																		+ arrays.Bz[i * arrays.bz_pitch + j - 1] + arrays.Bz[(i-1) * arrays.bz_pitch + j - 1]) / 4.0;

			B_mag[tidy * gridDim.x * blockDim.x + tidx] += sqrt(arrays.B[i * arrays.b_pitch + j * arrays.B_size_z + 2] * arrays.B[i * arrays.b_pitch + j * arrays.B_size_z + 2]);
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
	total_error = 0;

	/**
	 * @brief 
	 * setting up for comparison mode
	 */
	FILE *comp_file;			// comparison file pointer
	char comp_file_name[1024];
	char *buffer;				// buffer used to read comparison file
	size_t bufsize = 1024;
	if (comp_mode != 0) {
		sprintf(comp_file_name_base, "comp/comp_%d_%d%%s", specifics.X, specifics.Y);
		if (comp_mode == 1) {
			sprintf(comp_file_name, comp_file_name_base, "_mag.cmp");
			comp_file = fopen(comp_file_name, "r");
			if (comp_file == NULL) {
				printf("Unable to find comparison file %s\n", comp_file_name);
				exit(1);
			}
		}
	}

	printf("Running problem size %f x %f on a %d x %d grid.\n", specifics.lengthX, specifics.lengthY, specifics.X, specifics.Y);
	
	if (verbose) print_opts();
	
	allocate_arrays();

	dim3 blockShape = dim3(cuda_consts.block_x, cuda_consts.block_y);
	dim3 gridShape = dim3(cuda_consts.grid_x, cuda_consts.grid_y);
	int total_threads = gridShape.x * blockShape.x * gridShape.y * blockShape.y;
	
	problem_set_up<<<1,1>>>(arrays, specifics);
	
	// below vectors are needed to store E_mag and B_mag values in device and host
	double *E_mag_vec = (double *) calloc(total_threads, sizeof(double));
	double *B_mag_vec = (double *) calloc(total_threads, sizeof(double));
	double *d_E_mag_vec, *d_B_mag_vec;
	cudaMalloc(&d_E_mag_vec, total_threads * sizeof(double));
	cudaMalloc(&d_B_mag_vec, total_threads * sizeof(double));
	
	long e_pitch_host = arrays.E_size_y * arrays.E_size_z * sizeof(double);
	long b_pitch_host = arrays.B_size_y * arrays.B_size_z * sizeof(double);

	double t = 0.0;
	int i = 0;
	int comp_line_len = 0;
	long start, end;
    struct timeval timecheck;

    gettimeofday(&timecheck, NULL);
    start = (long)timecheck.tv_sec * 1000 + (long)timecheck.tv_usec / 1000;

	while (i < steps) {
		apply_boundary<<<gridShape, blockShape>>>(arrays);
		update_fields<<<gridShape, blockShape>>>(constants, specifics, arrays);
		t += specifics.dt;
		
		// Reading the line here (disregarding output frequency) to keep track of lines
		if (comp_mode == 1) {
			buffer = (char *) calloc(1024, sizeof(char));
			comp_line_len = getline(&buffer, &bufsize, comp_file);
		}

		if (i % output_freq == 0) {
			double E_mag = 0;
			double B_mag = 0;
			resolve_to_grid<<<gridShape, blockShape>>>(d_E_mag_vec, d_B_mag_vec, arrays);
			cudaMemcpy(E_mag_vec, d_E_mag_vec, total_threads * sizeof(double), cudaMemcpyDeviceToHost);
			cudaMemcpy(B_mag_vec, d_B_mag_vec, total_threads * sizeof(double), cudaMemcpyDeviceToHost);
			
			// summing over the vector to calculate total magnitude from all threads
			for (int u = 0; u < total_threads; u++) {
				E_mag += E_mag_vec[u];
				B_mag += B_mag_vec[u];
			}
			
			if (comp_mode == 1) {
				double mags[2] = { E_mag, B_mag };
				compare_line(comp_line_len, &buffer, mags);
			}
			printf("Step %8d, Time: %14.8e (dt: %14.8e), E magnitude: %14.8e, B magnitude: %14.8e\n", i, t, specifics.dt, E_mag, B_mag);

			if ((!no_output) && (enable_checkpoints)) {
				// used cudaMemcpy2D to adjust for pitch alignments when copying
				cudaMemcpy2D(&host_E[0][0][0], e_pitch_host, arrays.E, arrays.e_pitch * sizeof(double), e_pitch_host, arrays.E_size_x, cudaMemcpyDeviceToHost);
				cudaMemcpy2D(&host_B[0][0][0], b_pitch_host, arrays.B, arrays.b_pitch * sizeof(double), b_pitch_host, arrays.B_size_x, cudaMemcpyDeviceToHost);
				write_checkpoint(i);
			}
		}

		i++;
	}

	double E_mag, B_mag;
	resolve_to_grid<<<gridShape, blockShape>>>(d_E_mag_vec, d_B_mag_vec, arrays);
	cudaMemcpy(E_mag_vec, d_E_mag_vec, total_threads * sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(B_mag_vec, d_B_mag_vec, total_threads * sizeof(double), cudaMemcpyDeviceToHost);
	for (int u = 0; u < total_threads; u++) {
		E_mag += E_mag_vec[u];
		B_mag += B_mag_vec[u];
	}
	printf("Step %8d, Time: %14.8e (dt: %14.8e), E magnitude: %14.8e, B magnitude: %14.8e\n", i, t, specifics.dt, E_mag, B_mag);
	printf("Simulation complete.\n");

	gettimeofday(&timecheck, NULL);
    end = (long)timecheck.tv_sec * 1000 + (long)timecheck.tv_usec / 1000;
	
	printf("Elapsed time is %.4fs\n", (double)(end - start) / 1000.0);

	if (!no_output) {
		cudaMemcpy2D(&host_E[0][0][0], e_pitch_host, arrays.E, arrays.e_pitch * sizeof(double), e_pitch_host, arrays.E_size_x, cudaMemcpyDeviceToHost);
		cudaMemcpy2D(&host_B[0][0][0], b_pitch_host, arrays.B, arrays.b_pitch * sizeof(double), b_pitch_host, arrays.B_size_x, cudaMemcpyDeviceToHost);
		write_result();
	}
	
	if (comp_mode != 0)
		printf("Total error is %.15e\n", total_error);

	free_arrays();

	exit(0);
}


