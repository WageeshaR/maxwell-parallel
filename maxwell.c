#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <mpi.h>

#include "args.h"
#include "vtk.h"
#include "data.h"
#include "setup.h"

/**
 * @brief Update the magnetic and electric fields. The magnetic fields are updated for a half-time-step. The electric fields are updated for a full time-step.
 * 
 */
void update_fields(MPI_Datatype datatype1, MPI_Datatype datatype2, MPI_Datatype datatype3) {

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Sendrecv(Ey[0], 1, datatype2, left, 13, Ey[Ey_size_x-1], 1, datatype2, right, 13, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	for (int i = 0; i < Bz_size_x; i++) {
		for (int j = 0; j < Bz_size_y; j++) {
			Bz[i][j] = Bz[i][j] - (dt / dx) * (Ey[i+1][j] - Ey[i][j])
				                + (dt / dy) * (Ex[i][j+1] - Ex[i][j]);
		}
	}

	for (int i = 0; i < Ex_size_x; i++) {
		for (int j = 1; j < Ex_size_y-1; j++) {
			Ex[i][j] = Ex[i][j] + (dt / (dy * eps * mu)) * (Bz[i][j] - Bz[i][j-1]);
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Sendrecv(Bz[Bz_size_x-1], 1, datatype3, right, 13, Bz[Bz_size_x], 1, datatype3, left, 13, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	if (rank == 0) {
		for (int i = 1; i < Ey_size_x-1; i++) {
			for (int j = 0; j < Ey_size_y; j++) {
				Ey[i][j] = Ey[i][j] - (dt / (dx * eps * mu)) * (Bz[i][j] - Bz[i-1][j]);
			}
		}
	} else {
		for (int i = 0; i < Ey_size_x-1; i++) {
			for (int j = 0; j < Ey_size_y; j++) {
				if (i == 0)
					Ey[i][j] = Ey[i][j] - (dt / (dx * eps * mu)) * (Bz[i][j] - Bz[Bz_size_x][j]);
				else
					Ey[i][j] = Ey[i][j] - (dt / (dx * eps * mu)) * (Bz[i][j] - Bz[i-1][j]);
			}
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Sendrecv(Ex[Ex_size_x-1], 1, datatype1, right, 13, Ex[Ex_size_x], 1, datatype1, left, 13, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}

/**
 * @brief Apply boundary conditions
 * 
 */
void apply_boundary() {
	for (int i = 0; i < Ex_size_x; i++) {
		Ex[i][0] = -Ex[i][1];
		Ex[i][Ex_size_y-1] = -Ex[i][Ex_size_y-2];
	}

	for (int j = 0; j < Ey_size_y; j++) {
		if (rank == 0)
			Ey[0][j] = -Ey[1][j];
		if (rank == size -1)
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
	*E_mag = 0.0;
	*B_mag = 0.0;
	
	if (rank == 0) {
		for (int i = 1; i < E_size_x-1; i++) {
			for (int j = 1; j < E_size_y-1; j++) {
				E[i][j][0] = (Ex[i-1][j] + Ex[i][j]) / 2.0;
				E[i][j][1] = (Ey[i][j-1] + Ey[i][j]) / 2.0;

				*E_mag += sqrt((E[i][j][0] * E[i][j][0]) + (E[i][j][1] * E[i][j][1]));
			}
		}

		for (int i = 1; i < B_size_x-1; i++) {
			for (int j = 1; j < B_size_y-1; j++) {
				B[i][j][2] = (Bz[i-1][j] + Bz[i][j] + Bz[i][j-1] + Bz[i-1][j-1]) / 4.0;

				*B_mag += sqrt(B[i][j][2] * B[i][j][2]);
			}
		}
	} else {
		for (int i = 0; i < E_size_x-1; i++) {
			for (int j = 1; j < E_size_y-1; j++) {
				if (i == 0) {
					E[i][j][0] = (Ex[Bz_size_x][j] + Ex[i][j]) / 2.0;
					E[i][j][1] = (Ey[i][j-1] + Ey[i][j]) / 2.0;
				}
				else {
					E[i][j][0] = (Ex[i-1][j] + Ex[i][j]) / 2.0;
					E[i][j][1] = (Ey[i][j-1] + Ey[i][j]) / 2.0;
				}

				*E_mag += sqrt((E[i][j][0] * E[i][j][0]) + (E[i][j][1] * E[i][j][1]));
			}
		}

		for (int i = 0; i < B_size_x-1; i++) {
			for (int j = 1; j < B_size_y-1; j++) {
				if (i == 0)
					B[i][j][2] = (Bz[Bz_size_x][j] + Bz[i][j] + Bz[i][j-1] + Bz[Bz_size_x][j-1]) / 4.0;
				else
					B[i][j][2] = (Bz[i-1][j] + Bz[i][j] + Bz[i][j-1] + Bz[i-1][j-1]) / 4.0;

				*B_mag += sqrt(B[i][j][2] * B[i][j][2]);
			}
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
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	set_defaults();
	parse_args(argc, argv);
	setup();

	if (rank == 0)
		printf("Running problem size %f x %f on a %d x %d grid.\n", lengthX, lengthY, X*size, Y);
	
	if (verbose) print_opts();
	
	allocate_arrays();

	problem_set_up();

	// start at time 0
	double t = 0.0;
	int i = 0;
	double global_E_mag, global_B_mag;

	MPI_Datatype ex_column, ey_column, bz_column;
	MPI_Type_vector(1, Ex_size_y, Ex_size_y, MPI_DOUBLE, &ex_column);
	MPI_Type_commit(&ex_column);
	MPI_Type_vector(1, Ey_size_y, Ey_size_y, MPI_DOUBLE, &ey_column);
	MPI_Type_commit(&ey_column);
	MPI_Type_vector(1, Bz_size_y, Bz_size_y, MPI_DOUBLE, &bz_column);
	MPI_Type_commit(&bz_column);
	left = rank-1 < 0 ? MPI_PROC_NULL : rank-1;
	right = rank+1 >= size ? MPI_PROC_NULL: rank+1;
	
	MPI_Datatype global_3d_grid;
	MPI_Type_vector(E_size_x-1, E_size_z*E_size_y, E_size_z*E_size_y, MPI_DOUBLE, &global_3d_grid);
	MPI_Type_commit(&global_3d_grid);

	while (i < steps) {
		apply_boundary();
		update_fields(ex_column, ey_column, bz_column);
		t += dt;

		if (i % output_freq == 0) {
			double E_mag, B_mag;
			resolve_to_grid(&E_mag, &B_mag);
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Reduce(&E_mag, &global_E_mag, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
			MPI_Reduce(&B_mag, &global_B_mag, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
			if (rank == 0)
				printf("Step %8d, Time: %14.8e (dt: %14.8e), E magnitude: %14.8e, B magnitude: %14.8e\n", i, t, dt, global_E_mag, global_B_mag);
			
			MPI_Gather(E[0][0], 1, global_3d_grid, global_E[0][0], 1, global_3d_grid, 0, MPI_COMM_WORLD);
			MPI_Gather(B[0][0], 1, global_3d_grid, global_B[0][0], 1, global_3d_grid, 0, MPI_COMM_WORLD);
			if ((!no_output) && (enable_checkpoints) && rank == 0)
				write_checkpoint(i);
		}
		i++;
	}

	double E_mag, B_mag;
	resolve_to_grid(&E_mag, &B_mag);
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Reduce(&E_mag, &global_E_mag, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&B_mag, &global_B_mag, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	if (rank == 0) {
		printf("Step %8d, Time: %14.8e (dt: %14.8e), E magnitude: %14.8e, B magnitude: %14.8e\n", i, t, dt, global_E_mag, global_B_mag);
		printf("Simulation complete.\n");
	}

	MPI_Type_free(&ey_column);
	MPI_Gather(E[0][0], 1, global_3d_grid, global_E[0][0], 1, global_3d_grid, 0, MPI_COMM_WORLD);
	MPI_Gather(B[0][0], 1, global_3d_grid, global_B[0][0], 1, global_3d_grid, 0, MPI_COMM_WORLD);
	if (!no_output && rank == 0) 
		write_result();

	free_arrays();

	exit(0);
}


