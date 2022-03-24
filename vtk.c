#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <omp.h>

#include "vtk.h"
#include "data.h"
#include "setup.h"
#include "args.h"

char checkpoint_basename[1024];
char result_filename[1024];

/**
 * @brief Set the default basename for file output to out/maxwell
 *
 */
void set_default_base() {
    set_basename("maxwell");
}

/**
 * @brief Set the basename for file output
 *
 * @param base Basename string
 */
void set_basename(char *base) {
    checkpoint_basename[0] = '\0';
    result_filename[0] = '\0';
    sprintf(checkpoint_basename, "%s-%%d.vtk", base);
    sprintf(result_filename, "%s.vtk", base);
}

/**
 * @brief Get the basename for file output
 *
 * @return char* Basename string
 */
char *get_basename() {
    return checkpoint_basename;
}

/**
 * @brief Write a checkpoint VTK file (with the iteration number in the filename)
 *
 * @param iteration The current iteration number
 * @return int Return whether the write was successful
 */
int write_checkpoint(int iteration) {
    char filename[1024];
    sprintf(filename, checkpoint_basename, iteration);
    return write_vtk(filename);
}

/**
 * @brief Write the final output to a VTK file
 *
 * @return int Return whether the write was successful
 */
int write_result() {
    return write_vtk(result_filename);
}

/**
 * @brief Write a VTK file with the current state of the Electric and Magnetic Fields
 *
 * @param filename The filename to write out
 * @return int Return whether the write was successful
 */
int write_vtk(char* filename) {
    FILE *f, **files;
    char **file_names;
    int num_t = omp_num_threads;
    files = (FILE **) malloc(num_t * sizeof(FILE *));
    file_names = (char **) malloc(num_t * sizeof(char *));
    for (int i = 0; i < num_t; i++) {
        char f_name[1024];
        char num[16];
        sprintf(num, "%d", i);
        strcpy(f_name, filename);
        strcat(f_name, num);
        file_names[i] = (char *) malloc(1024 * sizeof(char));
        strcpy(file_names[i], f_name);
        files[i] = fopen(f_name, "w+");
        if (files[i] == NULL) {
            perror("Error opening temporary file");
            return -1;
        }
    }
    f = fopen(filename, "w");
    if (f == NULL) {
        perror("Error");
        return -1;
    }

    // Write the VTK header information
    fprintf(f, "# vtk DataFile Version 3.0\n");
    fprintf(f, "Karman Output\n");
    fprintf(f, "ASCII\n");
    fprintf(f, "DATASET RECTILINEAR_GRID\n");

    // Write out the grid information
    fprintf(f, "DIMENSIONS %d %d 1\n", (X+1), (Y+1));
    fprintf(f, "X_COORDINATES %d float\n", (X+1));
    for (int i = 0; i <= X; i++) fprintf(f, "   %.12e", (lengthX * ((double) i / (X+1))));
    fprintf(f, "\nY_COORDINATES %d float\n", (Y+1));
    for (int i = 0; i <= Y; i++) fprintf(f, "   %.12e", (lengthY * ((double) i / (Y+1))));
    fprintf(f, "\nZ_COORDINATES 1 float\n");
    fprintf(f, "  0.000000000000e+00");

    fprintf(f, "\nPOINT_DATA %d\n", ((X+1) * (Y+1)));

    if (enable_comparison == 1) {
        O_E = alloc_3d_array(E_size_x, E_size_y, E_size_z);
		O_B = alloc_3d_array(B_size_x, B_size_y, B_size_z);
        load_originals(filename);
    }

    // Write out the E and B vector fields
    fprintf(f, "VECTORS E_field float\n");

    #pragma omp parallel for schedule(static) collapse(2)
    for (int j = 0; j <= Y; j++) {
        for (int i = 0; i <= X; i++) {
            int t_num = omp_get_thread_num();
            fprintf(files[t_num], "  %.12e %.12e 0.000000000000e+00\n", E[i][j][0], E[i][j][1]);
            if (enable_comparison == 1)
                total_error += abs(O_E[i][j][0] - E[i][j][0]) + abs(O_E[i][j][1] - E[i][j][1]);
        }
    }

    f = freopen(filename, "a", f);
    for (int i = 0; i < num_t; i++) {
        files[i] = freopen(file_names[i], "r", files[i]);
        char buffer[256];
        while(fgets(buffer, 256, files[i]) != NULL) {
            fprintf(f, "%s", buffer);
        }
        files[i] = freopen(file_names[i], "w+", files[i]);
    }
    fprintf(f, "VECTORS B_field float\n");

    #pragma omp parallel for schedule(static) collapse(2)
    for (int j = 0; j <= Y; j++) {
        for (int i = 0; i <= X; i++) {
            int t_num = omp_get_thread_num();
            fprintf(files[t_num], "  0.000000000000e+00 0.000000000000e+00 %.12e\n", B[i][j][2]);
            if (enable_comparison == 1)
                total_error += abs(O_B[i][j][0] - B[i][j][0]) + abs(O_B[i][j][1] - B[i][j][1]);
        }
    }

    for (int i = 0; i < num_t; i++) {
        files[i] = freopen(file_names[i], "r", files[i]);
        char buffer[256];
        while(fgets(buffer, 256, files[i]) != NULL) {
            fprintf(f, "%s", buffer);
        }
        fclose(files[i]);
        remove(file_names[i]);
        free(file_names[i]);
    }
    
    free(files);
    free(file_names);
    fclose(f);

    if (enable_comparison == 1) {
        free_3d_array(O_E);
        free_3d_array(O_B);
    }
    return 0;
}
