#ifndef ARGS_H
#define ARGS_H

extern int verbose;
extern int input_problem;
extern int no_output;
extern int output_freq;
extern int enable_checkpoints;
extern int omp_num_threads;
extern int enable_comparison;

void parse_args(int argc, char *argv[]);
void print_opts();

#endif