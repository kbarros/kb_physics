/*
  alloc2d.h
  header file for alloc2d.c

  Erik Luijten - October 2000

  $Id: alloc2d.h,v 1.1 2002-07-15 10:51:19+02 luijten Exp luijten $
*/

double **allocate_double_matrix(int n_row, int n_col);
void free_double_matrix(double **mem_ptr);

int **allocate_int_matrix(int n_row, int n_col);
void free_int_matrix(int **mem_ptr);

