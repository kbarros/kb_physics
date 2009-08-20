/*
  alloc3di.h
  header file for alloc3d.c

  Erik Luijten - March 2000

  $Id: alloc3d.h,v 1.1 2002-07-15 10:51:10+02 luijten Exp luijten $
*/

int ***allocate_3d_int_matrix(int k, int l, int m);
void free_3d_int_matrix(int ***matrix);

short int ***allocate_3d_shortint_matrix(int k, int l, int m);
void free_3d_shortint_matrix(short int ***matrix);

double ***allocate_3d_double_matrix(int k, int l, int m);
void free_3d_double_matrix(double ***matrix);
