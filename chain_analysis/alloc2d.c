/*
  alloc2d.c

  Functions to dynamically allocate and free two-dimensional array.
  Works on 32-bit Intel, 64-bit Intel, 64-bit Alpha, 64-bit PPC.

  allocate_*_matrix() routines return NULL on allocation failure

  $Id: alloc2d.c,v 1.4 2006-08-24 12:28:07-05 luijten Exp luijten $
*/

static char rcsid[] =
	"$Id: alloc2d.c,v 1.4 2006-08-24 12:28:07-05 luijten Exp luijten $";

#include <stdio.h>
#include <stdlib.h>
#include "rcs.h"


#define PTR_SIZE sizeof(void *)

#define MAXINT 2147483647

double **allocate_double_matrix(int n_row, int n_col)
{

#undef  VARIABLE
#define VARIABLE double /* this is what is stored in the matrix */

	int i;
	VARIABLE **mem_ptr;
	char *tmp_ptr; /* we just need a byte-valued pointer */
	
	use_id(rcsid);
	
	/* allocate array of pointers to double */
	mem_ptr = malloc(n_row * PTR_SIZE);
	if (mem_ptr == NULL) return(NULL);
	
	/* allocate space for the actual matrix and store it temporarily
	   in tmp_ptr */
	if ((double)n_row*(double)n_col*(double)sizeof(VARIABLE) > MAXINT)
        {
/*              fprintf(stderr, "Overflow in alloc2d.\n"); */
                free((char *) mem_ptr);
                return(NULL);
        }
	tmp_ptr = malloc(n_row * n_col * sizeof(VARIABLE));
	if (tmp_ptr == NULL)
	{
		free((char *)mem_ptr);
		return(NULL);
	}
	
	/* fill mem_ptr with pointers to the rows of the matrix */
	for(i=0; i < n_row; i++)
	{
		mem_ptr[i] = (VARIABLE *) (tmp_ptr + i * n_col * sizeof(VARIABLE));
	}
	
/*	printf("Allocated memory for %d entries\n", n_row*n_col); */
	
	return(mem_ptr);
}

void free_double_matrix(double **mem_ptr)
{
	/* free actual matrix */
	free(mem_ptr[0]);
	/* free array of row_pointers */
	free((char *) mem_ptr);
}

int **allocate_int_matrix(int n_row, int n_col)
{

#undef  VARIABLE
#define VARIABLE int /* this is what is stored in the matrix */

	int i;
	VARIABLE **mem_ptr;
	char *tmp_ptr; /* we just need a byte-valued pointer */

	/* allocate array of pointers to double */
	mem_ptr = malloc(n_row * PTR_SIZE);
	if (mem_ptr == NULL) return(NULL);

	/* allocate space for the actual matrix and store it temporarily
	   in tmp_ptr */
	tmp_ptr = malloc(n_row * n_col * sizeof(VARIABLE));
	if (tmp_ptr == NULL) return(NULL);

	/* fill mem_ptr with pointers to the rows of the matrix */
	for(i=0; i < n_row; i++)
	{
		mem_ptr[i] = (VARIABLE *) (tmp_ptr + i * n_col * sizeof(VARIABLE));
	}

/*	printf("Allocated memory for %d entries\n", n_row*n_col); */

	return(mem_ptr);
}

void free_int_matrix(int **mem_ptr)
{
	/* free actual matrix */
	free(mem_ptr[0]);
	/* free array of row_pointers */
	free((char *) mem_ptr);
}

long int **allocate_longint_matrix(int n_row, int n_col)
{

#undef  VARIABLE
#define VARIABLE long int /* this is what is stored in the matrix */

	int i;
	VARIABLE **mem_ptr;
	char *tmp_ptr; /* we just need a byte-valued pointer */

	/* allocate array of pointers to double */
	mem_ptr = (VARIABLE **) malloc(n_row * PTR_SIZE);
	if (mem_ptr == NULL) return(NULL);

	/* allocate space for the actual matrix and store it temporarily
	   in tmp_ptr */
	tmp_ptr = (char *) malloc(n_row * n_col * sizeof(VARIABLE));
	if (tmp_ptr == NULL) return(NULL);

	/* fill mem_ptr with pointers to the rows of the matrix */
	for(i=0; i < n_row; i++)
	{
		mem_ptr[i] = (VARIABLE *) (tmp_ptr + i * n_col * sizeof(VARIABLE));
	}

/*	printf("Allocated memory for %d entries\n", n_row*n_col); */

	return(mem_ptr);
}


void free_longint_matrix(long int **mem_ptr)
{
	/* free actual matrix */
	free(mem_ptr[0]);
	/* free array of row_pointers */
	free((char *) mem_ptr);
}


char **allocate_char_matrix(int n_row, int n_col)
{

#undef  VARIABLE
#define VARIABLE char /* this is what is stored in the matrix */

	int i;
	VARIABLE **mem_ptr;
	char *tmp_ptr; /* we just need a byte-valued pointer */

	/* allocate array of pointers to double */
	mem_ptr = (VARIABLE **) malloc(n_row * PTR_SIZE);
	if (mem_ptr == NULL) return(NULL);

	/* allocate space for the actual matrix and store it temporarily
	   in tmp_ptr */
	tmp_ptr = (char *) malloc(n_row * n_col * sizeof(VARIABLE));
	if (tmp_ptr == NULL) return(NULL);

	/* fill mem_ptr with pointers to the rows of the matrix */
	for(i=0; i < n_row; i++)
	{
		mem_ptr[i] = (VARIABLE *) (tmp_ptr + i * n_col * sizeof(VARIABLE));
	}

/*	printf("Allocated memory for %d entries\n", n_row*n_col); */

	return(mem_ptr);
}

void free_char_matrix(char **mem_ptr)
{
	/* free actual matrix */
	free(mem_ptr[0]);
	/* free array of row_pointers */
	free((char *) mem_ptr);
}
