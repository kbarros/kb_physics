/*
  alloc3d.c

  functions to dynamically allocate and free a three-dimensional matrix,
  which is stored contiguously in memory

  allocate_3d_*_matrix functions return NULL on allocation failure

  Erik Luijten - March 2000

  $Id: alloc3d.c,v 1.3 2006-08-24 19:09:29-05 luijten Exp luijten $
*/

static char rcsid[] =
        "$Id: alloc3d.c,v 1.3 2006-08-24 19:09:29-05 luijten Exp luijten $";

#include <stdio.h>
#include <stdlib.h>
#include "rcs.h"

#define PTR_SIZE sizeof(void *)

int ***allocate_3d_int_matrix(int k, int l, int m)
{
#undef VARIABLE
#define VARIABLE int /* this is what is stored in the matrix */
	int i,j;
	VARIABLE ***mem_ptr;
	char *tmp_ptr, *plane_ptr; /* when setting the pointers, we access this in a byte-like fashion */

	use_id(rcsid);
	
	/* allocate the actual array */
	tmp_ptr = (char *) malloc( k * l * m *sizeof(VARIABLE) );
	if (tmp_ptr == NULL) return(NULL);
	
	/* we have k pointers to the planes of the 3D array */
	mem_ptr = (VARIABLE ***) malloc(k * PTR_SIZE);
	if (mem_ptr == NULL) return(NULL);
	
	/* we have a k*l plane of pointers to the 2nd columns (i.e., to the "m-values") */
	plane_ptr = (char *) malloc(k * l* PTR_SIZE);
	if (plane_ptr == NULL) return(NULL);
	
	/* now we have to set these pointers to the appropriate parts of the actual 3D array */
	for(i=0; i<k; i++)
	{
		/* set the k-values to the plane */
		mem_ptr[i] = (VARIABLE **) (plane_ptr + i * l * PTR_SIZE);
		/* set the pointers in the plane to the 2nd columns of the actual matrix */
		for(j=0; j<l; j++)
		{
/*			printf("(%d %d) pointing to %d\n", i,j,(i*l+j)*m); */
			mem_ptr[i][j] = (VARIABLE *) (tmp_ptr + ((i*l+j)*m) * sizeof(VARIABLE));
		}
	}

	printf("Alloc3d: Allocated memory for %d int entries.\n", k*l*m);

	return(mem_ptr);
}

void free_3d_int_matrix(int ***mem_ptr)
{
	/* free actual matrix */
	free(mem_ptr[0][0]);
	/* free the plane */
	free(mem_ptr[0]);
	/* free array of column pointers */
	free((char *) mem_ptr);
}

short int ***allocate_3d_shortint_matrix(int k, int l, int m)
{
#undef VARIABLE
#define VARIABLE short int /* this is what is stored in the matrix */
/* NOTE: the pointers themselves are "int"-sized, NOT short int! */
	int i,j;
	short int ***mem_ptr;
	char *tmp_ptr, *plane_ptr; /* when setting the pointers, we access this in a byte-like fashion */
	
	/* allocate the actual array */
	tmp_ptr = (char *) malloc( k * l * m *sizeof(VARIABLE) );
	if (tmp_ptr == NULL) return(NULL);
	
	/* we have k pointers to the planes of the 3D array */
	mem_ptr = (VARIABLE ***) malloc(k * PTR_SIZE);
	if (mem_ptr == NULL) return(NULL);
	
	/* we have a k*l plane of pointers to the 2nd columns (i.e., to the "m-values") */
	plane_ptr = (char *) malloc(k * l* PTR_SIZE);
	if (plane_ptr == NULL) return(NULL);
	
	/* now we have to set these pointers to the appropriate parts of the actual 3D array */
	for(i=0; i<k; i++)
	{
		/* set the k-values to the plane */
		mem_ptr[i] = (VARIABLE **) (plane_ptr + i * l * PTR_SIZE);
		/* set the pointers in the plane to the 2nd columns of the actual matrix */
		for(j=0; j<l; j++)
		{
/*			printf("(%d %d) pointing to %d\n", i,j,(i*l+j)*m); */
			mem_ptr[i][j] = (VARIABLE *) (tmp_ptr + ((i*l+j)*m) * sizeof(VARIABLE));
		}
	}

/*	printf("Alloc3d: Allocated memory for %d short int entries.\n", k*l*m); */

	return(mem_ptr);
}

void free_3d_shortint_matrix(short int ***mem_ptr)
{
	/* free actual matrix */
	free(mem_ptr[0][0]);
	/* free the plane */
	free(mem_ptr[0]);
	/* free array of column pointers */
	free((char *) mem_ptr);
}

double ***allocate_3d_double_matrix(int k, int l, int m)
{
#undef VARIABLE
#define VARIABLE double /* this is what is stored in the matrix */
	int i,j;
	VARIABLE ***mem_ptr;
	char *tmp_ptr, *plane_ptr; /* when setting the pointers, we access this in a byte-like fashion */
	
	/* allocate the actual array */
	tmp_ptr = (char *) malloc( k * l * m *sizeof(VARIABLE) );
	if (tmp_ptr == NULL) return(NULL);
	
	/* we have k pointers to the planes of the 3D array */
	mem_ptr = (VARIABLE ***) malloc(k * PTR_SIZE);
	if (mem_ptr == NULL) return(NULL);
	
	/* we have a k*l plane of pointers to the 2nd columns (i.e., to the "m-values") */
	plane_ptr = (char *) malloc(k * l* PTR_SIZE);
	if (plane_ptr == NULL) return(NULL);
	
	/* now we have to set these pointers to the appropriate parts of the actual 3D array */
	for(i=0; i<k; i++)
	{
		/* set the k-values to the plane */
		mem_ptr[i] = (VARIABLE **) (plane_ptr + i * l * PTR_SIZE);
		/* set the pointers in the plane to the 2nd columns of the actual matrix */
		for(j=0; j<l; j++)
		{
/*			printf("(%d %d) pointing to %d\n", i,j,(i*l+j)*m); */
			mem_ptr[i][j] = (VARIABLE *) (tmp_ptr + ((i*l+j)*m) * sizeof(VARIABLE));
		}
	}

/*	printf("Alloc3d: Allocated memory for %d double entries.\n", k*l*m); */

	return(mem_ptr);
}

void free_3d_double_matrix(double ***mem_ptr)
{
	/* free actual matrix */
	free(mem_ptr[0][0]);
	/* free the plane */
	free(mem_ptr[0]);
	/* free array of column pointers */
	free((char *) mem_ptr);
}
