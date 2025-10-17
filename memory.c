#include <stdio.h>
#include <stdlib.h>
#include "structs.h"
#include "nbrlist.h"

void alloc_vectors(struct Vectors *p_vectors, size_t num_part)
/* Allocate the arrays in 'vectors' needed to store information of all particles */
{
    p_vectors->type = (int *)malloc(num_part * sizeof(size_t));
    p_vectors->radius = (double *)malloc(num_part * sizeof(double));
    p_vectors->mass = (double *)malloc(num_part * sizeof(double));
    p_vectors->r = (struct Vec3D *)malloc(num_part * sizeof(struct Vec3D));
    p_vectors->dr = (struct Vec3D *)malloc(num_part * sizeof(struct Vec3D));
    p_vectors->v = (struct Vec3D *)malloc(num_part * sizeof(struct Vec3D));
    p_vectors->omega = (struct Vec3D *)malloc(num_part * sizeof(struct Vec3D));
    p_vectors->f = (struct Vec3D *)malloc(num_part * sizeof(struct Vec3D));
    p_vectors->T = (struct Vec3D *)malloc(num_part * sizeof(struct Vec3D));

    // initialise histogram arrays
    p_vectors->hist_vol_r = calloc(p_vectors->hist_num_bins_r, sizeof(double));
    p_vectors->hist_vol_z = calloc(p_vectors->hist_num_bins_z, sizeof(double));
}

void free_vectors(struct Vectors *p_vectors)
/* Free the arrays in 'vectors' */
{
    free(p_vectors->type);
    p_vectors->type = NULL;
    free(p_vectors->radius);
    p_vectors->radius = NULL;
    free(p_vectors->mass);
    p_vectors->mass = NULL;
    free(p_vectors->r);
    p_vectors->r = NULL;
    free(p_vectors->dr);
    p_vectors->dr = NULL;
    free(p_vectors->v);
    p_vectors->v = NULL;
    free(p_vectors->omega);
    p_vectors->omega = NULL;
    free(p_vectors->f);
    p_vectors->f = NULL;
    free(p_vectors->T);
    p_vectors->T = NULL;
    free(p_vectors->hist_vol_r);
    p_vectors->hist_vol_r = NULL;  
    free(p_vectors->hist_vol_z);
    p_vectors->hist_vol_z = NULL;
}

void alloc_memory(struct Parameters *p_parameters, struct Vectors *p_vectors, struct Nbrlist *p_nbrlist, struct Colllist *p_colllist)
/* Allocate all variables needed in the MD simulation */
{
    alloc_vectors(p_vectors, p_parameters->num_part);
    alloc_nbrlist(p_parameters, p_nbrlist);
    alloc_colllist(p_parameters, p_colllist);
}

void free_memory(struct Vectors *p_vectors, struct Nbrlist *p_nbrlist, struct Colllist *p_colllist)
/* Free the memory allocated by alloc_memory */
{
    free_vectors(p_vectors);
    free_nbrlist(p_nbrlist);
    free_colllist(p_colllist);
}
