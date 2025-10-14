#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "constants.h"
#include "structs.h"
#include "random.h"
#include "initialise.h"

// Initializes particle properties such as type, radius, and mass based on parameters.
void initialise_particles(struct Parameters *p_parameters, struct Vectors *p_vectors)
{
    int *type = p_vectors->type;
    double *radius = p_vectors->radius;
    double *mass = p_vectors->mass;
    double rho = (p_parameters->density);
    double R_max = p_parameters->R_max;
    double R_min = p_parameters->R_min;
    size_t num_part = p_parameters->num_part;
    for (size_t i = 0; i < num_part; i++)
    {
        type[i] = 0;
        radius[i] = R_min + (R_max-R_min)*generate_uniform_random();
        double V = PI*(4.0/3.0)*(radius[i]*radius[i]*radius[i]);
        mass[i] = V*rho;
    }
}

// Main initialization function that seeds the random generator, sets time to zero,
// and calls sub-functions to initialize particles, positions, and velocities.
void initialise(struct Parameters *p_parameters, struct Vectors *p_vectors)
{
    srand(SEED); // Positive integer as seed for random number generator
    p_vectors->time = 0.0; // Initialize the time to zero
    initialise_particles(p_parameters, p_vectors);
    initialise_positions(p_parameters, p_vectors);
    initialise_velocities(p_parameters, p_vectors);
    return;
}

// Initializes particle positions on a cubic lattice based on the box dimensions and particle count.
// Positions are evenly spaced, and the lattice size is adjusted to accommodate the largest particle.
void initialise_positions(struct Parameters *p_parameters, struct Vectors *p_vectors)
{
    struct Vec3D dr;
    struct Index3D n;
    double dl;
    size_t ipart;
    size_t num_part = p_parameters->num_part;
    double *R = p_vectors->radius;
    double R_max = 0;
    for (ipart =0; ipart<num_part; ++ipart)
        R_max = (R[ipart]>R_max? R[ipart]: R_max);
    //dl = pow(p_parameters->L.x * p_parameters->L.y * p_parameters->L.z / ((double)p_parameters->num_part), 1.0 / 3.0);
    dl = 2.1*R_max;
    n.i = (int)floor(p_parameters->L.x / dl);
    n.j = (int)floor(p_parameters->L.y / dl);
    n.k = (int) (p_parameters->num_part/(n.i*n.j)+1);
    dr.x = p_parameters->L.x / (double)n.i;
    dr.y = p_parameters->L.y / (double)n.j;
    dr.z = dl;
    ipart = 0;
    for (size_t i = 0; i < n.i; ++i)
        for (size_t j = 0; j < n.j; ++j)
            for (size_t k = 0; k < n.k; ++k, ++ipart)
            {
                if (ipart >= num_part)
                    break;
                p_vectors->r[ipart].x = (i + 0.5) * dr.x;
                p_vectors->r[ipart].y = (j + 0.5) * dr.y;
                p_vectors->r[ipart].z = (k + 0.5) * dr.z;
                //      p_vectors->r[ipart].x = p_parameters->L.x*generate_uniform_random();
                //      p_vectors->r[ipart].y = p_parameters->L.y*generate_uniform_random();
                //      p_vectors->r[ipart].z = p_parameters->L.z*generate_uniform_random();
            }
}

void initialise_velocities(struct Parameters *p_parameters, struct Vectors *p_vectors)
{
    struct Vec3D sumv = (struct Vec3D){
        0.0};

    for (size_t i = 0; i < p_parameters->num_part; i++)
    {
        double sqrttgm = sqrt(p_parameters->Tg / p_vectors->mass[i]);
        p_vectors->v[i].x = sqrttgm * gauss();
        p_vectors->v[i].y = sqrttgm * gauss();
        p_vectors->v[i].z = sqrttgm * gauss();
        sumv.x += p_vectors->v[i].x;
        sumv.y += p_vectors->v[i].y;
        sumv.z += p_vectors->v[i].z;
    }

    sumv.x /= ((double)(p_parameters->num_part)); /* remove average velocity */
    sumv.y /= ((double)(p_parameters->num_part)); /* so total momentum is zero */
    sumv.z /= ((double)(p_parameters->num_part));
    for (size_t i = 0; i < p_parameters->num_part; i++)
    {
        p_vectors->v[i].x -= sumv.x;
        p_vectors->v[i].y -= sumv.y;
        p_vectors->v[i].z -= sumv.z;
    }
    for (size_t i = 0; i < p_parameters->num_part; i++)
    {
        p_vectors->omega[i] = (struct Vec3D){0.0};
    }
}
