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
    initialise_positions_cylinder(p_parameters, p_vectors);
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

// This function initializes particle positions within a cylindrical container.
void initialise_positions_cylinder(struct Parameters *p_parameters, struct Vectors *p_vectors)
{
    size_t ipart = 0;
    size_t num_part = p_parameters->num_part;
    double *R = p_vectors->radius;
    double R_max = 0;
    for (size_t i = 0; i < num_part; ++i)
        R_max = (R[i] > R_max ? R[i] : R_max);

    double cx = 0.5 * p_parameters->L.x;
    double cy = 0.5 * p_parameters->L.y;
    double R_cyl = p_parameters->R_cyl;
    double dz = 2.1 * R_max;
    double dy = 2.1 * R_max;
    double dx = dy * sqrt(3.0) / 2.0; // hexagonal lattice spacing

    size_t nz = (size_t)(p_parameters->L.z / dz);
    size_t ny = (size_t)(2.0 * R_cyl / dy);
    size_t nx = (size_t)(2.0 * R_cyl / dx);

    for (size_t k = 0; k < nz; ++k) {
        double z = (k + 0.5) * dz;
        for (size_t j = 0; j < ny; ++j) {
            double y = cy - R_cyl + (j + 0.5) * dy;
            for (size_t i = 0; i < nx; ++i) {
                // Offset every other row for hexagonal packing
                double x = cx - R_cyl + (i + 0.5) * dx + ((j % 2) ? dx / 2.0 : 0.0);

                double dx_c = x - cx;
                double dy_c = y - cy;
                double dist_xy = sqrt(dx_c * dx_c + dy_c * dy_c);

                if (dist_xy + R_max <= R_cyl && z <= p_parameters->L.z) {
                    if (ipart >= num_part) break;
                    p_vectors->r[ipart].x = x;
                    p_vectors->r[ipart].y = y;
                    p_vectors->r[ipart].z = z;
                    ipart++;
                }
            }
        }
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
