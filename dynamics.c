#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "constants.h"
#include "structs.h"

// This function updates particle positions using their velocities.
// The positions are advanced by one full time step (dt), and displacement vectors
// (dr) for one time step are updated. The displacement since the last neighbor 
// list creation (stored in p_nbrlist->dr) is also updated for each particle.
void update_positions(struct Parameters *p_parameters, struct Nbrlist *p_nbrlist, struct Vectors *p_vectors)
{
    struct Vec3D dr_loc;
    struct Vec3D *r = p_vectors->r;     // Particle positions
    struct Vec3D *dr = p_vectors->dr;   // Displacement in one timestep
    struct Vec3D *v = p_vectors->v;     // Particle velocities
    struct DeltaR *dr_nbrlist = p_nbrlist->dr;  // Displacement since last neighbor list creation
    size_t num_part = p_parameters->num_part;
    double dt = p_parameters->dt;

    // Loop over all particles to update their positions
    for (size_t i = 0; i < num_part; i++)
    {
        dr_loc.x = v[i].x * dt;  // Compute displacement in x-direction for one timestep
        dr_loc.y = v[i].y * dt;  // Compute displacement in y-direction for one timestep
        dr_loc.z = v[i].z * dt;  // Compute displacement in z-direction for one timestep

        dr[i] = dr_loc;          // Store the displacement for this timestep
        r[i].x += dr_loc.x;      // Update position in x-direction
        r[i].y += dr_loc.y;      // Update position in y-direction
        r[i].z += dr_loc.z;      // Update position in z-direction

        // Update the displacement since last neighbor list creation
        dr_nbrlist[i].x += dr_loc.x;
        dr_nbrlist[i].y += dr_loc.y;
        dr_nbrlist[i].z += dr_loc.z;
        dr_nbrlist[i].sq = (dr_nbrlist[i].x) * (dr_nbrlist[i].x) + 
                           (dr_nbrlist[i].y) * (dr_nbrlist[i].y) + 
                           (dr_nbrlist[i].z) * (dr_nbrlist[i].z);  // Square of total displacement since last neighbor list creation
    }
}

void update_tangential_displacements(struct Parameters *p_parameters, struct Vectors *p_vectors,  struct Colllist *p_colllist)
{
/* Upon collision, if it is sticking, a tangential displacement builds up.
First the tangential displacement, vijt, is computed for pair ij.
This velocity is used to update the tangential displacement: tij(t+dt) = tij(t) + vijt(t+0.5*dt)*dt */
    size_t i,j,k;
    double fctr, dt = p_parameters->dt;
    const size_t num_nbrs = p_colllist->num_nbrs;
    struct DeltaR rij;
    struct Vec3D vij, vijn, vijt;
    struct Vec3D *v, *omega, *vw;
    struct Pair * nbr;
    struct DeltaR *tij;
    v = p_vectors->v;
    omega = p_vectors->omega;
    nbr = p_colllist->nbr;
    tij = p_colllist->tij;
    vw = p_colllist->vw;

    for (k=0; k < num_nbrs; k++)
    {
        // for each pair in collision list determine relative tangential velocity also taking into account rotation
        rij = nbr[k].rij;
        i = nbr[k].i;
        j = nbr[k].j;
        vij.x = v[i].x - v[j].x;
        vij.y = v[i].y - v[j].y;
        vij.z = v[i].z - v[j].z;
        fctr = (vij.x*rij.x+vij.y*rij.y+vij.z*rij.z)/rij.sq;
        vijn.x = fctr*rij.x;
        vijn.y = fctr*rij.y;
        vijn.z = fctr*rij.z;
        vijt.x = vij.x-vijn.x;
        vijt.y = vij.y-vijn.y;
        vijt.z = vij.z-vijn.z;
        vijt.x -= 0.5*((omega[i].y+omega[j].y)*rij.z-
                       (omega[i].z+omega[j].z)*rij.y);
        vijt.y -= 0.5*((omega[i].z+omega[j].z)*rij.x-
                       (omega[i].x+omega[j].x)*rij.z);
        vijt.z -= 0.5*((omega[i].x+omega[j].x)*rij.y-
                       (omega[i].y+omega[j].y)*rij.x);
        // for each pair in the collisionlist update relative tangential displacements
        tij[k].x += vijt.x*dt;
        tij[k].y += vijt.y*dt;
        tij[k].z += vijt.z*dt;
        tij[k].sq = tij[k].x *tij[k].x + tij[k].y*tij[k].y + tij[k].z*tij[k].z;
    }
    size_t num_w = p_colllist->num_w;
    size_t * indcs_w = p_colllist->indcs_w;
    struct DeltaR * riw = p_colllist->riw;
    struct DeltaR * tiw = p_colllist->tiw;
    for (j=0; j < num_w; ++j)
    {
        i = indcs_w[j];
        rij = riw[j];
        vij.x = v[i].x-vw[j].x;
        vij.y = v[i].y-vw[j].y;
        vij.z = v[i].z-vw[j].z;
        fctr = (vij.x*rij.x+vij.y*rij.y+vij.z*rij.z)/rij.sq;
        vijn.x = fctr*rij.x;
        vijn.y = fctr*rij.y;
        vijn.z = fctr*rij.z;
        vijt.x = vij.x-vijn.x;
        vijt.y = vij.y-vijn.y;
        vijt.z = vij.z-vijn.z;
        vijt.x -= (omega[i].y*rij.z-omega[i].z*rij.y);
        vijt.y -= (omega[i].z*rij.x-omega[i].x*rij.z);
        vijt.z -= (omega[i].x*rij.y-omega[i].y*rij.x);
        tiw[j].x += vijt.x*dt;
        tiw[j].y += vijt.y*dt;
        tiw[j].z += vijt.z*dt;
        tiw[j].sq = tiw[j].x *tiw[j].x + tiw[j].y*tiw[j].y + tiw[j].z*tiw[j].z;
    }
}

// This function updates particle velocities by half a time step using the current forces.
// The updated velocities are used in the velocity-Verlet integration scheme.
// The function also calculates and returns the kinetic energy of the system.
double update_velocities_half_dt(struct Parameters *p_parameters, struct Nbrlist *p_nbrlist, struct Vectors *p_vectors)
{
    double Ekin = 0.0, Ekin_rot = 0.0;
    double *R = p_vectors->radius;
    double *mass = p_vectors->mass;
    const double factor = 0.5 * p_parameters->dt;
    size_t num_part = p_parameters->num_part;
    struct Vec3D *v, *omg, *f, *T;
    v = p_vectors->v;
    omg = p_vectors->omega;
    f = p_vectors->f;
    T = p_vectors->T;
    for (size_t i = 0; i < p_parameters->num_part; i++)
    {
        double I = 0.4*mass[i]*R[i]*R[i];
        v[i].x += factor * f[i].x/mass[i];
        v[i].y += factor * f[i].y/mass[i];
        v[i].z += factor * f[i].z/mass[i];
        Ekin += mass[i]*(v[i].x * v[i].x + v[i].y * v[i].y + v[i].z * v[i].z);
        omg[i].x += factor*T[i].x/I;
        omg[i].y += factor*T[i].y/I;
        omg[i].z += factor*T[i].z/I;
        Ekin_rot += I*(omg[i].x*omg[i].x+omg[i].y*omg[i].y+omg[i].z*omg[i].z);
    }
    Ekin = 0.5* (Ekin + Ekin_rot);   
    return Ekin;
}

// This function applies periodic boundary conditions to ensure particles stay inside the simulation box.
// If a particle moves beyond the box, it is wrapped around to the opposite side.
void boundary_conditions(struct Parameters *p_parameters, struct Vectors *p_vectors)
{
    struct Vec3D invL;  // Inverse of the box size
    struct Vec3D *r = p_vectors->r;  // Particle positions
    struct Vec3D L = p_parameters->L;  // Box dimensions
    size_t num_part = p_parameters->num_part;  // Number of particles

    invL.x = 1.0 / L.x;
    invL.y = 1.0 / L.y;
    invL.z = 1.0 / L.z;

    // Loop over all particles and apply periodic boundary conditions
    for (size_t i = 0; i < num_part; i++)
    {
        r[i].x -= L.x * floor(r[i].x * invL.x);  // Apply periodic boundary in x-direction
        r[i].y -= L.y * floor(r[i].y * invL.y);  // Apply periodic boundary in y-direction
        r[i].z -= L.z * floor(r[i].z * invL.z);  // Apply periodic boundary in z-direction
    }
}
