#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "constants.h"
#include "structs.h"
#include "nbrlist.h"
#include "forces.h"

// Compute all forces on particles
// This function returns the total potential energy of the system.
double calculate_forces(struct Parameters *p_parameters, struct Colllist *p_colllist, struct Vectors *p_vectors)
{
    double Epot = 0.0;
    struct Vec3D *f = p_vectors->f;
    struct Vec3D *T = p_vectors->T;
    struct Vec3D g;
    g.x = p_parameters->g.x;
    g.y = p_parameters->g.y;
    g.z = p_parameters->g.z;
    double * mass = p_vectors->mass;
    // initialize the forces an torques with the gravitational force
    const size_t num_part = p_parameters->num_part;
    for (size_t i = 0; i < num_part; i++)
    {
        f[i].x = mass[i]*g.x; /*initialize forces to gravitational force*/
        f[i].y = mass[i]*g.y;
        f[i].z = mass[i]*g.z;
        T[i] = (struct Vec3D){0.0, 0.0, 0.0};
    }
    Epot += calculate_forces_pp(p_parameters, p_colllist, p_vectors);
    Epot += calculate_forces_pw(p_parameters, p_colllist, p_vectors);
    return Epot;
}

// Compute all forces on particles die to particle-particle contacts
// The function implement a soft-sphere model and used a collision list  
// This function returns the potential energy of (the concervative part of) these interactions
double calculate_forces_pp(struct Parameters *p_parameters, struct Colllist *p_colllist, struct Vectors *p_vectors)
{
    double Epot = 0.0;
    const double k_n_pp = p_parameters->k_n_pp;
    const double eta_n_pp = p_parameters->eta_n_pp;
    const double k_t_pp = p_parameters->k_t_pp;
    const double eta_t_pp = p_parameters->eta_t_pp;
    const double fric_pp = p_parameters->fric_pp;
    struct Pair *nbr = p_colllist->nbr;
    struct DeltaR *tijs = p_colllist->tij;
    const size_t num_nbrs = p_colllist->num_nbrs;
    const size_t num_part = p_parameters->num_part;
    double * R = p_vectors->radius;
    struct Vec3D *f = p_vectors->f;
    struct Vec3D *T = p_vectors->T;
    struct Vec3D *v = p_vectors->v;
    struct Vec3D *omega = p_vectors->omega;
    double * mass = p_vectors->mass;
    double inv_mass_ref = 1.0/p_parameters->mass_ref;
    for (size_t k = 0; k < num_nbrs; k++)
    {
        // for each pair in the neighbor list compute the pair forces
        struct DeltaR rij = nbr[k].rij;
        size_t i = nbr[k].i;
        size_t j = nbr[k].j;
        double mass_factor = sqrt(2.0*mass[i]*mass[j]/(mass[i]+mass[j])*inv_mass_ref);
        // normal spring force
        double r = sqrt(rij.sq);
        double overlap = R[i]+R[j] - r;
        double fr = k_n_pp * overlap / r;
        struct DeltaR dfn;
        dfn.x = fr * rij.x;
        dfn.y = fr * rij.y;
        dfn.z = fr * rij.z;
        Epot += 0.5 * k_n_pp * overlap * overlap;
        // normal dashpot force
        struct Vec3D vij, vijn;
        vij.x = v[i].x - v[j].x;
        vij.y = v[i].y - v[j].y;
        vij.z = v[i].z - v[j].z;
        double fctr = (vij.x * rij.x + vij.y * rij.y + vij.z * rij.z) / rij.sq;
        vijn.x = fctr * rij.x;
        vijn.y = fctr * rij.y;
        vijn.z = fctr * rij.z;
        fr = -mass_factor*eta_n_pp;
        dfn.x += fr * vijn.x;
        dfn.y += fr * vijn.y;
        dfn.z += fr * vijn.z;
        // tangential spring force
        struct DeltaR tij = tijs[k];
        fr = -k_t_pp;
        struct DeltaR dft;
        dft.x = fr * tij.x;
        dft.y = fr * tij.y;
        dft.z = fr * tij.z;
        // tangential dashpot force
        struct Vec3D vijt;
        vijt.x = vij.x - vijn.x;
        vijt.y = vij.y - vijn.y;
        vijt.z = vij.z - vijn.z;
        vijt.x -= 0.5 * ((omega[i].y + omega[j].y) * rij.z -
                         (omega[i].z + omega[j].z) * rij.y);
        vijt.y -= 0.5 * ((omega[i].z + omega[j].z) * rij.x -
                         (omega[i].x + omega[j].x) * rij.z);
        vijt.z -= 0.5 * ((omega[i].x + omega[j].x) * rij.y -
                         (omega[i].y + omega[j].y) * rij.x);
        fr = -mass_factor*eta_t_pp;
        dft.x += fr * vijt.x;
        dft.y += fr * vijt.y;
        dft.z += fr * vijt.z;

        //If tangential force is too large then sliding takes place
        dfn.sq = dfn.x * dfn.x + dfn.y * dfn.y + dfn.z * dfn.z;
        dft.sq = dft.x * dft.x + dft.y * dft.y + dft.z * dft.z;
        if (dft.sq >= fric_pp * fric_pp * dfn.sq) //sliding
        {
            fr = fric_pp * sqrt(dfn.sq / dft.sq);
            dft.x *= fr;
            dft.y *= fr;
            dft.z *= fr;
            // when sliding set the tangential displacement such that the sticking force (nearly) equals the sliding force
            tijs[k].x = -dft.x / k_t_pp;
            tijs[k].y = -dft.y / k_t_pp;
            tijs[k].z = -dft.z / k_t_pp;
        }
        else
            Epot += 0.5 * k_t_pp * tij.sq;
        struct Vec3D df;
        df.x = dfn.x + dft.x;
        df.y = dfn.y + dft.y;
        df.z = dfn.z + dft.z;
        struct Vec3D dT;
        dT.x = 0.5 * (rij.z * dft.y - rij.y * dft.z);
        dT.y = 0.5 * (rij.x * dft.z - rij.z * dft.x);
        dT.z = 0.5 * (rij.y * dft.x - rij.x * dft.y);
        f[i].x += df.x;
        f[i].y += df.y;
        f[i].z += df.z;
        f[j].x -= df.x;
        f[j].y -= df.y;
        f[j].z -= df.z;
        T[i].x += dT.x;
        T[i].y += dT.y;
        T[i].z += dT.z;
        T[j].x += dT.x;
        T[j].y += dT.y;
        T[j].z += dT.z;
    }
    return Epot; 
}

// Compute forces on particles due to particle-wall contacts
// The function implement a soft-sphere model and used a collision list
// This function returns the potential energy of (the concervative part of) these interactions
double calculate_forces_pw(struct Parameters *p_parameters, struct Colllist *p_colllist, struct Vectors *p_vectors)
{
    double Epot = 0.0;
    double *k_n_pw = p_parameters->k_n_pw;
    double *eta_n_pw = p_parameters->eta_n_pw;
    double *k_t_pw = p_parameters->k_t_pw;
    double *eta_t_pw = p_parameters->eta_t_pw;
    double *fric_pw = p_parameters->fric_pw;
    size_t num_w = p_colllist->num_w;
    size_t *indcs_w = p_colllist->indcs_w;
    unsigned int *wall_id = p_colllist->wall_id;
    struct DeltaR *riw = p_colllist->riw;
    struct DeltaR *tiw = p_colllist->tiw;
    struct Vec3D *vw = p_colllist->vw;
    double * R = p_vectors->radius;
    struct Vec3D *f = p_vectors->f;
    struct Vec3D *T = p_vectors->T;
    struct Vec3D *v = p_vectors->v;
    struct Vec3D *omega = p_vectors->omega;
    double * mass = p_vectors->mass;
    double inv_mass_ref = 1.0/p_parameters->mass_ref;
    for (size_t j = 0; j < num_w; ++j)
    {
        unsigned int w_id = wall_id[j];
        size_t i = indcs_w[j];
        struct DeltaR rij = riw[j];
        //normal elastic force
        double r = sqrt(rij.sq);
        double overlap = R[i] - r;
        double mass_factor = sqrt(mass[i]*inv_mass_ref);
        double fr = k_n_pw[w_id] * overlap / r;
        struct DeltaR dfn;
        dfn.x = fr * rij.x;
        dfn.y = fr * rij.y;
        dfn.z = fr * rij.z;
        Epot += 0.5 * k_n_pw[w_id] * overlap * overlap;
        // normal dashpot force
        struct Vec3D vij;
        vij.x = v[i].x - vw[j].x;
        vij.y = v[i].y - vw[j].y;
        vij.z = v[i].z - vw[j].z;
        double fctr = (vij.x * rij.x + vij.y * rij.y + vij.z * rij.z) / rij.sq;
        struct Vec3D vijn;
        vijn.x = fctr * rij.x;
        vijn.y = fctr * rij.y;
        vijn.z = fctr * rij.z;
        fr = -mass_factor*eta_n_pw[w_id];
        dfn.x += fr * vijn.x;
        dfn.y += fr * vijn.y;
        dfn.z += fr * vijn.z;

        struct DeltaR tij = tiw[j];
        fr = -k_t_pw[w_id];
        struct DeltaR dft;
        dft.x = fr * tij.x;
        dft.y = fr * tij.y;
        dft.z = fr * tij.z;
        struct Vec3D vijt;
        vijt.x = vij.x - vijn.x;
        vijt.y = vij.y - vijn.y;
        vijt.z = vij.z - vijn.z;
        vijt.x -= (omega[i].y * rij.z - omega[i].z * rij.y);
        vijt.y -= (omega[i].z * rij.x - omega[i].x * rij.z);
        vijt.z -= (omega[i].x * rij.y - omega[i].y * rij.x);
        fr = -mass_factor*eta_t_pw[w_id];
        dft.x += fr * vijt.x;
        dft.y += fr * vijt.y;
        dft.z += fr * vijt.z;

        dfn.sq = dfn.x * dfn.x + dfn.y * dfn.y + dfn.z * dfn.z;
        dft.sq = dft.x * dft.x + dft.y * dft.y + dft.z * dft.z;
        if (dft.sq >= fric_pw[w_id] * fric_pw[w_id] * dfn.sq)
        {
            //If tangential force is too large then sliding takes place:
            double fr = fric_pw[w_id] * sqrt(dfn.sq / dft.sq);
            dft.x *= fr;
            dft.y *= fr;
            dft.z *= fr;
            // when sliding set the tangential displacement such that the sticking force (nearly) equals the sliding force
            tiw[j].x = -dft.x / k_t_pw[w_id];
            tiw[j].y = -dft.y / k_t_pw[w_id];
            tiw[j].z = -dft.z / k_t_pw[w_id];
            tiw[j].sq = tiw[j].x *tiw[j].x + tiw[j].y*tiw[j].y + tiw[j].z*tiw[j].z;
        }
        else
            Epot += 0.5 * k_t_pw[w_id] * tij.sq;
        f[i].x += dfn.x + dft.x;
        f[i].y += dfn.y + dft.y;
        f[i].z += dfn.z + dft.z;
        T[i].x += (rij.z * dft.y - rij.y * dft.z);
        T[i].y += (rij.x * dft.z - rij.z * dft.x);
        T[i].z += (rij.y * dft.x - rij.x * dft.y);
    }
    return Epot; 
}