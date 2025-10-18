#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <stdbool.h>
#include "constants.h"
#include "structs.h"
#include "nbrlist.h"

/**
 * @file walls.c
 * Wall interaction helper functions. Currently only a bottom and top planar wall
 * are implemented. For the 2025 assignment you must implement an additional
 * cylindrical wall that confines the initial granular column (Task C1) and later
 * remove it to study the collapse (Task D1).
 *
 * The wall function signature:
 *   bool wall(struct Parameters *p_parameters, double radius, struct Vec3D *r,
 *             struct DeltaR *riw, struct Vec3D *vw)
 * returns true if the particle at position r (center) of given radius overlaps
 * with the wall. On true you must set:
 *   riw -> vector (and squared length) from particle centre to closest wall point
 *   vw  -> local wall velocity at that point (zero for static walls).
 *
 * (Task C1) Implement a function `cylindrical_wall` that detects overlap with
 * a vertical cylinder of radius R_cyl centered in the domain (suggest choose
 * centre at (0.5 Lx, 0.5 Ly)). Provide riw pointing outward (from wall point to particle?)
 * consistent with existing planar walls: currently riw = (particle - wallpoint).
 * Set vw = {0,0,0}. Register this wall via parameters->wall_function[...] in set_parameters.c.
 * After initial packing and before collapse (Task D1) remove the cylinder by
 * removing its function pointer from parameters->wall_function[...] in set_parameters.c.
 */

bool cylindrical_wall(struct Parameters *p_parameters, double radius, struct Vec3D *r, struct DeltaR *riw, struct Vec3D *vw)
{
    // Cylinder center at (0.5*Lx, 0.5*Ly), vertical axis along z
    double cx = 0.5 * p_parameters->L.x;
    double cy = 0.5 * p_parameters->L.y;
    double R_cyl = p_parameters->R_cyl; // Cylinder radius from parameters

    // Distance from particle center to cylinder axis in xy-plane
    double dx = r->x - cx;
    double dy = r->y - cy;
    double dist_xy = sqrt(dx * dx + dy * dy);

    // Overlap if particle edge is outside cylinder
    if (dist_xy + radius > R_cyl) {
        // Closest point on cylinder wall in xy-plane
        double wall_x = cx + dx * (R_cyl / dist_xy);
        double wall_y = cy + dy * (R_cyl / dist_xy);

        // Vector from wall point to particle center
        double riw_x = r->x - wall_x;
        double riw_y = r->y - wall_y;
        double riw_z = 0; // Cylinder is vertical, so z is unchanged
        double riw_sq = riw_x * riw_x + riw_y * riw_y;

        *riw = (struct DeltaR){riw_x, riw_y, riw_z, riw_sq};
        *vw = (struct Vec3D){0, 0, 0}; // Stationary wall

        return true;
    }
    return false;
}

bool bottom_wall(struct Parameters *p_parameters, double radius, struct Vec3D *r, struct DeltaR *riw, struct Vec3D *vw)
{
    if (r->z < radius) //check if there is interaction with the bottom wall
    {
        double d = r->z;
        *riw = (struct DeltaR){0, 0, d, d * d};
        *vw = (struct Vec3D){0};
        return true;
    }
    else
    {
        return false;
    }
}

bool top_wall(struct Parameters *p_parameters, double radius, struct Vec3D *r, struct DeltaR *riw, struct Vec3D *vw)
{
    if (r->z > p_parameters->L.z - radius) //check if there is interaction with the top wall
    {
        double d = r->z - p_parameters->L.z;
        *riw = (struct DeltaR){0, 0, d, d * d};
        *vw = (struct Vec3D){0};
        return true;
    }
    else
    {
        return false;
    }
}

bool check_remove_cylindrical_wall(struct Parameters *p_parameters, double Ekin, size_t step,
                                   struct Vectors *vectors, struct Nbrlist *nbrlist,
                                   struct Colllist *colllist)
{
    size_t settle_counter;
    size_t cyl_index = p_parameters->cyl_wall_index;
    bool cyl_removed = false;
    int n_part = p_parameters->num_part;
    int settle_pers_steps = p_parameters->settle_pers_steps;
    int num_walls = p_parameters->num_walls;
    double Ekin_tol = p_parameters->Ekin_tol;


    double Ekin_per_particle = Ekin / n_part;
        if (Ekin_per_particle < Ekin_tol) {
            settle_counter++;
        } else {
            settle_counter = 0;
        }
        if (settle_counter >= settle_pers_steps) {
            if (cyl_index >= 0 && cyl_index < num_walls) {
                if (num_walls > cyl_index) num_walls = cyl_index;
                update_colllist(p_parameters, vectors, nbrlist, colllist);
                cyl_removed = true;
                printf("Cylindrical wall removed automatically at step %lu (Ekin/part=%g)\n", (long unsigned)step, Ekin_per_particle);
            }
        }
    return cyl_removed;
}