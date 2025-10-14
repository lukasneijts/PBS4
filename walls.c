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