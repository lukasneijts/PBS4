/******************************************************************************/ 
/*                                                                            */
/*  A DEM simulation with linear soft-sphere contact model                    */
/*                                                                            */
/*	This code is part of the course "Particle-based Simulations"              */
/*  taught at Eindhoven University of Technology.                             */
/*  No part of this code may be reproduced without permission of the author:  */
/*  Dr. Ir. E.A.J.F. Peters                                                   */
/*                                                                            */
/*  Dr. Ir. J.T. Padding:    version 1.1, 30/1/2013                           */
/*  Jeroen Hofman:           version 1.2, 28/7/2015                           */
/*  Dr. Ir. E.A.J.F. Peters: version 6.1, 12/10/2025                          */
/******************************************************************************/ 

/**
 * 2025 PBS Assignment: Granular Column (DEM) – Required Code Extensions
 * ---------------------------------------------------------------
 * This student starter code must be extended according to Assignment 4.
 * The work is organized in four task groups (A–D). Below each item we point
 * to the file(s) where you implement the required functionality and metrics.
 * Keep instrumentation lightweight and modular so physical logic stays clear.
 *
 * A. VERIFICATION (restitution)
 *  A1 Single particle – wall rebound:
 *     - Reduce to num_part = 1 in set_parameters.c (temporary test mode).
 *     - Give initial downward velocity v0 and measure rebound velocity v1.
 *     - Compute e_n_pw = v1 / v0 and compare with tabulated bottom wall value.
 *     - (Optional) Automate a loop over several v0 and print CSV lines.
 *  A2 Two-particle head-on collision:
 *     - Set num_part = 2; place particles along z (or x) just apart.
 *     - Assign opposite velocities ±v0; measure post-collision speeds.
 *     - Compute e_n_pp; compare with specified particle–particle restitution.
 *
 * B. SLOPE STABILITY (rotated gravity)
 *  B1 Create a shallow packed bed.
 *  B2 Instead of tilting geometry, rotate gravity: set g = |g|*(sin θ, 0, -cos θ).
 *  B3 Scan θ until sustained particle motion occurs (define criterion: e.g. average |v| over bed > threshold for N steps).
 *  Output suggestion: theta_deg, mean_speed, coordination_number.
 *  Code locations: modify gravity in set_parameters.c (add parameter slope_angle_deg), optionally update each step if sweeping.
 *
 * C. CYLINDRICAL COLUMN INITIALIZATION
 *  C1 Cylindrical wall: Implement cylindrical_wall(...) in walls.c
 *      - Center at (0.5 Lx, 0.5 Ly); radius R_cyl chosen so aspect ratio H/R_cyl ≈ 0.8 (adjust number of particles / height).
 *      - Overlap condition: radial distance rc = sqrt((x-xc)^2 + (y-yc)^2). Overlaps when rc + particle_radius > R_cyl.
 *      - Populate riw = r - r_w (r_w is closest point on cylinder surface) and vw = 0.
 *      - Register via wall_function[...] and increment num_walls.
 *  C2 Packing inside cylinder (initialise.c):
 *      (i) Place particles initially on a loose lattice clipped to cylinder.
 *      (ii) Randomization: For example, temporarily set g = {0,0,0} and run until a random condition is met.
 *      (iii) Settling: restore gravity; run until kinetic energy per particle < tolerance.
 *  C3 Profiles: After settling, compute radial φ(r) and axial φ(z) solids-volume fraction profiles.
 *      - Bin radii (0..R_cyl) and heights (0..H_final); accumulate particle volumes portionally by full inclusion (simplified).
 *      - Write to files: dens_radial.dat (r_center, phi_r) and dens_axial.dat (z_center, phi_z).
 *
 * D. COLUMN COLLAPSE
 *  D1 Wall removal: At collapse start (t=0 of collapse phase) remove cylindrical wall.
 *  D2 Time series: record trajectories at higher frequency during early spreading, then revert.
 *  D3 Metrics:
 *      - h_max: max (z_i + R_i) after collapse.
 *      - R_base: max radial center distance plus its radius where particles rest (define threshold velocity).
 *      - Slope angle: fit tangent line to upper surface points (e.g. upper 10% of height) or compute arctan(h_max / R_base_eff).
 *  D4 Compare qualitatively with Lajeunesse et al. (scaling, shape trends). Keep code comments short; put discussion in report.
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "constants.h"
#include "structs.h"
#include "setparameters.h"
#include "initialise.h"
#include "nbrlist.h"
#include "forces.h"
#include "dynamics.h"
#include "memory.h"
#include "fileoutput.h"
#include "walls.h"
#include <stdbool.h>

/**
 * @brief main The main of the DEM code. After initialization, 
 * a velocity-Verlet scheme is executed for a specified number of time steps.
 * 
 * @return int 0 if successful
 */
int main(void)
{
    struct Vectors vectors;
    struct Parameters parameters;
    struct Nbrlist nbrlist;
    struct Colllist colllist;
    size_t step = 0;
    double Ekin, Epot;

    set_parameters(&parameters);
    alloc_memory(&parameters, &vectors, &nbrlist, &colllist);
    if(parameters.load_restart == 1)
    {
        load_restart(&parameters, &vectors);
    }
    else
        initialise(&parameters, &vectors);

    build_nbrlist(&parameters, &vectors, &nbrlist);
    update_colllist(&parameters, &vectors, &nbrlist, &colllist);
    Epot = calculate_forces(&parameters, &colllist, &vectors);
    record_trajectories_xyz(1, &parameters, &vectors);

    /* initialize profile accumulators (sample frequency uses parameters.num_dt_traj below) */
    profile_accumulators_init(&parameters);

    /* --- automatic settling detector / collapse handling --- */
    bool cyl_removed = false;
    int settle_counter = 0;
    bool use_manual_removal = (parameters.collapse_start_step >= 0);
    int cyl_index = parameters.cyl_wall_index;
    /* ----------------------------------------------------- */

    while (step < parameters.num_dt_steps) //start of the velocity-Verlet loop
    {

        step++;
        vectors.time += parameters.dt;

        Ekin = update_velocities_half_dt(&parameters, &nbrlist, &vectors);
        update_positions(&parameters, &nbrlist, &vectors);

        update_tangential_displacements(&parameters, &vectors, &colllist);
        boundary_conditions(&parameters, &vectors);
        update_nbrlist(&parameters, &vectors, &nbrlist);
        update_colllist(&parameters, &vectors, &nbrlist, &colllist);
        Epot = calculate_forces(&parameters, &colllist, &vectors);
        Ekin = update_velocities_half_dt(&parameters, &nbrlist, &vectors);

        /* --- detect settling and remove cylindrical wall when settled --- */
        if (!cyl_removed) {
            if (use_manual_removal) {
                if ((int)step == parameters.collapse_start_step) {
                    if (cyl_index >= 0 && cyl_index < parameters.num_walls) {
                        if (parameters.num_walls > cyl_index) parameters.num_walls = cyl_index;
                        update_colllist(&parameters, &vectors, &nbrlist, &colllist);
                        cyl_removed = true;
                        printf("Cylindrical wall manually removed at step %lu\n", (long unsigned)step);
                    }
                }
            } else {
                double Ekin_per_particle = Ekin / (double) parameters.num_part;
                if (Ekin_per_particle < parameters.Ekin_tol) {
                    settle_counter++;
                } else {
                    settle_counter = 0;
                }
                if (settle_counter >= parameters.settle_pers_steps) {
                    if (cyl_index >= 0 && cyl_index < parameters.num_walls) {
                        if (parameters.num_walls > cyl_index) parameters.num_walls = cyl_index;
                        update_colllist(&parameters, &vectors, &nbrlist, &colllist);
                        cyl_removed = true;
                        printf("Cylindrical wall removed automatically at step %lu (Ekin/part=%g)\n", (long unsigned)step, Ekin_per_particle);
                    }
                }
            }
        }
        /* --------------------------------------------------------------- */

       if (step%parameters.num_dt_printf ==0) printf("Step %lu, Time %g, Z %g, Epot %g, Ekin %g, Etot %g\n", (long unsigned) step, vectors.time,
               2.0*((double) colllist.num_nbrs)/((double) parameters.num_part),
               Epot, Ekin, Epot+Ekin); //Z is the coordination number

        if (step%parameters.num_dt_traj ==0) {
            record_trajectories_xyz(0,&parameters,&vectors);
            /* also sample profiles at the same frequency (averaging) */
            profile_accumulators_add_sample(&parameters, &vectors);
        }
        
        if (step%parameters.num_dt_restart == 0) save_restart(&parameters,&vectors); 
    }

    // write averaged profiles computed over all samples
    profile_accumulators_write_average(&parameters);

    // characterize final pile after collapse
    characterize_final_pile(&parameters, &vectors);

    compute_profiles(&parameters, &vectors);
    compute_profiles_center_based(&parameters, &vectors);
    save_restart(&parameters,&vectors);
    free_memory(&vectors, &nbrlist, &colllist);

    return 0;
}
