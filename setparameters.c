#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "constants.h"
#include "structs.h"
#include "walls.h"

void set_parameters(struct Parameters *p_parameters)
/* Set the parameters of this simulation */
{
  p_parameters->num_part = 2000;       // number of particles
  p_parameters->num_dt_steps = 10000;  // number of time steps
  p_parameters->density = 2500;        // mass density

  double R_min = 2.3e-3;             
  double R_max = R_min; 
  p_parameters->R_min = R_min;         // minimum particle radius
  p_parameters->R_max = R_max;         // maximum particle radius
  double kn = 1000;                    // spring stiffness of normal contact force
  double e_n_pp = 0.96, e_t_pp = 0.33, muf = 0.40; //restitution and friction coefficients for particle-particle interactions

  #define NUM_WALLS 3 // 2 walls are implements: bottom and top. Sides are periodic.
  p_parameters->num_walls = NUM_WALLS;           // number of walls in the system  
  p_parameters->wall_function[0] = bottom_wall;  // function used for wall 0                                                                           //velocity at bottom wall
  p_parameters->wall_function[1] = top_wall;     // function used for wall 1 
  p_parameters->wall_function[2] = cylindrical_wall;     // function used for wall 1 
  double e_n_pw[NUM_WALLS] = {0.96, 0.86, 0.86};  // normal restitution coefficients 
  double e_t_pw[NUM_WALLS] = {0.33, 0.33, 0.33};  // tangential restitution coefficients
  double muf_w[NUM_WALLS] = {0.40, 0.90, 0.15};   // friction coefficients 
  
  p_parameters->L = (struct Vec3D){4e1*R_min, 4e1*R_min, 1e2*R_max};                                                  //box size

  double g = 9.81;
  p_parameters->g = (struct Vec3D){0.0, 0.0, -g}; // gravitational acceleration                                                //gravitational acceleration vector
  p_parameters->num_dt_printf = 100;          // number of time steps between prints to screen
  p_parameters->num_dt_traj = 10;           //number of time steps between saves
  strcpy(p_parameters->filename_xyz, "trajectories");  //filename (without extension) for pdb file
  p_parameters->load_restart = 0;                      //if equal 1 restart file is loaded
  strcpy(p_parameters->restart_in_filename, "restart.dat");  //filename for loaded restart file
  p_parameters->num_dt_restart = 10000;                      // number of time steps between saves
  strcpy(p_parameters->restart_out_filename, "restart.dat"); //filename for saved restart file

  double mass_ref = p_parameters->density * (4.0 / 3.0) * PI * R_min * R_min * R_min; //mass_ref of a particle (later coefficients are corrected for real particle mass)
  double tcontact = sqrt(0.5 * mass_ref * (PI * PI + pow(log(e_n_pp), 2)) / kn); //computed contact time
  p_parameters->k_n_pp = kn;                                                     //normal elastic spring constant for particle-particle interactions
  p_parameters->eta_n_pp = -2.0 * log(e_n_pp) * (0.5 * mass_ref) / tcontact;     //normal dashpot damping constant for particle-particle interactions
  p_parameters->k_t_pp = (mass_ref / 7.0) * (PI * PI + pow(log(e_t_pp), 2)) / (tcontact * tcontact); //tangential elastic spring constant for particle-particle interactions
  p_parameters->eta_t_pp = -2.0 * log(e_t_pp) * (mass_ref / 7.0) / tcontact;     //tangential dashpot damping coeff. for particle-particle interactions
  p_parameters->fric_pp = muf;                                                   //friction coefficient for particle-particle interactions
  for(int i=0; i< p_parameters->num_walls; ++i)
  {
    p_parameters->k_n_pw[i] = mass_ref * (PI * PI + pow(log(e_n_pw[i]), 2)) / (tcontact * tcontact);               //normal elastic spring constant for particle-wall interactions
    p_parameters->eta_n_pw[i] = -2.0 * log(e_n_pw[i]) * (mass_ref) / tcontact;                                     //normal dashpot damping coeff. for particle-wall interactions
    p_parameters->k_t_pw[i] = (2.0 * mass_ref / 7.0) * (PI * PI + pow(log(e_t_pw[i]), 2)) / (tcontact * tcontact); //tangential elastic spring constant for particle-wall interactions
    p_parameters->eta_t_pw[i] = -2.0 * log(e_t_pw[i]) * (2.0 * mass_ref / 7.0) / tcontact;                         //tangential dashpot damping coeff. for particle-wall interactions
    p_parameters->fric_pw[i] = muf_w[i];
  }
  p_parameters->mass_ref = mass_ref; //mass_ref of a particle (later coefficients are corrected for real particle mass)
                                                                       
  double v_small = 1e-2 * R_max/tcontact;          //a velocity scale
  p_parameters->Tg = 0.5 * mass_ref * v_small * v_small; //here Tg denotes the granular temperature, an average kinetic energy used for initialization
  p_parameters->r_cut = 2.0 * R_max;               //cut-off distance for pair-par interactions
  p_parameters->dt = 0.05 * tcontact;              //integration time step
  p_parameters->r_shell = 0.2 * p_parameters->r_cut;             //shell thickness for neighbor list

  if (p_parameters->r_cut > p_parameters->L.x / 2.0)
      fprintf(stderr, "Warning! r_cut > Lx/2");
  if (p_parameters->r_cut > p_parameters->L.y / 2.0)
    fprintf(stderr, "Warning! r_cut > Ly/2");
  if (p_parameters->r_cut > p_parameters->L.z / 2.0)
    fprintf(stderr, "Warning! r_cut > Lz/2");

    
  // Parameters for cylindrical wall
  // p_parameters->H_R_ratio = 0.8;                        // height to radius ratio of cylindrical wall
  p_parameters->R_cyl = 0.5 * p_parameters->L.x;        // radius of cylindrical wall
  // p_parameters->L.x = 2.0 * p_parameters->R_cyl;         // box length in x direction
  // p_parameters->L.y = 2.0 * p_parameters->R_cyl;         // box length in y direction
  // p_parameters->L.z = p_parameters->H_R_ratio * p_parameters->R_cyl;     // height of cylindrical wall

}