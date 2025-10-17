#ifndef TYPES_MD_H_
#define TYPES_MD_H_
#include <stdbool.h>
#include "constants.h"

/* This header file contains definitions of struct types used in the molecular dynamics code */

/**
 * @brief Struct to store i, j, k indices of a 3D grid
 * 
 */
struct Index3D
{
    size_t i, j, k; //!< 3 indices: i,j, k
};

/**
 * @brief Struct to store x, y, and z component of a 3D vector.
 * 
 */
struct Vec3D
{
    double x, y, z; //!< Three three coordinates of a 3D vector
};

/**
 * @brief Struct to store a 3D vector and its square length. This is expecially useful for connecting vectors in e.g. neighbor lists.
 * 
 */
struct DeltaR
/* Structure to store a 3D vector and its square length. */
{
    double x, y, z; //!< x, y and z coordinates
    double sq;      //!< square length
};

/**
 * @brief Struct to store all parameters. These parameters are set by the function @ref set_parameters.
 * 
 */
struct Parameters
{
    size_t num_part;       //!< Number of particles
    size_t num_dt_steps;   //!< Number of time steps
    double dt;             //!< integration time step
    struct Vec3D L;        //!< Box sizes in 3 direction
    double Tg;             //!< Granular temperature. Can be used to initialize velocities. 1.5Tg is the average kinetic energy per particle.
    double density;        //!< Density of particles
    double mass_ref;       //!< Reference mass used for collision parameters
    double R_max;          //!< Maximum sphere radius
    double R_min;          //!< Minumum sphere radius 
    struct Vec3D g;        //!< gravitational acceleration vector
    double k_n_pp;         //!< normal elastic spring constant for particle-particle interactions
    double eta_n_pp;       //!< normal dashpot friction constant for particle-particle interactions
    double k_t_pp;         //!< tangential elastic spring constant for particle-particle interactions
    double eta_t_pp;       //!< tangential dashpot damping coeff. for particle-particle interactions
    double fric_pp;        //!< friction coefficient for particle-particle interactions
    unsigned int num_walls;//!< number of solid walls in the system
    bool (*wall_function[NUM_WALLS_MAX])(struct Parameters *, double, struct Vec3D *, struct DeltaR *, struct Vec3D *); //!< 10 function pointers that can be used to define walls
    double k_n_pw[NUM_WALLS_MAX];    //!< normal elastic spring constant for particle-wall interactions
    double eta_n_pw[NUM_WALLS_MAX];  //!< normal dashpot damping coeff. for particle-wall interactions
    double k_t_pw[NUM_WALLS_MAX];    //!< tangential elastic spring constant for particle-wall interactions
    double eta_t_pw[NUM_WALLS_MAX];  //!< tangential damping coeff. for particle-wall interactions
    double fric_pw[NUM_WALLS_MAX];   //!< friction coeff. for particle-wall interactions
    double r_cut;                    //!< Cut-off distance for LJ interaction
    double r_shell;                  //!< Shell thickness for neighbor list
    size_t num_dt_printf;            //!< Number of time steps between prints to screen
    size_t num_dt_traj;              //!< Number of time steps between trajectory saves
    char filename_xyz[1024];         //!< filename (without extension) for pdb file
    char load_restart;               //!< if equal 1 restart file is loaded
    size_t num_dt_restart;           //!< Number of time steps between saves of restart file
    char restart_in_filename[1024];  //!< filename for loaded restart file
    char restart_out_filename[1024]; //!< filename for saved restart file

    double H_R_ratio;                //!< height to radius ratio of cylindrical wall
    double R_cyl;                    //!< radius of cylindrical wall
    
};

/**
 * @brief Struct with pointers to all particle arrays relevant for a MD simulation
 * 
 */
struct Vectors
{
    double time;         //!< time stamp of the vectors
    int    *type;        //!< type
    double *mass;        //!< masses of particles
    double *radius;      //!< radii of particles
    struct Vec3D *r;     //!< positions
    struct Vec3D *dr;    //!< displacements
    struct Vec3D *v;     //!< velocities
    struct Vec3D *omega; //!< angular-velocity */
    struct Vec3D *f;     //!< forces
    struct Vec3D *T;     //!< torques
    double *hist_vol_r;    // accumulated particle volume per radial bin
    double *hist_vol_z;    // accumulated particle volume per axial bin
    size_t hist_num_bins_r;
    size_t hist_num_bins_z;
    size_t hist_samples;   // number of samples accumulated
};

/**
 * @brief Struct to store a pair of particles: its indices and connecting vector
 * 
 */
struct Pair
{
    size_t i, j;       //!< indices of the two particles forming a pair
    struct DeltaR rij; //!< The connecting vector between the pairs rij = r[i]-r[j] corrected for periodicity
};

/**
 * @brief Struct used to store a cell-linked-list
 * 
 */
struct Celllist
{
    size_t *head;                                  //!< head[icell] provides the head the list for cell icell
    size_t *list;                                  //!< list[i] provides the next particle index in the cell-linked-list. list[i]==SIZE_MAX encodes the end of the list.
    size_t *particle2cell;                         //!< provides the cell index for a particle
    size_t num_cells, num_cells_max, num_part_max; //!< number of cells used and number of cells and particles allocated for
    struct Index3D size_grid;                      //!< number of cells in each direction
};

/**
 * @brief Struct to store a neighbor list
 * 
 */
struct Nbrlist
{
    struct Celllist *p_celllist;   //!< pointer to celllist used to create the neighbor list
    size_t num_nbrs, num_nbrs_max; //!< number of neighbors and maximum number allocated
    struct Pair *nbr, *nbr_tmp;    //!< list of neighbor pairs
    struct DeltaR *dr;             //!< displacements particles with respect to nbrlist creation time
    size_t *nbr_cnt;               //!< counts number of neighbors of i with j<i. Used for sorting.
    //    size_t *nbr_cnt_tmp;  //!< counts number of neighbors of i with j<i. Used for sorting.
};

/**
 * @brief struct to store collision list
 * The collision list stores particle pairs that are currently in collision, 
 * i.e., have a distance closer then the sum of radii. 
 * This list is needed next to the neighbor list mostly to keep track of the tangential displacement.
 * 
 */
struct Colllist
{
    size_t num_nbrs;               //!< number of pairs in collision list
    struct Pair *nbr;              //!< pairs in collision list
    struct Pair *nbr_tmp;          //!< collision list for internal use
    struct DeltaR *tij;            //!< tangential displacements of pairs in collision list
    struct DeltaR *tij_tmp;        //!< tangential displacements for internal use
    size_t num_w;                  //!< number of collisions with wall
    size_t num_w_max;              //!< maximum number of array members allocated
    size_t *indcs_w;               //!< particle indices that experience a wall collision
    size_t *indcs_w_tmp;           //!< particle indices for internal use
    unsigned int *wall_id;         //!< ID of the wall with which the particle collides
    unsigned int *wall_id_tmp;     //!< wall ID array for internal use
    struct DeltaR *riw;            //!< vectors pointing from particle center wall riw = ri-rw
    struct DeltaR *tiw;            //!< tangential displacement vector for wall collision
    struct DeltaR *tiw_tmp;        //!<  array with tangential displacements for internal use
    struct Vec3D *vw;              //!< local velocity of wall at collision point
};

/**
 * @brief Struct to store data for a histogram
 * 
 */
struct Histogram
{
    size_t nbins;               //!< Number of bins
    double min, max;            //!< Range in speed units
    double bin_width;           //!< Width of the bins
    size_t *counts;             //!< Counts per bin
    double *bin_centers;        //!< Center bins
    double total_counts;        //!< Total samples added
    double delta;               //!< Change in velocity
};

#endif /* TYPES_MD_H_ */
