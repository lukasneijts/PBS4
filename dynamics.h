#ifndef DYNAMICS_H_
#define DYNAMICS_H_

/**
 * @brief Update particle positions by using velocities.
 * @param[in] p_parameters member: dt
 * @param[out] p_nbrlist used members: dr
 * @param[in,out] p_vectors members r, dr, v
 */
void update_positions(struct Parameters *p_parameters, struct Nbrlist *p_nbrlist, struct Vectors *p_vectors);

/**
 * @brief Update tangential displacements of particle pairs in collision
 * 
 * @param[in] p_parameters 
 * @param[in] p_vectors members v, omega 
 * @param[in,out] p_colllist members num_w, indcs_w, riw, tiw
 */
void update_tangential_displacements(struct Parameters *p_parameters, struct Vectors *p_vectors,  struct Colllist *p_colllist);

/**
 * @brief Update velocities for half a time step using forces.
 * @param[in] p_parameters used members: mass, dt
 * @param[in] p_nbrlist
 * @param[in, out] p_vectors used members: v, f
 * @return double kinetic energy
 */
double update_velocities_half_dt(struct Parameters *p_parameters, struct Nbrlist *p_nbrlist, struct Vectors *p_vectors);

/**
 * @brief Apply boundary conditions: particles folded back in periodic box.
 * @param[in] p_parameters used members: L
 * @param[in, out] p_vectors used members: r, dr
 */
void boundary_conditions(struct Parameters *p_parameters, struct Vectors *p_vectors);



#endif /* DYNAMICS_H_ */