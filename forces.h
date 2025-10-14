#ifndef FORCES_H_
#define FORCES_H_

/**
 * @brief Calculate forces and torques on particles
 * @param p_parameters
 * @param p_colllist
 * @param[out] p_vectors used members
 * @return double potential energy
 */
double calculate_forces(struct Parameters *p_parameters, struct Colllist *p_colllist, struct Vectors *p_vectors);

/**
 * @brief Calculate particle-particle forces and torques on particles
 * @param p_parameters
 * @param p_colllist
 * @param[out] p_vectors used members
 * @return double potential energy
 */
double calculate_forces_pp(struct Parameters *p_parameters, struct Colllist *p_colllist, struct Vectors *p_vectors);

/**
 * @briefCalculate particle-wall forces and torques on particles
 * @param p_parameters
 * @param p_colllist used members: num_nbrs, nbr, tij
 * @param[out] p_vectors used members: f
 * @return double potential energy
 */
double calculate_forces_pw(struct Parameters *p_parameters, struct Colllist *p_colllist, struct Vectors *p_vectors);


#endif /* FORCES_H_ */
