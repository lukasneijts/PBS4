#ifndef INITIALISE_H_
#define INITIALISE_H_

/**
 * @brief Initialise calls initialise_particles, initialise_positions and initialise_velocities
 * 
 * @param p_parameters parameters used for initialization
 * @param p_vectors used members: r, v
 * @see initialize_particles, initialise_positions, initialise_velocities
 */
void initialise(struct Parameters *p_parameters, struct Vectors *p_vectors);

/**
 * @brief Initialises properties of particles
 * @param p_parameters
 * @param p_vectors used members: type, radius, mass
 * @see initialize
 */
void initialise_particles(struct Parameters *p_parameters, struct Vectors *p_vectors);


/**
 * @brief Initialises positions on a cubic lattice
 * 
 * @param p_parameters used members: L
 * @param p_vectors used members: r
 * @see initialize
 */
void initialise_positions(struct Parameters *p_parameters, struct Vectors *p_vectors);

/**
 * @brief Initialises positions inside a cylinder
 * 
 * @param p_parameters used members: L, R_cyl
 * @param p_vectors used members: r, radius
 * @see initialize
 */
void initialise_positions_cylinder(struct Parameters *p_parameters, struct Vectors *p_vectors);

/**
 * @brief Initialises velocities according to a Gaussian distribution
 * 
 * @param p_parameters used members: kT, mass
 * @param p_vectors used members: r
 * @see initialize
 */
void initialise_velocities(struct Parameters *p_parameters, struct Vectors *p_vectors);

#endif /* INITIALISE_H_ */
