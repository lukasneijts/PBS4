#ifndef FILEOUTPUT_H_
#define FILEOUTPUT_H_

/**
 * @brief Output particle positions to xyz file
 * @param reset 1: write new file or overwrite existing file, 0: append data
 * @param p_parameters used members: filename_xyz
 * @param p_vectors used members: r
 * @param time time stamp
 */
void record_trajectories_xyz(int reset, struct Parameters * p_parameters, struct Vectors * p_vectors);

/**
 * @brief Save a restart file
 * 
 * @param p_parameters 
 * @param p_vectors 
 */
void save_restart(struct Parameters *p_parameters, struct Vectors *p_vectors);

/**
 * @brief Load a restart file
 * 
 * @param p_parameters 
 * @param p_vectors 
 */
void load_restart(struct Parameters *p_parameters, struct Vectors *p_vectors);

#endif /* FILEOUTPUT_H_ */
