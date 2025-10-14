#ifndef WALLS_H_
#define WALLS_H_

#include <stddef.h>
#include <stdbool.h>

/**
 * @brief Model for wall located at z=0.
 * This function can be used to model soft-sphere interactions with the bottom wall. 
 * 
 * @param[in] p_parameters member used: v_wall
 * @param[in] radius of a particle
 * @param[in] r position vector of a particle
 * @param[out] riw position of particle minus closest point on wall rw: riw = r-rw
 * @param[out] vw velocity of wall at point rw. Current implementation vw = p_parameters->v_wall
 * @return bool, true if particle and wall overlap false otherwise
 */
bool bottom_wall(struct Parameters *p_parameters, double radius, struct Vec3D *r, struct DeltaR *riw, struct Vec3D *vw);

/**
 * @brief Model for wall located at z=L.z.
 * This function can be used to model soft-sphere interactions with the top wall. 
 * 
 * @param[in] p_parameters member used: v_wall
 * @param[in] radius of a particle
 * @param[in] r position vector of a particle
 * @param[out] riw position of particle minus closest point on wall rw: riw = r-rw
 * @param[out] vw velocity of wall at point rw. Current implementation vw = 0.
 * @return bool, true if particle and wall overlap false otherwise
 */
bool top_wall(struct Parameters *p_parameters, double radius, struct Vec3D *r, struct DeltaR *riw, struct Vec3D *vw);

#endif  /* WALLS_H_ */