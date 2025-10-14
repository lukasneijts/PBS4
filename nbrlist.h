#ifndef NBRLIST_H_
#define NBRLIST_H_

#include <stddef.h>

/**
 * @brief Allocate arrays needed to store the cell-linked-list data
 *
 * @param p_parameters
 * @param p_cellist
 */
void alloc_celllist(struct Parameters *p_parameters, struct Celllist *p_celllist);

/** 
 * @brief Free arrays used for the cell-linked-list
 *
 * @param p_cellist
 */
void free_celllist(struct Celllist *p_cellist);

/**
 * @brief Build the cell-linked-list
 * 
 * @param p_parameters 
 * @param p_vectors 
 * @param p_celllist 
 */
void build_celllist(struct Parameters *p_parameters, struct Vectors * p_vectors, struct Celllist *p_celllist);

/**
 * @brief Allocate arrays needed to store the neighbor list 
 * 
 * @param p_parameters 
 * @param p_nbrlist 
 */
void alloc_nbrlist(struct Parameters *p_parameters, struct Nbrlist *p_nbrlist);

/**
 * @brief Free arrays use to store the neighbor list
 * 
 * @param p_nbrlist 
 */
void free_nbrlist(struct Nbrlist *p_nbrlist);

/**
 * @brief Build the neighbor list
 * 
 * @param p_parameters used members: rcut, rshell
 * @param p_vectors used members: r
 * @param p_nbrlist pointer to neighbor list
 */
void build_nbrlist(struct Parameters *p_parameters, struct Vectors *p_vectors, struct Nbrlist *p_nbrlist);

int cmp_sort_nbr(const void * p1, const void * p2);

/**
 * @brief Update the neighbor list
 * Checks if the neigbor lists needs to be rebuild be compairing the maximum displacement with rcut+rshell. 
 * If so build_nbrlist is called. If not, it updates positions and squared-distances of all pairs
 * 
 * @param p_parameters 
 * @param p_vectors 
 * @param p_nbrlist 
 * @return int Returns 1 if nbrlist is rebuild and 0 if it is only updated.
 */
int update_nbrlist(struct Parameters *p_parameters, struct Vectors *p_vectors, struct Nbrlist *p_nbrlist);

/**
 * @brief Allocate memory for the collision list
 * 
 * @param p_parameters 
 * @param p_colllist 
 */
void alloc_colllist(struct Parameters *p_parameters, struct Colllist *p_colllist);

/**
 * @brief Update the collision list from (the updated) neighbor list.
 * Uses the neigbor list to see if particle pairs are still in collision. If so the pairs are kept in the list. 
 * For newly detected collisions it adds a new entry in the collisoin list 
 * @param p_parameters 
 * @param p_vectors 
 * @param p_nbrlist 
 * @param p_colllist 
 */
void update_colllist(struct Parameters *p_parameters, struct Vectors *p_vectors, struct Nbrlist *p_nbrlist,  struct Colllist* p_colllist);

/**
 * @brief Free the memory allocated for the collision list
 * 
 * @param p_collist 
 */
void free_colllist(struct Colllist *p_collist);

#endif /* NBRLIST */
