#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <stdlib.h>
#include "constants.h"
#include "structs.h"
#include "nbrlist.h"
#include "walls.h"

void alloc_celllist(struct Parameters *p_parameters, struct Celllist *p_celllist)
/* Allocate arrays needed to store the cell-linked-list data*/
{
    struct Index3D size_grid;
    size_t Mtot;
    const double rlist = p_parameters->r_cut + p_parameters->r_shell;
    const size_t num_part = p_parameters->num_part;
    size_grid.i = floor(p_parameters->L.x / rlist);
    size_grid.j = floor(p_parameters->L.y / rlist);
    size_grid.k = floor(p_parameters->L.z / rlist);
    p_celllist->size_grid.i = size_grid.i;
    p_celllist->size_grid.j = size_grid.j;
    p_celllist->size_grid.k = size_grid.k;
    Mtot = size_grid.i * size_grid.j * size_grid.k;
    p_celllist->head = (size_t *)malloc(Mtot * sizeof(size_t));
    p_celllist->num_cells_max = Mtot;
    p_celllist->num_cells = Mtot;
    p_celllist->num_part_max = num_part;
    p_celllist->particle2cell = (size_t *)malloc(num_part * sizeof(size_t));
    p_celllist->list = (size_t *)malloc(num_part * sizeof(size_t));
}

void free_celllist(struct Celllist *p_celllist)
/* Free arrays used for the cell-linked-list */
{
    free(p_celllist->head);
    p_celllist->head = NULL;
    free(p_celllist->particle2cell);
    p_celllist->particle2cell = NULL;
    free(p_celllist->list);
    p_celllist->list = NULL;
}

void build_celllist(struct Parameters *p_parameters, struct Vectors *p_vectors, struct Celllist *p_celllist)
/* Build the cell-linked-list */
{
    struct Index3D size_grid, indx;
    size_t Mtot, icell;
    const double rlist = p_parameters->r_cut + p_parameters->r_shell;
    struct Vec3D mL;
    const size_t num_part = p_parameters->num_part;
    size_t *particle2cell, *head, *celllist;
    struct Vec3D *r;

    size_grid.i = floor(p_parameters->L.x / rlist);
    size_grid.j = floor(p_parameters->L.y / rlist);
    size_grid.k = floor(p_parameters->L.z / rlist);
    p_celllist->size_grid.i = size_grid.i;
    p_celllist->size_grid.j = size_grid.j;
    p_celllist->size_grid.k = size_grid.k;
    Mtot = size_grid.i * size_grid.j * size_grid.k;
    if (Mtot > p_celllist->num_cells_max)
    {
        p_celllist->head = (size_t *)realloc(p_celllist->head, Mtot * sizeof(size_t));
        p_celllist->num_cells_max = Mtot;
    }
    if (num_part > p_celllist->num_part_max)
    {
        p_celllist->particle2cell = (size_t *)realloc(p_celllist->particle2cell, num_part * sizeof(size_t));
        p_celllist->list = (size_t *)realloc(p_celllist->list, num_part * sizeof(size_t));
        p_celllist->num_part_max = num_part;
    }
    mL.x = ((double)size_grid.i) / p_parameters->L.x;
    mL.y = ((double)size_grid.j) / p_parameters->L.y;
    mL.z = ((double)size_grid.k) / p_parameters->L.z;
    head = p_celllist->head;
    for (icell = 0; icell < Mtot; ++icell)
        head[icell] = SIZE_MAX;
    particle2cell = p_celllist->particle2cell;
    celllist = p_celllist->list;
    r = p_vectors->r;
    for (size_t i = (num_part - 1); i != SIZE_MAX; --i)
    // Note that within a cell particles will be ordered in a descending way in the cell-linked list
    {
        indx.i = floor(r[i].x * mL.x);
        indx.j = floor(r[i].y * mL.y);
        indx.k = floor(r[i].z * mL.z);
        icell = indx.i + size_grid.i * (indx.j + indx.k * size_grid.j);
        particle2cell[i] = icell;
        celllist[i] = head[icell];
        head[icell] = i;
    }
}

void alloc_nbrlist(struct Parameters *p_parameters, struct Nbrlist *p_nbrlist)
/* Allocate arrays needed to store the neighbor list */
{
    size_t num_part = p_parameters->num_part;
    double Ndouble = (double)num_part;
    double rlist = p_parameters->r_cut + p_parameters->r_shell;
    // Estimate the number of neighbors using the average density
    double Nnbr = 0.5 * Ndouble * (4.0 / 3.0 * PI * rlist * rlist * rlist) / (p_parameters->L.x * p_parameters->L.y * p_parameters->L.z);
    size_t num_nbrs_max = ceil(0.6 * Ndouble * (Nnbr + 2.0 * sqrt(Nnbr)));
    p_nbrlist->p_celllist = (struct Celllist *)malloc(sizeof(struct Celllist));
    alloc_celllist(p_parameters, p_nbrlist->p_celllist);
    p_nbrlist->num_nbrs_max = num_nbrs_max;
    p_nbrlist->nbr = (struct Pair *)malloc(num_nbrs_max * sizeof(struct Pair));
    p_nbrlist->nbr_tmp = (struct Pair *)malloc(num_nbrs_max * sizeof(struct Pair));
    p_nbrlist->dr = (struct DeltaR *)malloc(p_parameters->num_part * sizeof(struct DeltaR));
    p_nbrlist->nbr_cnt = (size_t *)malloc((num_part) * sizeof(size_t));
    //    p_nbrlist->nbr_cnt_tmp = (size_t *)malloc(num_part * sizeof(size_t));
}

void free_nbrlist(struct Nbrlist *p_nbrlist)
/* Free arrays use to store the neighbor list */
{
    free_celllist(p_nbrlist->p_celllist);
    free(p_nbrlist->p_celllist);
    p_nbrlist->p_celllist = NULL;
    free(p_nbrlist->nbr_cnt);
    p_nbrlist->nbr_cnt = NULL;
    //    free(p_nbrlist->nbr_cnt_tmp);
    //    p_nbrlist->nbr_cnt_tmp = NULL;
    free(p_nbrlist->nbr_tmp);
    p_nbrlist->nbr_tmp = NULL;
    free(p_nbrlist->nbr);
    p_nbrlist->nbr = NULL;
    free(p_nbrlist->dr);
    p_nbrlist->dr = NULL;
}

void build_nbrlist(struct Parameters *p_parameters, struct Vectors *p_vectors, struct Nbrlist *p_nbrlist)
/* Build the neighbor list */
{
    struct Index3D indx, indx_nbr;
    size_t icell, inbr;
    size_t num_nbrs, num_nbrs_max;
    const double rlist = p_parameters->r_cut + p_parameters->r_shell; /* the radius for inclusion in the list is r_cut + r_shell */
    const double rlist_sq = rlist * rlist;
    struct Vec3D ri;
    struct Vec3D *r = p_vectors->r;
    struct DeltaR rij;
    const struct DeltaR dr = {0.0, 0.0, 0.0, 0.0};
    size_t *head, *particle2cell, *celllist;
    const int nbr_indcs[13][3] = {{0, 0, 1}, {0, 1, -1}, {0, 1, 0}, {0, 1, 1}, {1, -1, -1}, {1, -1, 0}, {1, -1, 1}, {1, 0, -1}, {1, 0, 0}, {1, 0, 1}, {1, 1, -1}, {1, 1, 0}, {1, 1, 1}};
    size_t num_part = p_parameters->num_part;

    // First build a cell-linked-list
    build_celllist(p_parameters, p_vectors, p_nbrlist->p_celllist);

    /*  Use the cell-linked-list to build a neighbor list.
        PairNBs are included to the neighbor list of their distance is less than r_cut+r_shell. */
    struct Index3D size_grid = p_nbrlist->p_celllist->size_grid;
    num_nbrs_max = p_nbrlist->num_nbrs_max;
    num_nbrs = 0;
    struct Pair *nbr = p_nbrlist->nbr;
    size_t *nbr_cnt = p_nbrlist->nbr_cnt;
    particle2cell = p_nbrlist->p_celllist->particle2cell;
    head = p_nbrlist->p_celllist->head;
    celllist = p_nbrlist->p_celllist->list;
    for (size_t i = 0; i < num_part; ++i)
        nbr_cnt[i] = 0;
    for (size_t i = 0; i < num_part; ++i)
    {
        // find neigbors of particle i in its own cell
        ri = r[i];
        for (size_t j = celllist[i]; j != SIZE_MAX; j = celllist[j]) // note that j < i
        {
            rij.x = ri.x - r[j].x;
            rij.y = ri.y - r[j].y;
            rij.z = ri.z - r[j].z;
            rij.sq = rij.x * rij.x + rij.y * rij.y + rij.z * rij.z;
            if (rij.sq < rlist_sq)
            {
                if (num_nbrs >= num_nbrs_max)
                {
                    num_nbrs_max += 5 * p_parameters->num_part;
                    nbr = (struct Pair *)realloc(nbr, num_nbrs_max * sizeof(struct Pair));
                    p_nbrlist->nbr = nbr;
                    p_nbrlist->num_nbrs_max = num_nbrs_max;
                }
                if (j > i) //is expected to true unless new particles have been inserted
                {
                    nbr[num_nbrs].i = i;
                    nbr[num_nbrs].j = j;
                    ++(nbr_cnt[i]);
                }
                else //swap particles i and j
                {
                    nbr[num_nbrs].i = j;
                    nbr[num_nbrs].j = i;
                    rij.x = -rij.x;
                    rij.y = -rij.y;
                    rij.z = -rij.z;
                    ++(nbr_cnt[j]);
                }
                nbr[num_nbrs].rij = rij;
                num_nbrs++;
            }
        }

        // next find neighbors of particle i in 1 of its 13 neighboring cells
        icell = particle2cell[i];
        indx.i = icell % size_grid.i;
        icell = icell / size_grid.i;
        indx.j = icell % size_grid.j;
        indx.k = icell / size_grid.j;
        for (size_t k = 0; k < 13; ++k)
        {
            ri = r[i];
            indx_nbr.i = (indx.i + nbr_indcs[k][0]);
            indx_nbr.j = (indx.j + nbr_indcs[k][1]);
            indx_nbr.k = (indx.k + nbr_indcs[k][2]);
            // The if-statements below implement periodic boundary conditions
            if (indx_nbr.i == SIZE_MAX)
            {
                ri.x += p_parameters->L.x;
                indx_nbr.i = size_grid.i - 1;
            }
            else if (indx_nbr.i >= size_grid.i)
            {
                ri.x -= p_parameters->L.x;
                indx_nbr.i = 0;
            }
            if (indx_nbr.j == SIZE_MAX)
            {
                ri.y += p_parameters->L.y;
                indx_nbr.j = size_grid.j - 1;
            }
            else if (indx_nbr.j >= size_grid.j)
            {
                ri.y -= p_parameters->L.y;
                indx_nbr.j = 0;
            }
            if (indx_nbr.k == SIZE_MAX)
            {
                ri.z += p_parameters->L.z;
                indx_nbr.k = size_grid.k - 1;
            }
            else if (indx_nbr.k >= size_grid.k)
            {
                ri.z -= p_parameters->L.z;
                indx_nbr.k = 0;
            }
            inbr = indx_nbr.i + size_grid.i * (indx_nbr.j + indx_nbr.k * size_grid.j);
            for (size_t j = head[inbr]; j != SIZE_MAX; j = celllist[j])
            {
                rij.x = ri.x - r[j].x;
                rij.y = ri.y - r[j].y;
                rij.z = ri.z - r[j].z;
                rij.sq = rij.x * rij.x + rij.y * rij.y + rij.z * rij.z;
                if (rij.sq < rlist_sq)
                {
                    if (num_nbrs >= num_nbrs_max)
                    {
                        num_nbrs_max += 5 * p_parameters->num_part;
                        nbr = (struct Pair *)realloc(nbr, num_nbrs_max * sizeof(struct Pair));
                        p_nbrlist->nbr = nbr;
                        p_nbrlist->num_nbrs_max = num_nbrs_max;
                    }
                    if (j > i)
                    {
                        nbr[num_nbrs].i = i;
                        nbr[num_nbrs].j = j;
                        ++(nbr_cnt[i]);
                    }
                    else //swap particles i and j
                    {
                        nbr[num_nbrs].i = j;
                        nbr[num_nbrs].j = i;
                        rij.x = -rij.x;
                        rij.y = -rij.y;
                        rij.z = -rij.z;
                        ++(nbr_cnt[j]);
                    }
                    nbr[num_nbrs].rij = rij;
                    num_nbrs++;
                }
            }
        }
    }

    // The neighbor list is fully sorted. The entries are such that nbr.i < nbr.j. nbr.i is ascending and for fixed i, j is ascending.
    size_t cnt_sum = 0;
    for (size_t i = 0; i < num_part; ++i)
    {
        size_t tmp = nbr_cnt[i];
        nbr_cnt[i] = cnt_sum;
        cnt_sum += tmp;
    }
    p_nbrlist->nbr_tmp = (struct Pair *)realloc(p_nbrlist->nbr_tmp, num_nbrs_max * sizeof(struct Pair));
    struct Pair *nbr_tmp = p_nbrlist->nbr_tmp;
    for (size_t k = 0; k < num_nbrs; ++k)
    {
        size_t i = (nbr[k].i < nbr[k].j ? nbr[k].i : nbr[k].j);
        nbr_tmp[nbr_cnt[i]++] = nbr[k];
    }
    p_nbrlist->nbr = nbr_tmp;
    p_nbrlist->nbr_tmp = nbr;
    nbr = p_nbrlist->nbr;
    nbr_tmp = p_nbrlist->nbr_tmp;
    size_t k = 0;
    for (size_t i = 0; i < num_part; k = nbr_cnt[i], i++)
        qsort(nbr + k, nbr_cnt[i] - k, sizeof(struct Pair), cmp_sort_nbr);
    p_nbrlist->num_nbrs = num_nbrs;
    p_nbrlist->dr = (struct DeltaR *)realloc(p_nbrlist->dr, num_part * sizeof(struct DeltaR));
    for (size_t i = 0; i < num_part; ++i) /*initialize particle displacements (with respect to creation time) to zero */
        p_nbrlist->dr[i] = dr;
}

int cmp_sort_nbr(const void *p1, const void *p2)
{
    const struct Pair *p_nbr1 = p1;
    const struct Pair *p_nbr2 = p2;
    if (p_nbr1->j < p_nbr2->j)
        return -1;
    else
        return (p_nbr1->j - p_nbr2->j);
}

int update_nbrlist(struct Parameters *p_parameters, struct Vectors *p_vectors, struct Nbrlist *p_nbrlist)
/* Update the connecting vectors of the neighbor list and if needed rebuild it.*/
{
    const double dr_sq_max = 0.25 * (p_parameters->r_shell * p_parameters->r_shell);
    int isRebuild = 0;
    struct DeltaR rij;
    struct Pair *nbr = p_nbrlist->nbr;
    struct Vec3D *dr = p_vectors->dr;
    // The neighbor list needs to be rebuild if one of the particles has displaced more then 0.5*r_shell
    for (size_t i = 0; i < p_parameters->num_part; ++i)
        if ((p_nbrlist->dr[i].sq) > dr_sq_max)
        {
            isRebuild = 1;
        }
    if (isRebuild) // rebuild neighbor list
        build_nbrlist(p_parameters, p_vectors, p_nbrlist);
    else // If no rebuild is needed, update the values of the connecting vectors
    {
        for (size_t k = 0; k < p_nbrlist->num_nbrs; ++k)
        {
            size_t i = nbr[k].i;
            size_t j = nbr[k].j;
            rij = nbr[k].rij;
            rij.x += (dr[i].x - dr[j].x);
            rij.y += (dr[i].y - dr[j].y);
            rij.z += (dr[i].z - dr[j].z);
            rij.sq = rij.x * rij.x + rij.y * rij.y + rij.z * rij.z;
            p_nbrlist->nbr[k].rij = rij;
        }
    }
    return isRebuild;
}

void alloc_colllist(struct Parameters *p_parameters, struct Colllist *p_colllist)
{
    p_colllist->num_nbrs = 0;
    p_colllist->nbr = (struct Pair *)malloc(0);
    p_colllist->nbr_tmp = (struct Pair *)malloc(0);
    p_colllist->tij = (struct DeltaR *)malloc(0);
    p_colllist->tij_tmp = (struct DeltaR *)malloc(0);
    p_colllist->num_w = 0;
    size_t num_w_max = 0;
    p_colllist->num_w_max = num_w_max;
    p_colllist->indcs_w = (size_t *)malloc(num_w_max * sizeof(size_t));
    p_colllist->indcs_w_tmp = (size_t *)malloc(num_w_max * sizeof(size_t));
    p_colllist->wall_id = (unsigned int *)malloc(num_w_max * sizeof(unsigned int));
    p_colllist->wall_id_tmp = (unsigned int *)malloc(num_w_max * sizeof(unsigned int));
    p_colllist->riw = (struct DeltaR *)malloc(num_w_max * sizeof(struct DeltaR));
    p_colllist->tiw = (struct DeltaR *)malloc(num_w_max * sizeof(struct DeltaR));
    p_colllist->tiw_tmp = (struct DeltaR *)malloc(num_w_max * sizeof(struct DeltaR));
    p_colllist->vw = (struct Vec3D *)malloc(num_w_max * sizeof(struct Vec3D));
}

void update_colllist(struct Parameters *p_parameters, struct Vectors *p_vectors, struct Nbrlist *p_nbrlist, struct Colllist *p_colllist)
/* The collision list is distilled from the neighbor list.
   Besides this information it stores the tangential displacement vector.
   For pairs already in collision the tangential displacement is copied from the old collision list.
   For new collision pairs the tangential displacement is set to zero.*/
{
    struct Pair *nbr_coll_old, *nbr_coll;
    struct Pair *nbr = p_nbrlist->nbr;

    size_t num_nbrs = p_nbrlist->num_nbrs;
    nbr_coll_old = p_colllist->nbr;
    nbr_coll = (struct Pair *)realloc(p_colllist->nbr_tmp, num_nbrs * sizeof(struct Pair));
    p_colllist->nbr_tmp = nbr_coll_old;

//    const double r_overlap_sq = 4.0 * p_parameters->R * p_parameters->R;
    double * R = p_vectors->radius;
    size_t m = 0;
    for (size_t k = 0; k < num_nbrs; k++)
    {
        // filter each pair in the neighbor list. Include in the collision list only of there is overlap
        struct DeltaR rij = nbr[k].rij;
        size_t i = nbr[k].i;
        size_t j = nbr[k].j;
        double sumR = R[i]+R[j];
        if (rij.sq < (sumR*sumR)) /*pair distance overlap*/
        {
            nbr_coll[m] = nbr[k];
            ++m;
        }
    }
    num_nbrs = m;
    nbr_coll = (struct Pair *)realloc(nbr_coll, num_nbrs * sizeof(struct Pair));
    p_colllist->nbr = nbr_coll;
    size_t num_nbrs_old = p_colllist->num_nbrs;
    p_colllist->num_nbrs = num_nbrs;

    struct DeltaR t0 = {0.0, 0.0, 0.0, 0.0};
    struct DeltaR *tij_old = p_colllist->tij;
    struct DeltaR *tij = (struct DeltaR *)realloc(p_colllist->tij_tmp, num_nbrs * sizeof(struct DeltaR));
    p_colllist->tij_tmp = tij_old;
    p_colllist->tij = tij;
    for (m = 0; m < num_nbrs; ++m)
        tij[m] = t0;

    size_t k;
    for (m = 0, k = 0; m < num_nbrs && k < num_nbrs_old;)
        if (nbr_coll[m].i < nbr_coll_old[k].i)
            ++m;
        else if (nbr_coll_old[k].i < nbr_coll[m].i)
            ++k;
        else // nbr_coll[m].i == nbr_coll_old[k].i
        {
            if (nbr_coll[m].j < nbr_coll_old[k].j)
                ++m;
            else if (nbr_coll_old[k].j < nbr_coll[m].j)
                ++k;
            else //collision pair is in old collision list
            {
                tij[m] = tij_old[k];
                ++m;
                ++k;
            }
        }

    /* Determine the particles in collision with a wall */
    size_t num_part = p_parameters->num_part;
    size_t num_w_old = p_colllist->num_w;
    size_t *indcs_w_old = p_colllist->indcs_w;
    size_t *indcs_w = p_colllist->indcs_w_tmp;
    unsigned int *wall_id_old = p_colllist->wall_id;
    unsigned int *wall_id = p_colllist->wall_id_tmp;
    struct DeltaR *riw = p_colllist->riw;
    struct Vec3D *vw = p_colllist->vw;
    /*     indcs_w = (size_t *)realloc(indcs_w, num_part * sizeof(size_t));
    wall_id = (unsigned int *)realloc(wall_id, num_part * sizeof(unsigned int));
    riw = (struct DeltaR *)realloc(riw, num_part * sizeof(struct DeltaR));
    vw = (struct Vec3D *)realloc(vw, num_part * sizeof(struct Vec3D)); */
    struct Vec3D *r = p_vectors->r;
    unsigned int num_walls = p_parameters->num_walls;
    size_t num_w_max = p_colllist->num_w_max;
    k = 0;
    for (size_t i = 0; i < num_part; i++)
    {
        for (unsigned int j = 0; j < num_walls; ++j)
        {
            struct DeltaR riw_loc;
            struct Vec3D vw_loc;
            if (p_parameters->wall_function[j](p_parameters, R[i], &r[i], &riw_loc, &vw_loc))
            {
                if (k >= num_w_max)
                {
                    num_w_max += 10;
                    indcs_w = (size_t *)realloc(indcs_w, num_w_max * sizeof(size_t));
                    wall_id = (unsigned int *)realloc(wall_id, num_w_max * sizeof(unsigned int));
                    riw = (struct DeltaR *)realloc(riw, num_w_max * sizeof(struct DeltaR));
                    vw = (struct Vec3D *)realloc(vw, num_w_max * sizeof(struct Vec3D));
                }
                indcs_w[k] = i;
                wall_id[k] = j;
                riw[k] = riw_loc;
                vw[k] = vw_loc;
                ++k;
            }

        }
            
    }
    size_t num_w = k;
    /*     indcs_w = (size_t *)realloc(indcs_w, num_w * sizeof(size_t));
    wall_id = (unsigned int *)realloc(wall_id, num_w * sizeof(unsigned int)); */
    /*    if (num_w>0)
        printf("check\n");*/
    /*     riw = (struct DeltaR *)realloc(riw, num_w * sizeof(struct DeltaR));
    vw = (struct Vec3D *)realloc(vw, num_w * sizeof(struct Vec3D)); */
    struct DeltaR *tiw_old = p_colllist->tiw;
    struct DeltaR *tiw = p_colllist->tiw_tmp;
    tiw = (struct DeltaR *)realloc(tiw, num_w_max * sizeof(struct DeltaR));
    for (size_t k = 0; k < num_w; ++k)
        tiw[k] = (struct DeltaR){0};
    for (m = 0, k = 0; m < num_w && k < num_w_old;)
        if (indcs_w[m] < indcs_w_old[k])
            ++m;
        else if (indcs_w_old[k] < indcs_w[m])
            ++k;
        else // indcs_w[m] == indcs_w_old[k]
        {
            if (wall_id[m] < wall_id_old[k])
                ++m;
            else if (wall_id_old[k] < wall_id[m])
                ++k;
            else //collision pair is in old collision list
            {
                tiw[m] = tiw_old[k];
                ++m;
                ++k;
            }
        }
    p_colllist->num_w = num_w;
    p_colllist->num_w_max = num_w_max;
    p_colllist->indcs_w = indcs_w;
    p_colllist->wall_id = wall_id;
    p_colllist->riw = riw;
    p_colllist->tiw = tiw;
    p_colllist->vw = vw;
    p_colllist->indcs_w_tmp = (size_t *)realloc(indcs_w_old, num_w_max * sizeof(size_t));
    p_colllist->wall_id_tmp = (unsigned int *)realloc(wall_id_old, num_w_max * sizeof(unsigned int));
    p_colllist->tiw_tmp = (struct DeltaR *)realloc(tiw_old, num_w_max * sizeof(struct DeltaR));
}

void free_colllist(struct Colllist *p_colllist)
{
    free(p_colllist->nbr);
    free(p_colllist->nbr_tmp);
    free(p_colllist->tij);
    free(p_colllist->tij_tmp);
    free(p_colllist->indcs_w);
    free(p_colllist->indcs_w_tmp);
    free(p_colllist->wall_id);
    free(p_colllist->wall_id_tmp);
    free(p_colllist->riw);
    free(p_colllist->tiw);
    free(p_colllist->tiw_tmp);
    free(p_colllist->vw);
}