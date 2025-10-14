#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "constants.h"
#include "memory.h"
#include "structs.h"

void record_trajectories_xyz(int reset, struct Parameters *p_parameters, struct Vectors *p_vectors)
/*  Write the particle positions to a xyz file
    The filename (without extension) is given by p_parameters->filename_xyz.
    If reset = 1 the data is written to the file deleting data it possibly contained.
    If reset = 0 the data is appended. */
{
  FILE *fp_traj;
  char filename[1024];

  snprintf(filename, 1024, "%s%s", p_parameters->filename_xyz, ".xyz");
  if (reset == 1)
  {
    fp_traj = fopen(filename, "w");
  }
  else
  {
    fp_traj = fopen(filename, "a");
  }

  fprintf(fp_traj, "%lu\n", p_parameters->num_part);
  fprintf(fp_traj, "time = %f\n", p_vectors->time);
  struct Vec3D *r = p_vectors->r;
  struct Vec3D *v = p_vectors->v;
  double * R = p_vectors->radius;
  for (size_t i = 0; i < p_parameters->num_part; i++)
  {
      fprintf(fp_traj, "  C        %10.5f %10.5f %10.5f %10.5f %10.5f\n", r[i].x, r[i].y, r[i].z, 
          R[i], sqrt(v[i].x*v[i].x + v[i].y*v[i].y +v[i].z*v[i].z));
  }

  fclose(fp_traj);
}

void save_restart(struct Parameters *p_parameters, struct Vectors *p_vectors)
/* save arrays in vectors to binary file */
{
  FILE* p_file = fopen( p_parameters->restart_out_filename, "wb");
  size_t num_part = p_parameters->num_part;
  size_t sz = num_part*sizeof(struct Vec3D);
  fwrite(&p_vectors->time, sizeof(double), 1, p_file);
  fwrite(&num_part, sizeof(size_t), 1, p_file);
  fwrite(p_vectors->radius, num_part*sizeof(double), 1, p_file);
  fwrite(p_vectors->mass, num_part*sizeof(double), 1, p_file);
  fwrite(p_vectors->r, sz, 1, p_file);
  fwrite(p_vectors->v ,sz, 1, p_file);
  fwrite(p_vectors->omega, sz, 1, p_file);
  fwrite(p_vectors->f, sz, 1, p_file);
  fwrite(p_vectors->T, sz, 1, p_file);
  fclose(p_file);
}

void load_restart(struct Parameters *p_parameters, struct Vectors *p_vectors)
/* load arrays in vectors from binary file */
{
  FILE* p_file = fopen( p_parameters->restart_in_filename, "rb" );
  size_t num_part;
  fread(&p_vectors->time, sizeof(double), 1, p_file);
  fread(&num_part, sizeof(size_t), 1, p_file);
  size_t sz = num_part*sizeof(struct Vec3D);
  alloc_vectors(p_vectors,num_part);
  p_parameters->num_part = num_part;
  fread(p_vectors->radius, num_part*sizeof(double), 1, p_file);
  fread(p_vectors->mass, num_part*sizeof(double), 1, p_file);
  fread(p_vectors->r, sz, 1, p_file);
  fread(p_vectors->v, sz, 1, p_file);
  fread(p_vectors->omega, sz, 1, p_file);
  fread(p_vectors->f, sz, 1, p_file);
  fread(p_vectors->T, sz, 1, p_file);
  fclose(p_file);
}