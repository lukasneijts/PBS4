#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "constants.h"
#include "memory.h"
#include "structs.h"

/* --- profile accumulator state (internal to this file) --- */
static size_t acc_num_bins_r = 50;
static size_t acc_num_bins_z = 50;
static double *acc_phi_r_sum = NULL;
static double *acc_phi_z_sum = NULL;
static double *acc_vol_bin_r_tot = NULL;
static double *acc_vol_bin_z_tot = NULL;
static size_t acc_sample_count = 0;

/* Initialize accumulators (call once before time loop) */
void profile_accumulators_init(struct Parameters *p_parameters)
{
    double R_cyl = p_parameters->R_cyl;
    double H = p_parameters->L.z;
    double dr = R_cyl / (double)acc_num_bins_r;
    double dz = H / (double)acc_num_bins_z;

    acc_phi_r_sum = calloc(acc_num_bins_r, sizeof(double));
    acc_phi_z_sum = calloc(acc_num_bins_z, sizeof(double));
    acc_vol_bin_r_tot = calloc(acc_num_bins_r, sizeof(double));
    acc_vol_bin_z_tot = calloc(acc_num_bins_z, sizeof(double));
    acc_sample_count = 0;

    if (!acc_phi_r_sum || !acc_phi_z_sum || !acc_vol_bin_r_tot || !acc_vol_bin_z_tot) {
        fprintf(stderr, "Error: failed to allocate profile accumulators\n");
        free(acc_phi_r_sum); free(acc_phi_z_sum); free(acc_vol_bin_r_tot); free(acc_vol_bin_z_tot);
        acc_phi_r_sum = acc_phi_z_sum = acc_vol_bin_r_tot = acc_vol_bin_z_tot = NULL;
        return;
    }

    for (size_t ir = 0; ir < acc_num_bins_r; ++ir) {
        double r1 = ir * dr;
        double r2 = (ir + 1) * dr;
        acc_vol_bin_r_tot[ir] = PI * (r2*r2 - r1*r1) * H;
    }
    for (size_t iz = 0; iz < acc_num_bins_z; ++iz) {
        double z1 = iz * dz;
        double z2 = (iz + 1) * dz;
        acc_vol_bin_z_tot[iz] = PI * R_cyl * R_cyl * (z2 - z1);
    }
}

/* Add one sample (call periodically during simulation) */
void profile_accumulators_add_sample(struct Parameters *p_parameters, struct Vectors *p_vectors)
{
    if (!acc_phi_r_sum || !acc_phi_z_sum || !acc_vol_bin_r_tot || !acc_vol_bin_z_tot) return;

    size_t num_part = p_parameters->num_part;
    double R_cyl = p_parameters->R_cyl;
    double H = p_parameters->L.z;
    double dr = R_cyl / (double)acc_num_bins_r;
    double dz = H / (double)acc_num_bins_z;

    double *vol_bin_r = calloc(acc_num_bins_r, sizeof(double));
    double *vol_bin_z = calloc(acc_num_bins_z, sizeof(double));
    if (!vol_bin_r || !vol_bin_z) { free(vol_bin_r); free(vol_bin_z); return; }

    for (size_t i = 0; i < num_part; ++i) {
        double x = p_vectors->r[i].x;
        double y = p_vectors->r[i].y;
        double z = p_vectors->r[i].z;
        double r_part = sqrt((x - 0.5 * p_parameters->L.x) * (x - 0.5 * p_parameters->L.x) +
                             (y - 0.5 * p_parameters->L.y) * (y - 0.5 * p_parameters->L.y));
        double R = p_vectors->radius[i];
        double vol = (4.0/3.0) * PI * R * R * R;

        int ir = (int)(r_part / dr);
        if (ir >= 0 && ir < (int)acc_num_bins_r)
            vol_bin_r[ir] += vol;

        int iz = (int)(z / dz);
        if (iz >= 0 && iz < (int)acc_num_bins_z)
            vol_bin_z[iz] += vol;
    }

    /* accumulate per-bin volume fraction */
    for (size_t ir = 0; ir < acc_num_bins_r; ++ir) {
        if (acc_vol_bin_r_tot[ir] > 0.0) {
            double phi = vol_bin_r[ir] / acc_vol_bin_r_tot[ir];
            acc_phi_r_sum[ir] += phi;
        }
    }
    for (size_t iz = 0; iz < acc_num_bins_z; ++iz) {
        if (acc_vol_bin_z_tot[iz] > 0.0) {
            double phi = vol_bin_z[iz] / acc_vol_bin_z_tot[iz];
            acc_phi_z_sum[iz] += phi;
        }
    }

    acc_sample_count++;

    free(vol_bin_r);
    free(vol_bin_z);
}

/* Write averaged profiles (call once after simulation ends) */
void profile_accumulators_write_average(struct Parameters *p_parameters)
{
    if (acc_sample_count == 0) {
        fprintf(stderr, "No profile samples to write (acc_sample_count == 0)\n");
        return;
    }
    /* radial */
    FILE *fr = fopen("data/dens_radial_avg.csv", "w");
    if (fr) {
        fprintf(fr, "radius,volume_fraction_avg\n");
        for (size_t ir = 0; ir < acc_num_bins_r; ++ir) {
            double r_center = ((ir + 0.5) * (p_parameters->R_cyl / (double)acc_num_bins_r))/p_parameters->R_cyl;
            double phi_avg = acc_phi_r_sum[ir] / (double)acc_sample_count;
            fprintf(fr, "%g,%g\n", r_center, phi_avg);
        }
        fclose(fr);
    } else {
        fprintf(stderr, "Error: cannot open data/dens_radial_avg.csv for writing\n");
    }

    /* axial */
    FILE *fz = fopen("data/dens_axial_avg.csv", "w");
    if (fz) {
        fprintf(fz, "height,volume_fraction_avg\n");
        for (size_t iz = 0; iz < acc_num_bins_z; ++iz) {
            double z_center = ((iz + 0.5) * (p_parameters->L.z / (double)acc_num_bins_z))/p_parameters->L.z;
            double phi_avg = acc_phi_z_sum[iz] / (double)acc_sample_count;
            fprintf(fz, "%g,%g\n", z_center, phi_avg);
        }
        fclose(fz);
    } else {
        fprintf(stderr, "Error: cannot open data/dens_axial_avg.csv for writing\n");
    }

    /* cleanup */
    free(acc_phi_r_sum); acc_phi_r_sum = NULL;
    free(acc_phi_z_sum); acc_phi_z_sum = NULL;
    free(acc_vol_bin_r_tot); acc_vol_bin_r_tot = NULL;
    free(acc_vol_bin_z_tot); acc_vol_bin_z_tot = NULL;
    acc_sample_count = 0;
}

void compute_profiles(struct Parameters *p_parameters, struct Vectors *p_vectors)
{
    size_t num_part = p_parameters->num_part;
    double R_cyl = p_parameters->R_cyl;
    double H = p_parameters->L.z;
    size_t num_bins_r = 25;
    size_t num_bins_z = 25;
    double dr = R_cyl / num_bins_r;
    double dz = H / num_bins_z;

    double *vol_bin_r = calloc(num_bins_r, sizeof(double));
    double *vol_bin_z = calloc(num_bins_z, sizeof(double));
    double *vol_bin_r_tot = calloc(num_bins_r, sizeof(double));
    double *vol_bin_z_tot = calloc(num_bins_z, sizeof(double));

    for (size_t i = 0; i < num_part; ++i) {
        double x = p_vectors->r[i].x;
        double y = p_vectors->r[i].y;
        double z = p_vectors->r[i].z;
        double r_part = sqrt((x - 0.5 * p_parameters->L.x) * (x - 0.5 * p_parameters->L.x) +
                             (y - 0.5 * p_parameters->L.y) * (y - 0.5 * p_parameters->L.y));
        double R = p_vectors->radius[i];
        double vol = (4.0/3.0) * PI * R * R * R;

        int ir = (int)(r_part / dr);
        if (ir >= 0 && ir < num_bins_r)
            vol_bin_r[ir] += vol;

        int iz = (int)(z / dz);
        if (iz >= 0 && iz < num_bins_z)
            vol_bin_z[iz] += vol;
    }

    for (size_t ir = 0; ir < num_bins_r; ++ir) {
        double r1 = ir * dr;
        double r2 = (ir + 1) * dr;
        vol_bin_r_tot[ir] = PI * (r2*r2 - r1*r1) * H;
    }
    for (size_t iz = 0; iz < num_bins_z; ++iz) {
        double z1 = iz * dz;
        double z2 = (iz + 1) * dz;
        vol_bin_z_tot[iz] = PI * R_cyl * R_cyl * (z2 - z1);
    }

    // Write radial profile as CSV
    FILE *fr = fopen("data/dens_radial.csv", "w");
    fprintf(fr, "radius,volume_fraction\n");
    for (size_t ir = 0; ir < num_bins_r; ++ir) {
        double r_center = ((ir + 0.5) * dr)/p_parameters->R_cyl;
        double phi_r = vol_bin_r[ir] / vol_bin_r_tot[ir];
        fprintf(fr, "%g,%g\n", r_center, phi_r);
    }
    fclose(fr);

    // Write axial profile as CSV
    FILE *fz = fopen("data/dens_axial.csv", "w");
    fprintf(fz, "height,volume_fraction\n");
    for (size_t iz = 0; iz < num_bins_z; ++iz) {
        double z_center = ((iz + 0.5) * dz)/p_parameters->L.z;
        double phi_z = vol_bin_z[iz] / vol_bin_z_tot[iz];
        fprintf(fz, "%g,%g\n", z_center, phi_z);
    }
    fclose(fz);

    free(vol_bin_r);
    free(vol_bin_z);
    free(vol_bin_r_tot);
    free(vol_bin_z_tot);
}

/* Compute and write a center-count based radial profile.
   This assigns each particle to the annulus containing its center and
   converts counts to a volume fraction using the average particle volume. */
void compute_profiles_center_based(struct Parameters *p_parameters, struct Vectors *p_vectors)
{
    size_t num_part = p_parameters->num_part;
    double R_cyl = p_parameters->R_cyl;
    double H = p_parameters->L.z;
    size_t num_bins_r = 25;
    double dr = R_cyl / (double)num_bins_r;

    size_t *count_bin = calloc(num_bins_r, sizeof(size_t));
    if (!count_bin) return;

    /* accumulate counts and total solid volume */
    double total_solid_vol = 0.0;
    for (size_t i = 0; i < num_part; ++i) {
        double R = p_vectors->radius[i];
        total_solid_vol += (4.0/3.0) * PI * R * R * R;
        double x = p_vectors->r[i].x;
        double y = p_vectors->r[i].y;
        double r_part = sqrt((x - 0.5 * p_parameters->L.x) * (x - 0.5 * p_parameters->L.x) +
                             (y - 0.5 * p_parameters->L.y) * (y - 0.5 * p_parameters->L.y));
        int ir = (int)(r_part / dr);
        if (ir >= 0 && (size_t)ir < num_bins_r) count_bin[ir]++;
    }

    double avg_particle_vol = (num_part > 0) ? (total_solid_vol / (double)num_part) : 0.0;

    /* prepare geometric bin volumes (annulus * height) */
    double *vol_bin_geom = calloc(num_bins_r, sizeof(double));
    if (!vol_bin_geom) { free(count_bin); return; }
    for (size_t ir = 0; ir < num_bins_r; ++ir) {
        double r1 = ir * dr;
        double r2 = (ir + 1) * dr;
        vol_bin_geom[ir] = PI * (r2*r2 - r1*r1) * H;
    }

    /* write CSV (one file) */
    FILE *fc = fopen("data/dens_radial_center.csv", "w");
    if (fc) {
        fprintf(fc, "radius,count,volume_fraction_center\n");
        for (size_t ir = 0; ir < num_bins_r; ++ir) {
            double r_center = (ir + 0.5) * dr;
            double phi_center = 0.0;
            if (vol_bin_geom[ir] > 0.0) {
                phi_center = (double)count_bin[ir] * avg_particle_vol / vol_bin_geom[ir];
            }
            fprintf(fc, "%g,%zu,%g\n", r_center, count_bin[ir], phi_center);
        }
        fclose(fc);
    } else {
        fprintf(stderr, "Error: cannot open data/dens_radial_center.csv for writing\n");
    }

    free(count_bin);
    free(vol_bin_geom);
}

void characterize_final_pile(struct Parameters *parameters, struct Vectors *vectors)
{
    double cx = 0.5 * parameters->L.x;
    double cy = 0.5 * parameters->L.y;
    double h_max = -INFINITY;
    double r_max = 0.0;
    bool reset_file = parameters->reset_final_pile;

    for (int i = 0; i < parameters->num_part; ++i) {
        struct Vec3D ri = vectors->r[i];
        double dx = ri.x - cx;
        double dy = ri.y - cy;
        double r_xy = sqrt(dx*dx + dy*dy);
        double top_z = ri.z + parameters->R_max; /* conservative top using R_max */

        if (top_z > h_max) h_max = top_z;
        if (r_xy > r_max) r_max = r_xy;
    }

    double R_base = r_max + parameters->R_max; /* include particle radius at edge */
    double slope_rad = atan2(h_max, R_base);
    double slope_deg = slope_rad * 180.0 / PI;

    printf("\nFINAL PILE CHARACTERISATION\n");
    printf("  h_max   = %g (m)\n", h_max);
    printf("  R_base  = %g (m)\n", R_base);
    printf("  slope   = %g deg (%g rad)\n\n", slope_deg, slope_rad);

    if (reset_file) {
        FILE *fp = fopen("data/final_pile_characterisation.csv", "w");
        if (!fp) {
            fprintf(stderr, "Error: cannot open data/final_pile_characterisation.csv for writing\n");
        } 
        fprintf(fp, "h_max,R_base,slope_deg,slope_rad,num_particles,time_steps,radius\n");
        fprintf(fp, "%g,%g,%g,%g,%zu,%d,%g\n", h_max, R_base, slope_deg, slope_rad, (size_t)parameters->num_part, parameters->num_dt_steps, parameters->R_cyl);

        fclose(fp);
    } else {
        FILE *fp = fopen("data/final_pile_characterisation.csv", "a");
        if (!fp) {
            fprintf(stderr, "Error: cannot open data/final_pile_characterisation.csv for appending\n");
        } 
        fprintf(fp, "%g,%g,%g,%g,%zu,%d,%g\n", h_max, R_base, slope_deg, slope_rad, (size_t)parameters->num_part, parameters->num_dt_steps, parameters->R_cyl);

        fclose(fp);
    }
}

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