#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include<iostream>
#include<cmath>
#include<fstream>
#include<sstream>
#include<cstdlib>
#include<time.h>
#include<vector>
#include "data.h"
#include <gsl/gsl_vector.h>

using namespace std;

inline double sq(double x){
  return x*x;
}

inline double pow_3(double x){
  return x*x*x;
}

inline double pow_4(double x){
  return x*x*x*x;
}

inline double pow_5(double x){
  return x*x*x*x*x;
}

inline int rand_unif_int(int x_min, int x_max){
  return x_min + rand()%(x_max-x_min+1);
}

void read_input();

void read_potential();

void lattice(double mat[][4], int n, double a, int el, char lattice_t);

// void file_xyz(string El1, string El2, double mat[][4], ofstream &file_out);

double dist(const int part1, const int part2, double mat[][4]);

double E_rep(const double r_ij, const int el_i, const int el_j);

double E_att2(const double r_ij, const int el_i, const int el_j);

void calc_par_en();

void neigh_list_in(double mat[][4], int n_nl[]);

// void neigh_list_update_auto(double mat_old[], double mat[], int neigh_list[], double energy[], double energy_att2_atom[], double energy_rep_atom[], int n_nl[], double energy_rep[], double energy_att2[]);

double energy_check(double mat[][4]);

double energy_nl(double mat[][4], int neigh_list[][500], int n_nl[]);

void energy_in(double mat[], double energy[], int neigh_list[], double energy_att2_atom[], double energy_rep_atom[], double energy_rep[], double energy_att2[]);

double energy_tot(double energy[]);

void energy_update_sing_somm(double mat[], int atom, int neigh_list[], double energy[], double energy_att2_atom[], double energy_rep_atom[], double energy_rep[], double energy_att2[]);

void set_impurity(double mat[][4], int atom);

double energy_obj_imp(const gsl_vector *pos, void* params);

double energy_obj(double x, void* params);

double energy_obj_ca(const gsl_vector *x, void* params);

void relax(double mat[][4], int max_iter, double m_expected, double min, double max, double m_yp);

void relax_all(double mat[][4], double eps, double first_step, int iter);

void relax_ca(double mat[][4]);

string getLastLineFromFile(ifstream &filename);

string format_val(double v, int n=3);
#endif
