#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include<iostream>
#include<cmath>
#include<fstream>
#include<sstream>
#include<cstdlib>
#include<time.h>
#include<vector>

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

void read_input();

void read_potential();

void read_lattice(double mat[], string file_lattice);

void lattice(double mat[], int n);

void lattice_ordering(double mat[]);
  
void file_xyz(string El1, string El2, double mat[], ofstream &file_out, int &step, double &E);

double dist(const int part1, const int part2, double mat[][4]);

double E_rep(const double r_ij, const int el_i, const int el_j);

double E_att2(const double r_ij, const int el_i, const int el_j);

void calc_par_en();

void neigh_list_in(double mat[], int neigh_list[], int n_nl[]);

void neigh_list_update_auto(double mat_old[], double mat[], int neigh_list[], double energy[], double energy_att2_atom[], double energy_rep_atom[], int n_nl[], double energy_rep[], double energy_att2[]);

double energy_check(double mat[]);

void energy_in(double mat[], double energy[], int neigh_list[], double energy_att2_atom[], double energy_rep_atom[], double energy_rep[], double energy_att2[]);

double energy_tot(double energy[]);

void energy_update_sing_somm(double mat[], int atom, int neigh_list[], double energy[], double energy_att2_atom[], double energy_rep_atom[], double energy_rep[], double energy_att2[]);

//void energy_null(int atom);

void shake(double mat[], double &E_tot, int neigh_list[], ofstream &file_out, double energy[], double energy_att2_atom[], double energy_rep_atom[], double energy_rep[], double energy_att2[]);

void exchange_random(double mat[], double &E_tot, int neigh_list[], ofstream &file_out, double energy[], double energy_att2_atom[], double energy_rep_atom[], double energy_rep[], double energy_att2[]);

void exchange_weighted(double mat[], double &E_tot, int neigh_list[], ofstream &file_out, double energy[], double energy_att2_atom[], double energy_rep_atom[], double energy_rep[], double energy_att2[]);

void mc_volume(double mat[], double &E_tot, int neigh_list[], double max_delta_V, ofstream &file_out, double energy[], double energy_att2_atom[], double energy_rep_atom[], double energy_rep[], double energy_att2[]);

void MC(double mat[], ofstream &file_images_out, ofstream &file_data_out,double energy[], double energy_att2_atom[], double energy_rep_atom[], int n_nl[], vector<double> &data, vector<double> &volume, int neigh_list[], double energy_rep[], double energy_att2[], int &MC_step_run);

void read_file_xyz(ifstream &file_in, double mat[][4]);

void write_file_out(int step, double &E, ofstream &file_data_out);

double average(vector<double> &data);

double VAR(vector<double> &data);

void convergence(vector<double> &data, vector<double> &volume, vector<vector<double>> &prop, double comp_el1, string name_file_conv, ofstream &file_prop, int &check_conv, double eps_min = 1e-8, double eps_max = 5e-11, int block_size = 100);

void change_composition(double mat[], double perc_f_el);

string getLastLineFromFile(ifstream &filename);

string format_val(double v, int n=3);
#endif
