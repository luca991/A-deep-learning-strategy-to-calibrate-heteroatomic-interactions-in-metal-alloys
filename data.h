#ifndef DATA_H
#define DATA_H

#include<vector>
#include<string>
using namespace std;

namespace data{
  extern double k_B;
  extern double e;
  extern double N_A;
  extern int N, n, pbc, n0, n1, file_xyz_yn;
  extern int length_nl;
  extern vector<double> perc_first_el;
  extern int n_shake;
  extern int n_shake_rej;
  extern int n_shake_tot;
  extern int n_shake_rej_tot;
  extern int n_vol_rej;
  extern int n_vol;
  extern int n_exch_rand_rej;
  extern int n_exch_rand;
  extern int n_exch_weigh_rej;
  extern int n_exch_weigh;
  extern string pot_file_name, file_out_data;
  extern char lattice_box, disp_in;
  extern string file_out_images; 
  extern int MC_steps, MC_step_imm;
  extern double T, R_max_shake, delta_NL;
  extern int seed;
  extern vector<int> MC_moves;
  extern vector<string> elements, file_par_origin;
  extern double **p;
  extern double **q;
  extern double **A;
  extern double **xi; 
  extern vector<double> E_coh;
  extern vector<double> r_at;
  extern vector<double> mass;
  extern double **cut_min;
  extern double **cut_max;
  extern vector<vector<double>> r0;
  extern vector<double> a0;
  extern vector<vector<double>> V0;               //Å3
  extern double *L;                        //Å
  extern double **x5;
  extern double **x4;
  extern double **x3;
  extern double **a5;
  extern double **a4;
  extern double **a3;
  //  extern int **neigh_list;
  extern int **neigh_list_inv;
  /*  extern double **energy_rep;
      extern double **energy_att2;*/
  extern vector<vector<int>> neigh_list;
  extern vector<vector<double>> energy_att2;
  extern vector<vector<double>> energy_rep;
  extern double *energy_rep_atom;
  extern double *energy_att2_atom;
  extern double *energy;
  extern int count_up_nei;
  extern int *atom_el0;
  extern int *atom_el1;
  //  extern vector<int> n_nl;
  extern int *n_nl;
  extern vector<vector<double>> mat_fix;
  extern vector<string> lattice_type;
  extern vector<int> n_cells_hcp;
}
#endif
