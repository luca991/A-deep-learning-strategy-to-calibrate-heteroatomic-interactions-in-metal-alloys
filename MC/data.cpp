#include<vector>
#include<string>

using namespace std;

namespace data{
  double k_B = 8.617333262e-5;      // eV/K
  double e = 1.602176634e-19;       // Coulomb    eV = J/e
  double N_A = 6.02214076e23;       // 1/mol
  int N, n, pbc, n0, n1, file_xyz_yn;
  int length_nl;
  vector<double> perc_first_el;
  int n_shake = 0;
  int n_shake_rej = 0;
  int n_shake_tot = 0;
  int n_shake_rej_tot = 0;
  int n_vol_rej = 0;
  int n_vol = 0;
  int n_exch_rand_rej = 0;
  int n_exch_rand = 0;
  int n_exch_weigh_rej = 0;
  int n_exch_weigh = 0;
  string pot_file_name, file_out_data;
  char lattice_box, disp_in;
  string file_out_images; 
  int MC_steps, MC_step_imm;
  double T, R_max_shake, delta_NL;
  int seed;
  vector<int> MC_moves;
  vector<string> elements, file_par_origin;
  double **p;
  double **q;
  double **A;
  double **xi; 
  vector<double> E_coh;                      //eV
  vector<double> r_at;                       //Å
  vector<double> mass;                       //amu
  double **cut_min;                          //Å
  double **cut_max;                          //Å
  vector<vector<double>> r0;                         //Å (a0/sqrt(2))
  vector<double> a0;                         //Å
  vector<vector<double>> V0;                 //Å3
  double *L;
  double **x5;
  double **x4;
  double **x3;
  double **a5;
  double **a4;
  double **a3;
  //  int **neigh_list;
  int **neigh_list_inv;
  /*  double **energy_rep;
      double **energy_att2;*/
  vector<vector<int>> neigh_list;
  vector<vector<double>> energy_att2;
  vector<vector<double>> energy_rep;
  double *energy_rep_atom;
  double *energy_att2_atom;
  double *energy;
  int count_up_nei = 0;
  int *atom_el0;
  int *atom_el1;
  //  vector<int> n_nl;
  int *n_nl;
  vector<int> n_cells_hcp;
}
