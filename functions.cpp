#include <numeric>
#include <iomanip>
#include <algorithm>
#include <iterator>
#include <cstring>
#include "functions.h"
#include "data.h"
#include <gsl/gsl_multifit.h>
#include <omp.h>
#include <unistd.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_min.h>

void read_input(){
  fstream data_in;
  string data_in_name = "leggi.in";
  data_in.open(data_in_name);
  if(!data_in.is_open()){
    cerr << "Error opening 'leggi.in'" <<endl;
    exit(0);
  }else{
    cout<<data_in_name<<" opened successfully"<<endl<<endl;
  }
  string line;
  int temp;
  int N_line = 0;
  while(getline(data_in, line))
    N_line++;
  if(N_line!=27){
    cout<<"The input file 'leggi.in' is not in the correct format"<<endl;
    exit(0);
  }
  data_in.clear();
  data_in.seekg(0, ios::beg);
  for(int i = 0; i < 26; i++){
    getline(data_in, line);
    istringstream iss(line);
    if(i == 1) {
      iss >> data::pot_file_name;
    } else if(i == 2) {
      iss >> data::lattice_box;
    } else if(i == 3) {
      iss >> data::n;
      if(data::n==0){
        cerr<<"Incorrect box size"<<endl;
        exit(0);
      }      
    } else if(i == 4) {
      int temp;
      for(int j = 0; j<3; j++){
      iss >> temp;
      if(temp==0){
        cerr<<"Incorrect box size fo hcp"<<endl;
        exit(0);
      }   
      data::n_cells_hcp.push_back(temp);
      }   
    } else if(i == 5) {
      double a;
      while(iss >> a){
	      data::perc_first_el.push_back(a/100);
      }
    } else if(i == 5) {
      iss >> data::disp_in;
    } else if(i == 9) {
      iss >> data::MC_steps;
    } else if(i == 10) {
      iss >> data::T;
    } else if(i == 11) {
      iss >> data::pbc;
    }  else if(i == 12) {
      iss >> data::seed;
    } else if(i == 13) {
      iss >> data::R_max_shake;
    } else if(i == 14) {
      iss >> data::delta_NL;
    } else if(i == 17) {
      iss>>temp;
      data::MC_moves.push_back(temp);
    } else if(i == 18) {
      iss>>temp;
      data::MC_moves.push_back(temp);
    } else if(i == 19) {
      iss>>temp;
      data::MC_moves.push_back(temp);
    }  else if(i == 20) {
      iss>>temp;
      data::MC_moves.push_back(temp);
    } else if(i == 23) {
      iss >> data::file_out_images;
    } else if(i == 24) {
      iss >> data::file_out_data;
    } else if(i == 25){
      char temp;
      iss >> temp;
      if(temp == 'y')
	data::file_xyz_yn = 1;
      else if(temp == 'n')
	data::file_xyz_yn = 0;
      else{
	cout<<"It has not been specified whether you want the simulation image or not"<<endl;
	exit(0);
      }
    } else if(i == 26) {
      iss >> data::MC_step_imm;
    }
  }
  data_in.close();
}

void read_potential(){
  ifstream data_pot;
  data_pot.open(data::pot_file_name);
  if(!data_pot.is_open()){
    cout << "Error opening '"<<data::pot_file_name<<"'" <<endl;
    exit(0);
  }else{
    cout<<"'"<<data::pot_file_name<<"' opened successfully"<<endl<<endl;
  }
  string line, temp;
  getline(data_pot, line);
  getline(data_pot, line);
  getline(data_pot, line);
  istringstream iss(line);
  //while(iss >> temp){
  for(int i = 0; i<4; i++){
    if(i<2){
      string temp;
      iss>>temp;
      data::elements.push_back(temp);
    }else{
      string temp;
      iss>>temp;
      data::lattice_type.push_back(temp);
    }
  }
  data::cut_min = new double* [data::elements.size()];
  data::cut_max = new double* [data::elements.size()];
  data::p = new double* [data::elements.size()];
  data::q = new double* [data::elements.size()];
  data::A = new double* [data::elements.size()];
  data::xi = new double* [data::elements.size()];
  for(size_t i = 0; i<data::elements.size(); i++){
    data::cut_min[i] = new double[data::elements.size()];
    data::cut_max[i] = new double[data::elements.size()];
    data::p[i] = new double [data::elements.size()];
    data::q[i] = new double [data::elements.size()];
    data::A[i] = new double [data::elements.size()];
    data::xi[i] = new double [data::elements.size()];
  }
  getline(data_pot, line);
  iss.clear();
  iss.str(line);
  int j = 0;
  while(iss>>temp){
    if(j == 3){
      if(temp == "),"){
	cout<<"Error: potential parameter for the first elements not read correctly!"<<endl;
	exit(0);
      }
    }else if(j == 7){
      if(temp == "),"){
	cout<<"Error: potential parameter for the second elements not read correctly!"<<endl;
	exit(0);
      }
    }
    j++;
  }
  getline(data_pot, line);
  
  getline(data_pot, line);
  iss.clear();
  iss.str(line);
  int i = 0;
  while(iss>>temp){
    if(i<2){
      data::p[i][i] = stod(temp);
    }else if(i == 2){
      data::p[0][1] = stod(temp);
      data::p[1][0] = stod(temp);
    }
    i++;
  }
  getline(data_pot, line);
  iss.clear();
  iss.str(line);
  i = 0;
  while(iss>>temp){
    if(i<2){
      data::q[i][i] = stod(temp);
    }else if(i == 2){
      data::q[0][1] = stod(temp);
      data::q[1][0] = stod(temp);
    }
    i++;
  }

  getline(data_pot, line);
  iss.clear();
  iss.str(line);
  i = 0;
  while(iss>>temp){
    if(i<2){
      data::A[i][i] = stod(temp);
    }else if(i == 2){
      data::A[0][1] = stod(temp);
      data::A[1][0] = stod(temp);
    }
    i++;
  }

  getline(data_pot, line);
  iss.clear();
  iss.str(line);
  i = 0;
  while(iss>>temp){
    if(i<2){
      data::xi[i][i] = stod(temp);
    }else if(i == 2){
      data::xi[0][1] = stod(temp);
      data::xi[1][0] = stod(temp);
    }
    i++;
  }
  
  getline(data_pot, line);
  getline(data_pot, line);
  getline(data_pot, line);
  iss.clear();
  iss.str(line);
  size_t k = 0;
  do{
    iss>>temp;
    data::E_coh.push_back(stod(temp));
    k++;
  }while(k<data::elements.size());
  
  getline(data_pot, line);
  iss.clear();
  iss.str(line);
  k = 0;
  do{
    iss>>temp;
    data::r_at.push_back(stod(temp));
    k++;
  }while(k<data::elements.size()+1);
  getline(data_pot, line);
  iss.clear();
  iss.str(line);
  k = 0;
  do{
    iss>>temp;
    data::mass.push_back(stod(temp));
    k++;
  }while(k<data::elements.size());
  getline(data_pot, line);
  getline(data_pot, line);
  getline(data_pot, line);
  iss.clear();
  iss.str(line);
  i = 0;
  do{
    iss>>temp;
    data::cut_min[0][0] = stod(temp);
    i++;
    iss>>temp;
    data::cut_max[0][0] = stod(temp);
    i++;
  }while(i < 2);

  getline(data_pot, line);
  iss.clear();
  iss.str(line);
  i = 0;
  do{
    iss>>temp;
    data::cut_min[1][1] = stod(temp);
    i++;
    iss>>temp;
    data::cut_max[1][1] = stod(temp);
    i++;
  }while(i < 2);

  getline(data_pot, line);
  iss.clear();
  iss.str(line);
  i = 0;
  do{
    iss>>temp;
    data::cut_min[0][1] = stod(temp);
    data::cut_min[1][0] = stod(temp);
    i++;
    iss>>temp;
    data::cut_max[0][1] = stod(temp);
    data::cut_max[1][0] = stod(temp);
    i++;
  }while(i < 2);
  data_pot.close();

  data::r0.resize(data::elements.size(), vector<double>(data::elements.size()));
  for(size_t i = 0; i<data::elements.size(); i++){
    data::a0.push_back(2 * data::r_at[i]/(sqrt(2)));
    data::r0[i][i] = (data::r_at[i]);
  }
  data::r0[0][1] = data::r_at[2];
  data::r0[1][0] = data::r_at[2];

  data_pot.close();
}

void lattice(double mat[][4], int n, double a, int el, char lattice_t){
  if(lattice_t == 'f'){
    data::L = new double [3];
    for(int i = 0; i<3; i++){
      data::L[i] = a*n;
    }
    int ii = 0;
    for(int ix = 0; ix<data::n; ix++){
      for(int iy = 0; iy<data::n; iy++){
	for(int iz = 0; iz<data::n; iz++){
	  mat[ii][0] = el;
	  mat[ii][1] = ix*a;
	  mat[ii][2] = iy*a;
	  mat[ii][3] = iz*a;
	  ii++;

	  mat[ii][0] = el;
	  mat[ii][1] = (ix+0.5)*a;
	  mat[ii][2] = (iy+0.5)*a;
	  mat[ii][3] = iz*a;
	  ii++;

	  mat[ii][0] = el;
	  mat[ii][1] = (ix+0.5)*a;
	  mat[ii][2] = iy*a;
	  mat[ii][3] = (iz+0.5)*a;
	  ii++;

	  mat[ii][0] = el;
	  mat[ii][1] = ix*a;
	  mat[ii][2] = (iy+0.5)*a;
	  mat[ii][3] = (iz+0.5)*a;
	  ii++;
	}
      }
    }
  }else if(lattice_t == 'h'){
    data::L = new double [3];
    data::L[0] = a * data::n_cells_hcp[0];
    data::L[1] = sqrt(3) * a * data::n_cells_hcp[1];
    data::L[2] = 2 * sqrt(2./3) * a * data::n_cells_hcp[2];
		
    int ii = 0;
    for(int ix = 0; ix<data::n_cells_hcp[0]; ix++){
      for(int iy = 0; iy<data::n_cells_hcp[1]; iy++){
	for(int iz = 0; iz<2*data::n_cells_hcp[2]; iz++){
	  mat[ii][0] = el;
	  mat[ii][1] = ix*a + 0.5*a * (iz%2);
	  mat[ii][2] = iy*a*sqrt(3) + a*sqrt(1./12) * (iz%2);
	  mat[ii][3] = iz*a*sqrt(2./3);
	  ii++;

	  mat[ii][0] = el;
	  mat[ii][1] = (ix+0.5)*a + 0.5*a * (iz%2);
	  mat[ii][2] = (iy+0.5)*a*sqrt(3) + a*sqrt(1./12) * (iz%2);
	  mat[ii][3] = iz*a*sqrt(2./3);
	  ii++;
	}
      }
    }
  }else if(lattice_t == 'b'){
    data::L = new double [3];
    data::L[0] = a * data::n;
    data::L[1] = a * data::n;
    data::L[2] = a * data::n;
		
    int ii = 0;
    for(int ix = 0; ix<data::n; ix++){
      for(int iy = 0; iy<data::n; iy++){
	for(int iz = 0; iz<data::n; iz++){
	  mat[ii][0] = el;
	  mat[ii][1] = ix*a;
	  mat[ii][2] = iy*a;
	  mat[ii][3] = iz*a;
	  ii++;

	  mat[ii][0] = el;
	  mat[ii][1] = (ix+0.5)*a;
	  mat[ii][2] = (iy+0.5)*a;
	  mat[ii][3] = (iz+0.5)*a;
	  ii++;
	}
      }
    }
  }
}


double dist(const int part1, const int part2, double mat[][4]){
  double dx = mat[part1][1] - mat[part2][1];
  double dy = mat[part1][2] - mat[part2][2];
  double dz = mat[part1][3] - mat[part2][3];
  if (data::pbc == 1) {
    dx -= round(dx/data::L[0]) * data::L[0];
    dy -= round(dy/data::L[1]) * data::L[1];
    dz -= round(dz/data::L[2]) * data::L[2];
  }
  return sq(dx) + sq(dy) + sq(dz);
};

double E_rep(const double r_ij, const int el_i, const int el_j){
  double E_rep_val;
  if(r_ij > sq(data::cut_max[el_i][el_j])){
    E_rep_val = 0;
  }else if(r_ij > sq(data::cut_min[el_i][el_j])){
    E_rep_val = data::a5[el_i][el_j] * pow_5(sqrt(r_ij) - data::cut_max[el_i][el_j]) + data::a4[el_i][el_j] * pow_4(sqrt(r_ij) - data::cut_max[el_i][el_j]) + data::a3[el_i][el_j] * pow_3(sqrt(r_ij) - data::cut_max[el_i][el_j]);
  }else{
    E_rep_val = data::A[el_i][el_j] * exp(-data::p[el_i][el_j] * ((sqrt(r_ij)/data::r0[el_i][el_j]) - 1));
  }
  return E_rep_val;
}

double E_att2(const double r_ij, const int el_i, const int el_j){
  double E_att2_val;
  if(r_ij > sq(data::cut_max[el_i][el_j])){
    E_att2_val = 0;
  }else if(r_ij > sq(data::cut_min[el_i][el_j])){
    E_att2_val = sq(data::x5[el_i][el_j] * pow_5(sqrt(r_ij) - data::cut_max[el_i][el_j]) + data::x4[el_i][el_j] * pow_4(sqrt(r_ij) - data::cut_max[el_i][el_j]) + data::x3[el_i][el_j] * pow_3(sqrt(r_ij) - data::cut_max[el_i][el_j]));
  }else{
    E_att2_val = sq(data::xi[el_i][el_j]) * exp(-2 * data::q[el_i][el_j] * (sqrt(r_ij) / data::r0[el_i][el_j] - 1));
  }
  return E_att2_val;
}

void calc_par_en(){
  data::x5 = new double* [data::elements.size()];
  data::x4 = new double* [data::elements.size()];
  data::x3 = new double* [data::elements.size()];
  data::a5 = new double* [data::elements.size()];
  data::a4 = new double* [data::elements.size()];
  data::a3 = new double* [data::elements.size()];
  for(size_t i = 0; i<data::elements.size(); i++){
    data::x5[i] = new double[data::elements.size()];
    data::x4[i] = new double[data::elements.size()];
    data::x3[i] = new double[data::elements.size()];
    data::a5[i] = new double[data::elements.size()];
    data::a4[i] = new double[data::elements.size()];
    data::a3[i] = new double[data::elements.size()];
  }

  for(size_t i = 0; i<data::elements.size(); i++){
    for(size_t j = i; j<data::elements.size(); j++){
      double E_rep_0 = data::A[i][j] * exp(-data::p[i][j] * (data::cut_min[i][j]/data::r0[i][j] - 1));
      double E_rep_1 = - (data::p[i][j]/data::r0[i][j]) * data::A[i][j] * exp(-data::p[i][j] * (data::cut_min[i][j]/data::r0[i][j] - 1));
      double E_rep_2 = sq(data::p[i][j]/(data::r0[i][j])) * data::A[i][j] * exp(-data::p[i][j] * (data::cut_min[i][j]/data::r0[i][j] - 1));
      double E_att_0 = data::xi[i][j] * exp(-data::q[i][j] * (data::cut_min[i][j]/data::r0[i][j] - 1));
      double E_att_1 = - (data::q[i][j]/data::r0[i][j]) * data::xi[i][j] * exp(-data::q[i][j] * (data::cut_min[i][j]/data::r0[i][j] - 1));
      double E_att_2 = sq(data::q[i][j]/data::r0[i][j]) * data::xi[i][j] * exp(-data::q[i][j] * (data::cut_min[i][j]/data::r0[i][j] - 1)); 
      data::a5[i][j] = (12 * E_rep_0 - 6 * E_rep_1 * (data::cut_min[i][j] - data::cut_max[i][j]) + E_rep_2 * sq(data::cut_min[i][j] - data::cut_max[i][j]))/(2 * pow_5(data::cut_min[i][j] - data::cut_max[i][j]));
      data::a4[i][j] = (-15 * E_rep_0 + 7 * E_rep_1 * (data::cut_min[i][j] - data::cut_max[i][j]) - E_rep_2 * sq(data::cut_min[i][j] - data::cut_max[i][j]))/(pow_4(data::cut_min[i][j] - data::cut_max[i][j]));
      data::a3[i][j] = (20 * E_rep_0 - 8 * E_rep_1 * (data::cut_min[i][j] - data::cut_max[i][j]) + E_rep_2 * sq(data::cut_min[i][j] - data::cut_max[i][j]))/(2 * pow_3(data::cut_min[i][j] - data::cut_max[i][j]));

      data::x5[i][j] = (12 * E_att_0 - 6 * E_att_1 * (data::cut_min[i][j] - data::cut_max[i][j]) + E_att_2 * sq(data::cut_min[i][j] - data::cut_max[i][j]))/(2 * pow_5(data::cut_min[i][j] - data::cut_max[i][j]));
      data::x4[i][j] = (-15 * E_att_0 + 7 * E_att_1 * (data::cut_min[i][j] - data::cut_max[i][j]) - E_att_2 * sq(data::cut_min[i][j] - data::cut_max[i][j]))/(pow_4(data::cut_min[i][j] - data::cut_max[i][j]));
      data::x3[i][j] = (20 * E_att_0 - 8 * E_att_1 * (data::cut_min[i][j] - data::cut_max[i][j]) + E_att_2 * sq(data::cut_min[i][j] - data::cut_max[i][j]))/(2 * pow_3(data::cut_min[i][j] - data::cut_max[i][j]));
    }
  }
  data::a5[1][0] = data::a5[0][1];
  data::a4[1][0] = data::a4[0][1];
  data::a3[1][0] = data::a3[0][1];
  data::x5[1][0] = data::x5[0][1];
  data::x4[1][0] = data::x4[0][1];
  data::x3[1][0] = data::x3[0][1];
}

void neigh_list_in(double mat[][4], int n_nl[]){
  int length_nl = 500;
  double r_cut = max(data::cut_max[0][0], data::cut_max[1][1]);
  for(int i = 0; i<data::N; i++){
    for(int j = 0; j<length_nl; j++){
      data::neigh_list[i][j] = 0;
      //     data::neigh_list_inv[i][j] = 0;
    }
    n_nl[i] = 0;
  }

  //In this function we have renumbered the indices of the particles: instead of starting from 0 to N-1, the enumeration starts from 1 to N. Use this numeration in 'neigh_mat'
  double d;
  for(int i = 0; i<data::N-1; i++){
    for(int j = i+1; j<data::N; j++){
      d = dist(i, j, mat);
      if(d<sq(r_cut + data::delta_NL)){
	n_nl[i] = n_nl[i] + 1;
	n_nl[j] = n_nl[j] + 1;
	data::neigh_list[i][n_nl[i]] = j + 1;
	data::neigh_list[j][n_nl[j]] = i + 1;
      }
    }
  }
}

double energy_check(double mat[][4]){
  double E_temp1, E_temp2, d;
  double E_tot = 0;

  for(int i = 0; i<data::N; i++){
    E_temp1 = 0;
    E_temp2 = 0;
    for(int j = 0; j<data::N; j++){
      if(i!=j){
	d = dist(i, j, mat);
	E_temp1 += E_rep(d, mat[i][0], mat[j][0]);
	E_temp2 += E_att2(d, mat[i][0], mat[j][0]);
      }
    }
    E_tot += E_temp1 - sqrt(E_temp2);
  }
  return E_tot;
}

void set_impurity(double mat[][4], int atom){
  if(mat[atom][0]==1)
    mat[atom][0] = 0;
  else if(mat[atom][0] == 0)
    mat[atom][0] = 1;
}

double energy_obj_imp(const gsl_vector *pos, void* params) {
  double pos_pos[data::N][4];
  int *el = (int *)params;
  int k = 0;
  for (int i = 0; i < data::N; i++) {
    for (int j = 0; j < 4; j++) {
      if(j==0){
	pos_pos[i][j] = el[i];
      }else{
	pos_pos[i][j] = gsl_vector_get(pos, k);
	k++;
      }
    }
  }
  
  double E_temp1, E_temp2, d;
  double E_tot = 0;
  
  int i = 0;
  do{
    E_temp1 = 0;
    E_temp2 = 0;
    for(int j = 1; j<500; j++){
      int k = data::neigh_list[i][j];
      if(k==0){
	goto jump;
      }
      k = k-1;
      //			if(i!=j){
      d = dist(i, k, pos_pos);
      E_temp1 += E_rep(d, pos_pos[i][0], pos_pos[k][0]);
      E_temp2 += E_att2(d, pos_pos[i][0], pos_pos[k][0]);
      //			}
    }
  jump:
    E_tot += E_temp1 - sqrt(E_temp2);
    i++;
  }while(i<data::N);
  //  cout<<"En  min : "<<E_tot<<endl;
  
  return E_tot;
}

double energy_obj(double x, void* params) {
  cout<<"x: "<<x<<endl;
  double pos_pos[data::N][4];
  for (int i = 0; i < data::N; i++) {
    for (int j = 0; j < 4; j++) {
      if(j==0){
	pos_pos[i][j] = data::mat_fix[i][j];
      }else{
	pos_pos[i][j] = data::mat_fix[i][j] * x;
      }
    }
  }
  
  for(int i = 0; i<3; i++){
    data::L[i] = data::L[i]*x;
  }

  double E_temp1, E_temp2, d;
  double E_tot = 0;
  
  for(int i = 0; i<data::N; i++){
    E_temp1 = 0;
    E_temp2 = 0;
    for(int j = 0; j<data::N; j++){
      if(i!=j){
	d = dist(i, j, pos_pos);
	E_temp1 += E_rep(d, pos_pos[i][0], pos_pos[j][0]);
	E_temp2 += E_att2(d, pos_pos[i][0], pos_pos[j][0]);
      }
    }
    E_tot += E_temp1 - sqrt(E_temp2);
  }
  cout<<"En  min : "<<E_tot<<endl;
  	
  for(int i = 0; i<3; i++){
    data::L[i] = data::L[i]/x;
  }
  return E_tot;
}

double energy_obj_ca(const gsl_vector *x, void* params){
  double pos_pos[data::N][4];
  for (int i = 0; i < data::N; i++) {
    for (int j = 0; j < 4; j++) {
      if(j==0){
	pos_pos[i][j] = data::mat_fix[i][j];
      }else if(j == 3){
	pos_pos[i][j] = data::mat_fix[i][j] * gsl_vector_get(x, 1);
      }else{
	pos_pos[i][j] = data::mat_fix[i][j] * gsl_vector_get(x, 0);
      }
    }
  }
  for(int i = 0; i<2; i++){
    data::L[i] = data::L[i]*gsl_vector_get(x, 0);
  }
  data::L[2] = data::L[2]*gsl_vector_get(x, 1);

  double E_temp1, E_temp2, d;
  double E_tot = 0;

  for(int i = 0; i<data::N; i++){
    E_temp1 = 0;
    E_temp2 = 0;
    for(int j = 0; j<data::N; j++){
      if(i!=j){
	d = dist(i, j, pos_pos);
	E_temp1 += E_rep(d, pos_pos[i][0], pos_pos[j][0]);
	E_temp2 += E_att2(d, pos_pos[i][0], pos_pos[j][0]);
      }
    }
    E_tot += E_temp1 - sqrt(E_temp2);
  }
  cout<<"En  min : "<<E_tot<<endl;
  for(int i = 0; i<2; i++){
    data::L[i] = data::L[i]/gsl_vector_get(x, 0);
  }
  data::L[2] = data::L[2]/gsl_vector_get(x, 1);
  return E_tot;
}

void relax(double mat[][4], int max_iter, double m_expected, double min, double max, double m_yp){
  int status;
  int iter = 0;
  const gsl_min_fminimizer_type *T;
  gsl_min_fminimizer *s;
  gsl_function F;
  F.function = &energy_obj;
  F.params = 0;

  //T = gsl_min_fminimizer_brent;
  T = gsl_min_fminimizer_goldensection;
  //T = gsl_min_fminimizer_quad_golden;
  s = gsl_min_fminimizer_alloc (T);
  gsl_min_fminimizer_set (s, &F, m_yp, min, max);

  printf ("using %s method\n", gsl_min_fminimizer_name (s));

  printf ("%5s [%9s, %9s] %9s %10s %9s\n", "iter", "lower", "upper", "min", "err", "err(est)");

  printf ("%5d [%.7f, %.7f] %.7f %+.7f %.7f\n", iter, min, max, m_yp, m_yp - m_expected, max - min);

  do{
    iter++;
    status = gsl_min_fminimizer_iterate (s);

    m_yp = gsl_min_fminimizer_x_minimum (s);
    min = gsl_min_fminimizer_x_lower (s);
    max = gsl_min_fminimizer_x_upper (s);
    status = gsl_min_test_interval (min, max, 0.0000001, 0.0);
    if (status == GSL_SUCCESS)
      printf ("Converged:\n");

    printf ("%5d [%.7f, %.7f] " "%.7f %+.7f %.7f\n", iter, min, max, m_yp, m_yp - m_expected, max - min);
  }while (status == GSL_CONTINUE && iter < max_iter);
  gsl_min_fminimizer_free (s);
  for(int i = 0; i<data::N; i++){
    for(int j = 0; j<4; j++){
      if(j==0)
	mat[i][j] = data::mat_fix[i][j];
      else
	mat[i][j] = data::mat_fix[i][j]*m_yp;
    }
  }
  for(int i=0; i<3; i++){
    data::L[i]*=m_yp;
  }
}

void relax_all(double mat[][4], double eps, double first_step, int iter){
  size_t iter_m = 0;
  int status = 0;
  double size = 10;

  const gsl_multimin_fminimizer_type *T_m = gsl_multimin_fminimizer_nmsimplex2;
  //const gsl_multimin_fminimizer_type *T_m = gsl_multimin_fminimizer_nmsimplex2rand;
  gsl_multimin_fminimizer *s_m = NULL;
	
  //Atomic species in parameter vector
  int par_m[data::N];
  for(int i = 0; i<data::N; i++){
    par_m[i] = mat[i][0];
  }
	
  //Atomic position in gsl_vector
  gsl_vector *pos, *ss;
  pos = gsl_vector_alloc(3*data::N);
  int k = 0;
  for(int i = 0; i<data::N; i++){
    for(int j = 1; j<4; j++){
      gsl_vector_set(pos, k, mat[i][j]);
      k++;
    }
  }
  
  gsl_multimin_function minex_func;
  minex_func.n = data::N*3;
  minex_func.f = energy_obj_imp;
  minex_func.params = par_m;
	
  ss = gsl_vector_alloc (data::N*3);
  gsl_vector_set_all (ss, first_step);
	
  s_m = gsl_multimin_fminimizer_alloc (T_m, data::N*3);
  gsl_multimin_fminimizer_set (s_m, &minex_func, pos, ss);
  
  do{
    iter_m++;
    status = gsl_multimin_fminimizer_iterate(s_m);
    
    if (status)
      break;
    size = gsl_multimin_fminimizer_size (s_m);
    status = gsl_multimin_test_size (size, eps);
    
    if (status == GSL_SUCCESS){
      printf ("converged to minimum at\n");
    }
    if(iter_m%100 == 0){
      printf ("%5d %10.3e %10.3e f() = %7.8f size = %.3f\n", iter_m, gsl_vector_get (s_m->x, 0), gsl_vector_get (s_m->x, 1), s_m->fval, size);
    }
  }  while (status == GSL_CONTINUE && iter_m < iter);
  
  const gsl_vector *optimized_x = gsl_multimin_fminimizer_x(s_m);
  k = 0;
  for (int i = 0; i < data::N; i++) {
    for (int j = 0; j < 4; j++) {
      if(j==0){
	      mat[i][j] = par_m[i];
      }else{
        mat[i][j] = gsl_vector_get(optimized_x, k);
        k++;
      }
    }
  }
  gsl_vector_free(pos);
  gsl_vector_free(ss);
  gsl_multimin_fminimizer_free (s_m);
}

void relax_ca(double mat[][4]){
  gsl_vector *scale, *ss;
  scale = gsl_vector_alloc(2);
  gsl_vector_set(scale, 0, 1);
  gsl_vector_set(scale, 1, 1);

  size_t iter_m = 0;
  int status = 0;
  double size = 0;

  const gsl_multimin_fminimizer_type *T_m = gsl_multimin_fminimizer_nmsimplex2;
  //const gsl_multimin_fminimizer_type *T_m = gsl_multimin_fminimizer_nmsimplex2rand;
  gsl_multimin_fminimizer *s_m = NULL;

  gsl_multimin_function minex_func;
  minex_func.n = 2;
  minex_func.f = energy_obj_ca;
  minex_func.params = 0;


  ss = gsl_vector_alloc (2);
  gsl_vector_set_all (ss, 0.001);
  
  s_m = gsl_multimin_fminimizer_alloc (T_m, 2);
  gsl_multimin_fminimizer_set (s_m, &minex_func, scale, ss);
  
  do{
    iter_m++;
    status = gsl_multimin_fminimizer_iterate(s_m);
    
    if (status)
      break;
    size = gsl_multimin_fminimizer_size (s_m);
    status = gsl_multimin_test_size (size, 1e-4);
    
    if (status == GSL_SUCCESS){
      printf ("converged to minimum at\n");
    }
    if(iter_m%1000 == 0)
      printf ("%5d %10.3e %10.3e f() = %7.8f size = %.3f\n", iter_m, gsl_vector_get (s_m->x, 0), gsl_vector_get (s_m->x, 1), s_m->fval, size);
  }
  while (status == GSL_CONTINUE && iter_m < 1e4);
  
  const gsl_vector *optimized_x = gsl_multimin_fminimizer_x(s_m);
	
  for (int i = 0; i < data::N; i++) {
    for (int j = 0; j < 4; j++) {
      if(j==0){
	mat[i][j] = data::mat_fix[i][j];
      }else if(j==3){
	mat[i][j] = data::mat_fix[i][j] * gsl_vector_get(optimized_x, 1);
      }else{
	mat[i][j] = data::mat_fix[i][j] * gsl_vector_get(optimized_x, 0);
      }
    }
  }
  for(int i=0; i<2; i++){
    data::L[i]*=gsl_vector_get(optimized_x, 0);
  }
  data::L[2]*=gsl_vector_get(optimized_x, 1);
  //gsl_vector_free(pos);
  //gsl_vector_free(ss);
  //gsl_multimin_fminimizer_free (s_m);
}

string getLastLineFromFile(ifstream &filename) {
  std::string lastLine, line;
  while (std::getline(filename, line)) {
    if (!line.empty()) {
      lastLine = line;
    }
  }
    
  return lastLine;
};

string format_val(double v, int n){
    ostringstream oss;
    oss << fixed << setprecision(n) << v;
    return oss.str();
}
