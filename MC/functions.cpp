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

void read_input(){
  fstream data_in;
  // Opens the file "leggi.in", counts the lines and checks that it is written correctly
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
  // Read the file and save the data
  for(int i = 0; i < 27; i++){
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
    } else if(i == 4){
      int temp;
      for(int j = 0; j<3; j++){
	iss>>temp;
	if(temp ==0){
	  cerr<<"Incorrect box size fo hcp"<<endl;
	  exit(0);
	}
	data::n_cells_hcp.push_back(temp);
      }
    }else if(i == 5) {
      double a;
      while(iss >> a){
	data::perc_first_el.push_back(a/100);
      }
    } else if(i == 6) {
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
      iss>>temp;
      data::MC_step_imm = temp;
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
  int k = 0;
  do{
    string temp;
    iss >> temp;
    data::elements.push_back(temp);
    k++;
  }while(k<2);
  data::cut_min = new double* [data::elements.size()];
  data::cut_max = new double* [data::elements.size()];
  data::p = new double* [data::elements.size()];
  data::q = new double* [data::elements.size()];
  data::A = new double* [data::elements.size()];
  data::xi = new double* [data::elements.size()];
  for(int i = 0; i<data::elements.size(); i++){
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
  i = 0;
  do{
    iss>>temp;
    data::E_coh.push_back(stod(temp));
    i++;
  }while(i<data::elements.size());
	
  getline(data_pot, line);
  iss.clear();
  iss.str(line);
  i = 0;
  do{
    iss>>temp;
    data::r_at.push_back(stod(temp));
    i++;
  }while(i<data::elements.size()+1);

  getline(data_pot, line);
  iss.clear();
  iss.str(line);
  i = 0;
  do{
    iss>>temp;
    data::mass.push_back(stod(temp));
    i++;
  }while(i<data::elements.size());
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
  for(int i = 0; i<data::elements.size(); i++){
    //data::a0.push_back(4 * data::r_at[i]/(sqrt(2)));
    //data::r0[i][i] = (2 * data::r_at[i]);
    data::a0.push_back(2 * data::r_at[i]/(sqrt(2)));
    data::r0[i][i] = (data::r_at[i]);
  }
  data::r0[0][1] = data::r_at[2];
  data::r0[1][0] = data::r_at[2];

  data_pot.close();
}

void lattice(double mat[], int n){
  switch(data::lattice_box){
  case 'f':
    mat[4 * 0 + 0] = 0;
    mat[4 * 0 + 1] = 0;
    mat[4 * 0 + 2] = 0;
    mat[4 * 0 + 3] = 0;
    mat[4 * 1 + 0] = 0;
    mat[4 * 1 + 1] = data::a0[0]*0.5;
    mat[4 * 1 + 2] = 0;
    mat[4 * 1 + 3] = data::a0[0]*0.5;
    mat[4 * 2 + 0] = 0;
    mat[4 * 2 + 1] = data::a0[0]*0.5;
    mat[4 * 2 + 2] = data::a0[0]*0.5;
    mat[4 * 2 + 3] = 0;
    mat[4 * 3 + 0] = 0;
    mat[4 * 3 + 1] = 0;
    mat[4 * 3 + 2] = data::a0[0]*0.5;
    mat[4 * 3 + 3] = data::a0[0]*0.5;
    
    int i = 1;
    int s = 4;
    double a0;
    int count_a = 0;
    int count_b = 0;
    a0 = data::a0[0];
    data::L = new double [3];
    data::L[0] = a0 * n;
    data::L[1] = a0 * n;
    data::L[2] = a0 * n;
    
    //creation along x
    do{
      for(int j=0; j<4; j++){
	mat[4 * s + 0] = mat[4 * j + 0];
	mat[4 * s + 1] = mat[4 * j + 1] + (a0 * i);
	mat[4 * s + 2] = mat[4 * j + 2];
	mat[4 * s + 3] = mat[4 * j + 3];
	s++;
      }
      i++;
    }while(i<n);

    //creation along y
    int k = 0;
    int nx = 4*data::n;
    do{
      i=1;
      do{
        for(int j=k; j<k+4; j++){
	  mat[4 * s + 0] = mat[4 * j + 0];
	  mat[4 * s + 1] = mat[4 * j + 1];
	  mat[4 * s + 2] = mat[4 * j + 2] + (a0 * i);
	  mat[4 * s + 3] = mat[4 * j + 3];
	  s++;
	}
        i++;
      }while(i<n);
      k+=4;
    }while(k<n*4);

    //creation along z
    i = 1;
    int nxy = data::n * data::n * 4;
    do{
      for(int j=0; j<nxy; j++){
	mat[4 * s + 0] = mat[4 * j + 0];
	mat[4 * s + 1] = mat[4 * j + 1];
	mat[4 * s + 2] = mat[4 * j + 2];
	mat[4 * s + 3] = mat[4 * j + 3] + (a0 * i);
	s++;
      }
      i++;
    }while(i<n);
    break;
  }
}

void read_lattice(double mat[], string file_lattice){
  ifstream file;
  file.open(file_lattice);
  if(!file.is_open()){
    cout << "Error opening '"<<file_lattice<<"'" <<endl;
    exit(0);
  }else{
    cout<<"'"<<file_lattice<<"' opened successfully"<<endl<<endl;
  }
  string line;
  getline(file, line);
  getline(file, line);
  istringstream iss(line);
  vector<string>temp;
  string temp_s;
  while(iss>>temp_s){
    temp.push_back(temp_s);
  }
  data::L = new double [3];
  data::L[0] = stod(temp[1]);
  data::L[1] = stod(temp[5]);
  data::L[2] = stod(temp[9]);
  
  data::atom_el0 = new int [data::n0];
  data::atom_el1 = new int [data::n1];
  int El1 = 0;
  int El2 = 0;
  for(int i = 0; i<data::N; i++){
    getline(file, line);
    iss.clear();
    iss.str(line);
    int j = 0;
    while(iss>>temp_s){
      if(temp_s == data::elements[0]){
	mat[i*4+j] = 0;
	data::atom_el0[El1] = i;
	El1++;
	j++;
      }else if(temp_s == data::elements[1]){
	mat[i*4+j] = 1;
	data::atom_el1[El2] = i;
	El2++;
	j++;
      }else{
	mat[i*4+j] = stod(temp_s);
	j++;
      }
    }
  }
  cout<<El1<<endl;
  cout<<data::n0<<endl;
  if(El1 != data::n0){
    cerr<<"Error in lattice composition!"<<endl;
    exit(0);
  }
}

void lattice_ordering(double mat[]){
  double dist_min = sqrt(sq(data::r0[0][0])/2) * (1+0.01);
  switch(data::disp_in){
  case 'o':
    {
      int i = 0;
      int m = 0;
      int j = 0;
      int k = 0;
      int l = 0;
      if(data::perc_first_el[0] == 1){
	break;
      }
      do{
	if(mat[4 * i + 3] <= j * dist_min && mat[4 * i + 3] > (j-1) * dist_min){
	  if(mat[4 * i + 2] <= k * dist_min && mat[4 * i + 2] > (k-1) * dist_min){
	    if(j%2 == 0){
	      if(k%2==0){
		if(mat[4 * i + 1] <= l * (2*dist_min) && mat[4 * i + 1] > (l-1) * (2*dist_min)){
		  mat[4 * i + 0] = 1;
		  m++;
		  l++;
		}
	      }else{
		if(mat[4 * i + 1] <= (l+1) * (2*dist_min) && mat[4 * i + 1] > (l) * (2*dist_min)){
		  mat[4 * i + 0] = 1;
		  m++;
		  l++;
		}
	      }
	    }else{
	      if(k%2==1){
		if(mat[4 * i + 1] <= l * (2*dist_min) && mat[4 * i + 1] > (l-1) * (2*dist_min)){
		  mat[4 * i + 0] = 1;
		  m++;
		  l++;
		}
	      }else{
		if(mat[4 * i + 1] <= (l+1) * (2*dist_min) && mat[4 * i + 1] > (l) * (2*dist_min)){
		  mat[4 * i + 0] = 1;
		  m++;
		  l++;
		}
	      }
	    }
	    
	    if(m%data::n==0){
	      l=0;
	      i = 0;
	      k++;
	    }
	    if(m%(data::n * data::n *2) == 0){
	      j++;
	      k = 0;
	    }
	  }
	}
	i++;
      }while(m<data::n1);
      data::atom_el0 = new int [data::n0];
      data::atom_el1 = new int [data::n1];
      int count0 = 0;
      int count1 = 0;
      for(int j = 0; j<data::N; j++){
	if(mat[4 * j + 0] == 0){
	  data::atom_el0[count0] = j;
	  count0++;
	}else if(mat[4 * j + 0] == 1){
	  data::atom_el1[count1] = j;
	  count1++;
	}
      }
      break;
    }
  case 'r':
    {
      data::atom_el0 = new int [data::n0];
      data::atom_el1 = new int [data::n1];

      if(data::perc_first_el[0] == 1){
	/*for(int i = 0 ; i<data::N; i++){
	  mat[4*i] = 1;
	  }*/
	break;
      }
      int i = 0;
      srand(data::seed);
      do{
	int atom = rand() % data::N;
	if(mat[4 * atom + 0] != 1){
	  mat[4 * atom + 0] = 1;
	  i++;
	}
      }while(i<data::n1);

      int count0 = 0;
      int count1 = 0;
      for(int j = 0; j<data::N; j++){
	if(mat[4 * j + 0] == 0){
	  data::atom_el0[count0] = j;
	  count0++;
	}else if(mat[4 * j + 0] == 1){
	  data::atom_el1[count1] = j;
	  count1++;
	}
      }
      break;
    }
  }
}

void file_xyz(string El1, string El2, double mat[], ofstream &file_out, int &step, double &E){
  file_out<<"\t"<<data::N<<endl;
  //  file_out<<El1<<"\t"<<El2<<"\t"<<step<<"\t"<<E<<endl;
  file_out<<"Lattice=\" "<<data::L[0]<<" 0.0 0.0 0.0 "<<data::L[1]<<" 0.0 0.0 0.0 "<<data::L[2]<<" \" Properties=species:S:1:pos:R:3"<<endl;
  for(int i=0;i<data::N;i++){
    for(int j=0; j<4; j++){
      if(j==0){
	if(mat[4 * i + j]==0){
	  file_out<<El1<<"\t";
	}else if(mat[4 * i + j]==1){
	  file_out<<El2<<"\t";
	}else{
	  file_out<<"ERROR!!!!"<<endl;
	  break;
	}
      }else{
	file_out<<mat[4 * i + j]<<"\t";
      }
    }
    file_out<<endl;  
  }
}

double dist(const int part1, const int part2, double mat[]){
  double dx = mat[4 * part1 + 1] - mat[4 * part2 + 1];
  double dy = mat[4 * part1 + 2] - mat[4 * part2 + 2];
  double dz = mat[4 * part1 + 3] - mat[4 * part2 + 3];
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
  for(int i = 0; i<data::elements.size(); i++){
    data::x5[i] = new double[data::elements.size()];
    data::x4[i] = new double[data::elements.size()];
    data::x3[i] = new double[data::elements.size()];
    data::a5[i] = new double[data::elements.size()];
    data::a4[i] = new double[data::elements.size()];
    data::a3[i] = new double[data::elements.size()];
  }

  for(int i = 0; i<data::elements.size(); i++){
    for(int j = i; j<data::elements.size(); j++){
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

void neigh_list_in(double mat[], int neigh_list[], int n_nl[]){
  double r_cut = max(data::cut_max[0][0], data::cut_max[1][1]);
  for(int i = 0; i<data::N; i++){
    for(int j = 0; j<data::length_nl; j++){
      neigh_list[data::length_nl * i + j] = 0;
      //     data::neigh_list_inv[i][j] = 0;
    }
    n_nl[i] = 0;
  }
  for(int i = 0; i<data::N; i++){
    for(int j = 0; j<data::N; j++){
      data::neigh_list_inv[i][j] = 0;
    }
  }

  //In this function we have renumbered the indices of the particles: instead of starting from 0 to N-1, the enumeration starts from 1 to N. Use this numeration in 'neigh_mat'
  double d;
  for(int i = 0; i<data::N-1; i++){
    for(int j = i+1; j<data::N; j++){
      d = dist(i, j, mat);
      if(d<sq(r_cut + data::delta_NL)){
	n_nl[i] = n_nl[i] + 1;
	n_nl[j] = n_nl[j] + 1;
	neigh_list[data::length_nl * i + n_nl[i]] = j + 1;
	neigh_list[data::length_nl * j + n_nl[j]] = i + 1;
	data::neigh_list_inv[i][j] = n_nl[i];
	data::neigh_list_inv[j][i] = n_nl[j];
      }
    }
  }
}

void neigh_list_update_auto(double mat_old[], double mat[], int neigh_list[], double energy[], double energy_att2_atom[], double energy_rep_atom[], int n_nl[], double energy_rep[], double energy_att2[]){
  double max_displ = 0;
  double displ = 0;
  double sum2 = 0;
  int k;
  for(int i = 0; i<data::N; i++){
    double dr[3];
    sum2 = 0;
    k = 0;
    for(int j = 1; j<4; j++){
      dr[k] = mat_old[4 * i + j] - mat[4 * i + j];
      k++;
    }
    for(int j = 0; j<3; j++){
      sum2 += sq(dr[j]);
    }
    displ = (sum2);
    if(displ > max_displ){
      max_displ = displ;
    }
  }

  //  max_displ = sqrt(max_displ);
  // if(max_displ>0){
  if(max_displ>sq(0.2 * data::delta_NL)){
    for(int i = 0; i<data::N * 4; i++){
      mat_old[i] = mat[i];
    }
    //    cout<<"neigh_list_update_auto"<<endl;
    neigh_list_in(mat, neigh_list, n_nl);
    energy_in(mat, energy, neigh_list, energy_att2_atom, energy_rep_atom, energy_rep, energy_att2);
    data::count_up_nei++;
  }
}

double energy_check(double mat[]){
  double E_temp1, E_temp2, E_temp3, d;
  double E_tot = 0;
  E_temp3 = 0;
  for(int i = 0; i<data::N; i++){
    E_temp1 = 0;
    E_temp2 = 0;
    for(int j = 0; j<data::N; j++){
      if(i!=j){
	d = dist(i, j, mat);
	E_temp1 += E_rep(d, mat[4 * i + 0], mat[4 * j + 0]);
	E_temp2 += E_att2(d, mat[4 * i + 0], mat[4 * j + 0]);
      }
    }
    E_tot += E_temp1 - sqrt(E_temp2);
  }
  return E_tot;
}

void energy_in(double mat[], double energy[], int neigh_list[], double energy_att2_atom[], double energy_rep_atom[], double energy_rep[], double energy_att2[]){
  for(int i = 0; i<data::N; i++){
    for(int j = 0; j<data::length_nl; j++){
      energy_rep[data::length_nl * i + j] = 0;
      energy_att2[data::length_nl * i + j] = 0;
    }
    energy_rep_atom[i] = 0;
    energy_att2_atom[i] = 0;
  }
  
  int i = 0;
  int k = 0;
  //  const int len_nl = 150;
  double temp_rep, temp_att2, d;
  do{
    for(int j = 1; j<data::length_nl; j++){
      k = neigh_list[data::length_nl * i + j]-1;
      if(k+1 == 0){
	goto jump;
      }
      d = dist(i, k, mat);
      temp_rep = E_rep(d, mat[4 * i + 0], mat[4 * k + 0]);
      temp_att2 = E_att2(d, mat[4 * i + 0], mat[4 * k + 0]);
      
      energy_rep[data::length_nl * i + j] = temp_rep;
      energy_att2[data::length_nl * i + j] = temp_att2;
    
      energy_rep_atom[i] += temp_rep;
      energy_att2_atom[i] += temp_att2;
    }
  jump:
    energy[i] = energy_rep_atom[i] - sqrt(energy_att2_atom[i]);
    i++;
  }while(i<data::N);
}

double energy_tot(double energy[]){
  double E_tot = 0;
  
  for(int i = 0; i<data::N; i++){
    //    data::energy[i] = data::energy_rep_atom[i] - sqrt(data::energy_att2_atom[i]);
    E_tot += energy[i];
  }
  return E_tot;
}


void energy_update_sing_somm(double mat[], int atom, int neigh_list[], double energy[], double energy_att2_atom[], double energy_rep_atom[], double energy_rep[], double energy_att2[]){
  //  const int len_nl = 150;
  int k, k_inv;
  double E_temp_rep, E_temp_att2, d;
  energy_rep_atom[atom] = 0;
  energy_att2_atom[atom] = 0;
  
  for(int i = 1; i<data::length_nl; i++){
    k = neigh_list[data::length_nl * atom + i]-1;
    //   cout<<i<<endl;
    if(k+1 == 0){
      goto jump0;
      //      break;
    }
    k_inv = data::neigh_list_inv[k][atom];
    energy_rep_atom[k] -= energy_rep[data::length_nl * atom + i];
    energy_att2_atom[k] -= energy_att2[data::length_nl * atom + i];
      
    d = dist(atom, k, mat);
    E_temp_rep = E_rep(d, mat[4 * atom + 0], mat[4 * k + 0]);
    E_temp_att2 = E_att2(d, mat[4 * atom + 0], mat[4 * k + 0]);
    energy_rep_atom[k] += E_temp_rep;
    energy_att2_atom[k] += E_temp_att2;
    energy_rep_atom[atom] += E_temp_rep;
    energy_att2_atom[atom] += E_temp_att2;

    energy_rep[data::length_nl * atom + i] = E_temp_rep;
    energy_att2[data::length_nl * atom + i] = E_temp_att2;
    energy_rep[data::length_nl * k + k_inv] = E_temp_rep;
    energy_att2[data::length_nl * k + k_inv] = E_temp_att2;
    
    energy[k] = energy_rep_atom[k] - sqrt(energy_att2_atom[k]);
  }
 jump0:
  energy[atom] = energy_rep_atom[atom] - sqrt(energy_att2_atom[atom]);
}

void shake(double mat[], double &E_tot, int neigh_list[], ofstream &file_out, double energy[], double energy_att2_atom[], double energy_rep_atom[], double energy_rep[], double energy_att2[]){
  int count = 0;
  double r_shake, costheta_shake, theta_shake, phi_shake, r, E_old;
  int choice;
  do{
    choice = rand()%data::N;
    r_shake = ((double)rand()/RAND_MAX) * data::R_max_shake;
    theta_shake = ((double)rand()/RAND_MAX) * M_PI *2;
    phi_shake = ((double)rand()/RAND_MAX) * M_PI;
    
    mat[4 * choice + 1] += (r_shake * sin(phi_shake) * cos(theta_shake));
    mat[4 * choice + 2] += (r_shake * sin(phi_shake) * sin(theta_shake));
    mat[4 * choice + 3] += (r_shake * cos(phi_shake));
    energy_update_sing_somm(mat, choice, neigh_list, energy, energy_att2_atom, energy_rep_atom, energy_rep, energy_att2);
   
    E_old = E_tot;
    E_tot = energy_tot(energy);
    r = (double)rand()/RAND_MAX;
   
    if(r > exp(-(E_tot-E_old)/(data::k_B * data::T))){
      mat[4 * choice + 1] -= (r_shake * sin(phi_shake) * cos(theta_shake));
      mat[4 * choice + 2] -= (r_shake * sin(phi_shake) * sin(theta_shake));
      mat[4 * choice + 3] -= (r_shake * cos(phi_shake));
      energy_update_sing_somm(mat, choice, neigh_list, energy, energy_att2_atom, energy_rep_atom, energy_rep, energy_att2);
      E_tot = E_old;
      data::n_shake_rej++;
      data::n_shake_rej_tot++;
      //     cout<<"Move shake refused!"<<endl;
    }
    //    else{
    //      cout<<"Move shake accepted!"<<endl;
    //    }
    count++;
    data::n_shake++;
    data::n_shake_tot++;
  }while(count<data::N);
  
  //  if(data::n_shake%(data::N*50000)==0){
  if(data::n_shake_tot%(data::N*5000)==0){
    double acc_rate = 1.0 - ((double)data::n_shake_rej_tot)/((double)data::n_shake_tot);
    data::n_shake_rej = 0;
    data::n_shake = 0;
    if(acc_rate < 0.3){
      data::R_max_shake = data::R_max_shake * 0.9;
    }else if(acc_rate > 0.6){
      data::R_max_shake = data::R_max_shake * 1.1;
      if(data::R_max_shake>=(data::delta_NL*1/3)){
	data::R_max_shake /=1.1;
      }
    }
    cout<<"R_max: "<<data::R_max_shake<<endl;
    cout<<"acc rate: "<<acc_rate<<endl;
  }
}

void exchange_random(double mat[], double &E_tot, int neigh_list[], ofstream &file_out, double energy[], double energy_att2_atom[], double energy_rep_atom[], double energy_rep[], double energy_att2[]){
  //  cout<<"Exchange random"<<endl;
  int atom_0, atom_1, temp_1, temp_0;
  double r0[3];
  double r1[3];
  int count = 0;
  double r, E_old;
  do{
    temp_0 = rand()%data::n0;
    temp_1 = rand()%data::n1;

    atom_0 = data::atom_el0[temp_0];
    atom_1 = data::atom_el1[temp_1];
 
    mat[4 * atom_0 + 0] = 1;
    data::atom_el0[temp_0] = atom_1;
    
    mat[4 * atom_1 + 0] = 0;
    data::atom_el1[temp_1] = atom_0;
    
    energy_update_sing_somm(mat, atom_0, neigh_list, energy, energy_att2_atom, energy_rep_atom, energy_rep, energy_att2);
    
    //    mat[4 * atom_1 + 0] = 0;
    //    data::atom_el1[temp_1] = atom_0;

    energy_update_sing_somm(mat, atom_1, neigh_list, energy, energy_att2_atom, energy_rep_atom, energy_rep, energy_att2);
   
    E_old = E_tot;
    E_tot = energy_tot(energy);

    r = (double)rand()/RAND_MAX;
    
    if(r > exp(-(E_tot - E_old)/(data::k_B * data::T))){
      mat[4 * atom_0 + 0] = 0;
      data::atom_el0[temp_0] = atom_0;
      mat[4 * atom_1 + 0] = 1;
      data::atom_el1[temp_1] = atom_1;
      
      energy_update_sing_somm(mat, atom_0, neigh_list, energy, energy_att2_atom, energy_rep_atom, energy_rep, energy_att2);

      //mat[4 * atom_1 + 0] = 1;
      //      data::atom_el1[temp_1] = atom_1;
      energy_update_sing_somm(mat, atom_1, neigh_list, energy, energy_att2_atom, energy_rep_atom, energy_rep, energy_att2);

      E_tot = E_old;
      data::n_exch_rand_rej++;
      //      cout<<"Exchange refused!"<<endl;
    }
    //  else{
    //    cout<<"Exchange accepted!"<<endl;
    //  }
    data::n_exch_rand++;
    count++;
  }while(count<data::N);
}
 

void exchange_weighted(double mat[], double &E_tot, int neigh_list[], ofstream &file_out, double energy[], double energy_att2_atom[], double energy_rep_atom[], double energy_rep[], double energy_att2[]){
  //  cout<<"Exchange weighted "<<endl;
  int atom_0, atom_1, temp_0, temp_1;
  int l = 0;
  double W_old, coiche_j, E_old, W_new, r, temp;
  int nu = 2;
  do{
    temp_0 = rand()%data::n0;
    atom_0 = data::atom_el0[temp_0];
    
    vector<double> dist_exch_old;
    for(int j = 0; j<data::n1; j++){
      dist_exch_old.push_back(sqrt(dist(atom_0, data::atom_el1[j], mat)));
    }
    vector<double> weight_old;

    for(int j = 0; j<dist_exch_old.size(); j++){
      weight_old.push_back(1/pow(dist_exch_old[j], nu));
    }
    
    W_old = accumulate(weight_old.begin(), weight_old.end(),  0.0);
    vector<double> weight_rel_old;
    for(int j = 0; j<weight_old.size(); j++){
      weight_rel_old.push_back(weight_old[j]/W_old);
    }
    
    coiche_j = (double)rand()/RAND_MAX;
    temp_1 = 0;
    temp = 0;
    do{
      temp += weight_rel_old[temp_1];
      temp_1++;
    }while(temp<coiche_j);
    temp_1 = temp_1-1;
    atom_1 = data::atom_el1[temp_1];

    mat[4 * atom_0 + 0] = 1;
    data::atom_el0[temp_0] = atom_1;
    energy_update_sing_somm(mat, atom_0, neigh_list, energy, energy_att2_atom, energy_rep_atom, energy_rep, energy_att2);
     
    mat[4 * atom_1 + 0] = 0;
    data::atom_el1[temp_1] = atom_0;
    energy_update_sing_somm(mat, atom_1, neigh_list, energy, energy_att2_atom, energy_rep_atom, energy_rep, energy_att2);
    
    E_old = E_tot;
    E_tot = energy_tot(energy);

    vector<double> dist_exch_new;
    
    for(int j = 0; j<data::n1; j++){
      dist_exch_new.push_back(sqrt(dist(atom_1, data::atom_el1[j], mat)));
    }
    vector<double> weight_new;
    for(int j = 0; j<dist_exch_new.size(); j++){
      weight_new.push_back(1/pow(dist_exch_new[j], nu));
    }
    
    W_new = accumulate(weight_new.begin(), weight_new.end(),  0.0);
  
    r = (double)rand()/RAND_MAX;
    
    if(r > (W_old/W_new) * exp(-(E_tot - E_old)/(data::k_B * data::T))){
      mat[4 * atom_0 + 0] = 0;
      data::atom_el0[temp_0] = atom_0;
      energy_update_sing_somm(mat, atom_0, neigh_list, energy, energy_att2_atom, energy_rep_atom, energy_rep, energy_att2);
 
      mat[4 * atom_1 + 0] = 1;
      data::atom_el1[temp_1] = atom_1;
      energy_update_sing_somm(mat, atom_1, neigh_list, energy, energy_att2_atom, energy_rep_atom, energy_rep, energy_att2);

      E_tot = E_old;
      data::n_exch_weigh_rej++;
      //      cout<<"Exchange refused!"<<endl;
    }
    //  else{
    //    cout<<"Exchange accepted!"<<endl;
    //  }
    data::n_exch_weigh++;
    l++;
  }while(l<data::N);
}

void mc_volume(double mat[], double &E_tot, int neigh_list[], double max_delta_V, ofstream &file_out, double energy[], double energy_att2_atom[], double energy_rep_atom[], double energy_rep[], double energy_att2[]){
  //  cout<<"Change of size"<<endl;
  double old_length_x = data::L[0];
  double old_length_y = data::L[1];
  double old_length_z = data::L[2];
  double V_old = data::L[0] * data::L[1] * data::L[2];
  double delta_V = ((double)rand()/RAND_MAX)*2*max_delta_V - max_delta_V;
  double V_new = V_old + delta_V;
  double shift = pow(V_new, 1./3);
  double shift_ratio = pow(V_new/V_old, 1./3);
  for(int k = 0; k<3; k++){
    data::L[k] *= shift_ratio;
  }

  for(int j = 0; j<data::N; j++){
    for(int k = 1; k<4; k++){
      mat[4 * j + k] = mat[4 * j + k] * shift_ratio;
    }
  }

  energy_in(mat, energy, neigh_list, energy_att2_atom, energy_rep_atom, energy_rep, energy_att2);
  double E_old = E_tot;
  E_tot = energy_tot(energy);
   
  double r = (double)rand()/RAND_MAX;
  
  if(r > pow(V_new/V_old, data::N) * exp(-(E_tot-E_old)/(data::k_B * data::T))){
    for(int j = 0; j<data::N; j++){
      for(int k = 1; k<4; k++){
	mat[4 * j + k] = mat[4 * j + k] / shift_ratio;
      }
    }
    
    data::L[0] = old_length_x;
    data::L[1] = old_length_y;
    data::L[2] = old_length_z;
    
    energy_in(mat, energy, neigh_list, energy_att2_atom, energy_rep_atom, energy_rep, energy_att2);
    E_tot = E_old;
    data::n_vol_rej++;
    //    cout<<"Move volume refused!"<<endl;
  }
  //   else{
  //     cout<<"Move volume accepted!"<<endl;
  // }
  data::n_vol++;
}

void MC(double mat[], ofstream &file_images_out, ofstream &file_data_out, double energy[], double energy_att2_atom[], double energy_rep_atom[], int n_nl[], vector<double> &data, vector<double> &volume, int neigh_list[], double energy_rep[], double energy_att2[], int &MC_step_run){
  energy_in(mat, energy, neigh_list, energy_att2_atom, energy_rep_atom, energy_rep, energy_att2);
  double E_tot_in = energy_tot(energy);
  //  cout<<"Initial total energy: "<<setprecision(12)<<E_tot_in<<" eV"<<endl<<endl;
  double E_tot = E_tot_in;
  vector<double> modes(4);
  data.push_back(E_tot);
  int tot_move = 0;
  for(int i = 0; i<data::MC_moves.size(); i++){
    modes[i] = data::MC_moves[i];
  }
	
  if(data::n0 == 0 || data::n1 == 0){
    modes[1] = 0;
    modes[2] = 0;
  }else{
    modes[1] = data::MC_moves[1];
    modes[2] = data::MC_moves[2];
  }
  if(modes[0]==0){
    int prpopor = 100-modes[1] - modes[2];
    modes[3] = prpopor;
  }

  double choice_move;
  double E_check;
  double eps = 1e-5;
  double box_old[data::N * 4];
	
  for(int i = 0; i<data::N; i++){
    for(int j = 0; j<4; j++){
      box_old[4 * i + j] = mat[4 * i + j];
    }
  }

  for(int i = MC_step_run; i<data::MC_steps; i++){
    if(i>0){
      neigh_list_update_auto(box_old, mat, neigh_list, energy, energy_att2_atom, energy_rep_atom, n_nl, energy_rep, energy_att2);
    }
		
    if(i%1000 == 0){
      E_check = energy_check(mat);
      cout<<"MC step: "<<i<<endl;
      cout<<"Energy check!"<<endl;
      cout<<"E_check: "<<E_check<<endl;
      cout<<"E_tot: "<<E_tot<<endl<<endl;
      if(E_check != E_tot){
	if((E_check+eps) < E_tot || (E_check-eps) > E_tot){
	  cout<<"Error in energy calculation!!"<<endl;
	  exit(0);
	}else if(isnan(E_tot)){
	  cout<<"Total energy is not valid!\n";
	  exit(0);
	}
      }
    }
    choice_move = (double)rand()/RAND_MAX;
    if(choice_move < (double)modes[1]/100){
      //Exchange random move
      exchange_random(mat, E_tot, neigh_list, file_images_out, energy, energy_att2_atom, energy_rep_atom, energy_rep, energy_att2);
    }else if(choice_move >=  (double)modes[1]/100 && choice_move < ((double)modes[1]/100 + (double)modes[2]/100)){
      //Exchange weighted move
      exchange_weighted(mat, E_tot, neigh_list, file_images_out, energy, energy_att2_atom, energy_rep_atom, energy_rep, energy_att2);
    }else if(choice_move >=  ((double)modes[1]/100 + (double)modes[2]/100) && choice_move < ((double)modes[1]/100 + (double)modes[2]/100 + (double)modes[3]/100)){
      // Volume move
      double max_V_change = 0.1 * data::N;
      mc_volume(mat, E_tot, neigh_list, max_V_change, file_images_out, energy, energy_att2_atom, energy_rep_atom, energy_rep, energy_att2);
    }else{
      // Shake move
      shake(mat, E_tot, neigh_list, file_images_out, energy, energy_att2_atom, energy_rep_atom, energy_rep, energy_att2);
    }
    if(data::file_xyz_yn == 1){
      if(i%data::MC_step_imm == 0 || i == data::MC_steps-1){
	file_xyz(data::elements[0], data::elements[1], mat, file_images_out, i, E_tot);
      }
    }
    
    if(i%(data::MC_step_imm/2) == 0 || i == data::MC_steps-1){
      write_file_out((i), E_tot, file_data_out);
    }
		
    data.push_back(E_tot);
    volume.push_back(data::L[0] * data::L[1] * data::L[2]);
  }
  MC_step_run = data::MC_steps;
}
  
void write_file_out(int step, double &E, ofstream &file_data_out){
  file_data_out<<step<<"\t"<<setprecision(12)<<E/data::N<<"\t"<<data::L[0] *data::L[1] * data::L[2] <<endl;
}

double average(vector<double> &data){
  int el =  data.size();
  double av = 0;
  for(double i : data)
    av += i;
  return av / el;
}

double VAR(vector<double> &data){
  double av = average(data);

  // Calculate sum of squared differences
  double sumSquaredDifferences = 0.0;
  for (const double& value : data) {
    double difference = value - av;
    sumSquaredDifferences += difference * difference;
  }

  double variance = sumSquaredDifferences / (data.size());
  return variance;
}

void convergence(vector<double> &data, vector<double> &volume, vector<vector<double>> &prop, double comp_el1, string name_file_conv, ofstream &file_prop, int &check_conv, double eps_min, double eps_max, int block_size){
  //Set parameter
  int tot_val = data::MC_steps + 1;
  int l = tot_val/block_size;
  vector<double> data_conv;
 
  //Block
  vector<double> av_block;
  int cont = 0;
  int o = 0;
  double av_block_temp;
  do{
    av_block_temp = 0;
    for(int i = o; i<(o+block_size); i++){
      if((o+block_size)>tot_val) break;
      av_block_temp += data[i]/data::N;
    }
    av_block.push_back(av_block_temp/block_size);
    o = o+block_size;
    cont++;
  }while(cont<l);
  
  //STD for the blocks
  vector<double> STD_block, x;
  cont = 0;
  o = 0;
  do{
    double std_block_temp = 0;
    for(int i = o; i<(o+block_size); i++){
      if((o+block_size)>tot_val) break;
      std_block_temp += sq((data[i]/data::N-av_block[cont]));
    }
    STD_block.push_back(sqrt(std_block_temp/(block_size-1)));
    x.push_back(o+block_size/2);
    o = o+block_size;
    cont++;
  }while(cont<l);
  
  //Fitting
  double inter;
  int point_conv = 0;
  double ang1 = 1;
  int initial_step = 50;
		
  int n;
  for(int r = initial_step; r<(l-60); r = r+10){
    n = l - r;
    gsl_matrix *X = gsl_matrix_alloc(n, 2);
    gsl_vector *Y = gsl_vector_alloc(n);
    gsl_vector *W = gsl_vector_alloc(n);
    gsl_vector *c = gsl_vector_alloc(2);
    gsl_matrix *cov = gsl_matrix_alloc(2, 2);

    double chisq;

    for (size_t i = 0; i < n; ++i) {
      gsl_matrix_set(X, i, 0, 1.0);
      gsl_matrix_set(X, i, 1, x[i+r]);
      gsl_vector_set(Y, i, av_block[i+r]);
      gsl_vector_set(W, i, 1.0 / (STD_block[i+r] * STD_block[i+r]));  // Inversesquare of the error as weight
    }
		
    // Perform the fit
    gsl_multifit_linear_workspace *work = gsl_multifit_linear_alloc(n, 2);
    gsl_multifit_wlinear(X, W, Y, c, cov, &chisq, work);
    gsl_multifit_linear_free(work);
			
    inter = gsl_vector_get(c, 1);
    if(abs(inter)<eps_min){
      if(ang1 == 1){
	ang1 = inter;
      }
      if(abs(inter)<=abs(ang1)){
	//Update the start convergence point
	ang1 = inter;
	point_conv = r*block_size;
      }
      if(chisq/n<1 && abs(ang1)<eps_max)
	break;
    }
  }
  
  //Write the data in a file, starting from the convergence point
  vector<double> vol_conv;
  vector<double> vol_conv_2;
  if(point_conv !=0){
    ofstream file_conv;
    file_conv.open(name_file_conv);
    file_conv<<"MC_step\tEn/N\tV(A^3)"<<endl;
    for(int i = (point_conv + 10*block_size); i<tot_val-1; i++){
      //  for(int i = 0; i<tot_val; i++){
      file_conv<<i<<"\t"<<setprecision(12)<<data[i]/data::N<<"\t"<<volume[i]<<endl;
      data_conv.push_back(data[i]);
      vol_conv.push_back(volume[i]);
      vol_conv_2.push_back(sq(volume[i]));
    }
    check_conv = 1;
    file_conv.close();
  }else{
    cout<<"Doesn't converge!!"<<endl;
    check_conv = 0;
    return;
  }
  
  //Correlation time
  tot_val = data_conv.size();
  int tot_val_break = 10000;
	
  //Total average
  double mean_tot = 0;
  double mean_tot_2 = 0;
  for(int i = 0; i<tot_val; i++){
    mean_tot += data_conv[i];
    mean_tot_2 += sq(data_conv[i]);
  }
  mean_tot /= tot_val;
  mean_tot_2 /= tot_val;
  
  double rho_0 = 0;
  int tau = 0;
  double E_0 = 0;
  double E_t = 0;
  double rho = 0;
  do{
    for(int t = 0; t<tot_val_break; t++){
      E_0 = 0;
      E_t = 0;
      rho = 0;
      for(int tp = 0; tp<tot_val-t; tp++){
	E_0 += data_conv[tp];
	E_t += data_conv[tp + t];
      }
      E_0 = E_0 / (tot_val - t);
      E_t = E_t / (tot_val - t);
			
      if(t == 0){
	for(int tp = 0; tp<tot_val-t; tp++){
	  rho_0 += (data_conv[tp] - E_0) * (data_conv[tp+t] - E_t); 
	}
	rho_0 = rho_0 / (tot_val - t);
      }else{
	for(int tp = 0; tp<tot_val-t; tp++){
	  rho += (data_conv[tp] - E_0) * (data_conv[tp+t] - E_t); 
	}
	rho = rho /(rho_0 * (tot_val - t));
	if(rho<2*exp(-1.0) && tau == 0){
	  tau = (int)t;
	  //	cout<<"tau: "<<tau<<endl;
	}
      }
    }
    //  }while(tau>exp(-1.0));
  }while(tau == 0);
  tau = 3*tau;
  if(tot_val/tau<20){
    check_conv = 0;
    return;
  }
  /*To check*/
  const double atomic_mass_comp = comp_el1 * data::mass[0] + (1-comp_el1)*data::mass[1];
  //const double conv_fact = data::N_A * data::e /atomic_mass_comp;    //eV/atom -> J/g
  const double conv_fact = data::N_A * data::e; //J/mol
  //
  
  //Blocking method
  int M_b = tot_val / tau;
  vector<double> E_av_b(M_b);
  vector<double> E_av2_b(M_b);
  cont = 0;
  int tempo = 0;
  do{
    for(int i = cont; i<cont+tau; i++){
      E_av_b[tempo] += data_conv[i];
      E_av2_b[tempo] += data_conv[i] * data_conv[i];
    }
    cont += tau;
    tempo++;
  }while(cont<(M_b)*tau);
	
  for(int i = 0; i<M_b; i++){
    E_av_b[i] /= (double)tau; 
    E_av2_b[i] /= tau;
  }
  
  double E_av = 0;
  double E_av_2 = 0;
  for(int i = 0; i<M_b; i++){
    E_av += E_av_b[i];
    E_av_2 += E_av2_b[i];
  }
  E_av /= M_b;
  E_av_2 /= M_b;
  double V_av = average(vol_conv);
  double V_av_2 = average(vol_conv_2);
  vector<double> latt_par_conv;
  for(int i = 0; i<vol_conv.size(); i++){
    latt_par_conv.push_back(pow(vol_conv[i], 1./3)/data::n);
  }
	
  double sig_E = sqrt((1/((double)M_b-1)) * VAR(E_av_b));
  vector<double> data_conv_2;
  for(int i = 0; i<data_conv.size(); i++){
    data_conv_2.push_back(sq(data_conv[i]));
  }
  double E_mean = average(data_conv);
  double E_2_mean = average(data_conv_2);
  double latt_par_mean = average(latt_par_conv);
  double compr = (V_av_2 - sq(V_av))/(V_av*data::T*data::k_B);
  double delta_E_2 = E_2_mean - sq(E_mean);
  double h_c = (delta_E_2/(sq(data::T)*data::k_B) + 1.5*data::N*data::k_B)/data::N;
  h_c *= conv_fact;
  file_prop<<data::T<<"\t"<<comp_el1*100<<"\t"<<E_av<<"\t"<<sig_E<<"\t"<<h_c<<"\t"<<h_c/atomic_mass_comp<< "\t"<<V_av<<"\t"<<latt_par_mean<<"\t"<<compr<<endl;
  prop.push_back({comp_el1*100, E_av, sig_E, h_c/atomic_mass_comp, h_c});
}

void change_composition(double mat[], double perc_f_el){
  int n0 = data::N*perc_f_el;
  int n0_old = data::n0;
  if(n0 != data::n0){
    int atom, atom_v;
    //for(int i = 0; i<(n0-data::n0); i++){
    int cont_new0 = 0;
    do{
      atom_v = rand()%data::n1;
      atom = data::atom_el1[atom_v];
      if(mat[4 * atom + 0] == 1){
	cont_new0++;
      }
      mat[4 * atom + 0] = 0;
    }while(cont_new0 < (n0-n0_old));
    data::n0 = n0;
    data::n1 = data::N - n0;
    delete[] data::atom_el0;
    delete[] data::atom_el1;
    data::atom_el0 = new int[data::n0];
    data::atom_el1 = new int[data::n1];
    
    int cont0 = 0;
    int cont1 = 0;
    for(int i = 0; i<data::N; i++){
      if(mat[4 * i + 0] == 0){
	data::atom_el0[cont0] = i;
	cont0++;
      }else if(mat[4 * i + 0] == 1){
	data::atom_el1[cont1] = i;
	cont1++;
      }
    }
  }
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