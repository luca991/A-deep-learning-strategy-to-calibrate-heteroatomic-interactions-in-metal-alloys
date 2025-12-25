#include "data.h"
#include "functions.h"
#include <iomanip>
#include <cmath>
using namespace std;

inline string getLastLineFromFile(ifstream &filename) {
    std::string lastLine, line;
    while (std::getline(filename, line)) {
        if (!line.empty()) {
            lastLine = line;
        }
    }
    
    return lastLine;
};

int main(int argc, char* argv[]){
  ifstream file_check;
  
  srand(time(0));
  
  cerr<<"impure started run_"<<argv[1]<<endl;
  
  cout<<setprecision(14);

  //Initialization
  read_input();
  read_potential();
  calc_par_en();
  srand(data::seed);

  //Open data output file
  ofstream file_out_en;
  file_out_en.open("impurity_energies_"+(data::elements[0])+(data::elements[1])+".out");
  file_out_en<<"Energy of an impurity in a pure bulk, expressed in eV\n1"+data::elements[1]+"in"+data::elements[0]+"\t1"+data::elements[0]+"in"+data::elements[1]+"\n";
  ofstream file_out_fig;
  char latt_t;
  double energy, latt_par, e_pure_0, e_pure_1, e_imp_0, e_imp_1, e_impurity_1, e_impurity_0;
  int imp;
  
  //First system - Pure El0
  cout<<"Pure "<<data::elements[0]<<": \n";
  if(data::lattice_type[0] == "fcc"){
    data::N = pow(data::n, 3)*4;
    latt_t = 'f';
    latt_par = data::r0[0][0] * sqrt(2);
  }else if(data::lattice_type[0] == "hcp"){
    data::N = 4*data::n_cells_hcp[0]*data::n_cells_hcp[1]*data::n_cells_hcp[2];
    latt_t = 'h'; 
    latt_par = data::r0[0][0];
  }else if(data::lattice_type[0] == "bcc"){
    data::N = 2*pow(data::n, 3);
    latt_t = 'b'; 
    latt_par = data::r0[0][0]*(2./sqrt(3));
  }

  data::length_nl = data::N + data::N/10;
  double mat[data::N][4];
  data::neigh_list.resize(data::N, vector<int>(500));
  int n_nl[data::N];

  lattice(mat, data::n, latt_par, 0, latt_t);
  energy = energy_check(mat);
  cout<<"Initial energy: "<<energy<<endl;
  
  data::mat_fix.resize(data::N, vector<double>(4));
  for(int i = 0; i<data::N; i++){
    for(int j = 0; j<4; j++){
      data::mat_fix[i][j] = mat[i][j];
    }
  }

  if(latt_t == 'h'){
    relax_ca(mat);
  }else{
    relax(mat, 100, 1, 1e-10, 2, 1);
  }
  e_pure_0 = energy_check(mat);
  cout<<"E pure sys 1 min: "<<e_pure_0<<endl;
  cout<<"E atomic pure sys 1 min: "<<e_pure_0/data::N<<endl;
  neigh_list_in(mat, n_nl);

  //Impure El0
  cout<<endl<<"Pure "<<data::elements[0]<<" with one "<<data::elements[1]<<": \n";
  data::neigh_list.resize(data::N, vector<int>(500));
  imp = rand_unif_int(0, data::N);
  set_impurity(mat, imp);
  neigh_list_in(mat, n_nl);
  relax_all(mat, 5e-3, 0.5, 2e5);
  e_imp_0 = energy_check(mat);
  cout<<"E min imp sys 1: "<<e_imp_0<<endl;

  cout<<"First system completed\n";
  
  //Second system - Pure El1
  cout<<endl<<"Pure "<<data::elements[1]<<": \n";
  if(data::lattice_type[1] == "fcc"){
    data::N = pow(data::n, 3)*4;
    latt_t = 'f';
    latt_par = data::r0[1][1] * sqrt(2);
  }else if(data::lattice_type[1] == "hcp"){
    data::N = 4*data::n_cells_hcp[0]*data::n_cells_hcp[1]*data::n_cells_hcp[2];
    latt_t = 'h'; 
    latt_par = data::r0[1][1];
  }else if(data::lattice_type[1] == "bcc"){
    data::N = 2*pow(data::n, 3);
    latt_t = 'b'; 
    latt_par = data::r0[1][1]*(2./sqrt(3));
  }

  data::length_nl = data::N + data::N/10;
  double mat1[data::N][4];
  lattice(mat1, data::n, latt_par, 1, latt_t);

  energy = energy_check(mat1);
  cout<<"Initial energy: "<<energy<<endl;
	
  data::mat_fix.resize(data::N, vector<double>(4));
  for(int i = 0; i<data::N; i++){
    for(int j = 0; j<4; j++){
      data::mat_fix[i][j] = mat1[i][j];
    }
  }

  if(latt_t == 'h'){
    relax_ca(mat1);
  }else{
    relax(mat1, 100, 1, 1e-10, 2, 1);
  }
  e_pure_1 = energy_check(mat1);
  cout<<"E pure sys 1 min: "<<e_pure_1<<endl;
  cout<<"E atomic pure sys 1 min: "<<e_pure_1/data::N<<endl;

  //Impure El1
  imp = rand_unif_int(0, data::N);
  set_impurity(mat1, imp);
  data::neigh_list.resize(data::N, vector<int>(500));
  neigh_list_in(mat1, n_nl);
  relax_all(mat1, 5e-3, 0.5, 2e5);
  e_imp_1 = energy_check(mat1);

  cout<<"E imp sys min: "<<e_imp_1<<endl;
	
  e_impurity_1 = e_imp_0 - ((double)(data::N-1)/(double)data::N)*e_pure_0 - (e_pure_1/(double)data::N);
  e_impurity_0 = e_imp_1 - ((double)(data::N-1)/(double)data::N)*e_pure_1 - (e_pure_0/(double)data::N);
  file_out_en<<e_impurity_1<<"\t"<<e_impurity_0<<endl;

  cout<<"Second system completed\n";
  file_out_en.close();
  cerr<<"impure finished run_"<<argv[1]<<endl;
  return 0;
}
