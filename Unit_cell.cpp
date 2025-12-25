#include<iomanip>

#include"functions.h"
#include"data.h"

using namespace std;

int main(int argc, char* argv[]){
  if(argc != 2){
    cerr<<"Incorrect number of arguments. Expected: number of the simulation\n";
    return 1;
  }
  cerr<<"MC started run_"<<argv[1]<<endl;
  ifstream file_check;
  
  // Check that the initial lattices are ready
  try{
    file_check.open("qq");
    if(!file_check.is_open()){
      throw runtime_error("Unable to open 'qq'");
    }
    string last_line;
    last_line = getLastLineFromFile(file_check);
    file_check.close();
    if(last_line != "initial lattices created"){
      cout<<"initial lattices was not created!\n";
      exit(0);
    }
  }catch(const exception &e){
    cerr<<"Error: "<<e.what()<<endl;
    return 1;
  }
  
  //Read input file and initialization of quantities
  cout<<setprecision(12);

  // Read 'leggi.in'
  read_input();
  srand(data::seed);

  vector<vector<double>> prop;
  vector<double> data;//[data::MC_steps*10+1];
  vector<double> volume;
  ofstream cell;
  ofstream file_out_data;
  ofstream file_out_data_conv;
  string name_file_out_data_conv;
  string file_in_latt;
  int check_conv;
  int MC_steps = data::MC_steps;
	
  // Reading the "potential.in" file
  read_potential();
  
  // Calculate the parameters for the cut_off function
  calc_par_en();
  
  //Loop over the compositions
  ofstream file_prop("properties.out");
  if(!file_prop.is_open()){
    cerr << "Error opening 'properties.out'" <<endl;
    exit(0);
  }else{
    cout<<"'properties.out' is opened successfully"<<endl<<endl;
  }
  file_prop<<"T\tcompo El1\tE mean\tsig E\t h_c (J/(mol K))\th_c (J/(g K))\tV_mean(A^3)\tlatt_par\tisotermal_compr\n";
  
  for(double u:data::perc_first_el){
    cout<<"%"<<data::elements[0]<<": "<<u*100<<endl;
    data.clear();
    volume.clear();
    
    //Definition the name of the files
    if(u*100<1 && u!=0){
      if(data::file_xyz_yn == 1){
	      cell.open(data::file_out_images+"_0"+to_string((int)(u*1000))+".xyz");
      }
      file_out_data.open(data::file_out_data+"_0"+to_string((int)(u*1000))+".out");
      name_file_out_data_conv = data::file_out_data+"_conv_0"+to_string((int)(u*1000))+".out";
    }else if (u*100>99 && u!=1){
      if(data::file_xyz_yn == 1){
	      cell.open(data::file_out_images+"_"+to_string((int)(u*1000))+".xyz");
      }
      file_out_data.open(data::file_out_data+"_"+to_string((int)(u*1000))+".out");
      name_file_out_data_conv =data::file_out_data+"_conv_"+to_string((int)(u*1000))+".out";
    }else{
      if(data::file_xyz_yn == 1){
	      cell.open(data::file_out_images+"_"+to_string((int)(u*100))+".xyz");
      }
      int u_t = (double)u*100;
      file_in_latt = "cell_in_"+to_string(u_t)+".xyz";
      ifstream file_temp;
      file_temp.open(file_in_latt);
      file_temp>>data::N;
      file_temp.close();
      cout<<data::N<<endl;
      data::n0 = data::N * u;
      data::n1 = data::N - data::n0;
      cout<<data::n0<<"\t"<<data::n1<<endl;
      file_out_data.open(data::file_out_data+"_"+to_string(u_t)+".out");
      name_file_out_data_conv = data::file_out_data+"_conv_"+to_string(u_t)+".out";
    }

    data::length_nl = data::N + data::N/10;
    int *neigh_list = new int [data::N * data::length_nl];
    double *energy_rep = new double [data::N * data::length_nl];
    double *energy_att2 = new double[data::N * data::length_nl];
    data::neigh_list_inv = new int* [data::N];
    for(int i = 0; i<data::N; i++){
      data::neigh_list_inv[i] = new int[data::N];
    }
    
    double box[data::N * 4];
    double energy[data::N];
    double energy_att2_atom[data::N];
    double energy_rep_atom[data::N];
    int n_nl[data::N];
    
    // Reads the file with the crosshairs
    read_lattice(box, file_in_latt);
    file_out_data<<"MC_step\tEn/N\tVol (Å^3)"<<endl;
    int step_in = 0;
    cout<<(int)(data::N*u)<<endl;
    // Check that the system is correct
    if(data::n0 != (int)(data::N*u)){
      change_composition(box, u);
    }
    // Creation of the neighbour list
    neigh_list_in(box, neigh_list, n_nl);

    // Initialization of energy matrices
    energy_in(box, energy, neigh_list, energy_att2_atom, energy_rep_atom, energy_rep, energy_att2);
     
    double E_tot = energy_tot(energy);
     
    data::n_shake_rej_tot = 0;
    data::n_exch_rand_rej = 0;
    data::n_exch_weigh_rej = 0;
    data::n_vol_rej = 0;
    data::n_shake_tot = 0;
    data::n_exch_rand = 0;
    data::n_exch_weigh = 0;
    data::n_vol = 0;
    data::count_up_nei = 0;
     
    cout<<setprecision(12);
    if(u == 0){
      data::MC_steps = MC_steps/2;
    }else if(u == 1){
      data::MC_steps = MC_steps/2;
    }else{
      data::MC_steps = MC_steps;
    }
    cout<<data::MC_steps<<endl;
    check_conv = 0;
    int MC_step_run = 0;

    // Performs Monte Carlo for a maximum of 1000000 MCsteps
    do{
      MC(box, cell, file_out_data, energy, energy_att2_atom, energy_rep_atom, n_nl, data, volume, neigh_list, energy_rep, energy_att2, MC_step_run);

      // Check the convergence
      convergence(data, volume, prop, u, name_file_out_data_conv, file_prop, check_conv);
      data::MC_steps +=1e5;
      if(data::MC_steps > 1e6){
	cerr<<"Data does not converge!"<<endl;
	exit(0);
      }
    }while(check_conv != 1);
     
    cout<<"count_update_nei: "<<data::count_up_nei<<endl;
    double perc_shake_acc = (1 - (double)data::n_shake_rej_tot/(double)data::n_shake_tot) * 100;
    double perc_exch_rand_acc = (1 - (double)data::n_exch_rand_rej/(double)data::n_exch_rand) * 100;
    double perc_exch_weigh_acc = (1 - (double)data::n_exch_weigh_rej/(double)data::n_exch_weigh) * 100;
    double perc_vol_acc = (1 - (double)data::n_vol_rej/(double)data::n_vol) * 100;
    cout<<endl;
    cout<<setprecision(4);
    printf("Shake moves called: %d\n", data::n_shake_tot/data::N);
    cout<<"Rate of shake moves accepted: "<<perc_shake_acc<<"%"<<endl;
    printf("\nExchange random moves called: %d\n", data::n_exch_rand/data::N);
    cout<<"Rate of exchange random moves accepted: "<<perc_exch_rand_acc<<"%"<<endl;
    printf("\nExchange weighted moves called: %d\n", data::n_exch_weigh/data::N);
    cout<<"Rate of exchange weighted moves accepted: "<<perc_exch_weigh_acc<<"%"<<endl;
    printf("\nVolume moves called: %d\n", data::n_vol);
    cout<<"Rate of volume resize moves accepted: "<<perc_vol_acc<<"%"<<endl<<endl;
    cell.close();
    file_out_data.close();
    file_out_data_conv.close();
  }
  
  cout<< "MC finished!"<<endl;
  cerr<<"MC finished run_"<<argv[1]<<endl;
  return 0;
}
