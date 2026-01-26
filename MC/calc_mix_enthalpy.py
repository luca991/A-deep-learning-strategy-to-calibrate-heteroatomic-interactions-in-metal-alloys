#This code calculates the mixing enthalpy using the pure energies calculated with my MC code in tha same session
import sys
import subprocess
import numpy as np
import math
from scipy.interpolate import interp1d
import scipy.optimize as opt
from scipy.constants import Boltzmann, eV

kB = Boltzmann/eV

def fit_func(x, a, b):
	return a + b*x


def delete_last_lines(file_path, num_lines):
	with open(file_path, 'r') as file:
		lines = file.readlines()[:-num_lines]

	with open(file_path, 'w') as file:
		file.writelines(lines)

def copy_file(original_file, new_file):
	with open(original_file, 'r') as original:
		with open(new_file, 'w') as new:
			new.write(original.read())

#Main
if len(sys.argv)>1:
	tempNameRun = sys.argv[1]
else:
	print("The simulation name was not entered during launch\n")
	sys.stderr.write("The simulation name was not entered during launch\n")
	exit()


file_leggi = open("leggi.in", "r")
cont = 0
T = 0
n = 0
n_hcp = []
name_file_pot = 0
compo = []
for x in file_leggi:
	cont += 1
	if(cont == 2):
		x = x.split()
		name_file_pot = x[0]
	if(cont == 4):
		x = x.split()
		n = int(x[0])
	if(cont == 5):
		x = x.split()
		n_hcp.append(int(x[0]))
		n_hcp.append(int(x[1]))
		n_hcp.append(int(x[2]))
	if(cont == 6):
		x = x.split()
		for j in x:
			try:
				compo.append(float(j))
			except ValueError:
				pass
	elif(cont == 11):
		x = x.split()
		T = float(x[0])
file_leggi.close()
print(T)
print(compo)


name_file_El = []
latt_type= []

try:
	with open(name_file_pot, 'r') as file_pot:
		print(name_file_pot, " opened correctly")
		cont = 0
		for x in file_pot:
			cont += 1
			if(cont == 3):
				x = x.split()
				latt_type.append(x[2])
				latt_type.append(x[3])
except FileNotFoundError:
	print("File ", name_file_pot, " not found")
	sys.stderr.write(tempNameRun + ": File " + name_file_pot + " not found\n")
	exit()
except IOError:
	print("Error in opening ", name_file_pot)
	sys.stderr.write(tempNameRun + ": Error in opening " + name_file_pot + "\n")
	exit()


H_sim_0 = 0
H_sim_100=0
H_sim, delta_H_sim, hc, conc, hc_g, V, latt_par_0, latt_par, compre = [], [], [], [], [], [], [], [], []
try:
	with open('properties.out', 'r') as file_prop:
		print("File 'properties.out' opened correctly")
		cont = 0
		for x in file_prop:
			cont +=1
			if(cont >= 2):
				x = x.split()
				if(float(x[1])==0):
					H_sim_0 = float(x[2])
				elif(float(x[1])==100):
					H_sim_100=float(x[2])
				H_sim.append(float(x[2]))
				delta_H_sim.append(float(x[3]))
				hc.append(float(x[4]))
				hc_g.append(float(x[5]))
				V.append(float(x[6]))
				latt_par.append(float(x[7]))
				compre.append(float(x[8]))
except FileNotFoundError:
	print("File 'properties.dat' not found")
	sys.stderr.write(tempNameRun + ": File 'properties.dat' not found\n")
	exit()
except IOError:
	print("Error in opening 'properties.dat'")
	sys.stderr.write(tempNameRun + ": Error in opening 'properties.dat'\n")
	exit()

H_mix, delta_H_mix = [], []
N = 0
cont = 0
print(compo)

print(compo)
for x in compo:
	compos = float(x)
	if(float(x) < 50):
		if(latt_type[1] == 'fcc'):
			N = 4*n**3
			H_mix.append(H_sim[cont]/N - ((compos/100) * (H_sim_100/N) + ((100-compos)/100)*(H_sim_0/N)))
		elif(latt_type[1] == 'hcp'):
			N = 4*n_hcp[0]*n_hcp[1]*n_hcp[2]
			H_mix.append(H_sim[cont]/N - ((compos/100) * (H_sim_100/N) + ((100-compos)/100)*(H_sim_0/N)))
		elif(latt_type[1] == 'bcc'):
			N = 2*n**3
			H_mix.append(H_sim[cont]/N - ((compos/100) * (H_sim_100/N) + ((100-compos)/100)*(H_sim_0/N)))
	else:
		print(x)
		if(latt_type[0] == 'fcc'):
			N = 4*n**3
			H_mix.append(H_sim[cont]/N - ((compos/100) * (H_sim_100/N) + ((100-compos)/100)*(H_sim_0/N)))
		elif(latt_type[0] == 'hcp'):
			N = 4*n_hcp[0]*n_hcp[1]*n_hcp[2]
			H_mix.append(H_sim[cont]/N - ((compos/100) * (H_sim_100/N) + ((100-compos)/100)*(H_sim_0/N)))
		elif(latt_type[0] == 'bcc'):
			N = 2*n**3
			H_mix.append(H_sim[cont]/N - ((compos/100) * (H_sim_100/N) + ((100-compos)/100)*(H_sim_0/N)))
	cont = cont+1


try:
	with open('H_mix.dat', 'a') as file_out:
		print("File 'properties_fin.dat' opened correctly")
		for i in range(0, len(compo)):
			file_out.write(str(T) + "\t" + str(compo[i]) + "\t" + str(H_mix[i]) +"\n")
		file_out.close()
		sys.stderr.write(tempNameRun + ": Finished\n")
except FileNotFoundError:
	print("File 'properties_fin.dat' not found")
	sys.stderr.write("File 'properties_fin.dat' not found\n")
	exit()
except IOError:
	print("Error in opening'properties_fin.dat'")
	sys.stderr.write("Error in opening'properties_fin.dat'\n")
	exit()
