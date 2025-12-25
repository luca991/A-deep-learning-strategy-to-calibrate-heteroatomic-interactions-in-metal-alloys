import sys
import subprocess
import numpy as np
import math
from scipy.interpolate import interp1d
import scipy.optimize as opt

def fit_func(x, a, b):
	return a + b*x

#Main
if len(sys.argv)>1:
	tempNameRun = sys.argv[1]
else:
	print("The simulation name was not entered during launch\n")
	sys.stderr.write("The simulation name was not entered during launch\n")
	exit()

try:
	with open('qq', 'r') as tmp:
		print('qq', "opened correctly")
		lastLine = None
		for line in tmp:
			lastLine = line
except FileNotFoundError:
	print("File 'qq' not found")
	sys.stderr.write(tempNameRun + ": File 'qq' not found\n")
	exit()
except IOError:
	print("Error in opening 'qq'")
	sys.stderr.write(tempNameRun + ": Error in opening 'qq'\n")
	exit()
tmp.close()
lastLine = lastLine.strip()

if(lastLine == 'Error in energy calculation!!'):
	print("The H_mix cannot be calculated because an error occurred in the energy calculation!\n")
	sys.stderr.write(tempNameRun + ": The H_mix cannot be calculated because an error occurred in the energy calculation!\n")
	exit()
elif(lastLine == "Doesn't converge!!"):
	print("The H_mix cannot be calculated because the system has not reached the convergence!\n")
	sys.stderr.write(tempNameRun + ": The H_mix cannot be calculated because the system has not reached the convergence!\n")
	exit()
elif(lastLine == "Total energy is not valid!"):
	print("The H_mix cannot be calculated because the the total energy is not valid!\n")
	sys.stderr.write(tempNameRun + ": The H_mix cannot be calculated because the the total energy is not valid!\n")
	exit()
elif(lastLine == "Error: potential parameter for the first elements not read correctly!"):
	print("Error: potential file for the first elements not selected!\n")
	sys.stderr.write(tempNameRun + ": Error: potential file for the first elements not selected!\n")
	exit()
elif(lastLine == "Error: potential parameter for the second elements not read correctly!"):
	print("Error: potential file for the second elements not selected!\n")
	sys.stderr.write(tempNameRun + ": Error: potential file for the second elements not selected!\n")
	exit()
	
file_leggi = open("leggi.in", "r")
cont = 0
T = 0
n = 0
name_file_pot = 0
for x in file_leggi:
	cont += 1
	if(cont == 2):
		x = x.split()
		name_file_pot = x[0]
	if(cont == 4):
		x = x.split()
		n = int(x[0])
	elif(cont == 10):
		x = x.split()
		T = float(x[0])
file_leggi.close()

N = n**3*4

name_file_El = []
A, p, q, xi, r0, rc00, rc11, rc01 = [], [], [], [], [], [], [], []
name_sim = 0
try:
	with open(name_file_pot, 'r') as file_pot:
		print(name_file_pot, " opened correctly")
		cont = 0
		for x in file_pot:
			cont += 1
			if(cont == 4):
				x = x.split()
				name_file_El.append(x[3])
				name_file_El.append(x[7])
				name_sim = x[11]
			if(cont == 6):
				x = x.split()
				p.append(float(x[0]))
				p.append(float(x[1]))
				p.append(float(x[2]))
			if(cont == 7):
				x = x.split()
				q.append(float(x[0]))
				q.append(float(x[1]))
				q.append(float(x[2]))
			if(cont == 8):
				x = x.split()
				A.append(float(x[0]))
				A.append(float(x[1]))
				A.append(float(x[2]))
			if(cont == 9):
				x = x.split()
				xi.append(float(x[0]))
				xi.append(float(x[1]))
				xi.append(float(x[2]))
			if(cont == 13):
				x = x.split()
				r0.append(float(x[0]))
				r0.append(float(x[1]))
			if(cont == 17):
				x = x.split()
				rc00.append(float(x[0]))
				rc00.append(float(x[1]))
			if(cont == 18):
				x = x.split()
				rc11.append(float(x[0]))
				rc11.append(float(x[1]))
			if(cont == 19):
				x = x.split()
				rc01.append(float(x[0]))
				rc01.append(float(x[1]))
except FileNotFoundError:
	print("File ", name_file_pot, " not found")
	sys.stderr.write(tempNameRun + ": File " + name_file_pot + " not found\n")
	exit()
except IOError:
	print("Error in opening ", name_file_pot)
	sys.stderr.write(tempNameRun + ": Error in opening " + name_file_pot + "\n")
	exit()
print("Parameter source: ", name_file_El)
T_inf = np.zeros(2)
H_inf = np.zeros(2)
T_sup = np.zeros(2)
H_sup = np.zeros(2)
pos_in_sortedMt = np.zeros(2)
err_H_pure = np.zeros(2)
j = 0
for file in name_file_El:
	matrix_temp = []
	try:
		with open('../pure_metals/thermo_files/thermo_'+file+'.dat') as file_El:
			print("File 'thermo_", file,".dat' opened correctly")
			cont = 0
			for x in file_El:
				cont +=1
				if(cont>2):
					x = x.split()
					matrix_temp.append([float(x[1]), float(x[5])])
	except FileNotFoundError:
		print("File thermo_", file, " not found")
		sys.stderr.write(tempNameRun + ": File thermo_" + file + " not found\n")
		exit()
	except IOError:
		print("Error in opening thermo_", file)
		sys.stderr.write(tempNameRun + ": Error in opening thermo_"+ file + "\n")
		exit()
	matrix = np.array(matrix_temp)
	sorted_matrix = matrix[matrix[:, 0].argsort()]
	#print(sorted_matrix)
	cont = 0
	for k in range(0, len(sorted_matrix)):
		i = sorted_matrix[k]
		if(i[0]<T):
			T_inf[j] = float(i[0])
			H_inf[j] = float(i[1])
			pos_in_sortedMt[j] = k
			cont += 1
		else:
			break
	if(T_inf[j] == 0):
		print("This temeprature range is not in the database!")
		sys.stderr.write(tempNameRun + ": This temeprature range is not in the database!\n")
		exit()
	pos_temp = int(pos_in_sortedMt[j])
	T_sup[j] = sorted_matrix[cont, 0]
	H_sup[j] = sorted_matrix[cont, 1]
	if(pos_in_sortedMt[j]<=10):
		x = sorted_matrix[:(pos_temp+10), 0]
		y = sorted_matrix[:(pos_temp+10), 1]
	elif(pos_in_sortedMt[j]>=(len(sorted_matrix)-10)):
		x = sorted_matrix[(pos_temp-10):, 0]
		y = sorted_matrix[(pos_temp-10):, 1]
	else:
		x = sorted_matrix[(pos_temp-10):(pos_temp+10), 0]
		y = sorted_matrix[(pos_temp-10):(pos_temp+10), 1]
	optPar, pcov = opt.curve_fit(fit_func, x, y)
	resi = y - (optPar[1]*x + optPar[0])
	err_H_pure[j] = np.std(resi)
	j += 1

H_fin_sing = np.zeros(2)
for i in range (0, 2):
	interp = np.zeros((1000, 2))
	x = np.zeros(2)
	y = np.zeros(2)
	x[0] = T_inf[i]
	x[1] = T_sup[i]
	y[0] = H_inf[i]
	y[1] = H_sup[i]
	linear_interp = interp1d(x, y, kind='linear', fill_value='extrapolate')
	interp[:, 0] = np.linspace(x[0], x[1], 1000)
	interp[:, 1] = linear_interp(interp[:, 0])
	for temp1 in interp:
		if temp1[0]<=T:
			H_fin_sing[i] = temp1[1]

H_sim = 0
delta_H_sim = 0
hc = 0
conc = 0

try:
	with open('properties.out', 'r') as file_prop:
		print("File 'properties.dat' opened correctly")
		cont = 0
		for x in file_prop:
			cont +=1
			if(cont>1):
				x = x.split()
				conc = float(x[1])
				H_sim = float(x[2])
				delta_H_sim = float(x[3])
				hc = float(x[4])
except FileNotFoundError:
	print("File 'properties.dat' not found")
	sys.stderr.write(tempNameRun + ": File 'properties.dat' not found\n")
	exit()
except IOError:
	print("Error in opening 'properties.dat'")
	sys.stderr.write(tempNameRun + ": Error in opening 'properties.dat'\n")
	exit()

H_mix = H_sim/N - ((conc/100) * H_fin_sing[0] + ((100-conc)/100)*H_fin_sing[1])
delta_H_mix = math.sqrt(((1-conc/100)*err_H_pure[1])**2 + ((conc/100)*err_H_pure[0])**2 + (delta_H_sim/N)**2)
try:
	with open('../properties_fin.dat', 'a') as file_out:
		print("File 'properties_fin.dat' opened correctly")
		file_out.write("\n" + name_sim + "\t" + str(p[0]) + "\t" + str(p[1]) + "\t" + str(p[2]) + "\t" + str(q[0]) + "\t" + str(q[1]) + "\t" + str(q[2]) + "\t" + str(A[0]) + "\t" + str(A[1]) + "\t" + str(A[2]) + "\t" + str(xi[0]) + "\t" + str(xi[1]) + "\t" + str(xi[2]) + "\t" + str(r0[0]) + "\t" + str(r0[1]) + "\t" + str(rc00[0]) + "\t" + str(rc11[0]) + "\t" + str(rc01[0]) + "\t" + str(rc00[1]) + "\t" + str(rc11[1]) + "\t" + str(rc01[1]) + "\t" + str(conc) + "\t" + str(T) + "\t" + str(H_mix) + "\t" + str(delta_H_mix) + "\t" + str(hc))
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
