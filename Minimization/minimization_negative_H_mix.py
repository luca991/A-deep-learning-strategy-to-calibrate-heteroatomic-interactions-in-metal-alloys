import pandas as pd
import math
import random as rnd
from collections import OrderedDict
from tensorflow.keras.models import load_model
from sklearn.preprocessing import LabelEncoder
from scipy.optimize import minimize
import tensorflow as tf
import os
import sys

def get_min(name, df):
	row = df[df['name'] == name]
	if not row.empty:
		min_value = row['min'].values[0]
		return float(min_value)
	else:
		return None
	
def get_max(name, df):
	row = df[df['name'] == name]
	if not row.empty:
		max_value = row['max'].values[0]
		return float(max_value)
	else:
		return None

def remove_nth_el(original_list, n):
	return [original_list[i] for i in range(len(original_list)) if i!=n]

#Input: 'input' is a 1-line dataframe and 'data_in_tot' is the entire database with only the colimns selected in 'sel_col'  
def minimization(input, resc_data):
	#Transforming the input dataframe into a list
	input_list = input.iloc[0].tolist()
	print(input_list)

	#Defining two different lists: 'fix' for fixed values dunring minimization (pure parameters, physical quantities) and 'var' list for values that may vary during minimization (mixed interaction parameters: a, p, q, xi, r0)
	fix = input_list[:20]
	var = input_list[20:]
	print(fix)
	print(var)
	
	#Deifing the minimization bounds for p, q, r0 between the pures using these steps:
	#  - Select the min and the max value from the pure ones
	#  - Use an if loop to identify whether the min comes from El0 or El1 and perform the following steps:
	#    - Rescale the min and max values respectively depending on the metal they correspond to to obtain the 'real' values
	#    - Rescale the 'real' max and min values with respect to the values of the corresponding parameter of the mixed interaction. This last step must be done in order to define the correct range in which these parameters can vary during minimization
	r_min = min(fix[8], fix[9])
	r_max = max(fix[8], fix[9])
	#print(r_min, "\t", r_max)

	if(r_min == fix[8]):
		r_min_or = r_min * (get_max('r0_00', resc_data) - get_min('r0_00', resc_data)) + get_min('r0_00', resc_data)
		r_max_or = r_max * (get_max('r0_11', resc_data) - get_min('r0_11', resc_data)) + get_min('r0_11', resc_data)
		r_min = (r_min_or - get_min('r0_01', resc_data)) /(get_max('r0_01', resc_data) - get_min('r0_01', resc_data))
		r_max = (r_max_or - get_min('r0_01', resc_data)) /(get_max('r0_01', resc_data) - get_min('r0_01', resc_data))
	elif(r_min == fix[9]):
		r_min_or = r_min * (get_max('r0_11', resc_data) - get_min('r0_11', resc_data)) + get_min('r0_11', resc_data)
		r_max_or = r_max * (get_max('r0_00', resc_data) - get_min('r0_00', resc_data)) + get_min('r0_00', resc_data)
		r_min = (r_min_or - get_min('r0_01', resc_data)) /(get_max('r0_01', resc_data) - get_min('r0_01', resc_data))
		r_max = (r_max_or - get_min('r0_01', resc_data)) /(get_max('r0_01', resc_data) - get_min('r0_01', resc_data))
	
	p_min = min(fix[6], fix[7])
	p_max = max(fix[6], fix[7])
	#print(p_min, "\t", p_max)
	
	if(p_min == fix[6]):
		p_min_or = p_min * (get_max('p_00', resc_data) - get_min('p_00', resc_data)) + get_min('p_00', resc_data)
		p_max_or = p_max * (get_max('p_11', resc_data) - get_min('p_11', resc_data)) + get_min('p_11', resc_data)
		p_min = (p_min_or - get_min('p_01', resc_data)) /(get_max('p_01', resc_data) - get_min('p_01', resc_data))
		p_max = (p_max_or - get_min('p_01', resc_data)) /(get_max('p_01', resc_data) - get_min('p_01', resc_data))
	if(p_min == fix[7]):
		p_min_or = p_min * (get_max('p_11', resc_data) - get_min('p_11', resc_data)) + get_min('p_11', resc_data)
		p_max_or = p_max * (get_max('p_00', resc_data) - get_min('p_00', resc_data)) + get_min('p_00', resc_data)
		p_min = (p_min_or - get_min('p_01', resc_data)) /(get_max('p_01', resc_data) - get_min('p_01', resc_data))
		p_max = (p_max_or - get_min('p_01', resc_data)) /(get_max('p_01', resc_data) - get_min('p_01', resc_data))
	
	q_min = min(fix[2], fix[3])
	q_max = max(fix[2], fix[3])
	if(q_min == fix[2]):
		q_min_or = q_min * (get_max('q_00', resc_data) - get_min('q_00', resc_data)) + get_min('q_00', resc_data)
		q_max_or = q_max * (get_max('q_11', resc_data) - get_min('q_11', resc_data)) + get_min('q_11', resc_data)
		q_min = (q_min_or - get_min('q_01', resc_data)) /(get_max('q_01', resc_data) - get_min('q_01', resc_data))
		q_max = (q_max_or - get_min('q_01', resc_data)) /(get_max('q_01', resc_data) - get_min('q_01', resc_data))
	elif(q_min == fix[3]):
		q_min_or = q_min * (get_max('q_11', resc_data) - get_min('q_11', resc_data)) + get_min('q_11', resc_data)
		q_max_or = q_max * (get_max('q_00', resc_data) - get_min('q_00', resc_data)) + get_min('q_00', resc_data)
		q_min = (q_min_or - get_min('q_01', resc_data)) /(get_max('q_01', resc_data) - get_min('q_01', resc_data))
		q_max = (q_max_or - get_min('q_01', resc_data)) /(get_max('q_01', resc_data) - get_min('q_01', resc_data))

	#Perform the minimization of the 'obj_function' in term of the vector 'var'
	res = minimize(obj_function, var, args=(fix), method='Nelder-Mead', tol = 1e-4, bounds=((r_min, r_max), (q_min, q_max), (p_min, p_max), (0,1), (0,1)))
	var_new = res.x
	val = res.fun
	print(res.message)
	return var_new, val

def obj_function(var, fix):
    var = var.tolist()

    #To also vary the cut-off radii of the mixed interaction as the r0 of the mixed interaction varies during the minimization we carry out the following steps, using the data in the vector saved as r0:
    #  - Rescale r0 of the mixed interaction in order to get value between the min and the max in the database -> unscaled
    #  - Calculate the cut-off radii with the role used to create the database 
    #  - Reascale the new cut-off radii respect to the min and the max values of the cut-off radii of the mixed interaction in the database 
        
    print('H_mix_5 prediction')
    fix_tempo = remove_nth_el(fix, 11)
    y_true = fix[11]
    fix_var = fix_tempo[:10] + [fix_tempo[10]] + var[:] + fix_tempo[17:19]
    input_tensor = tf.convert_to_tensor([fix_var], dtype=tf.float32)  
    out = Hmix5(input_tensor, training=False)
    y_pred = out.numpy().flatten()[0]
    loss_1 = abs(y_pred-y_true)/y_true
    print(y_true, "\t", y_pred)
    print("difference: ", loss_1, "\n")

    print('H_mix_95 prediction')
    fix_tempo = remove_nth_el(fix, 12)
    y_true = fix[12]
    fix_var = fix_tempo[:10] + [fix_tempo[10]] + var[:] + fix_tempo[17:19]
    input_tensor = tf.convert_to_tensor([fix_var], dtype=tf.float32)  
    out = Hmix95(input_tensor, training=False)
    y_pred = out.numpy().flatten()[0]
    loss_2 = abs(y_pred-y_true)/y_true
    print(y_true, "\t", y_pred)
    print("difference: ", loss_2, "\n")

    print('latt_5 prediction')
    fix_tempo = remove_nth_el(fix, 15)
    y_true = fix[15]
    fix_var = fix_tempo[:10] + [fix_tempo[14]] + var[:] + fix_tempo[17:19]
    x_data = pd.DataFrame([fix_var], columns=["A_00", "A_11", "q_00", "q_11", "xi_00", "xi_11", "p_00", "p_11", "r0_00", "r0_11", "T", "r0_01", "q_01", "p_01", "xi_01","A_01", "latt_type_El0_en", "latt_type_El1_en"])
    # y_pred = latt_5.predict(x_data, verbose=0)
    input_tensor = tf.convert_to_tensor([fix_var], dtype=tf.float32)  
    out = latt_5(input_tensor, training=False)
    y_pred = out.numpy().flatten()[0]
    loss_3 = abs(y_pred-y_true)/y_true
    print(y_true, "\t", y_pred)
    print("difference: ", loss_3, "\n")

    print('latt_95 prediction')
    fix_tempo = remove_nth_el(fix, 16)
    y_true = fix[16]
    fix_var = fix_tempo[:10] + [fix_tempo[14]] + var[:] + fix_tempo[17:19]
    input_tensor = tf.convert_to_tensor([fix_var], dtype=tf.float32)  
    out = latt_95(input_tensor, training=False)
    y_pred = out.numpy().flatten()[0]
    loss_4 = abs(y_pred-y_true)/y_true
    print(y_true, "\t", y_pred)
    print("difference: ", loss_4, "\n")

    print('latt_50 prediction')
    fix_tempo = remove_nth_el(fix, 17)
    y_true = fix[17]
    fix_var = fix_tempo[:10] + [fix_tempo[14]] + var[:] + fix_tempo[17:19]
    input_tensor = tf.convert_to_tensor([fix_var], dtype=tf.float32)  
    out = latt_50(input_tensor, training=False)
    y_pred = out.numpy().flatten()[0]
    loss_5 = abs(y_pred-y_true)/y_true
    print(y_true, "\t", y_pred)
    print("difference: ", loss_5, "\n")


    print('H_mix_50 prediction')
    fix_tempo = remove_nth_el(fix, 13)
    y_true = fix[13]
    print(y_true)
    fix_var = fix_tempo[:10] + [fix_tempo[10]] + var[:] + fix_tempo[17:19]
    input_tensor = tf.convert_to_tensor([fix_var], dtype=tf.float32)  
    out = Hmix50(input_tensor, training=False)
    y_pred = out.numpy().flatten()[0]
    loss_6 = abs(y_pred-y_true)/y_true
    print(y_true, "\t", y_pred)
    print("difference: ", loss_6, "\n")

    loss = (((2*(loss_1)**2 + 2*(loss_2)**2 + loss_3**2 + loss_4**2 + loss_5**2 + 2*(loss_6)**2)/9)**0.5)
    print('r0_q_p_xi_A_loss:\t', var[0], '\t', var[1], '\t', var[2], '\t', var[3], '\t', var[4], '\t', loss)
    
    print('Loss: ', loss, "\n")

    return loss

### MAIN ###
alloy = sys.argv[1]

El0 = alloy[:2]
El1 = alloy[2:]
name_alloy = El0+'-'+El1
print(El0+'\t'+El1+'\t'+name_alloy)

#Define the number of indipendent minimizations
n_m = 15

# Import NN

Hmix5 = load_model("H_mix_5.h5")      

Hmix95 = load_model("H_mix_95.h5")   

latt_5= load_model("latt_par_5.h5")  

latt_95= load_model("latt_par_95.h5") 

Hmix50 = load_model("H_mix_50.h5")   

latt_50= load_model("latt_par_50.h5") 

label_encoder = LabelEncoder()


quanties=['A_00', 'A_11', 'q_00', 'q_11', 'xi_00', 'xi_11', 'p_00', 'p_11', 'r0_00', 'r0_11', 'T', 'H_mix_5', 'H_mix_95', 'H_mix_50', 'T_latt_par', "latt_par_5", "latt_par_95", "latt_par_50", "LS_0", "LS_1", 'r0_01', 'q_01', 'p_01', 'xi_01', 'A_01']

# Quantities for which you minimize
to_skip = ['A_01', 'p_01', 'q_01', 'xi_01', 'r0_01']

to_fill = [c for c in quanties if c not in to_skip]

# Quantities which, if they are not presenti in 'to_skip', are assigned the mean value
av = ['A_01', 'p_01', 'q_01', 'xi_01', 'r0_01']

to_the_av = [c for c in av if c not in to_skip]

x = pd.DataFrame(columns=quanties)

resc_data = pd.read_csv('min_max_values.dat', sep='\t', header=None, names=['name', 'min', 'max'])

pure=pd.read_csv('Pure_elements.dat', sep='\t')
pure_col = pure.columns.to_list()

byn = pd.read_csv('Byn_prop_resc.dat', sep='\t')
byn_col = byn.columns.to_list()

x_col = x.columns.to_list()
print("x_col:\n", x_col)

# Value assigment
#The order for x and y must be respected:
# latt_type_El0_en & latt_type_El1_en must be:
# 		0 if the lattice is fcc
# 		1 if the lattice is hcp

colonne_mappa = {
    "p": "p",
    "q": "q",
    "A": "a",
    "xi": "xi",
    "r0": "r0",
    "LS": "LS"
}

valori_riga = {}
for i in x_col:
	#print(i)

	if '_00' in i:
		elemento = El0
		chiave_base = i.replace('_00', '')
	elif '_11' in i:
		elemento = El1
		chiave_base = i.replace('_11', '')
	elif '_01' in i:
		elemento = 'Null'
		chiave_base = i.replace('_01', '')
	elif i.endswith('_0'):
		elemento = El0
		chiave_base = i.replace('_0', '')
	elif i.endswith('_1'):
		elemento = El1
		chiave_base = i.replace('_1', '')

	colonna = colonne_mappa.get(chiave_base)

	if elemento!='Null':
		if colonna:
			valore = pure.loc[pure['El'] == elemento, colonna].values[0]
			if colonna=='LS':
				if valore=='fcc':
					valori_riga[i] = 0
				elif valore=='hcp':
					valori_riga[i] = 1
			else: 
				valori_riga[i] = valore
	else:
		if colonna:
			valore=byn.loc[byn["El0-El1"] == name_alloy, colonna].values[0]
			valori_riga[i] = valore

	if 'rc1' in i:
		if '_00' in i:
			valori_riga[i] = valori_riga['r0_00']*(3**0.5)
		elif '_11' in i:
			valori_riga[i] = valori_riga['r0_11']*(3**0.5)
		elif '_01' in i:
			valori_riga[i] = valori_riga['r0_01']*(3**0.5)
	elif 'rc2' in i:
		if '_00' in i:
			valori_riga[i] = valori_riga['r0_00']*(4**0.5)
		elif '_11' in i:
			valori_riga[i] = valori_riga['r0_11']*(4**0.5)
		elif '_01' in i:
			valori_riga[i] = valori_riga['r0_01']*(4**0.5)

	if i in byn_col:
		valore=byn.loc[byn["El0-El1"] == name_alloy, i].values[0]
		valori_riga[i] = valore

for i in to_the_av:
	print(i)
	
	if "p" in i:
		valori_riga[i] = ((pure.loc[pure['El']==El0, 'p'].values[0]+pure.loc[pure['El']==El1, 'p'].values[0])/2)
	if "q" in i:
		valori_riga[i] = ((pure.loc[pure['El']==El0, 'q'].values[0]+pure.loc[pure['El']==El1, 'q'].values[0])/2)
	if "A" in i:
		valori_riga[i] = ((pure.loc[pure['El']==El0, 'A'].values[0]+pure.loc[pure['El']==El1, 'A'].values[0])/2)
	if "xi" in i:
		valori_riga[i] = ((pure.loc[pure['El']==El0, 'xi'].values[0]+pure.loc[pure['El']==El1, 'xi'].values[0])/2)
	if "r0" in i:
		valori_riga[i] = ((pure.loc[pure['El']==El0, 'r0'].values[0]+pure.loc[pure['El']==El1, 'r0'].values[0])/2)
	
valori_riga = pd.DataFrame([valori_riga])
print(valori_riga)

# Define the ranges for the extraction of initial values for p, q and r, the extremes being given by the values of the pure 
r_min = min(valori_riga["r0_00"][0], valori_riga["r0_11"][0])
r_max = max(valori_riga["r0_00"][0], valori_riga["r0_11"][0])
q_min = min(valori_riga["q_00"][0], valori_riga["q_11"][0])
q_max = max(valori_riga["q_00"][0], valori_riga["q_11"][0])
p_min = min(valori_riga["p_00"][0], valori_riga["p_11"][0])
p_max = max(valori_riga["p_00"][0], valori_riga["p_11"][0])

print(r_min, "\t", r_max)
print(q_min, "\t", q_max)
print(p_min, "\t", p_max)

print(valori_riga)

mse = []
x_pred_norm_min = []
for i in range(0, n_m):
    print("\nIteration: ", i)
    data_in = valori_riga.iloc[[0]]
    print("\nData_in_0:\n", valori_riga)
    print("\nData_in:\n", data_in)

    #Random extraction of the starting value of the parameters to be minimized
    data_in.loc[0, "r0_01"] = float(rnd.uniform(r_min, r_max))

    data_in.loc[0, "q_01"] = float(rnd.uniform(q_min, q_max))

    data_in.loc[0, "p_01"] = float(rnd.uniform(p_min, p_max))

    data_in["xi_01"] =rnd.uniform(float(resc_data.loc[resc_data['name']=='xi_00'].iloc[0,1]), float(resc_data.loc[resc_data['name']=='xi_11'].iloc[0,2]))
    data_in["A_01"] = rnd.uniform(float(resc_data.loc[resc_data['name']=='A_00'].iloc[0,1]), float(resc_data.loc[resc_data['name']=='A_11'].iloc[0,2]))  
    print("\nData_in:\n", data_in)

    #Normalisation of all quantities (defined in "sel_col"), of the simulation chosen in "id" for minimization, with respect to the same quantities in the entire database
    data_in_norm = data_in
    print("\nData_in_norm:\n", data_in_norm)
    print("\nData_in:\n", data_in)
    for index, row in resc_data.iterrows():
        name = row['name']
        #print(name)
        min_val = row['min']
        max_val = row['max']
        if name in data_in_norm.columns:
            data_in_norm[name] = (float(data_in_norm[name]) - float(min_val))/(float(max_val)-float(min_val))
            if name == 'T':
                print(name)
                data_in_norm['T_latt_par'] = (float(data_in_norm['T_latt_par']) - float(min_val))/(float(max_val)-float(min_val))
    print("\nData_in_norm:\n", data_in_norm)
    print("\nData_in:\n", data_in)
    print("\nData_in_0:\n", valori_riga)
    print(data_in["r0_01"], "\t", data_in["q_01"], "\t", data_in["p_01"], "\t", data_in["A_01"], "\t", data_in["xi_01"])

    #Minimization is performed
    x_pred_norm, toll = minimization(data_in_norm, resc_data)
    mse.append(toll)
    x_pred_norm_min.append(x_pred_norm)
    print(data_in_norm)
    print(valori_riga)


mse_min = 1e20
x_pred_norm_min_tab = []
x_pred_l = []
par_fin = []
for i in range(0, len(mse)):
	x_pred_norm_min_tab.append(x_pred_norm_min[i].tolist())
	x_pred_norm_min_df = pd.DataFrame([x_pred_norm_min_tab[i]], columns=["r0_01", "q_01", "p_01","xi_01","A_01"])
	#x_pred_resc = data[["r0_01", "q_01", "p_01","xi_01","A_01"]]
	#print(data_in[["r0_01", "q_01", "p_01","xi_01","A_01"]])
	x_pred = x_pred_norm_min_df
	for index, row in resc_data.iterrows():
		name = row['name']
		min_val = row['min']
		max_val = row['max']
		if name in x_pred.columns:
			x_pred[name] = x_pred[name] *(float(max_val)-float(min_val)) + float(min_val)	
	x_pred_l.append(x_pred.iloc[0].tolist())
	print(mse[i])
	print(x_pred_l[i], "\n")
	if(mse[i] < mse_min):
		mse_min = mse[i]
		par_fin = x_pred_l[i]

#for i in range(0, len(f)):
print(mse_min)
#print(data_in_tot.loc[0, ["r0_01", "q_01", "p_01","xi_01","A_01"]])
print(par_fin)
