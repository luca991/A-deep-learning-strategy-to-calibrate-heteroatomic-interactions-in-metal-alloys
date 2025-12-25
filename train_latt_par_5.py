import pandas as pd
import tensorflow as tf
import numpy as np
from tensorflow.keras.models import Sequential, load_model
from tensorflow.keras.layers import Dense, Dropout, CategoryEncoding
from tensorflow.keras.callbacks import CSVLogger, ModelCheckpoint
from sklearn.preprocessing import MinMaxScaler, StandardScaler, LabelEncoder
import matplotlib.pyplot as plt
import keras
import random

#Seed generation
generated_seed = random.randint(0, 2**32-1)
seed = generated_seed
#seed = int()
file_out = open("NN_params.dat", "w")
file_out.write(f"NN seed: {seed}\n")
file_out.close()

#Set seed
np.random.seed(seed)
random.seed(seed)
tf.random.set_seed(seed)

#Open database
data = pd.read_csv('data_ML_article.dat', sep='\t')
resc_data = pd.read_csv("min_max_values.dat", sep='\t', header=None, names=['name', 'min', 'max'])
print(data)

#Selection data
x = data[["A_00", "A_11", "q_00", "q_11", "xi_00", "xi_11", "p_00", "p_11", "r0_00", "r0_11", "T", "r0_01", "q_01", "p_01", "xi_01","A_01"]]
y = data[["latt_par_5"]]
label_encoder = LabelEncoder()
x["latt_type_El0_en"] = label_encoder.fit_transform(data["LS_0"])
x["latt_type_El1_en"] = label_encoder.fit_transform(data["LS_1"])
print(x)

for index, row in resc_data.iterrows():
	name = row['name']
	min_val = row['min']
	max_val = row['max']
	if name in x.columns:
		x[name] = (x[name] - float(min_val))/(float(max_val)-float(min_val))
	if name in y.columns:
		y[name] = (y[name] - float(min_val))/(float(max_val)-float(min_val))
normalized_x = x
normalized_y = y

input_shape = (normalized_x.shape[1],)


normalized_x_train = normalized_x[:37500]
normalized_y_train = normalized_y[:37500]

#Design NN
model = Sequential()
model.add(Dense(256, activation='relu', input_shape=input_shape))
model.add(Dense(128, activation='relu'))
model.add(Dense(64, activation='relu'))
model.add(Dense(32, activation='relu'))
model.add(Dense(16, activation='relu'))
model.add(Dense(8, activation='relu'))
model.add(Dense(4, activation='relu'))
model.add(Dense(normalized_y.shape[1], activation='sigmoid'))

optimizer = keras.optimizers.Adam()

model.compile(optimizer=optimizer, loss=keras.losses.MeanAbsoluteError(), metrics=['accuracy'])

check_point = ModelCheckpoint(
    filepath='model_test_1_best_{epoch:04d}.h5',
    monitor='val_loss',
    verbose=1,
    mode='min',
    save_best_only=True,
    save_freq='epoch'
)

model.build()
model.summary()
csv_logger = CSVLogger('training_1.log')

history = model.fit(normalized_x_train, normalized_y_train, epochs=5000, verbose=2, validation_split=0.15, callbacks=[check_point, csv_logger])
history_dict = history.history

#Save NN in final configuration
model.save('model_test_1_fin.h5')
