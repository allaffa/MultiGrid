#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  3 11:45:37 2020

@author: 7ml
"""

import numpy as np
import keras 
from keras.models import Sequential
from keras.layers import Dense
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt 


meshes = np.loadtxt('meshes.txt', skiprows=0)
#gradients = np.loadtxt('gradients.txt', skiprows=0)
solutions = np.loadtxt('solutions.txt', skiprows=0)

input_dim = solutions.shape[1]
output_dim = solutions.shape[1]
training_size = solutions.shape[0]

#Split between train and validaiton set
x_train, x_test, y_train, y_test = train_test_split(solutions, meshes, test_size=0.2, random_state=42)

#build neural network
model = Sequential()
model.add(Dense(100, input_dim=input_dim, activation='relu'))
model.add(Dense(100, activation='relu'))
model.add(Dense(100, activation='relu'))
model.add(Dense(100, activation='relu'))
model.add(Dense(output_dim, activation=None))

# compile the keras model
model.compile(loss='mse', optimizer='adam', metrics=['accuracy'])

# fit the keras model on the dataset
model.fit(x_train, y_train, epochs=1000, batch_size=50, validation_split = 0.2)

results = model.evaluate(x_test, y_test, batch_size=128)
print('test loss, test acc:', results)

preds = model.predict(x_test)