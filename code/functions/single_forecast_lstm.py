# Function to perform a single one-step-ahead forecast using an LSTM neural network.
# Code modified from https://github.com/MITEcology/JRSI-Cenci-Medeiros-Sugihara-Saavedra-2019

import numpy as np
import pandas as pd
import importlib, os, sys
from sklearn.metrics import mean_squared_error
np.random.seed(5)
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'

def single_forecast_lstm(training_data, n_test, cv):
    forecast = lstm_forecast(training_data, n_test, do_cv = cv)
    return forecast
