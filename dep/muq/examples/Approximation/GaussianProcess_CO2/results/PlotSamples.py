import h5py
import matplotlib.pyplot as plt
import numpy as np
# Read origional data
f1 = h5py.File('../data/MaunaLoaCO2.h5')
date1 = f1['/Weekly/Dates']
conc1 = f1['/Weekly/Concentrations']

# Read the prediction results
f2 = h5py.File('CO2_Prediction.h5')
date2 = f2['/Predict/Dates']
conc2_mu = f2['/Predict/Concentrations']
conc2_var = f2['/Predict/ConcentrationVariance']

# Define HDF attributes
f2['/Predict'].attrs['Source'] = 'Results from GaussianProcess_CO2.cpp'
f2['/Predict/Concentrations'].attrs['Units'] = 'CO2 mole fraction as parts per million (ppm)'
f2['/Predict/Dates'].attrs['Units'] = 'decimal year'

# Plot them together
plt.fill_between(date2[0,:], conc2_mu[0,:]+2.0*np.sqrt(conc2_var[0,:]), conc2_mu[0,:]-2.0*np.sqrt(conc2_var[0,:]),  facecolor='green', interpolate=True, label='Posterior $\pm 2\sigma$')
plt.plot(date2[0,:], conc2_mu[0,:],label='Posterior Mean')
plt.plot(date1, conc1,'.k',markersize=3, label='Observations')
#plt.legend()
plt.show()
