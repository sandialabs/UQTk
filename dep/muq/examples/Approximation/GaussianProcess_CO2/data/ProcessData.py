import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import h5py

# Read in raw datafile
csv_fileName = 'co2_weekly_mlo.csv'
df = pd.read_csv(csv_fileName, skiprows=49, usecols=[3,4], delim_whitespace=True, header=None)
df.columns = ['Date','Concentrations']

#plt.plot(df['Date'], df['Concentrations'])
#plt.show()

# Write cleaned data to h5 file
h5_filename = 'MaunaLoaCO2.h5'
f = h5py.File(h5_filename)
f['/Weekly/Dates'] = df['Date']
f['/Weekly/Concentrations'] = df['Concentrations']

# Define HDF attributes
f['/Weekly'].attrs['Mauna Loa Observatory Website URL'] = 'https://www.esrl.noaa.gov/gmd/ccgg/trends/data.html'
f['/Weekly'].attrs['Weekly Data Text File'] = 'ftp://aftp.cmdl.noaa.gov/products/trends/co2/co2_weekly_mlo.txt'
f['/Weekly'].attrs['Save Weekly Data'] = 'Save weekly data as CSV file'
f['/Weekly/Concentrations'].attrs['Units'] = 'CO2 mole fraction as parts per million (ppm)'
f['/Weekly/Dates'].attrs['Units'] = '(year,month,day) as a decimal'
