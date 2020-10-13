import h5py
import pandas as pd
import numpy as np

saveFile = h5py.File('LTA01.h5', header=0)
readFile = 'LTA01.csv'

data = pd.read_csv(readFile)
units = data.iloc[0,:]
data = data.iloc[1:,:]

strain = np.array(data['Strain (Gauge1)'],dtype=float)
stress = np.array(data['Engr. Stress'],dtype=float)

# Remove the dip typically used to measure Young's Modulus
newStress = []
newStrain = []

prevStress = stress[0]-1e-1
for i in range(stress.shape[0]):
    if((stress[i] > prevStress) & (strain[i]<0.22)):
        newStress.append(stress[i])
        newStrain.append(strain[i])
        prevStress = stress[i]

newStress = np.array(newStress)
newStrain = np.array(newStrain)


saveFile['/Strain'] = newStrain
saveFile['/Stress'] = newStress

saveFile['/'].attrs['Source'] = 'This data is a processed version of data provided with the 2017 Sandia Fracture Challenge (SFC3) in file ' + readFile
saveFile['/Strain'].attrs['Units'] = str(units['Strain (Gauge1)'])
saveFile['/Stress'].attrs['Units'] = str(units['Engr. Stress'])

