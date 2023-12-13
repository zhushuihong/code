import numpy as np
import math
import MDAnalysis as mda
from MDAnalysis.analysis.base import (AnalysisBase, AnalysisFromFunction, analysis_class)
import MDAnalysis.transformations
import pandas as pd
import matplotlib.pyplot as plt

dataname = "system.data"
dumpname = "shear_ep0.4t1.4.dcd"
output_file = "shear_ep0.4t1.4.png"
csv_file = "end_to_end_distance_ep0.4t1.4.csv"

# calculate the end to end distance
def end_end_distance(batomgroup, eatomgroup, polymer_num):
    #coordinates change for each frame
    bcoordinate  = batomgroup.positions
    ecoordinate  = eatomgroup.positions
    

    # get squared distance between two group
    di_sq = (bcoordinate - ecoordinate)**2

    # sum the unweighted positons
    sq = np.sum(di_sq, axis=1)  
    sq_x = di_sq[:,0]
    sq_y = di_sq[:,1]
    sq_z = di_sq[:,2]
    
    # make into array
    sq_di = np.array([sq, sq_x, sq_y, sq_z])
    
    # weight positions
    per_di = np.sum(sq_di, axis=1)/polymer_num
    # square root and return
    return np.sqrt(per_di)

# Define constants
natoms = 10000
nbeads = 100
npoly = natoms/nbeads

# read file 
u = mda.Universe(dataname, dumpname)
#solveing pbc
#ag = u.atoms
#transform = mda.transformations.unwrap(ag)
#u.trajectory.add_transformations(transform)


begin_atoms = u.atoms[0:10000:100]
end_atoms = u.atoms[99:10000:100]

print(begin_atoms)
print(begin_atoms.indices)
print(end_atoms)
print(end_atoms.indices)

# write my analysis
"""
rog = AnalysisFromFunction(end_end_distance, u.trajectory, begin_atoms, end_atoms, npoly)
rog.run()

print(rog.results.timeseries.shape)

for col, label in zip(rog.results['timeseries'].T, labels):
    plt.plot(col, lw=2, label=label)

plt.ylabel('Radius of gyration (ai)')
plt.xlabel('Frame')
plt.legend(loc="upper right")
plt.savefig("end_to_end_distance")
"""
labels = ['all', 'x-axis', 'y-axis', 'z-axis']
rog1 = AnalysisFromFunction(end_end_distance, u.trajectory, begin_atoms, end_atoms, npoly)

time = np.arange(0,1481.6,3.2)
rog1.run(start=0, stop=926, step=2)
#print("This is result")
#print(rog1.results.timeseries.T)

columns = ['end to end distance (all)','end to end distance (x-axis)', \
                   'end to end distance (y-axis)', 'end to end distance (z-axis)' ]

ll = pd.DataFrame(rog1.results.timeseries, columns=columns)
ll['Time(ps)'] = time
ll.to_csv(csv_file)

for col, label in zip(rog1.results.timeseries.T, labels):
    plt.plot(time, col*3, lw=2, label=label)
plt.ylabel('End to end distance (nm)')
plt.xlabel('Time(ns)')
plt.legend(loc="upper right")
plt.savefig(output_file)
    

