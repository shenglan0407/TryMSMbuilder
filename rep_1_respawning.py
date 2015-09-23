#---------------------------------------------------------------------
# Author: Shenglan Qiao <shenglan@stanford.edu>
# Contributors:
# Copyright (c) 2015, Stanford University
# All rights reserved.
# 
# This script takes simulation data in replicate 1 of the b2AR-alprenolol binding 
# simulation and demonstrates how to determine the starting conformations if the next
# epoch of simulations. 

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------
import numpy as np
import pickle
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import utilities as util

import mdtraj as md
from mdtraj.geometry import distance as md_dist

from msmbuilder.featurizer import RawPositionsFeaturizer
from msmbuilder.cluster import KCenters
from msmbuilder.cluster import KMedoids
from msmbuilder.msm import MarkovStateModel
#-----------------------------------------------------------------------------
# matplotlib settings
#-----------------------------------------------------------------------------
mpl.rcParams['font.size'] = 24.0
#-----------------------------------------------------------------------------
# Code
#-----------------------------------------------------------------------------

sim_path = '/scratch/PI/rondror/MD_simulations/amber/b2AR_ligand_binding/alprenolol/ten_ligands/production/ten_ligands/1/'
sim_file = 'b2AR_ALP_Prod1to9_skip1_reimaged.nc'
topology = '/scratch/PI/rondror/MD_simulations/amber/b2AR_ligand_binding/alprenolol/ten_ligands/system.psf'

traj = md.load(sim_path+sim_file,top = topology, stride = 10)
print traj