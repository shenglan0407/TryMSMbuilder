# Author: Shenglan Qiao <shenglan@stanford.edu>
# Contributors:
# Copyright (c) 2015, Stanford University
# All rights reserved.

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------


import numpy as np
import os

import mdtraj as md
from msmbuilder.featurizer import RawPositionsFeaturizer


#-----------------------------------------------------------------------------
# Code
#-----------------------------------------------------------------------------


def load_trajs(run_dir, parent_dir,top_dir,load_stride = None):
    """Loads trajectories from the same run with the same condition 
    into a list using mdtraj.load_dcd
    
    Keyword arguments
    -----------------
    load_stride : int, default = None 
    mdtraj to read every stride-th frame
    
    Arguments
    ---------
    run_dir :  str
    the directory where all files for a single simulation run lives
    parent_dir : str
    the directory where all simulations are stored
    
    Returns
    -------
    List, mdtraj objects that are all from the same simulation run, in order
    """
    
    trajs = []
    for this_file in os.listdir(parent_dir+run_dir):
        if this_file.endswith('.dcd'):
            topology = top_dir + '/pnas2011a-' + this_file.split('-')[1] \
            + '-0-no-water-no-lipid.pdb'
        
            if load_stride == None:
                this_traj = md.load_dcd(parent_dir + run_dir + '/' + this_file
                                    , top = topology
                                    )
            else:
                this_traj = md.load_dcd(parent_dir + run_dir + '/' + this_file
                                    , top = topology
                                    , stride = load_stride
                                    )
            trajs.append(this_traj)
    
    return trajs

def average_center_seq(sequence):
    """
    """
    
    center_sequence = []
    for c_pos in sequence:
        this_x = np.mean(c_pos[::3])
        this_y = np.mean(c_pos[1::3])
        this_z = np.mean(c_pos[2::3])
        center_sequence.append([this_x,this_y,this_z])
    return np.array(center_sequence)
    
def calc_time_step(times_file_dir,stride = None): 
    '''Returns the time step between each loaded frame in ns
    '''
    times_file = open(times_file_dir)
    first_line = times_file.readlines()[1].split(',')
    sim_time_step = int(first_line[2]) #ps, this is the time step used  in simulations
    if stride == None:
        return sim_time_step*1e-3
    else:
        return sim_time_step*1e-3*stride
        
def featurize_RawPos(indices,trajs,average=False):
    '''takes indices of atoms in ligand and return sequences of average position
    each sequence is n by m where n is the number of frames in the simulations and m = 3 for now
    trajs is a list of mdtraj objects. Right now the are actually all from the same simulation run.
    That's why the sequences for the same ligand are stitched together.
    '''
    
    sequences = []
    for ii in indices:
        RawPosFeaturizer = RawPositionsFeaturizer(ii)
        this_seq = []
        for this_traj in trajs:
            this_RawPosSeq = RawPosFeaturizer.partial_transform(this_traj)
            if average:
                this_seq.extend(average_center_seq(this_RawPosSeq))
            else:
                this_seq.extend(this_RawPosSeq)
    
        this_seq = np.array(this_seq)
        sequences.append(this_seq)
    
    return sequences

def convert_to_pdb(filepath,savepath):
    crds = np.loadtxt(filepath)
    outfile = open(savepath,'wb')
    atom_num = crds.shape[0]

    outfile.write('%s     %s \n' % ('TITLE', 'cluster centers'))
    outfile.write('%s        %d \n' % ('MODEL' , 1))
    for ii in range(atom_num):
        outfile.write("%6s%5.2d %4s %s %s  %s   %7.3f %7.3f %7.3f %s %5s %12s  \n" % ('HETATM', ii,'C','CLU','I',str(ii)+'D'
        ,crds[ii,0]*10,crds[ii,1]*10,crds[ii,2]*10,'1.00','0.00','C')) 
    #the factor of 10 multiplied to the coordinate is a conversion from nm to angstrom, to match the .mae file provided by DesRes
    
    outfile.write('%s \n' % 'ENDMDL')
    outfile.write('%s    \n' % 'END')
    
    print('Done converting! \nConverted .pdb file lives in %s.' % savepath)