
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 15 23:49:36 2017

@author: casper
"""
import matlab.engine
eng = matlab.engine.start_matlab()
eng.cd('../ddCRP-master/')
from loadData import *
datasetPath     = '../../ZebrafishRecordings/photo_omr_blob_etc_data/20150120_2_1_photo_OMR_prey_blob_blue_cy74_6d_20150120_220356/Registered/'
# datasetPath      = '../../ZebrafishRecordings/sample_data_code_spontaneous/20140805_1_3_cy74_6d_spon_20140805_170009/'
# datasetPath     = '../../ZebrafishRecordings/photo_omr_blob_etc_data/20150122_2_2_cy74_6d_photo_OMR_prey_blob_blue_20150122_190028/Registered/'
#    datasetPath     = '../ZebrafishRecordings/photo_omr_blob_etc_data/20150122_2_2_cy74_6d_photo_OMR_prey_blob_blue_20150122_190028/Registered/'
dataName        = datasetPath.split('/')[-3] # assuming it has registered subfolder
                                      # otherwise edit this
dataPath        = '../Data/'

loadThis        = dataPath + dataName + '.hdf5'


frameTurn        = loadFrameTurn(datasetPath + 'frame_turn.mat')
conditionsAll    = frameTurn[16,:].round()
subsets          = {'omr': [[12, 16]], 'rest': [[99]]}
conditions       = seperateConditions(conditionsAll, subsets)
behavior         = frameTurn[14, :] #INDEX into here to get the behavior for that section
coordinates      = loadCoordinates(datasetPath)[0]
samplingRate     = fromfile(datasetPath + 'Stack_frequency.txt', sep = '\n')[0]
#tmp loadfile [change this]
#%%
print(dataName)
voxSize = [5,5,4]
try:
    mp, A   = createSuperVox(coordinates, voxSize, createAdjacency = 1)
    data = loadSuperVox(mp, conditions, loadThis)
    data = array(data['omr'][0])
    nT   = len(data.T)
    data = signal.detrend(data, bp = arange(0,nT,nT//5))
    coordinates = array(list(mp.keys()))
except FileNotFoundError:
    loadConditions(conditions, dataName)
    mp   = createSuperVox(coordinates,voxSize )[0]
    data = loadSuperVox(mp, conditions, loadThis)

#%%


opts = {'steps': 5, 'hyp': {'a0' : 2., 'b0' : 1., 'mu0' : 0., 'kappa0' : 1.} }
opts = {'steps': 20, 'hyp': {'a0' : 2., 'b0' : 1., 'mu0' : 0., 'kappa0' : 1.} }


# %%

import matlab.engine
print('Starting engine')
eng = matlab.engine.start_matlab()
eng.cd('../ddCRP-master/')

#%%
print('Creating datasets')
toMatlab =  matlab.double(data.T.tolist())
# %%
A = [matlab.int32(i) for i in A] # otherwise the transfer crashes
print('Running algorithm')
toMatlab =  matlab.double(data.T.tolist())
# %%
A = [matlab.int32(i) for i in A] # otherwise the transfer crashes

out     = eng.PMC_ddCRP_NG(toMatlab, A, opts)
