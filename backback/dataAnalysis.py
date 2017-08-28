
# -*- coding: utf-8 -*-
"""
Created on Thu May 18 16:17:51 2017

@author: caspe
"""
import sys
sys.path.append('../Packages/OASIS/')
import oasis
from zbf import *
from loadData import *
from sklearn import decomposition
from visualize import fastImshow
from scipy import stats
import os
#%%

datasetPath     = '../../ZebrafishRecordings/photo_omr_blob_etc_data/20150120_2_1_photo_OMR_prey_blob_blue_cy74_6d_20150120_220356/Registered/'
# datasetPath      = '../../ZebrafishRecordings/sample_data_code_spontaneous/20140805_1_3_cy74_6d_spon_20140805_170009/'
# datasetPath     = '../../ZebrafishRecordings/photo_omr_blob_etc_data/20150122_2_2_cy74_6d_photo_OMR_prey_blob_blue_20150122_190028/Registered/'
#    datasetPath     = '../ZebrafishRecordings/photo_omr_blob_etc_data/20150122_2_2_cy74_6d_photo_OMR_prey_blob_blue_20150122_190028/Registered/'

# datasetPath     = '../../ZebrafishRecordings/test1/test/Registered/'
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

try:
    samplingRate     = fromfile(datasetPath + 'Stack_frequency.txt', sep = '\n')[0]
except FileNotFoundError:
    print('Stack frequency not found, skipping')

#tmp loadfile [change this]
#%%
print(dataName)
try:
    mp   = createSuperVox(coordinates, [1,1,1])[0]
    data = loadSuperVox(mp, conditions, loadThis)
    coordinates = array(list(mp.keys()))
except:
    loadConditions(conditions, dataName)
#    mp   = createSuperVox(coordinates, [4,4,2])[0]
    data = loadSuperVox(mp, conditions, loadThis)




def detrend(data):
    for key, datasets in data.items():
        print('detrending {0}'.format(key))
        for idx, dataset in enumerate(datasets):
            data[key][idx] = signal.detrend(dataset)

    return data
data = detrend(data)

def createQt(object, cfg):
    '''
    I don't fully understand how to correctly create QT object
    This seems to be a work around
    '''
    import sys
    from PyQt5.QtCore import QCoreApplication
    app = QCoreApplication.instance()
#    app = QCoreApplication(sys.argv)
    if app is None:
        app = QApplication(sys.argv)
        w = object(**cfg)
    else:
        print('QApplication instance already exists: %s' % str(app))
        w = object(**cfg)
    sys.exit(app.exec_())
#%%
omr0 = data['omr'][0]
from sklearn.decomposition import PCA
pca = PCA(10)
fit  = pca.fit_transform(omr0)
r    = pca.inverse_transform(fit).T
cfg = {'coordinates': coordinates, 'data': r}
createQt(ZebraFishViewer, cfg)
fig, ax = subplots(); ax.plot(pca.explained_variance_ratio_)
ax.plot(diff(pca.explained_variance_ratio_))
show()

# from functions import deconvolve
# tmpdata = data['omr'][0]
# d = []
# for idx, i in enumerate(tqdm(tmpdata)):
#     i = float64(i)
#     d.append(deconvolve(i))

#%%

#cfg = {'coordinates': coordinates, 'data': data['omr'][0]}
# cfg = {'coordinates': coordinates, 'data': tt}
# createQt(ZebraFishViewer, cfg)


#w = fastImshow(float32(data['rest'][0][::5,:]))

#%%
# pca = decomposition.PCA(n_components=5)
# tmp  = []
# fig = figure()
# idx =  1
# for key, subsets in data.items():
#
#     for subset in subsets:
#         ax = fig.add_subplot(2,2, idx)
#         nT = len(subset.T)
#         subset = signal.detrend(subset, bp = arange(0, nT, nT//4))
#         fit = pca.fit_transform(subset.T)
#         print('Explained variance')
#         print(pca.explained_variance_ratio_)
#         tmp.append(fit)
#         ax.plot(fit)
#         idx +=1



#%%
def corrWithComponent(data, components):
    '''
    Input:
        components : nd array (n components, n features)
        data       : nd array (cell x  time)
    '''
    r = []
    for di in tqdm(data):
        r.append([stats.spearmanr(di, component)[0] for component in components])
    r = array(r).T
    return r

def pcaOverTime(data, nWin, nOverlap):
    '''
    Returns pca windows over time with overlap
    '''

    from sklearn.decomposition import PCA
    pca         = PCA()
    comps       = []
    remaining   = mod(len(data.T), nWin)
    add         = zeros((data.shape[0], nWin - remaining))
    data        = hstack((data, add)) #zeropad data to prevent error
    windows = range(0, len(data.T), nOverlap)
    for window in tqdm(windows):
        if window+nWin > len(data.T):
            break
        tdata = data[:, window: window + nWin]
        comps.append(pca.fit_transform(tdata.T))
    return array(comps).T

def corrPCAOverTime(comps, threshold = .5):
    '''
    comps = nComp  x nComp x nTime
    returns the correlations of 1 time step back with specified threshold
    '''

    from scipy import stats
    nT = comps.shape[-1]
    r  = zeros(comps.shape)
    for t in range(1, nT):
        c1 = comps[..., t]
        c2 = comps[..., t - 1]
        for idx, i in enumerate(c1):
            for jdx, j in enumerate(c2):
                if idx <= jdx:
                    cor = stats.spearmanr(i, j)[0]
                    if cor > threshold:
                        r[idx, jdx, t - 1] = cor
    return r

comps = pcaOverTime(data['rest'][0], 50, 10)
r     = corrPCAOverTime(comps)
#fig, ax = subplots()
#ax.plot(comps[:,2,:])
#%%
# motor = {'Motor': frameTurn[14, conditions['omr'][0]], 'Conditions': conditionsAll[conditions['omr'][0]]}

# cfg = {'coordinates': coordinates, 'data' : detrend['omr'][0], 'motorAndBehavior' : motor}

# createQt(ZebraFishViewer, cfg)
#%%
def dumpData(file, data):
    '''
    This is a tmp function that can be used for large data dumps
    To prevent it to load it using the slow process of h5py, be sure to remove it
    later or switching datasets
    '''
    print('starting data dump')
    with File(file) as f:
        for key, values in data.items():
            print('Group {0}'.format(key))
            f.create_group(key)
            for subset, value in enumerate(values):
                f[key].create_dataset(str(subset), data = value)
#dumpData('tmp', data)
