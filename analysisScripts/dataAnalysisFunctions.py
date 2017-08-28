
# -*- coding: utf-8 -*-
"""
Created on Thu May 18 16:17:51 2017

@author: caspe
"""
import sys
sys.path.append('../Packages/OASIS/')

from zbf import *
from loadData import *
from sklearn import decomposition
from visualize import fastImshow
from scipy import stats
import os
from functions import deconvolve

def removeBaseline(data, windowSize = 300, prct = 15):
    nT  = len(data.T)
    out = zeros(data.shape)
    for i in tqdm(range(0, nT, windowSize)):
        if i - windowSize < 0:
            start = 0
        else:
            start = i - windowSize
        if i + windowSize > nT:
            stop = nT
        else:
            stop = i + windowSize
        tmp                = data[:, start:stop]
        p                  = percentile(tmp, prct, 1)
        tmp                = p[:, None] * ones((tmp.shape)) # element wise product
        out[:, start:stop] = tmp 
    out = data - out
#    out = stats.zscore(data - out, 1)
    return out


    
def mp_deconvolved(x):
    '''
    Returns spike estimates and deconvoled data '''
    return array(deconvolve(x)[:2]).T
def mp_run(x):
    import multiprocess as mp
    pool = mp.Pool(mp.cpu_count()-1) #TODO: will crash for 1 systems
    x = pool.map(mp_deconvolved, float64(x)) # should return ordered list
    pool.close()
    return array(x, dtype = float32)
# %%
    
# %%
    
def plot3():
    fig, ax = subplots(subplot_kw = {'projection':'3d'})
    return fig, ax
#%% 
def detrend(data):
    for key, datasets in data.items():
        print('detrending {0}'.format(key))
        for idx, dataset in enumerate(datasets):
#            dataset = signal.detrend(dataset, bp = arange(0, len(dataset.T), len(dataset.T)//10))
            dataset = signal.detrend(dataset)
            data[key][idx] = dataset
#            data[key][idx] = stats.zscore(dats)
    return data

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

def pcaOverTime(data, nWin, nOverlap, nComps = 5):
    '''
    Returns pca windows over time with overlap
    '''

    from sklearn.decomposition import PCA
    pca         = PCA(nComps)

    remaining   = mod(len(data.T), nWin)
    add         = zeros((data.shape[0], nWin - remaining))
    data        = hstack((data, add)) #zeropad data to prevent error
    windows = range(0, len(data.T), nOverlap)
    comps = zeros((len(windows)-1), dtype = object)
    d   = zeros((data.shape[0], nWin, len(windows)))
    for idx, window in enumerate(tqdm(windows)):
        if window+nWin > len(data.T):
            break
        pca         = PCA(nComps)
        tdata = data[:, window: window + nWin]
        pca.fit(tdata)
        comps[...,idx] = pca
        d[...,idx]  = tdata

    return comps, d
def cleanConditions(conditions):
    ''' Removes odd numbers in conditions numbers
    '''
    for i in range(2, len(conditions)):
        middle = conditions[i-1]
        start = conditions[i-2]
        stop  = conditions[i]
#        print(start, middle, stop)
        # if center has odd value change it to neighbor
        if middle != stop and middle != start:
            conditions[i-1] = start
    return conditions
        
def corrPCAOverTime(comps, threshold = .5):
    '''
    comps = nComp  x nComp x nTime
    returns the correlations of 1 time step back with specified threshold
    '''
    print(comps)
    comps = array([i.components_.T for i in comps]).T
    print(comps.shape)
    from scipy import stats
    nC, nW, nT = comps.shape
    r  = zeros((nC, nC, nT))
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

#%%
def corrMap(data, theta = .5):
    rmap = {}
    for cellIdx, celli in enumerate(data):
        print('cell {0}'.format(cellIdx))
        for cellJdx, cellj in enumerate(tqdm(data)):
            if cellIdx < cellJdx:
                if cellIdx not in rmap.keys():
                    rmap[cellIdx] = []
                else:
                    r = stats.spearmanr(celli, cellj)[0]
                    if r is not nan and abs(r) > theta:
                        rmap[cellIdx].append((cellJdx, r))
    return rmap




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
