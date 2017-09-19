# -*- coding: utf-8 -*-
"""
Created on Tue Aug 15 18:06:46 2017

@author: Casper
"""

# %%
from dataAnalysisFunctions import *
from scipy import ndimage
import sklearn
assert 0 
# %%
#import dask.array as darray
# Setup paths
dataPathRoot    = '../Data/'
rawPathRoot     = '/media/casper/External/ZebrafishRecordings/'
#rawPathRoot     = '../../ZebrafishRecordings/'
#dataName        = '20140805_1_3_cy74_6d_spon_20140805_170009' # name of the folder
dataName        = '20150120_2_1_photo_OMR_prey_blob_blue_cy74_6d_20150120_220356' # part of the omr set
#dataName        = '20140827_1_7_GAD_H2B_6d_spon_20140827_121957' # part of the spontanaeous set
#dataName        = '20150120_2_1_photo_OMR_prey_blob_blue_cy74_6d_20150120_220356'
rawFilePath     = findFolder(dataName, rawPathRoot)           # will find the root
print(rawFilePath)
# load conditions
step = 1; reduce = 0
# conditions are not always present
try:
    frameTurn        = loadFrameTurn(rawFilePath + 'frame_turn.mat')
    # for spontaneous sets the direction is flipped
    if len(frameTurn) > len(frameTurn.T):
        frameTurn = frameTurn.T
    data             = loadData(dataName, rawPathRoot, dataPathRoot) [::step, :]
    coordinates      = loadCoordinates(rawFilePath)[0] [::step, :]
    if reduce :
        voxSize = [8,8,8]
        mp, A = createSuperVox(coordinates, voxSize, createAdjacency = 0)
        c = array(list(mp.keys()))
        coordinates = array(list(mp.keys()))
    
        subdata = zeros((coordinates.shape[0], data.shape[1]))
        for i, idx in enumerate(tqdm(mp.values())):
            subdata[i,:] = data[idx,:].mean(0)
        data = subdata
    samplingRate     = fromfile(rawFilePath + 'Stack_frequency.txt', sep = '\n')[0]
    conditionsAll    = frameTurn[16,:].round()            # condition numbers
    idx              = where(conditionsAll < 100)[0]
    conditionsAll    = conditionsAll[idx]
    conditionsAll    = cleanConditions(conditionsAll)       # remove odd jumps in the data
    subsets          = {'omr': [[12, 16]], 'rest': [[99]], 'moving dot' : [[22, 24]], 'blob': [[29, 33]], 'photo':[[1, 6]]}  # condition numbers
    conditions       = seperateConditions(conditionsAll, subsets) # [::step, :]
    data             = data[:, idx]
    frameTurn        = frameTurn[:, idx]
#    data = data[:, conditions['omr'][0]]
except:
    pass
# %%
tmp = '/media/casper/External/ZebrafishRecordings/sample_data_code_spontaneous'
c = []
d = [] 

for i in os.walk(tmp):
    print(i[0])
    try:
        j = loadCoordinates(i[0] + '/')[0]
        c.append(j)
    except:
        continue
# %%
j = sklearn.neighbors.NearestNeighbors(1).fit(c[0])
p = {i: [] for i in arange(len(c[0]))}
for ci in c:
    _, n = j.kneighbors(ci)
    for idx, ni in zip(ci, n):
        p[ni[0]].append(idx)
pp = []
for i in arange(len(c[0])):
    pp.append(mean(p[i], axis = 0))
pp = array(pp)
# %% meeting wednesday 27 - 7
data  = removeBaseline(data[:, :100])
#zdata = stats.zscore(data, 1)

# %%
voxSize = [12,12,12]
mp, A = createSuperVox(coordinates, voxSize, createAdjacency = 0)
c = array(list(mp.keys()))
print(c.shape)
# %%
subdata = zeros((c.shape[0], data.shape[1]))
for i, idx in enumerate(tqdm(mp.values())):
    subdata[i,:] = data[idx,:].mean(0)
# %%
f = lambda x, a,b,c: a * exp(-b*x) + c
theta = 1e-5
h = []
for ai in tqdm(a):
    try:
        p   = scipy.optimize.curve_fit(f, bins, ai)[0]
        idx = where(diff(f(bins, *p), 2) <= theta)[0][0]
        h.append(bins[idx])
    except:
        h.append(0)
h = array(h)
hh = (h-h.min())/(h.max() - h.min())
# %%
#cfg = {'coordinates' : c, 'data' : {'scatter main': {'data':subdata}}}
cfg = {'data' : {'scatter main': {'data': data, 'coordinates': coordinates, 'cluster' : False}}}
w = createQt(cfg)
# %%
p = []
for i in mp.values():
    points = h[i,:]
    p.append(median(points, axis = 0))
p = array(p)
print(p.shape)
# %%
def singleCor(cellIdx):
    global zdata, bins, coordinates
    coordinate = coordinates[cellIdx, :]
    cellData = zdata[cellIdx,:]
#    use = arange(len(zdata))
    dist =  linalg.norm(coordinates - coordinate, axis = 1) # 2 norm
    idx  = digitize(dist, bins) - 1
    r    = zdata.dot(cellData) / (len(cellData) - 1)
#    r = array([stats.spearmanr(zi, cellData)[0] for zi in zdata])
    
    storage = zeros((len(bins)))
    u = unique(idx)
    if -1 in u:
        print('there is a false god')
    for ui in u:
        j = where(idx == ui)[0] # find the index that belong to this bin
        if len(j) > 0 :
            if cellIdx in j:
                j = delete(j, where(j == cellIdx)[0])
#                print(j, cellIdx)
            storage[ui] = abs(r[j]).mean() if len(j) > 0 else 0
    return storage
        
def totalCor(idx):
    from multiprocess import cpu_count, Pool, Array
    with Pool(cpu_count() - 1 ) as p:
        cors = p.map(singleCor, idx)
            
    return array(cors)
#bins = z
#tmp = totalCor(arange(len(zdata)))
# %%
from joblib import delayed, Parallel
bins = arange(0, 500, 5)
i = 2
from functools import partial
#tmp = array([singleCor(d, zdata, coordinates[:, :], z) for d in tqdm(arange(len(zdata)))])    
#[singleCor(i, zdata, coordinates, z) for i in range(10)]
#names = ['photo', 'moving dot', 'blob', 'omr']
names = []
corrDistance = {}
for key, values in tqdm(conditions.items()):
    if key not in names:
        if key not in corrDistance:
            corrDistance[key] = []
        for value in values:
            zdata = stats.zscore(data[:, value], axis = 1)
            corrDistance[key].append(array(Parallel(n_jobs  = 3)(delayed(singleCor)(i) for i in arange(len(zdata)))))
#            corrDistance[key].append(array([singleCor(d) for d in tqdm(arange(len(zdata)))]))

# %%
tmp = []
for key, values in corrDistance.items():
    print(key)
    [tmp.append(i[:,2]) for i in values]
tmp = array(tmp).T

# %%
cfg = {'coordinates': coordinates, 'data': {'scatter' : {'data': tmp}}}
createQt(cfg)