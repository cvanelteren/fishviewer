#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 23 13:19:50 2017

@author: casper
"""
from dataAnalysisFunctions import *
#%%

#datasetPath     = '../../ZebrafishRecordings/photo_omr_blob_etc_data/20150120_2_1_photo_OMR_prey_blob_blue_cy74_6d_20150120_220356/Registered/'
# datasetPath      = '../../ZebrafishRecordings/sample_data_code_spontaneous/20140805_1_3_cy74_6d_spon_20140805_170009/'
datasetPath     = '../../ZebrafishRecordings/photo_omr_blob_etc_data/20150122_2_2_cy74_6d_photo_OMR_prey_blob_blue_20150122_190028/Registered/'
#    datasetPath     = '../ZebrafishRecordings/photo_omr_blob_etc_data/20150122_2_2_cy74_6d_photo_OMR_prey_blob_blue_20150122_190028/Registered/'

# datasetPath     = '../../ZebrafishRecordings/test1/test/Registered/'
dataName        = datasetPath.split('/')[-3] # assuming it has registered subfolder
                                      # otherwise edit this
dataPath        = '../Data/'

loadThis        = dataPath + dataName + '.hdf5'

# load conditions
frameTurn        = loadFrameTurn(datasetPath + 'frame_turn.mat')
conditionsAll    = frameTurn[16,:].round()              # condition numbers 
subsets          = {'omr': [[12, 16]], 'rest': [[99]]}  # condition numbers
conditions       = seperateConditions(conditionsAll, subsets)
behavior         = frameTurn[14, :] #INDEX into here to get the behavior for that section
coordinates      = loadCoordinates(datasetPath)[0]

# load sampling rate
try:
    samplingRate     = fromfile(datasetPath + 'Stack_frequency.txt', sep = '\n')[0] 
except FileNotFoundError:
    print('Stack frequency not found, skipping')
    
#%%

addSubsetToDataset(loadThis, conditions)
test = loadSubsetOfData(loadThis, 'omr', 0)
assert 0 
#tmp loadfile [change this]
#%% h5py or create it by loading the .mat
voxSize = [1,1,1]
try:
    mp, A   = createSuperVox(coordinates, voxSize, createAdjacency = False)
    data = loadSuperVox(mp, conditions, loadThis, condition = 'Raw')
    coordinates = array(list(mp.keys()))
except:
    loadConditions(conditions, dataName)
#    mp   = createSuperVox(coordinates, [4,4,2])[0]
    data = loadSuperVox(mp, conditions, loadThis, voxSize, condition = 'Raw')

cond = 'rest'
w    = 0
analyzeThis = signal.detrend(data[cond][w])
#%% binarize left right forward (i.e. bind together info)
tmp = [0,1,3] # left, right, forward
binBehavior = zeros(analyzeThis.shape[-1])
for i,j in enumerate(tmp):
    idx = where(abs(frameTurn[j, conditions[cond][w]]) > 0)[0]
    binBehavior[idx] = i + 1
#%%  deconvolve data
print('Deconvolving data')
deconvolved = zeros((*analyzeThis.shape, 2))
for idx, di in enumerate(tqdm(analyzeThis)):
    di =  deconvolve(float64(di)) #somereason error raises when not float64
    for j in range(2):
        deconvolved[idx, :, j] = di[j]
# %% PCA windows over time (overlap)
comps, de = pcaOverTime(deconvolved[...,0], 50, 10, 20)
comps = comps[:-4]
rr = corrPCAOverTime(comps)
#%% plot duration of component
fig, ax = subplots(9,2, sharex = 'all')
tax = fig.add_subplot(1,1,1)
tax.set_xlabel('Window', labelpad = 20)
tax.set_ylabel('Correlation (r)', labelpad = 20)
fig.suptitle('Window correlation per component\n win length = 50, overlap = 10, threshold = .5')
tax.grid('off')
tax.set_xticks([])
tax.set_yticks([])
tax.patch.set_alpha(0)

ax = ax.flatten()
for i in range(len(ax)):
    ax[i].plot(rr[i,i,:])
    ax[i].set_title('Component {0}'.format(i))
subplots_adjust(hspace=.5)
#%% reconstruct component i
i = 1
rs  = []
for c, d in zip(comps, de.T):
    d = d.T
    r  = (d - c.mean_[i]).dot(c.components_[i,:].T)
    r  = r[:,None].dot(c.components_[[i],:]) + c.mean_[i]
    rs.append(r.T)
rs = array(rs).T;
tmp = {'scatter {0}'.format(i) : {'data': rs[...,i::2].reshape((rs.shape[0],-1)) } for i in range(2)}
cfg = {'coordinates': coordinates, 'data': tmp}
createQt(cfg)

m = rs.mean(-1)
mm = m.mean(0)
fig, ax = subplots()
ax.errorbar(range(len(mm)), mm, yerr = m.std(0))

#%% plot condition data
#i = 2
#j = where(rr[i,i, :] > 0)[0] * 10
#i = range(180, 260)

tcfg = {'scatter':{'data': deconvolved[...,0]}, 'line': {'data': binBehavior}}
cfg = {'coordinates': coordinates, 'data': tcfg}
createQt(cfg)

#%%
cond = 'omr'
w = 0

tmp = signal.detrend(data[cond][w])
i = conditions[cond][w]

j = frameTurn[9, i]
ll = []
for i in tqdm(tmp):
    ll.append(stats.spearmanr(i, j)[0])

fig, ax = subplots();
ax.hist(ll, bins = linspace(-1,1,21))


#%% 
with File(loadThis) as f:
    for i in f :print(i)
#    gr = f['omr']
#    f.__delitem__('rest')
#    f.__delitem__('omr')
#    f.__delattr__('rest')
#    f.__delattr__('omr')
#    for i in f:print(i)
#with File(loadThis) as f:
#    for key, item in data.items():
#        try:
#            f.create_group(key)
#        except:
#            pass
#        for idx, subset in enumerate(item):
#            f[key].create_dataset(str(idx), data = subset)