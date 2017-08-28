#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 23 13:19:50 2017

@author: casper
"""

# %%
from dataAnalysisFunctions import *
from scipy import ndimage
assert 0 
# %%
#import dask.array as darray
# Setup paths
dataPathRoot    = '../Data/'
rawPathRoot     = '/media/casper/External/ZebrafishRecordings/'
rawPathRoot     = '../../ZebrafishRecordings/'
#dataName        = '20140805_1_3_cy74_6d_spon_20140805_170009' # name of the folder
dataName        = '20150122_2_2_cy74_6d_photo_OMR_prey_blob_blue_20150122_190028' # part of the omr set
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
        voxSize = [20,20,10]
        mp, A = createSuperVox(coordinates, voxSize, createAdjacency = 0)
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
    subsets          = {'omr': [[12, 16]], 'rest': [[99]]}  # condition numbers
    conditions       = seperateConditions(conditionsAll, subsets) # [::step, :]
    data             = data[:, idx]
    frameTurn        = frameTurn[:, idx]
#    data = data[:, conditions['omr'][0]]
except:
    pass


#i = where(logical_and(conditionsAll >= 12, conditionsAll <= 16))[0]
#data = data[:, i]
sparsity = lambda x : (1 - nanmean(x, 1)**2 / nanvar(x,1))
#data = data[:,1000:]
# %% meeting wednesday 27 - 7
data  = removeBaseline(data)
zdata = stats.zscore(data, 1)
# %%
def singleCor(cellData, coordinate,  zdata, coordinates, bins):
#    cellData, coordinate = x
    dist =  linalg.norm(coordinates - coordinate, axis = 1) # 2 norm
    idx  = digitize(dist, bins) - 1

    r    = zdata.dot(cellData) / (len(cellData) - 1)
    storage = zeros((len(bins)))
    u = unique(idx)
    if -1 in u:
        print('there is a false god')
    for ui in u:
        j = where(idx == ui)[0] # find the index that belong to this bin
        if len(j) > 0 :
            storage[ui] = abs(r[j]).mean()
    return storage
        
def totalCor(zdata, coordinates, bins = arange(0, 500, 5)):
    from multiprocess import cpu_count, Pool
    from functools import partial
    func = partial(singleCor, zdata = zdata, coordinates = coordinates, bins = bins)
    with Pool(cpu_count() - 1 ) as p:
        cors = p.map(func, zip(zdata, coordinates))
            
    return array(cors)
z = arange(0, 500, 5)
i = 2
from functools import partial
tmp = array([singleCor(d, c, zh, coordinates[::2, :], z) for d,c in tqdm(zip(zh[:,:], coordinates[::2,:]))])      

# %%

b = frameTurn[0, 1000:]
b = ndimage.filters.gaussian_filter(b,3)
# %%
ica = decomposition.FastICA(5)
ica.fit(h)
r = zdata[i,:].dot(stats.zscore(ica.components_).T) / (len(zdata.T) - 1)
rr = argmax(r, 1)
x = linspace(0, 1, len(ica.components_))

fig, ax = subplots(); ax.imshow(ica.components_, aspect = 'auto')
# %%
rr = array([argmax(ri) for ri in r])
# %%
i = where(tmp[:,1] > .25)[0]

r = zdata[i,:].dot(stats.zscore(b)) / (len(b) - 1)
# %%
#j = where(rr == 29)[0]
fig, ax = plot3(); g = ax.scatter(*coordinates[i,:].T, c = r, cmap = 'RdBu'); cb =  colorbar(g, ax = ax)
cb.set_label('Pearson correlation (r)')
ax.set_title('Correlation with behavior vector')
setp(ax, 'xlabel', 'x', 'ylabel', 'y', 'zlabel', 'z')
#%%
h = zdata[i,:]
pca = decomposition.PCA(10)
pca.fit(h)
r = h.dot(stats.zscore(pca.components_.T)) / (len(h.T) - 1)
rr = argmax(r, 1)
x = linspace(0, 1, len(ica.components_))

fig, ax = subplots(); ax.imshow(ica.components_, aspect = 'auto')
# %%
jj = []
import statsmodels.api as sm
gc =  sm.tsa.stattools.grangercausalitytests
for i, x in enumerate(h):
    for j, y in enumerate(h):
        if i < j:
            res = gc(vstack((x,y)).T, 1, verbose = False)
            if res[1][0]['lrtest'][1] < .05:
                jj.append((i,j))
                
# %%
idx = []
for i in tqdm(a):
    ii = gradient(gradient(i))
    b = [allclose(j, 0, 0,.001) for j in ii]
    idx.append(i[where(b)[0][0]])
idx = array(idx)
# %%
#cfg = {'coordinates': c[i,:], 'data': {'scatter' : {'data' : k[i,:], 'xlabel' : 'time[step]'}}}
#cfg = {'coordinates': coordinates[:,:], 'data': {'scatter' : {'data' : stats.zscore(idx[:,None], axis = 0)}}}
#cfg = {'coordinates': coordinates[i,:], 'data': {'scatter' : {'data' :x[rr[:,None]]}}, 'cmap' : 'hsl'}
#cfg = {'coordinates': coordinates[i,:], 'data': {'scatter' : {'data' : r}}}
pcfg = {'coordinates': coordinates[:,:], 'data' : {'scatter main' : {'data' : h, 'xlabel' :'Time[step]','cluster': 0}}}
#pcfg = {'coordinates': coordinates[i[j],:], 'data' : {'scatter main' : {'data' : data[i[j],:],'cluster': False, 'xlabel' :'Time[step]'},
#        'line behavior': {'data': b}}, 'cmap' :'viridis'}
w =createQt(pcfg)
# %%
y = zeros((len(data[i,:]), len(unique(ap.labels_))))
for u, yi in zip(unique(ap.labels_), y.T):
    j = where(ap.labels_ == u)[0]
    yi[j] = 1
    
# %%
fig, ax = subplots();
#ax.plot(data.mean(0), label  = 'raw')
ax.plot(detrendData.mean(0), label = 'baseline remove')
ax.plot(decData.mean(0)[...,0], label = 'friedrich')
ax.plot(decData.mean(0)[...,1])
#ax.plot(stats.zscore(tmp,1).mean(0), label = 'linear detrend')
ax.legend()
ax.set_xlabel('Time[step]')
savefig('../Figures/comparison')

#fig, ax = subplots();
#ax.plot(ddata.mean(0), label = 'raw deconvolved')
#ax.plot(drdata.mean(0), label = 'baseline janelia deconvolved')
# %%
f, px = signal.periodogram(decData[...,0].mean(0), samplingRate)
ff, ppx = signal.periodogram(detrendData.mean(0), samplingRate)
fig, ax = subplots(2,1, sharex = 'all')
ax[0].plot(f, px)
ax[1].plot(ff, ppx)
ax[0].set_title('Friedrich FR')
ax[1].set_title('Detrended FR')
xlabel('Frequency [Hz]')
[axi.set_ylabel('Power(A**2)') for axi in ax]
#%%
z = arange(0, 500, 5)
zz = arange(0, len(z))
h = []
for idx, c in enumerate(tqdm(coordinates)):
    d = linalg.norm(coordinates - c, axis = 1)
    idx = digitize(d, z)
    for zi in zz:
        h.append(where(idx == zi)[0])
    
# %%
for i in [detrendData, decData[...,0], decData[...,1]]:
    print(sparsity(i))
# %%
#i = where(logical_and(sp > .5, sp < .7))[0]

tmp = data.mean(1)
i = where(tmp[:,2] > 0.2)[0]; print(i.shape)

h = data[i,:]

# %%
w =createQt({'coordinates': coordinates[:, :], 'data':{'scatter':{'data':data}}})
# %%
rr = array([])
for idx, values in enumerate(mp.values()):
    rr = hstack((rr,r[values, idx]))
rr = array(rr)
# %%

# %%
def createData(fileName, names):
    with File(fileName) as f:
        for name in names:
            f.create_dataset(name, data = globals()[name])
names = 'coordinates samplingRate data frameTurn z tmp'.split(' ')
createData(dataName + '.hdf5', names)

# %%
with File('20150122_2_2_cy74_6d_photo_OMR_prey_blob_blue_20150122_190028.h5') as f:
    tmp = []
    gr = f['Conditions']
    for i in gr:
        tmp.append(gr[i].value)
print(tmp[:10])
        
# %%
with File(g) as f:
    for i in f:
        print(f[i].value.shape)

# %%
from scipy import ndimage
#zscore = lambda x: (x - x.mean(-1)[:,None]) / x.std(-1)[:,None]

#tdata  = array([ndimage.filters.gaussian_filter(float32(i), 1) for i in tqdm(data)])

#tb = frameTurn[0, conditions[c][cc]]
tb = frameTurn[0, :]
#tb = frameTurn[:,0]
tb = ndimage.filters.gaussian_filter(tb, 3)
tb = stats.zscore(tb)
r  = array([stats.spearmanr(x,tb)[0] for x in tqdm(data)])

sortR = sort(r)

#i = .005
i = 0.005
theta1 = sortR[int(1-i * len(r))]
theta2 = sortR[int(i * len(r))]
idx = where(r >= theta1)[0]
jdx = where(r <= theta2)[0]
kdx = hstack((idx, jdx))
#plot(tb)
#w = createQt({'coordinates': coordinates[i, :], 'data':{'scatter':{'data': data[i, :]}}})

# %%
        
def calcA(x, coordinates, maxDist, maxNum = 200):
    '''
    Should be used with multi A
    '''
    i, x = x
    dist = linalg.norm((coordinates - x), axis = 1)
    idx = where(logical_and(dist <= maxDist, dist > 0))[0]
    if len(idx) > maxNum:
        jdx = argsort(dist[idx])
        idx = idx[jdx[:maxNum]]
#        d = 1/dist[idx]
#        d = d / sum(d)
#        print(sum(d))
#        i = random.choice(range(len(idx)), maxNum, False, d)
#        idx = idx[i]
#        print(len(i))
#        hist(dist[idx[i]])
#        assert 0
    return (idx, dist[idx])

def multiA(fileName, coordinates, maxDist, maxNum):
    '''
    Multi threading adjacency computation, @threshold indicates the max distance to use
    returns adjacency @A which consists of the neighbour index and the distance
    '''
    from multiprocessing import cpu_count, Pool
    pool = Pool(cpu_count())
    from functools import partial
    func = partial(calcA, coordinates = coordinates, maxDist = maxDist, maxNum = maxNum) # allow partial parse
    # parallel write to file, but for every return
    A = pool.map(func, enumerate(coordinates))
#    pool.join()
    pool.close()
    return array(A)
#    with File(fileName, 'w') as f:
#        gr = f.create_group('A')
#        res = tqdm(pool.imap(func, enumerate(coordinates)))
##        print()
#        for n, d in res:
#            gr.create_dataset(str(n), data = d)

maxDist = 400
maxNeighbor = 500
#A = multiA(dataName, coordinates, maxDist)
A = array([calcA((i,j), coordinates, maxDist, maxNeighbor) for i,j in enumerate(tqdm(coordinates))])
#createQt({'coordinates': coordinates[a, :], 'data':{'scatter':{'data': b[:,None] }}})

# %% 
import statsmodels.api as sm
gc = sm.tsa.stattools.grangercausalitytests
tmp = []
for i, xi in enumerate(tqdm(yy)):
    for j, xj in enumerate(yy):
        if i < j:
            test = gc(vstack((xi,xj)).T,1, verbose = False)
            if test[1][0]['lrtest'][1] < .05:
                tmp.append((i,j))
# %%
def singleCellStat(celli, A, data):
    '''
    Computes the correlation between @celli and neighbours in A
    returns correlation of @celli and neighbours @A[celli]
    '''
    localCor = []
    try:
        ai = A[celli, ...] # idx first column, dist 2nd
        ai = array(ai, dtype = int)
        x = data[[celli], :]
        if len(ai) == 1:
            y = data[[ai], :]
        else:
            y = data[ai, :]
        r = x.dot(y.T) / (len(x.flatten()) - 1)
       
        
#        xi = x - x.mean(axis = -1)[:, None]
#        yi = y - y.mean(axis = -1)[:, None]
#        xi_std = x.std(axis = -1)[:, None]
#        yi_std = y.std(axis = -1)[:, None]
#        n = x.shape[1]
#        tmp = xi.dot(yi.T)
#        tmp2 = xi_std * yi_std * (n - 1)
#        r = tmp / tmp2.T
        localCor = r
        return localCor.T
    except:
        return localCor

def localCorrelation(A, x):
    '''
    Produces local correlation based on Adjacency
    A : adjacecny matrix containing tuple of @distance, and @idx of cells
    data : data, nCells x nTime
    '''
    #TODO: this does not seem to recruit all cores as @multiA

    import multiprocessing as mp
    from functools import partial
    func = partial(singleCellStat,  A = A, data = x)
    pool = mp.Pool(mp.cpu_count() - 1) #TODO: will crash for 1 systems
    lc = pool.map(func, range(len(A))) # should return ordered list
    pool.close()
    return array(lc)
#tmp = singleCellStat(0, A, data[::2,:])
#lc = localCorrelation(A[:,0,:], stats.zscore(detrendData, 1))

lc = array([singleCellStat(i, A[:,0, :], zdata) for i in tqdm(range(len(A)))]).squeeze()
#lc = array([singleCellStat(i, A[:,0, :], zdata) for i in tqdm(kdx)]).squeeze()
dist = A[:,1,:]


# %%

# %%
from multiprocess import Pool
from functools import partial

def tmp(i, x):
    linalg.norm(x - i)
def tmp2(x):
    func = partial(tmp, x = x)
    with Pool(3) as p:
        p.imap(func, x)


# %%
z = arange(0, 500, 10)
d = []

for i in tmp:
    i = array(i).squeeze()
    j = digitize(i[:,0], z)
    dd = zeros((len(z)+1))
    for u in unique(j):
        idx = where(j == u)[0]
        if len(idx) > 0:
            dd[u] = i[idx,1].mean()

    d.append(dd)
d = array(d)
print(d.shape)
        
# %%
c = []
h = []
for ii, i in enumerate(d):
    l = i.shape[0]
    if l not in c:
        c.append(l)
        h.append(iih)

# %%
def widthOfCorrelation(localCor, distances):
    out = zeros((len(localCor),1))
    for idx, (lc, disti) in enumerate(tqdm(zip(localCor, distances))):
        tmp = linspace(0, max(disti), 10)
        j   = digitize(disti, tmp)
        r  = []
        for ji in unique(j):
            i = where(j == ji)[0]
            r.append(mean(abs(lc[i])))
#        plot(gradient(r))
        tmp = gradient(r)
        i = argmax(tmp)
#        print(i)
        out[idx,:] = tmp[i]
    return out

idx = []
for ii, i in enumerate(lc):
    tmp = argmin(i)
    idx.append(tmp)
idx = array(idx)
d = array([d[i] for d, i in zip(dist, idx)])
#idx = array([where(abs(i) < .2)[0][0] for i in lc])
#tmp = widthOfCorrelation(lc, dist)
# %%
w = createQt({'coordinates': coordinates[:, :], 'data':{'scatter':{'data': data}}})
#w = createQt({'coordinates': coordinates[:, :], 'data':{'scatter':{'data': lc}}})
# %%
z = arange(0,maxDist, 10)
binnedDistance = digitize(A[:,1,:].flatten(), z)
bins = unique(binnedDistance)

localCors = lc.flatten()
r = []
s = []
for bi in bins:
    i = where(binnedDistance == bi)[0]
    r.append(mean(localCors[i]))
    s.append(localCors[i].std()*2)
    

#%%
fig, ax = subplots();
#ax.errorbar(bins, r, s)
ax.plot(bins, r)
setp(ax, 'xlabel','distance $\mu$m', 'ylabel', 'correlation (r)')
# %%
with File(dataName) as f:
    for i in f['A']:
        v = f['A/' + i].value
        if v != []:
            print(v)
# %%
lc = localCorrelation(A, signal.detrend(data))
#%%
def pcaRecon(data, comps, idx):
    mu = data.mean(0)
    r = (data - mu).dot(comps[idx,:].T)
    r = r.dot(comps[idx,:]) + mu
    return r
pca = decomposition.PCA(3)
pca.fit(tmp)

#%%
pca=decomposition.PCA(4)
pca.fit(tmp)
c = conditionsAll[conditions['omr'][0]]
c = (c-c.min()) / (c.max() - c.min())
#cmap= cm.viridis(linspace(0,1, len(c)))
cmap = cm.viridis(c)
#%matplotlib notebook
from mpl_toolkits.mplot3d import Axes3D
fig, ax = subplots(subplot_kw = {'projection':'3d'}) ; ax.scatter(*pca.components_[1:,:], color = cmap)
show()
#%%
r = pcaRecon(tmp, range(2,5))
# %% KDE plot of distance and local correlation
h = []
for j in tqdm(range(100)):
    num = random.randint(0, len(lc))
    i = lc[num]
    if i != [] and num not in h :
        if j == 0:
            tmp = i.T
        tmp = hstack((tmp, i.T))
    h.append(num)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
#tmp = vstack((tmp, tmp2))
#fig, ax = subplots(); ax.hist(z,z,weights = tmp)
#cfg ={'title' : 'Local correlation, max distance = {0}'.format(maxDist),
# 'xlabel' : 'Correlation  (r)' ,'ylabel' :'count'}
#setp(ax, **cfg)
#
#savefig('../Figures/Local Correlation')
#
#fig, ax = subplots()
#h = ax.imshow(s, extent = [zz[0], zz[-1], z[0], z[-1]], interpolation = None)
#cfg ={'title' : 'Local correlation, max distance = {0}'.format(maxDist),
# 'xlabel' : 'distance' ,'ylabel' :'r'}
#setp(ax, **cfg)
#savefig('../Figures/Local Correlation 2d')
##colorbar(h)
## %%
#import seaborn as sb
#sb.jointplot(tmp[1,:], tmp[0,:], kind = 'kde').set_axis_labels('Distance', 'Correlation (r)')
##xlabel('test')
##cfg ={'title' : 'Local correlation, max distance = {0}'.format(maxDist),
## 'xlabel' : 'distance' ,'ylabel' :'r'}
##setp(ax, **cfg)
#savefig('../Figures/Local Correlation 2d')
## %%
#for i in range(1): print(i)
#from scipy.io import savemat
#cfg ={'coordinates': coordinates, 'data' :data,  'A': A}
#savemat(dataName, cfg)
##createQt(cfg)
#
##%%
#createQt({'coordinates': coordinates, 'data':{'scatter':{'data': data}}})
##tmp = loadmat('../Data/' + dataName + '.mat')
## %%
#
def save2file(fileName, data, dataName):
    with File(fileName) as f :
        for cellIdx, dataRow in enumerate(tqdm(data)):
            f.create_dataset(dataName + '/{0}'.format(cellIdx), data = dataRow)
##save2file(dataName + '.h5', data, 'Raw')
#save2file(dataName + '.h5', A, 'Adjacency')
#save2file(dataName + '.h5', lc, 'Local correlation')
#save2file(dataName + '.h5', dfiltDat, 'Deconvolved')
#save2file(dataName + '.h5', A, 'Adjacency')
# addSubsetToDataset(loadThis, conditions)
# test = loadSubsetOfData(loadThis, 'omr', 0)
#tmp loadfile [change this]
# %% h5py or create it by loading the .mat
# try:
#
#try:
#    mp, A   = createSuperVox(coordinates, voxSize, createAdjacency = False)
#    cond = 'omr'
#    w    = 0
#    analyzeThis = loadSuperVox(mp, loadThis, cond + '/' + str(w))
##    analyzeThis = signal.detrend(analyzeThis)
#except:
#    addSubsetToDataset(loadThis, conditions)
#    analyzeThis = loadSuperVox(mp, loadThis, cond + '/' + str(w))

# addSubsetToDataset(loadThis, conditions)
# test = loadSubsetOfData(loadThis, 'omr', 0)
#tmp loadfile [change this]
# %% h5py or create it by loading the .mat
# try:
#
#try:
#    mp, A   = createSuperVox(coordinates, voxSize, createAdjacency = False)
#    cond = 'omr'
#    w    = 0
#    analyzeThis = loadSuperVox(mp, loadThis, cond + '/' + str(w))
##    analyzeThis = signal.detrend(analyzeThis)
#except:
#    addSubsetToDataset(loadThis, conditions)

# cond = 'rest'
# w    = 0
# analyzeThis = signal.detrend(data[cond][w])
# #%%
# left = frameTurn[9, conditions[cond][w]]
# corrs = []
# for cell in tqdm(stats.zscore(analyzeThis,1)):
#     corrs.append(stats.spearmanr(cell, left)[0])
# corrs = array(corrs)
# #%%
# cfg = {'scatter': {'data': analyzeThis}}
# cfg = {'coordinates' : coordinates, 'data' : cfg}
# createQt(cfg)
#
# # %% binarize left right forward (i.e. bind together info)
# tmp = [0,1,3] # left, right, forward
#
# binBehavior = zeros(analyzeThis.shape[-1])
# for i,j in enumerate(tmp):
#     idx = where(abs(frameTurn[j, conditions[cond][w]]) > 0)[0]
#     binBehavior[idx] = i + 1
# # %%  deconvolve data
# #print('Deconvolving data')
# #deconvolved = zeros((*analyzeThis.shape, 2))
# #for idx, di in enumerate(tqdm(analyzeThis)):
# #    di =  deconvolve(float64(di)) #somereason error raises when not float64
# #    for j in range(2):
# #        deconvolved[idx, :, j] = di[j]
#
# # %%

# def mp_deconvolved(x):
#     '''
#     Returns spike estimates and deconvoled data '''
#
#     return array(deconvolve(float64(x))[:2]).T
#
# def mp_run(data):
#     '''
#     multithreaded approach to deconvolution
#     '''
#     import multiprocessing as mp
#
#     pool = mp.Pool(mp.cpu_count()-1) #TODO: will crash for 1 systems
#     res = pool.map(mp_deconvolved, data) # should return ordered list
#     pool.close()
#     pool.join()

#     return array(res)
# print('Deconvolving data')
# deconvolved = mp_run(analyzeThis)
# assert 0
# # %% PCA windows over time (overlap)
# nWin = 50; nStep = 10; nComp = 20
# comps, de = pcaOverTime(deconvolved[...,0], nWin, nStep, nComp)
# comps = comps[:-4]
# rr = corrPCAOverTime(comps)
#
# #%%
# print(deconvoled.shape)
# # %% plot duration of component
# fig, ax = subplots(9,2, sharex = 'all')
# tax = fig.add_subplot(1,1,1)
# tax.set_xlabel('Window', labelpad = 20)
# tax.set_ylabel('Correlation (r)', labelpad = 20)
# fig.suptitle('Window correlation per component\n win length = {0}, overlap = {1}, threshold = {2}'.format(nWin, nWin - nStep, .5))
# tax.grid('off')
# tax.set_xticks([])
# tax.set_yticks([])
# tax.patch.set_alpha(0)
#
# ax = ax.flatten()
# for i in range(len(ax)):
#     ax[i].plot(rr[i,i,:])
#     ax[i].set_title('Component {0}'.format(i))
# subplots_adjust(hspace=.5)
# show()
#
# # %% reconstruct component i
#
# i = 1
# rs  = []
# for c, d in zip(comps, de.T):
#     d = d.T
#     r  = (d - c.mean_[i]).dot(c.components_[i,:].T)
#     r  = r[:,None].dot(c.components_[[i],:]) + c.mean_[i]
#     rs.append(r.T)
#
# r = 0
# d = 0
# c = 0
#
# rs = array(rs).T;
# tmp = {'scatter {0}'.format(i) : {'data': rs[...,i::2].reshape((rs.shape[0],-1)) } for i in range(2)}
# cfg = {'coordinates': coordinates, 'data': tmp}
# createQt(cfg)
#
# m = rs.mean(-1)
# mm = m.mean(0)
# fig, ax = subplots()
# ax.errorbar(range(len(mm)), mm, yerr = m.std(0))
#
# # %% plot condition data
# #i = 2
# #j = where(rr[i,i, :] > 0)[0] * 10
# #i = range(180, 260)
#
# tcfg = {'scatter':{'data': deconvolved[...,0]}, 'line': {'data': binBehavior}}
# cfg = {'coordinates': coordinates, 'data': tcfg}
# createQt(cfg)
#
# # %%
# cond = 'omr'
# w = 0
#
# tmp = signal.detrend(data[cond][w])
# i = conditions[cond][w]
#
# j = frameTurn[9, i]
# ll = []
# for i in tqdm(tmp):
#     ll.append(stats.spearmanr(i, j)[0])
#
# fig, ax = subplots();
# ax.hist(ll, bins = linspace(-1,1,21))
#
#
# # %%
# with File(loadThis) as f:
#     for i in f :print(i)
# #    gr = f['omr']
# #    f.__delitem__('rest')
# #    f.__delitem__('omr')
# #    f.__delattr__('rest')
# #    f.__delattr__('omr')
# #    for i in f:print(i)
# #with File(loadThis) as f:
# #    for key, item in data.items():
# #        try:
# #            f.create_group(key)
# #        except:
#            pass
#        for idx, subset in enumerate(item):
#            f[key].create_dataset(str(idx), data = subset)
#memn