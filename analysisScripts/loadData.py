

# -*- coding: utf-8 -*-
"""
Created on Fri May 12 17:24:32 2017

@author: Casper van Elteren
"""

#from scipy import loadmat
import matplotlib
matplotlib.use('Qt5Agg')

#from matplotlib.pyplot import *
from numpy import *
from h5py import  File
from tqdm import tqdm

from matplotlib.pyplot import *
from mpl_toolkits.mplot3d import Axes3D
from scipy.io import loadmat

from scipy import stats, signal

# %%
def loadCoordinates(datasetPath, xycorrection = .409, zcorrection = 4):
    '''
    Loads the coordinates from the mat file in a directory indicated in path
    The xy, coordiantes are divided by 4 as this is assumed by Ahrens [2012]; there is a measurement uncertainty of 4-8 micrometer, 4 is assumed

    '''

    ci = loadmat( datasetPath + 'cell_info.mat')['cell_info']           # load the cell_info struct; numpy puts the fields in dtype
#    print(ci.dtype)                                                     # print the fields of the struct
    N = len(ci.T)

    # locate memory
    # TO DO: make this dynamic by having a struct, to reduce loading large redundant arrays
    location = zeros( (N, 3) )
    area     = zeros( (N, 1) )
    inds     = zeros( (N, 1) )
    print('Obtaining coordinates from cells...')
    for cellIdx, cellInfo in enumerate(tqdm(ci.T)):
        x, y, z = hstack((cellInfo['center'][0][0].T, cellInfo['slice'][0].flatten())) # struct issue loading causes this uggly indexing
        location[cellIdx, :] = [x * xycorrection, y * xycorrection, z * zcorrection]
        area[cellIdx, :]     = cellInfo['area'][0][0][0]
        inds                 = cellInfo['inds']
    return location, area, inds

def findFolder(fileName, dataFolder):
    '''
    Walks the root filde structure and returns the project
    subdirectory (registed is appended)
    '''
    import os
    # list all subdirectories look for the file name
    for directory in os.walk(dataFolder):
        # first is always the full filePath
        if fileName in directory[0]:
            if 'Registered' in os.listdir(directory[0]):
                return directory[0] + '/' + 'Registered/'
            else:
                return directory[0] + '/'


def loadData(fileName,
            rawDataFolder = '../../ZebrafishRecordings/',
            dataStoreFolder = '../Data/', subset = 'Raw'):
    '''
    This function checks whether a h5 data structure is present in the
    @dataStoreFolder, otherwise it searches for a dataset in @rawDataFolder
    '''
    from scipy.io import loadmat
    from h5py import File
    from numpy import array, fromfile
    h5Name = fileName + '.h5'
    filePath = findFolder(fileName, rawDataFolder)
    if h5Name in dataStoreFolder:
        # check if the @subset exists
        if subset in f:
            print('h5 found loading data')
            with File(dataStoreFolder + h5Name) as f:
                data = array([i.value for i in f[subset]])
            return data
        else:
            print('Subset not found, check file contents')
    # load the binaries
    else:
        print('h5 not found loading binaries')
        # TODO: the subfolders show different types of stackf check what they are
        print(filePath)
        cell_resp_dim = loadmat(filePath + 'cell_resp_dim.mat')['cell_resp_dim'][0]
        data = fromfile(filePath + 'cell_resp_lowcut.stackf', dtype = float32).reshape(cell_resp_dim[::-1]) # data is float32
        return data.T # transposed form




# %%
def createAdjacency(pos, threshold = 5):
    """
    Returns adjacency matrix based on @threshold
    """
    dims = pos.max(axis = 0)
    A = []
    for cellIdx, coord in enumerate(tqdm(pos)):
        dist          = linalg.norm((coord - pos), axis = 1)
        idx            = where(logical_and(dist <= threshold, dist > 0))[0]
        A.append((idx + 1, dist[idx]))            # init list
    return A
# def createAdjacency(pos, file, threshold = 10):
#     '''
#     Pdist memery bullshitery
#     This function will create the adjacency matrix based on a distance threshold
#     It is intended to be used with the entire grid, and will circumvent large memory capacity used in pdist
#     by using loops (slower)
#     '''
#     import tqdm
#     grid = array([unique(i) for i in pos.T])
#            A[tuple(coord)] = la[0]            # init list
#     cellIdx = 0
#     with File(loadThis) as f:
#         f.__delitem__('Adjacency')
#         f.create_group('Adjacency')
#         d = zeros(pos.shape[0])
#         for cellIdx, coord in tqdm(enumerate(pos)):
#             dist          = linalg.norm((coord - pos), axis = 1)
#             la            = where(logical_and(dist <= threshold, dist > 0))
# #                    print(la);
#             d[cellIdx]    = len(la)
#
# #            A[tuple(coord)] = la[0]            # init list
#
#             gr.create_dataset(str(cellIdx), data = la)
#         cellIdx += 1
#
#
#     return d
# %%
def createSuperVox(coordinates,
                voxSize = [1, 1, 1],
                plot = False,
                createAdjacency = False,
                threshold = 1):
    '''
    This function creates super voxels given the desired voxel dimensions and coordinate structure
    Returns:
        Adjacency dictionary per voxel
        Dictionary mapping of voxel to original coordinates
    '''

    print('Creating super voxels')
    print('Voxel size {0}'.format(voxSize))
    newCoords   = coordinates - coordinates.min(axis = 0)   # make it start at zero index
    spacing     = []                                        # contains the bins per dimension
    digits      = []                                        # this will contain the bin number per cell
    print('Making grid')
    for dim, c in zip(voxSize, newCoords.T):
#        tmp = arange(0, c.max() + dim, dim)   # create new coordinate axis
        tmp = linspace(c.min(), c.max(), c.max() / dim)
        spacing.append(tmp)
        digits.append(digitize(c, tmp, right = True))  # obtain the binned indices

    digits  = array(digits).T
    spacing = array(spacing).T
    print('Dimensions')
    [print(i.min(), i.max()) for i in digits.T]
    # create dictionary mapping from super voxel coordinate to
    # cell number
    MAP     = {}
    count   = 0

    # creating dictionary map
    tmp = [tuple(row) for row in digits] # convert in to tuples
    tmp = array(list(set(tmp)))         # convert into array
    # add idx of coordinate to the tuple key
    for c in tmp:
        idx = where(all(digits == c,1))[0]
        MAP[tuple(c)] = idx
    A = []
    if createAdjacency:
        tmpMAP = array(list(MAP.keys()))
        print('Creating adjacecency matrix')
        for idx, coordinate in enumerate(tqdm(tmpMAP)):
            dist = linalg.norm(coordinate - tmpMAP, axis = 1)
            dist[idx] = 0 #nan issues in 1/dist (#TODO: not doing this?)
            idx = where(logical_and(dist > 0, dist <= threshold))[0] # TODO : change threshold  pick 1 if gridized
            # idx += 1 #zero index shananingans
            A.append((idx + 1, dist[idx]))                # add indices to list

    # show distribution of supervoxels memberships
    if plot:
        superVoxMember = [len(i) for i in MAP.values()]
        fig, ax = plt.subplots()
        ax.hist(superVoxMember, unique(superVoxMember))
        ax.set_xlabel('Membership size')
        ax.set_ylabel('Count')
        ax.set_title('Distribution of membership size of super voxels')
    return MAP, A

def splitConditions(data, conditions):
    '''
    Assumes the input is ncell x ntime and conditions is ntime long
    it will split the data per trial and time according to the uniques in conditions
    '''

    uConditions = unique(conditions)
    dataPerCondition = {}
#    dataPerConditionAndTrial = {}
    for idx, condition in enumerate(uConditions):
        marker = zeros((data.shape[1],))
        condIdx          = where(conditions == condition)[0]
        marker[condIdx]  = 1                                     # mark where condition holds
        trialMarker      = where(diff(marker))[0]                # supposed to be where the difference is 1 which is at the boundaries
        trialDur         = diff(trialMarker)[0]                  # assuming trial duration is constant, just take first marker
        nTrial           = len(where(marker)[0]) / trialDur      # you always count 1 edge less if edge includes starting point thus floor

        d = mod(nTrial * trialDur,10)
        # NOTE: ZEROPADDING IF NTRIALS IS NOT INTEGER!
        if d != 0:
            d = 10 - d # add the difference in zeros to the array
            print(data.shape[0])
            dd =  zeros((data.shape[0], d))

            print(dd.shape)
            tmp = hstack((data[:, condIdx],dd))

        else:
            tmp = data[:, condIdx]
        dataPerCondition[condition] = tmp.reshape(-1,trialDur, nTrial.round())

#        plt.plot(diff(trialMarker),'.')
#        print(trialMarker[1::2], trialMarker[0::2])
    return dataPerCondition

def loadSubset(MAP, conditions, file, Demean = True):
    '''
    Loads a subset of the data based on create super voxel
    '''

    subset = zeros((len(MAP.keys()), len(conditions)))
    with File(loadThis, 'r') as f:
        gr = f['Raw']
        for idx, (key, values) in enumerate(tqdm(MAP.items())):
            tmp = []
            for member in values:
#                print(member)
                tmp.append(gr[str(member)].value[conditions])

            subset[idx, :] = mean(array(tmp), axis = 0)
    return subset

#def loadSuperVox(MAP, file, subset):
#    '''
#    subset = subset path i.e. omr/0
#    map is created with createSuperVox
#    file is file path
#    '''
#    print(subset)
#    with File(file, 'r') as f:
#        gr = f[subset]
#        store = []
#        for key, values in tqdm(MAP.items()):
#            tmp = []
#            for member in values:
#                tmp.append(gr[str(member)].value)
#            tmp = array(tmp).mean(0) # supervoxel data
#            store.append(tmp)
#    return array(store)
# old versiont
def loadSuperVox(MAP, conditions, file, norm = True, condition = 'Deconvolved'):
    # TODO: make it work with new structure, i.e. h5p omr 0 tests
     subset = {}
     with File(file, 'r') as f:
         gr = f[condition]
         # obtain super voxel average
         for key, values in tqdm(MAP.items()):
             tmp = []
             for member in values:
                 tmp.append(gr[str(member)].value)
             tmp = array(tmp).mean(0) # supervoxel data
             # now sort according to conditions
             for condition, idxs in conditions.items():
                 for set, idx in enumerate(idxs):
                     # test if key already exists
                     try:
                         subset[condition][set].append(tmp[idx])
                     # else make list size of idxs
                     except KeyError:
                         subset[condition] = [[] for i in range(len(idxs))]
                         subset[condition][set].append(tmp[idx])
     for key, values in subset.items():
         for idx, value in enumerate(values):
 #            tmp = signal.detrend(value, 1)
 #            tmp = stats.zscore(value, 1)
             subset[key][idx] = array(value)
     return subset
def loadConditions(conditions, subsets, rawFilePath):
    '''
    Loads the conditions from the raw rawFilePath folder
    conditions : array, containing condition numbers
    subsets : dicts containing bounds of conditions
    rawFilePath : string pointing to the project folder
    '''
    idxOfConditions = {}
    for condition, lim in subsets.values():
        if condition not in idxOfConditions.keys():
            idxOfConditions[condition] = []

# def loadConditions(conditions, file, subset = 'Raw',
#                     core = '../../ZebrafishRecordings/',
#                     registered = '/Registered/'):
#     '''
#     Loads file from either a h5py or a binary file
#     Conditions : a dictionary where the values indicate the time steps that should be extracted
#     file    : name of datasets (not the directory)
#     subset  : if h5py exists loads a subset
#     Core : base directory where the raw binaries are stored
#
#     returns a dictionary per condition
#     '''
#     from scipy import stats, signal
#     import os
#
    # def addToConditions(conditions, data):
    #     dataPerConditions = {}
    #     for condition, idxs in conditions.items():
    #         # for all the indices of the subsets
    #         if condition not in dataPerConditions.keys():
    #             dataPerConditions[condition] = [[]] * len(idxs)
    #         for num, idx in enumerate(idxs):
    #             dataPerConditions[condition][num] = data[:, idx]
    #     return dataPerConditions
    # def findDataset(file):
    #     dirs = os.listdir(core)
    #     for dir in dirs:
    #         try:
    #             print(dir)
    #             if file in os.listdir(core + dir):
    #                 print('Datset found')
    #                 if 'Registered' in os.listdir(core + dir + '/' + file):
    #                     file = core + dir + '/' + file  + registered
    #                 else:
    #                     file = core + dir + '/' + file
    #                 return file
    #         except NotADirectoryError:
    #             print('Not a directory, skipping')
    #     return None
    # file = findDataset(file)
    #
    # if file != None:
    #     data  = createOrLoadHDF(file)
    #     data  = addToConditions(conditions, data)
    #     return data
    # else:
    #     print ('File not found')

def loadSubsetOfData(fileName, condition, subset):
    '''
    Loads specific subset from condition in in the dataset
    This function is aimed to be used after addSubsetToDataset is run
    '''
    from tqdm import tqdm
    print('Loading condition {0} set {1}'.format(condition, subset))
    with File(fileName) as f :
        gr   = f[condition][str(subset)]
        data = array([gr[member].value for member in tqdm(gr)])
    return data
def addSubsetToDataset(fileName, conditions, dataPath = '../Data/'):
    """
    Loads the raw set and subsets the data into different datasets
    Conditions is a dict with name as keys and array of subsets with indices of the raw subset
    """
    from os import listdir
    from tqdm import tqdm
    # TODO: clean up try except
    # obtain the raw dataset eiher from h6py or from the mat
    compute = True
    try:
        print('Loading from h5py')
        print(fileName)
        with File(fileName) as f:
            if all([f.__contains__(key) for key in conditions.keys()]):
                print('conditions e xists, skipping')
                compute = False
            else:
                raw = array([f['Raw'][member].value for member in tqdm(f['Raw'])])
    except:
        print('Creating h5py from binaries')
        raw = createOrLoadHDF(fileName)

    # if the conditions are not present add them to file
    if compute:
        # keep file active
        with File(fileName) as f:
            # loop through dictionariy
            for conditionName, subsets in conditions.items():
                # create the different stimulus conditions
                try:
                    print('Creating group {0}'.format(conditionName))
                    f.create_group(conditionName)
                except:
                    print('Group {0} exists'.format(conditionName))
                gr = f[conditionName] # Group
                for subset, timePoints in enumerate(subsets):
                    # check if group is already made
                    try:
                        gr.create_group(str(subset))
                    except:
                        continue
                    # make dataset per cell
                    for cellIdx, cellData in enumerate(tqdm(raw)):
                        try:
                            gr[str(subset)].create_dataset(str(cellIdx), data = cellData[timePoints])
                        except:
                            print('Dataset already exists')


def seperateConditions(conditions, subsets):
    '''
    returns a dictionary condtaining list of indices
    '''
    idxReturn = {}
    for key, indices in subsets.items():
        idxReturn[key] = []
        for idxSet in indices:
            # if only one number is present, just search this condition
            if len(idxSet) > 1:
                start, stop = idxSet
                tmp = where(logical_and(conditions >= start, conditions <= stop))[0]
#                    dif = where(diff(tmp)> 1)[0][0]
#                    tmp = [tmp[:dif], tmp[dif + 1:]]
            # rest condition has only 1 marker
            else:
                tmp = where(conditions == idxSet)[0]
            j = where(diff(tmp) > 1)[0]
            y = []
            j = hstack((j, len(tmp)))
            start = 0
            for stop in j:
                idxReturn[key].append(tmp[start:stop])
                start = stop
    return idxReturn

def createOrLoadHDF(datasetName):
    '''
    Loads from mat or from h5py if present
    returns raw dataset
    '''
    import os
    if 'Registered' in datasetName:
        idx = -3
    else:
        idx = -1
    fname = datasetName.split('/')[idx]
    print('Dataset name {0}'.format(fname))

    loadThis =  fname + '.hdf5'
    if loadThis in os.listdir('../Data/'):
        print('Datset found loading hdf5')
        with File('../Data/' + loadThis) as f:
            for i in f:print(i)
            gr = f['Raw']
            for idx, member in enumerate(tqdm(gr)):
                cellData = gr[member].value
                # append seems to cause temporary memory issue
                if idx is 0:
                    data = zeros((len(gr), len(cellData)))
                data[idx, :] = cellData
#                data.append(gr[member].value)

    else:
        print('HDF5 not found\nLoading binaries')
        datasetPath     = datasetName
        if 'Registered' in os.listdir(datasetPath):
             datasetPath = datasetPath + 'Registered'
#        assert 0
        from scipy.io import loadmat
        dims = loadmat(datasetPath + 'cell_resp_dim.mat')['cell_resp_dim'].flatten().tolist()
        data = fromfile(datasetPath + 'cell_resp_lowcut.stackf', dtype = float32).reshape(dims[::-1]).T
        print('Data size', data.shape)
        print('Creating HDF5')
        with File('../Data/' + loadThis, 'w') as f:
            f.create_group('Raw')
            for idx, cell in enumerate(tqdm(data)):
                f['Raw'].create_dataset(str(idx), data = cell)
    print('Data shape {0}'.format(data.shape))
    return data


def loadFrameTurn(file):
    from scipy.io import loadmat
    # if mat file try and loat with loadmat
    try:
        print('Attempting load with loadmat')
        return loadmat(file)['frame_turn']
    # check if frame_turn new
    except:
        tmp = file.split('.mat')
        tmp = tmp[0] + '_new.mat'
        print('Attempting load of {0}', tmp)
        # load with mat again else try h6py
        try:
            return loadmat(tmp)['frame_turn']
        except NotImplementedError:
            from h5py import File
            with File(tmp) as f:
                return f['frame_turn'].value
#    with File(file) as f:
#        frameTurn = f['frame_turn'].value
#    return frameTurn
# %%


# %%
opts = {'steps': 1, 'hyp': {'a0' : 2., 'b0' : 1., 'mu0' : 0., 'kappa0' : 1.} }

# %%
#print(MAPP)
#store = array(store); print(store.shape)
## %%
#behavior = frame_turn[14, conditions]
#from scipy import stats
#res = array([stats.pearsonr(i, behavior.flatten()) for i in store])
#
#import matlab.engine
#eng = matlab.engine.start_matlab()
#eng.cd('../ddCRP-master/')
## %%
#toMatlab =  matlab.double(data.T.tolist())
## %%
#A = [matlab.int32(i) for i in A] # otherwise the transfer crashes
#MAPP, samples = eng.PMC_ddCRP_NG(toMatlab, A, opts, nargout = 2)
#import matlab.engine
#eng = matlab.engine.start_matlab()
#eng.cd('../cell_view_package/')
# %%
#test = eng.test(A, nargout = 1); print(test) #test script for interfacing with matlab


if __name__ == '__main__':
    # test part of these functions
    datasetPath     = '../ZebrafishRecordings/photo_omr_blob_etc_data/20150120_2_1_photo_OMR_prey_blob_blue_cy74_6d_20150120_220356/Registered/'
    tmp = '20150122_2_2_cy74_6d_photo_OMR_prey_blob_blue_20150122_190028'
    data = loadData(tmp)
    print(data.shape)
