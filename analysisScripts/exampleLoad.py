# -*- coding: utf-8 -*-
"""
Created on Sun Jul 30 10:42:00 2017

@author: Casper
"""
from numpy import *
from tqdm import tqdm # status bar for loading, [not necessary but remove tqdm reference below]

fileName = '20150122_2_2_cy74_6d_photo_OMR_prey_blob_blue_20150122_190028.h5'

def loadCellDataSimple(fileName, groupName = 'Frame turn'):
    from h5py import File
    with File(fileName) as f:
        # for group/x where x is cell data
        try:
            return array([f[groupName + '/' + i].value for i in tqdm(f[groupName])])
        # for single vectors, i.e. conditions
        except TypeError:
            return array(f[groupName].value)
            
c = loadCellDataSimple(fileName)