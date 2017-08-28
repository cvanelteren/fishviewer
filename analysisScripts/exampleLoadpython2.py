# -*- coding: utf-8 -*-
"""
Created on Wed Aug  9 13:15:59 2017

@author: Casper
"""

from h5py import File
def exampleOld(fileName, condition):
    out = [] 
    with File(fileName) as f:
        gr = f[condition]
        for i in gr:
            out.append(gr[i].value)
    return out



g = '20150122_2_2_cy74_6d_photo_OMR_prey_blob_blue_20150122_190028.h5'
c = exampleOld(g, 'Conditions')

