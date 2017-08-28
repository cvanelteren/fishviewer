#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 22 11:48:22 2017

@author: casper
"""

from numpy import *
from vispy import use, test
use('pyqt5')
from vispy.plot import Fig
fig  = Fig()

a = random.rand(100,2)
fig[0,0].plot(a)


test('pyqt5')
