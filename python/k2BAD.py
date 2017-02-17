import math
import cmath
import numpy as np
import matplotlib.pyplot as plt
import os
import pdb

import brewer2mpl

import kali.k2


plt.ion()

agn_name = '220180147'
agn_campaign = 'c08'
home = os.environ['HOME']
outpath = os.path.join(home, 'Desktop', 'k2BADAnalysis')
if not os.path.isdir(outpath):
    os.mkdir(outpath)
maxColors = 12
primary = brewer2mpl.get_map('Paired', 'Qualitative', maxColors).hex_colors

starList = ['220180454',
            '220180018',
            '220179459',
            '229228541',
            '229228434',
            '220180838',
            '220181273',
            '220177656',
            '220180300',
            '220182981',
            '220181691',
            '220182985',
            '220179903',
            '220182091',
            '220178326',
            '220177750',
            '220179064',
            '229228948',
            '220181652']

try:
    agnLC_raw = kali.k2.k2LC(name=agn_name, campaign=agn_campaign, processing='raw', path=outpath)
    starLCList_raw = list()
    for star in starList:
        try:
            starLCList_raw.append(kali.k2.k2LC(name=star, campaign=agn_campaign, processing='raw',
                                               path=outpath))
        except ValueError:
            pass
except ValueError:
    agnLC_raw = None
    starLCList_raw = None

try:
    agnLC_mast = kali.k2.k2LC(name=agn_name, campaign=agn_campaign, processing='mast', path=outpath)
    starLCList_mast = list()
    for star in starList:
        try:
            starLCList_mast.append(kali.k2.k2LC(name=star, campaign=agn_campaign, processing='mast',
                                                path=outpath))
        except:
            pass
except ValueError:
    agnLC_mast = None
    starLCList_mast = None

try:
    agnLC_vj = kali.k2.k2LC(name=agn_name, campaign=agn_campaign, processing='vj', path=outpath)
    starLCList_vj = list()
    for star in starList:
        try:
            starLCList_vj.append(kali.k2.k2LC(name=star, campaign=agn_campaign, processing='vj',
                                              path=outpath))
        except ValueError:
            pass
except ValueError:
    agnLC_vj = None
    starLCList_vj = None

fig_raw = agnLC_raw.plot(fig=100, colory=r'#000000', labely='AGN: ' + agnLC_raw.name)
for num, starLC in enumerate(starLCList_raw):
    if num < maxColors:
        starLC.plot(fig=100, colory=primary[num], labely='Star: ' + starLC.name, alphay=0.2, clearFig=False)

fig_mast = agnLC_mast.plot(fig=200, colory=r'#000000', labely='AGN: ' + agnLC_mast.name)
for num, starLC in enumerate(starLCList_mast):
    if num < maxColors:
        starLC.plot(fig=200, colory=primary[num], labely='Star: ' + starLC.name, alphay=0.2, clearFig=False)

fig_vj = agnLC_vj.plot(fig=300, colory=r'#000000', labely='AGN: ' + agnLC_vj.name)
for num, starLC in enumerate(starLCList_vj):
    if num < maxColors:
        starLC.plot(fig=300, colory=primary[num], labely='Star: ' + starLC.name, alphay=0.2, clearFig=False)

pdb.set_trace()
