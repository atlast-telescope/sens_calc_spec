import numpy as nm
import pylab as pl
import matplotlib.patches as patches

# load loading data
tmp = nm.loadtxt('am_runs/JJA_25.dat')
f_GHz  = tmp[1:,0]
atm_tx = tmp[1:,2]
T_atm  = tmp[1:,3]

pl.ion()
pl.plot(f_GHz,atm_tx)

patches.Rectangle(
        (130.0, 0.0), 30.0, 1.0,
        hatch='\\',
        fill=False)