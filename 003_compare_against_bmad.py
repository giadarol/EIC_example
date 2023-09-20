from cpymad.madx import Madx
import xtrack as xt
import xdeps as xd

import numpy as np

mad = Madx()
mad.call("hsr_thick-275-10_20230817.seq")
mad.beam(particle="proton",pc=275,charge=-1)
mad.use(sequence='hsr_thick')
t0 = xd.Table(mad.twiss(sequence='hsr_thick',chrom=True,exact=True))

line=xt.Line.from_madx_sequence(mad.sequence['hsr_thick'], allow_thick=True,
                                deferred_expressions=True)
line.particle_ref = xt.Particles(mass0=xt.PROTON_MASS_EV, q0=-1, p0c=275e9)
line.config.XTRACK_USE_EXACT_DRIFTS = True
line.configure_bend_model(edge='full', core='full')
t1 = line.twiss(method='4d', freeze_longitudinal=True)

import pandas as pd
df = pd.read_csv('twiss-275-10.txt', sep=' ', skipinitialspace=True)

import matplotlib.pyplot as plt


plt.close('all')

# Comparison Xsuite BMAD
plt.figure(1)
sp1 = plt.subplot(3,1,1)
plt.plot(df.S, df.BETA11, label='BMAD')
plt.plot(t1.s, t1.betx, label='Xsuite', linestyle='--')
plt.legend()
plt.ylabel(r'$\beta_x$ [m]')

plt.subplot(3,1,2,sharex=sp1)
plt.plot(df.S, df.X, label='BMAD')
plt.plot(t1.s, t1.x, label='Xsuite', linestyle='--')
plt.ylabel(r'$x$ [m]')

plt.subplot(3,1,3,sharex=sp1)
plt.plot(df.S, df.DX, label='BMAD')
plt.plot(t1.s, t1.dx, label='Xsuite', linestyle='--')
plt.ylabel(r'$D_x$ [m]')
plt.xlabel('s [m]')
plt.suptitle('Comparison Xsuite vs BMAD')

# Comparison Mad-X BMAD
plt.figure(2)
sp1 = plt.subplot(3,1,1)
plt.plot(df.S, df.BETA11, label='BMAD')
plt.plot(t0.s, t0.betx, label='MAD-X', linestyle='--')
plt.legend()
plt.ylabel(r'$\beta_x$ [m]')

plt.subplot(3,1,2,sharex=sp1)
plt.plot(df.S, df.X, label='BMAD')
plt.plot(t0.s, t0.x, label='MAD-X', linestyle='--')
plt.ylabel(r'$x$ [m]')

plt.subplot(3,1,3,sharex=sp1)
plt.plot(df.S, df.DX, label='BMAD')
plt.plot(t0.s, t0.dx, label='MAD-X', linestyle='--')
plt.ylabel(r'$D_x$ [m]')
plt.xlabel('s [m]')
plt.suptitle('Comparison MAD-X vs BMAD')

plt.show()

