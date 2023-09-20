from cpymad.madx import Madx
import xtrack as xt
import xdeps as xd

import numpy as np

mad = Madx()
mad.call("hsr_thick-275-10_20230817.seq")
mad.beam(particle="proton",pc=275,charge=-1)
mad.use(sequence='hsr_thick')

line_thick=xt.Line.from_madx_sequence(mad.sequence['hsr_thick'], allow_thick=True,
                                deferred_expressions=True)
line_thick.particle_ref = xt.Particles(mass0=xt.PROTON_MASS_EV, q0=-1, p0c=275e9)
line_thick.config.XTRACK_USE_EXACT_DRIFTS = True
line_thick.configure_bend_model(edge='full', core='full')
tw_thick = line_thick.twiss(method='4d', freeze_longitudinal=True)

line = line_thick.copy()

Strategy = xt.slicing.Strategy
Teapot = xt.slicing.Teapot

slicing_strategies = [
    Strategy(slicing=None),  # Default
    Strategy(slicing=Teapot(1), element_type=xt.Sextupole),
    Strategy(slicing=Teapot(3), name='^yo.*dh.*'), # Arc bends
    Strategy(slicing=Teapot(3), name='^yi.*dh.*'), # Arc bends
    Strategy(slicing=Teapot(3), name='^bo.*dh.*'), # Arc bends
    Strategy(slicing=Teapot(3), name='^bi.*dh.*'), # Arc bends
    Strategy(slicing=Teapot(10), name='^yo.*qf.*'), # Arc quads
    Strategy(slicing=Teapot(10), name='^yi.*qf.*'), # Arc quads
    Strategy(slicing=Teapot(10), name='^bo.*qf.*'), # Arc quads
    Strategy(slicing=Teapot(10), name='^bi.*qf.*'), # Arc quads
    Strategy(slicing=Teapot(10), name='^yo.*qd.*'), # Arc quads
    Strategy(slicing=Teapot(10), name='^yi.*qd.*'), # Arc quads
    Strategy(slicing=Teapot(10), name='^bo.*qd.*'), # Arc quads
    Strategy(slicing=Teapot(10), name='^bi.*qd.*'), # Arc quads
    Strategy(slicing=Teapot(3), name='^yi.*tq.*'),
    Strategy(slicing=Teapot(3), name='^yo.*tq.*'),
    Strategy(slicing=Teapot(3), name='^bo.*tq.*'),
    Strategy(slicing=Teapot(3), name='^bi.*tq.*'),
    Strategy(slicing=Teapot(10), name='^qds.*'),
    Strategy(slicing=Teapot(10), name='^qus.*'),
    Strategy(slicing=Teapot(10), name='^warm_quad.*'),
    Strategy(slicing=None, name='^qff.*'), # Touchy, not sliced
    Strategy(slicing=None, name='^q1.*'),  # Touchy, not sliced
    Strategy(slicing=None, name='^q2.*'),  # Touchy, not sliced
    Strategy(slicing=None, name='^q3.*'),  # Touchy, not sliced
    Strategy(slicing=None, name='^q4.*'),  # Touchy, not sliced
    Strategy(slicing=None, name='^q5.*'),  # Touchy, not sliced
    Strategy(slicing=None, name='b2pr'),
    Strategy(slicing=None, name='bxds9m2'),
    Strategy(slicing=None, name='bxdsds04'),
    Strategy(slicing=None, name='bxds9m1'),
    Strategy(slicing=None, name='bxds02disp1'),
    Strategy(slicing=None, name='bxds01b'),
    Strategy(slicing=None, name='bxds01a'),
    Strategy(slicing=None, name='bxsp01'),
    Strategy(slicing=None, name='bxus9m1'),
    Strategy(slicing=None, name='bxus9m2'),
    Strategy(slicing=None, name='bxus9m3'),
    Strategy(slicing=None, name='dsw_ir10'),
    Strategy(slicing=None, name='dwarm_ir10'),
    Strategy(slicing=None, name='dwarm_ir12'),
    Strategy(slicing=None, name='dsw_ir12'),
    Strategy(slicing=None, name='d5'),
    Strategy(slicing=None, name='d5:0'),
    Strategy(slicing=None, name='dwarm3'),
    Strategy(slicing=None, name='dwarm4'),
    Strategy(slicing=None, name='h5_dh4'),
    Strategy(slicing=None, name='b2pf'),
    Strategy(slicing=None, name='b1apf'),
    Strategy(slicing=None, name='b1pf'),
    Strategy(slicing=None, name='b0pf'),
    Strategy(slicing=None, name='b0apf'),
]

line.discard_tracker()
line.slice_thick_elements(slicing_strategies=slicing_strategies)
line.build_tracker()

tw = line.twiss(method='4d')

# Plot beta beating
tw = line.twiss(method='4d', only_markers=True)
tw_thick = line_thick.twiss(method='4d', only_markers=True)

tt = line.get_table()
tt_thick_elems = tt.rows[(~(tt.element_type=='Drift')) & tt.isthick]
print(f'{len(tt_thick_elems)} thick elements:')
tt_thick_elems.show()

import matplotlib.pyplot as plt
plt.close('all')

# Comparison Thick / Thin
plt.figure(1, figsize=(6.4*1.5, 4.8))
sp1 = plt.subplot(3,1,1)
plt.plot(tw_thick.s, tw_thick.betx, label='Thick')
plt.plot(tw.s, tw.betx, label='Thin', linestyle='--')
plt.legend()
plt.ylabel(r'$\beta_x$ [m]')

plt.subplot(3,1,2,sharex=sp1)
plt.plot(tw_thick.s, tw_thick.x, label='Thick')
plt.plot(tw.s, tw.x, label='Thin', linestyle='--')
plt.ylabel(r'$x$ [m]')

plt.subplot(3,1,3,sharex=sp1)
plt.plot(tw_thick.s, tw_thick.dx, label='Thick')
plt.plot(tw.s, tw.dx, label='Thin', linestyle='--')
plt.ylabel(r'$D_x$ [m]')
plt.xlabel('s [m]')
plt.suptitle('Comparison Thick vs Thin')

plt.figure(2, figsize=(6.4*1.5, 4.8))
sp1 = plt.subplot(2,1,1)
plt.plot(tw.s, tw.betx/tw_thick.betx - 1)
plt.ylabel(r'$\Delta \beta_x / \beta_x$')
sp2 = plt.subplot(2,1,2, sharex=sp1)
plt.plot(tw.s, tw.bety/tw_thick.bety - 1)
plt.ylabel(r'$\Delta \beta_y / \beta_y$')
plt.xlabel('s [m]')
plt.show()
