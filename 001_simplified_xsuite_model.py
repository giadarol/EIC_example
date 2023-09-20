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
    Strategy(slicing=Teapot(1)),  # Default catch-all as in MAD-X
    Strategy(slicing=Teapot(10), element_type=xt.Bend),
    Strategy(slicing=Teapot(50), element_type=xt.CombinedFunctionMagnet),
    Strategy(slicing=None, element_type=xt.Quadrupole),
    Strategy(slicing=None, element_type=xt.Solenoid),
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