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

expanded = []
full = []
for nn in line.element_names:
    ee = line[nn]
    if (isinstance(ee, xt.Bend) and (nn.startswith('yo') or nn.startswith('yi')
        or nn.startswith('bo') or nn.startswith('bi'))):
        ee.model = 'expanded'
        expanded.append(nn)
    elif isinstance(ee, xt.Bend):
        full.append(nn)

line['b2pr'].model          = 'full'
line['bxds9m2'].model       = 'full'
line['bxdsds04'].model      = 'full'
line['bxds9m1'].model       = 'full'
line['bxds02disp1'].model   = 'full'
line['bxds01b'].model       = 'full'
line['bxds01a'].model       = 'full'
line['bxsp01'].model        = 'full'
line['bxus9m1'].model       = 'full'
line['bxus9m2'].model       = 'full'
line['bxus9m3'].model       = 'full'
line['dsw_ir10'].model      = 'full'
line['dwarm_ir10'].model    = 'full'
line['dwarm_ir12'].model    = 'full'
line['dsw_ir12'].model      = 'full'
line['d5'].model            = 'full'
line['d5:0'].model          = 'full'
line['dwarm3'].model        = 'full'
line['dwarm4'].model        = 'full'
line['h5_dh4'].model        = 'full'
line['b2pf'].model          = 'full'
line['b1apf'].model         = 'full'
line['b1pf'].model          = 'full'
line['b0apf'].model         = 'full'

t2 = line.twiss(method='4d', freeze_longitudinal=True)

for nn in full:
    if ':' in nn: continue
    print(f'{nn:<15}', f'k0l = {mad.sequence.hsr_thick.expanded_elements[nn].k0 * mad.sequence.hsr_thick.expanded_elements[nn].l:<25}'
              f'angle = {mad.sequence.hsr_thick.expanded_elements[nn].angle}')


import matplotlib.pyplot as plt
plt.close('all')
plt.figure(1)
sp1 = plt.subplot(2,1,1)
plt.plot(t0.s, t0.dx, label='madx')
plt.ylabel(r'$D_x$')
plt.title('MAD-X')
sp2 = plt.subplot(2,1,2, sharex=sp1)
plt.plot(t1.s, t1.dx, label='Xsuite')
plt.ylabel(r'$D_x$')
plt.xlabel('s [m]')
plt.title('Xsuite')
