from cpymad.madx import Madx
import xtrack as xt
import xdeps as xd

mad = Madx()
mad.call("hsr_thin_xsuite.seq")
mad.beam(particle="proton",pc=275,charge=-1)
mad.use(sequence='hsr_thin')
t0 = xd.Table(mad.twiss(sequence='hsr_thin',chrom=True,exact=True))

line=xt.Line.from_madx_sequence(mad.sequence['hsr_thin'],deferred_expressions=True)
line.particle_ref = xt.Particles(mass0=xt.PROTON_MASS_EV, q0=-1, p0c=275e9)
line.config.XTRACK_USE_EXACT_DRIFTS = True
line.configure_bend_model(edge='full')
t1 = line.twiss(method='4d', freeze_longitudinal=True)

# twini=t1.get_twiss_init(0)

# twini.particle_on_co.x=t0.x[0]
# twini.particle_on_co.px=t0.px[0]
# t1 = line.twiss(method='4d', freeze_longitudinal=True,twiss_init=twini)
line.config.XTRACK_USE_EXACT_DRIFTS = True
two = line.twiss(ele_start=line.element_names[0], ele_stop=len(line)-1,
                 twiss_init=xt.TwissInit(betx=1, bety=1, x=t0.x[0], px=t0.px[0],
                                         dx=t0.dx[0], dpx=t0.dpx[0]))

import matplotlib.pyplot as plt
plt.close('all')
plt.figure(1)
sp1 = plt.subplot(3,1,1)
plt.plot(t0.s,t0.betx,label='madx')
plt.plot(t1.s,t1.betx,label='xtrack')

plt.subplot(3,1,2,sharex=sp1)
plt.plot(t0.s,t0.x,label='madx')
plt.plot(t1.s,t1.x,label='xtrack')

plt.subplot(3,1,3,sharex=sp1)
plt.plot(t0.s,t0.dpx,label='madx')
# plt.plot(t1.s,t1.dpx,label='xtrack')
plt.plot(two.s,two.dpx,label='xtrack2')

plt.legend()

plt.show()
