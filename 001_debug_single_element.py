from cpymad.madx import Madx
import xtrack as xt
import xdeps as xd

mad = Madx()
dr_len = 1e-6
mad.input(f"""
ss: sequence, l={dr_len};
    bxds01b..1: multipole, at=0, lrad=1.5, knl={0.005223894981}, angle=0.005;
    dr: drift, at={dr_len / 2}, l={dr_len};
endsequence;
beam;
use, sequence=ss;

twiss, betx=1, bety=1;
""")

tmad = xd.Table(mad.table.twiss)

line = xt.Line.from_madx_sequence(mad.sequence.ss)
line.particle_ref = xt.Particles(gamma0=mad.sequence.ss.beam.gamma)
txs = line.twiss(ele_start='ss$start', ele_stop='ss$end',
                 twiss_init=xt.TwissInit(
                     betx=1, bety=1, x=tmad.x[0], px=tmad.px[0],
                     dx=tmad.dx[0], dpx=tmad.dpx[0]))
