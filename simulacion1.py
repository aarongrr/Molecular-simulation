import hoomd
import hoomd.hpmc
from hoomd import md
import numpy as np 
from mpmath import *
import matplotlib as plt
mp.dps = 25; mp.pretty = True

t=2
hoomd.context.initialize()
hoomd.init.create_lattice(unitcell=hoomd.lattice.unitcell(N=4,
                                                        a1=[t,0,0],
                                                        a2=[0,t,0],
                                                        a3=[0,0,t],
                                                        dimensions=3,
                                                        position=[[0,0,0],[0,-t/2,-t/2],
                                                        [-t/2,0,-t/2],[-t/2,-t/2,0]],
                                                        type_name=["A","A","A","A"],
                                                        mass=[1,1,1,1],
                                                        charge=[0,0,0,0],
                                                        diameter=[1.0,1.0,1.0,1.0],
                                                        moment_inertia=None,
                                                        orientation=None),n=7)


# Potencial de Leonard Jhones
nl = md.nlist.cell()
lj = md.pair.lj(r_cut=3.0, nlist=nl)
lj.pair_coeff.set('A', 'A', epsilon=2.0, sigma=1.0)

all = hoomd.group.all()
md.integrate.mode_standard(dt=0.0005)
hoomd.md.integrate.langevin(group=all, kT=1.2, seed=4)

hoomd.dump.gsd("trajectory.gsd", period=2000, group=all, overwrite=True)
# datos
hoomd.analyze.log(filename="kinectic.dat",quantities=['kinetic_energy'], period=200, overwrite=True)
hoomd.analyze.log(filename="pressure.dat",quantities=['pressure'], period=200, overwrite=True)
hoomd.analyze.log(filename="potential.dat",quantities=['potential_energy'], period=200, overwrite=True)
hoomd.analyze.log(filename="volume.dat",quantities=['volume'], period=200, overwrite=True)
nl.set_params(r_buff = 1.6)

hoomd.run(1e5)

