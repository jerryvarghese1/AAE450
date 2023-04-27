import pygmo as pg
from pykep import epoch, epoch_from_string
from pykep.planet import jpl_lp
from pykep.trajopt import mga, mga_1dsm
from pykep.util import load_spice_kernel
from pykep.planet import spice
import pykep as pk
import numpy as np
from scipy.optimize import root_scalar, minimize_scalar, root, minimize
import plotly.graph_objects as go
import pandas as pd
import plotly
import plotly.io as pio
from scipy.interpolate import interp1d

from mga_mod import mga_mod
from mga_1dsm import mga_1dsm_mod
from mga_mod_minr import mga_mod_minr
import spiceypy as sppy

from scipy.spatial.transform import Rotation as R
from scipy.signal import savgol_filter

import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("Qt5Agg")

load_spice_kernel("de721_full.bsp")
load_spice_kernel("../current/ura111.bsp")
load_spice_kernel("pck00010.tpc")

sppy.furnsh("de721_full.bsp")
sppy.furnsh("pck00010.tpc")

if __name__ == '__main__':
    point_planet = False

    if point_planet:
        e = spice('EARTH', 'SUN', 'ECLIPJ2000', 'NONE', pk.MU_SUN, 3.986e14, 10.0, 10.0)

        mer = spice('MERCURY', 'SUN', 'ECLIPJ2000', 'NONE', pk.MU_SUN, 2.2032e13, 10.0, 10.0)
        m = spice('MARS', 'SUN', 'ECLIPJ2000', 'NONE', pk.MU_SUN, 4.282e13, 10.0, 10.0)
        v = spice('VENUS', 'SUN', 'ECLIPJ2000', 'NONE', pk.MU_SUN, 3.24859e14, 10.0, 10.0)
        j = spice('JUPITER BARYCENTER', 'SUN', 'ECLIPJ2000', 'NONE', pk.MU_SUN, 1.266865e17, 10.0, 10.0)
        s = spice('SATURN BARYCENTER', 'SUN', 'ECLIPJ2000', 'NONE', pk.MU_SUN, 3.793e16, 10.0, 10.0)

        u = spice('URANUS BARYCENTER', 'SUN', 'ECLIPJ2000', 'NONE', pk.MU_SUN, 5.793e15, 10.0, 10.0)

    else:
        e = spice('EARTH', 'SUN', 'ECLIPJ2000', 'NONE', pk.MU_SUN, 3.986e14, 6371000, 1.5 * 6371000)

        mer = spice('MERCURY', 'SUN', 'ECLIPJ2000', 'NONE', pk.MU_SUN, 2.2032e13, 2439400, 2439400 + 100_000)
        m = spice('MARS', 'SUN', 'ECLIPJ2000', 'NONE', pk.MU_SUN, 4.282e13, 3389500, 1.1 * 3389500)
        v = spice('VENUS', 'SUN', 'ECLIPJ2000', 'NONE', pk.MU_SUN, 3.24859e14, 6052000, 1.1 * 6052000)
        j = spice('JUPITER BARYCENTER', 'SUN', 'ECLIPJ2000', 'NONE', pk.MU_SUN, 1.266865e17, 69911000, 8 * 69911000)
        s = spice('SATURN BARYCENTER', 'SUN', 'ECLIPJ2000', 'NONE', pk.MU_SUN, 3.793e16, 136775000, 6 * 136775000)

        u = spice('URANUS BARYCENTER', 'SUN', 'ECLIPJ2000', 'NONE', pk.MU_SUN, 5.793e15, 25362000, 25362000 + 100_000)

    r_uranus = 25_362_000
    mu_uranus = 5.793939e15
    seq = [e, j, u]

    start = epoch_from_string('2028-01-20 23:59:54.003')

    starts = start.mjd2000 + np.linspace(1, 365*2, 300)
    all_epochs = [pk.epoch(a_start) for a_start in starts]

    end = epoch_from_string('2038-01-20 23:59:54.003')
    dvs = []
    udps = []
    pops = []

    alt_LEO = 200_000

    free_v_LEO = 6600
    v_c = np.sqrt(e.mu_self / (e.radius + alt_LEO))
    v_hyp = v_c + free_v_LEO
    spec_e = .5 * v_hyp ** 2 - e.mu_self / (e.radius + alt_LEO)

    vinf_free = np.sqrt(2 * spec_e) / 1000
    C3 = vinf_free**2

    T = 30 * pk.DAY2SEC

    r_p_mult1 = 5
    r_p_mult2 = 1.3

    tof_inter = 12
    mission_tof = 17

    sma = (T / 2 / np.pi * np.sqrt(mu_uranus)) ** (2 / 3)
    ecc1 = 1 - r_p_mult1 * r_uranus / sma
    ecc2 = 1 - r_p_mult2 * r_uranus / sma

    total_dv = [0]*len(starts)

    for i in range(len(starts)):
        udp = mga_mod(
            seq=seq,
            vinf=vinf_free,
            t0=[starts[i], end],
            tof= tof_inter * 365.25,
            orbit_insertion=True,
            rp_target=r_uranus*r_p_mult2,
            e_target=ecc2,
            tof_encoding='eta',
            dv_launch_weight=1,
            max_revs = 4
        )

        prob = pg.problem(udp)
        # We solve it!!
        #uda = pg.xnes(gen=8000, force_bounds=True, xtol=1e-12, ftol=1e-12)
        uda = pg.mbh(algo=pg.scipy_optimize(), perturb=.9, stop=100)
        algo = pg.algorithm(uda)
        algo.set_verbosity(0)

        # We construct a random population of 20 individuals (the initial guess)
        pop = pg.population(prob, size=1)
        # We solve the problem
        pop = algo.evolve(pop)

        dv_struct = udp._compute_dvs(pop.champion_x)
        seg_list = list((5000*np.ones(len(dv_struct[3]))).astype(int))

        total_dv[i] = np.sum(dv_struct[1]) + dv_struct[2]
        print(i)

    total_dv = np.array(total_dv)