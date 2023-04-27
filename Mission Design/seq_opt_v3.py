import pygmo as pg
from pykep import epoch, epoch_from_string
from pykep.planet import jpl_lp
from pykep.trajopt import mga, mga_1dsm
from pykep.util import load_spice_kernel
from pykep.planet import spice
import pykep as pk
import numpy as np
from scipy.optimize import root_scalar, minimize_scalar
import plotly.graph_objects as go
import pandas as pd
import plotly
import plotly.io as pio
from tqdm import tqdm
from parfor import parfor

from mga_mod import mga_mod
from mga_1dsm import mga_1dsm_mod

pio.renderers.default = "browser"

load_spice_kernel("../current/de721_full.bsp")

def calc_dV_LEO(alt, lambert_arc, earth, pop):
    # Assume v_c os on ecliptic plane
    # alt in km
    v_c = np.sqrt(earth.mu_self/(earth.radius + alt*1000))

    v_inf_helio = lambert_arc.get_v1()
    earth_helio = np.array(earth.eph(pop.champion_x[0])[1])
    v_inf_earth = v_inf_helio - earth_helio

    spec_e = .5*np.linalg.norm(v_inf_earth)**2
    v_p = np.sqrt(2*(spec_e + earth.mu_self/(earth.radius + alt*1000)))

    return v_p - v_c

def plot_all(planet_rs, traj):
    layout = go.Layout(
        width=1024,
        height=1024,
        scene=dict(aspectmode='data')
    )

    fig = go.Figure(layout = layout, data=go.Scatter3d(
        x=traj[:, 0], y=traj[:, 1], z=traj[:, 2],
        line=dict(
            color='darkblue',
            width=1
        ),
        mode='lines'
    ))

    for planet_r in planet_rs:
        fig.add_scatter3d(x=planet_r[:, 0], y=planet_r[:, 1], z=planet_r[:, 2], mode='lines')

    fig.show()


def lambert_xyz(l, N=60, sol=0, units=1.0, color='b', legend=False, axes=None, alpha=1.):
    """
    ax = plot_lambert(l, N=60, sol=0, units='pykep.AU', legend='False', axes=None, alpha=1.)
    - axes:     3D axis object created using fig.add_subplot(projection='3d')
    - l:        pykep.lambert_problem object
    - N:        number of points to be plotted along one arc
    - sol:      solution to the Lambert's problem we want to plot (must be in 0..Nmax*2)
                where Nmax is the maximum number of revolutions for which there exist a solution.
    - units:    the length unit to be used in the plot
    - color:    matplotlib color to use to plot the line
    - legend:   when True it plots also the legend with info on the Lambert's solution chosen
    Plots a particular solution to a Lambert's problem
    Example::
      import pykep as pk
      import matplotlib.pyplot as plt
      fig = plt.figure()
      ax = fig.add_subplot(projection='3d')
      t1 = pk.epoch(0)
      t2 = pk.epoch(640)
      dt = (t2.mjd2000 - t1.mjd2000) * pk.DAY2SEC
      pl = pk.planet.jpl_lp('earth')
      pk.orbit_plots.plot_planet(pl, t0=t1, axes=ax, color='k')
      rE,vE = pl.eph(t1)
      pl = pk.planet.jpl_lp('mars')
      pk.orbit_plots.plot_planet(pl, t0=t2, axes=ax, color='r')
      rM, vM = pl.eph(t2)
      l = lambert_problem(rE,rM,dt,pk.MU_SUN)
      pk.orbit_plots.plot_lambert(l, ax=ax, color='b')
      pk.orbit_plots.plot_lambert(l, sol=1, axes=ax, color='g')
      pk.orbit_plots.plot_lambert(l, sol=2, axes=ax, color='g', legend = True)
      plt.show()
    """
    from pykep import propagate_lagrangian, AU
    import numpy as np
    import matplotlib.pylab as plt
    from mpl_toolkits.mplot3d import Axes3D

    if sol > l.get_Nmax() * 2:
        raise ValueError(
            "sol must be in 0 .. NMax*2 \n * Nmax is the maximum number of revolutions for which there exist a solution to the Lambert's problem \n * You can compute Nmax calling the get_Nmax() method of the lambert_problem object")

    # We extract the relevant information from the Lambert's problem
    r = l.get_r1()
    v = l.get_v1()[sol]
    T = l.get_tof()
    mu = l.get_mu()

    # We define the integration time ...
    dt = T / (N - 1)

    # ... and allocate the cartesian components for r
    x = np.array([0.0] * N)
    y = np.array([0.0] * N)
    z = np.array([0.0] * N)

    # We calculate the spacecraft position at each dt
    for i in range(N):
        x[i] = r[0] / units
        y[i] = r[1] / units
        z[i] = r[2] / units
        r, v = propagate_lagrangian(r, v, dt, mu)

    return np.vstack([x, y, z]).T


def conc_lamberts(ls, Ns):
    xyz_list = []
    for i, l in enumerate(ls):
        xyz_list.append(lambert_xyz(l, N=Ns[i]))

    return np.vstack(xyz_list)

def calc_planets(udp, pop, seg_list, seq):
    epoch = pk.epoch(pop.champion_x[0], julian_date_type='mjd2000')
    tofs = udp._decode_tofs(pop.champion_x)

    tof_list = np.linspace(pop.champion_x[0], pop.champion_x[0]+tofs[0], seg_list[0])
    epoch_list = [pk.epoch(time, julian_date_type='mjd2000') for time in tof_list]

    start_e = epoch_list[-1]
    start_tof = tof_list[-1]

    for i in range(len(seg_list) - 1):
        cur_tof_list = np.linspace(start_tof, start_tof+tofs[i+1], seg_list[i+1])
        cur_epoch = [pk.epoch(time, julian_date_type='mjd2000') for time in cur_tof_list]

        epoch_list += cur_epoch
        tof_list += cur_tof_list

        start_tof = cur_tof_list[-1]
    all_rs = []
    for j in range(len(seq)):
        cur_rs = np.zeros([len(epoch_list), 3])

        for i in range(len(epoch_list)):
            cur_rs[i] = seq[j].eph(epoch_list[i])[0]

        all_rs.append(cur_rs)

    return epoch_list, all_rs


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
    alt_LEO = 200_000

    free_v_LEO = 6600
    v_c = np.sqrt(e.mu_self/(e.radius + alt_LEO))
    v_hyp = v_c + free_v_LEO
    spec_e = .5*v_hyp**2 - e.mu_self/(e.radius + alt_LEO)

    vinf_free = np.sqrt(2*spec_e)/1000
    seqs1 = [[e, v, e, u],
           [e, v, e, e, u],
           [e, v, v, v, u],
           [e, v, e, m, u],
           [e, v, e, e, m, u],
           [e, v, e, v, e, u],
           [e, v, e, m, e, u],
           [e, v, e, v, e, e, u],
           [e, v, e, m, e, e, u],
           [e, j, u],
           [e, v, e, j, u],
           [e, v, v, v, j, u],
           [e, v, e, e, j, u],
           [e, v, e, v, e, j, u]]

    seqs2 = [[e, v, e, u],
            [e, v, e, e, u],
            [e, v, v, v, u],
            [e, v, e, m, u],
            [e, v, e, e, m, u],
            [e, v, e, m, e, u],
            [e, v, e, m, e, e, u],
            [e, v, e, j, u]]

    seqs3 = [[e, v, e, j, u],
            [e, v, v, v, j, u],
            [e, v, e, e, j, u],
            [e, v, e, v, e, j, u]]

    seqs = [[e, j, u]]

    obj = []
    dV_launch_LEO = []
    dV_launch_direct = []
    dv_deep_space = []
    dv_capture = []

    tof = np.array([12])

    T = 30 * pk.DAY2SEC

    r_p_mult2 = 1.3

    sma2 = (T / 2 / np.pi * np.sqrt(mu_uranus)) ** (2 / 3)
    ecc2 = 1 - r_p_mult2 * r_uranus / sma2

    start = epoch_from_string('2028-01-20 23:59:54.003')
    end = epoch_from_string('2038-01-20 23:59:54.003')
    for to in tqdm(tof):
        obj_tmp = []
        launch_tmp_LEO = []
        deep_tmp = []
        dv_cap_tmp = []
        for seq in tqdm(seqs1):
            udp = mga_1dsm_mod(
                seq=seq,
                vinf=[0, vinf_free],
                t0=[start, end],
                tof=to * 365.25,
                orbit_insertion=True,
                rp_target=r_uranus * r_p_mult2,
                e_target=ecc2,
                tof_encoding='eta',
                max_revs=4
            )

            prob = pg.problem(udp)
            # We solve it!!
            #uda = pg.xnes(gen=200, force_bounds=True, xtol=1e-12)
            uda = pg.mbh(algo=pg.scipy_optimize(), perturb=.9, stop=300)
            algo = pg.algorithm(uda)
            algo.set_verbosity(1)

            # We construct a random population of 20 individuals (the initial guess)
            pop = pg.population(prob, size=10)
            # We solve the problem
            pop = algo.evolve(pop)

            dv_struct = udp._compute_dvs(pop.champion_x)
            seg_list = list((5000*np.ones(len(dv_struct[3]))).astype(int))
            data = conc_lamberts(dv_struct[3], seg_list)
            #pd.DataFrame(data).to_csv('test.csv')

            epoch_lists, planet_rs = calc_planets(udp, pop, seg_list, list(set(seq)))

            all_data = pd.DataFrame(np.hstack([np.hstack(planet_rs), data]))

            #plot_all(planet_rs, data)

            tot_tof = np.sum(np.array(udp._decode_tofs(pop.champion_x)))
            mjd = pop.champion_x[0] + tot_tof

            planet_encounters = pop.champion_x[0] + np.cumsum(udp.eta2direct(pop.champion_x))

            obj_tmp.append(udp.fitness(pop.champion_x)[0])
            launch_tmp_LEO.append(calc_dV_LEO(200, dv_struct[3][0], e, pop))
            deep_tmp.append(np.sum(dv_struct[1]))

            tot_tof = np.sum(np.array(udp._decode_tofs(pop.champion_x)))
            mjd = pop.champion_x[0] + tot_tof

            v_uran_arr = np.array(list(u.eph(pk.epoch(mjd, "mjd2000"))[1]))
            v_sc_helio_arr = dv_struct[3][-1].get_v2()

            v_rel_uran = v_sc_helio_arr - v_uran_arr
            dv_cap_tmp.append(dv_struct[2])

        tmp = np.vstack([obj_tmp, launch_tmp_LEO, deep_tmp, dv_cap_tmp]).T

        obj.append(obj_tmp)
        dV_launch_LEO.append(launch_tmp_LEO)
        dv_deep_space.append(deep_tmp)
        dv_capture.append(np.linalg.norm(dv_cap_tmp))






