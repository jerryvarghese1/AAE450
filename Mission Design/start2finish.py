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

plt.style.use('dark_background')

pio.renderers.default = "browser"

load_spice_kernel("de721_full.bsp")
load_spice_kernel("../current/ura111.bsp")
load_spice_kernel("pck00010.tpc")

sppy.furnsh("de721_full.bsp")
sppy.furnsh("pck00010.tpc")

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

def plot_all(planet_rs, traj, planet_names, colors):

    layout = go.Layout(
        width=1024,
        height=1024,
        scene=dict(aspectmode='data'),
        template='plotly_dark'
    )

    fig = go.Figure(layout = layout, data=go.Scatter3d(
        x=traj[:, 0], y=traj[:, 1], z=traj[:, 2],
        line=dict(
            color='navajowhite',
            width=2
        ),
        mode='lines',
        name='Trajectory'
    ))

    for i, planet_r in enumerate(planet_rs):
        fig.add_scatter3d(x=planet_r[:, 0], y=planet_r[:, 1], z=planet_r[:, 2], mode='lines', name=planet_names[i], marker=go.Marker(color=colors[i]))

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

def kepler_xyz(r0, v0, tof, mu, N=60, units=1, color='b', label=None, axes=None):
    """
    ax = plot_kepler(r0, v0, tof, mu, N=60, units=1, color='b', label=None, axes=None):
    - axes:     3D axis object created using fig.add_subplot(projection='3d')
    - r0:       initial position (cartesian coordinates)
    - v0:       initial velocity (cartesian coordinates)
    - tof:      propagation time
    - mu:       gravitational parameter
    - N:	number of points to be plotted along one arc
    - units:	the length unit to be used in the plot
    - color:	matplotlib color to use to plot the line
    - label 	adds a label to the plotted arc.
    Plots the result of a keplerian propagation
    Example::
        import pykep as pk
        pi = 3.14
        pk.orbit_plots.plot_kepler(r0 = [1,0,0], v0 = [0,1,0], tof = pi/3, mu = 1)
    """

    from pykep import propagate_lagrangian
    import matplotlib.pylab as plt
    from mpl_toolkits.mplot3d import Axes3D
    from copy import deepcopy

    if axes is None:
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
    else:
        ax = axes

    # We define the integration time ...
    dt = tof / (N - 1)

    # ... and calculate the cartesian components for r
    x = [0.0] * N
    y = [0.0] * N
    z = [0.0] * N

    # We calculate the spacecraft position at each dt
    r = deepcopy(r0)
    v = deepcopy(v0)
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

    tof_list = np.linspace(pop.champion_x[0], pop.champion_x[0]+np.sum(tofs[0]), seg_list[0])
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

#def change_frame(v_inf, NAIF_ID, when, transpose):
#    if not(transpose):
#        DCM = sppy.tipbod("ECLIPJ2000", NAIF_ID, when)
#    else:

#    return (DCM@v_inf.reshape(-1, 1)).flatten()

def bplane_rp_ecc(v_inf_helio, mu, r_p, NAIF_ID, when, in_out):
    DCM = sppy.tipbod("ECLIPJ2000", NAIF_ID, when)

    v_inf_mag = np.linalg.norm(v_inf_helio)

    spec_e = .5*v_inf_mag**2
    sma = mu/(2*spec_e)
    ecc = r_p/sma + 1
    delta = 2*np.arcsin(1/ecc)
    theta = (np.pi - delta)/2

    if in_out=="in":
        A = R.from_euler('Z', -theta).as_matrix().T
    else:
        A = R.from_euler('Z', theta).as_matrix().T

    v_inf_hat = v_inf_helio/v_inf_mag

    rp_wrt_sun = r_p*A@v_inf_hat.reshape(-1, 1)

    rp_wrt_body = DCM@(rp_wrt_sun.reshape([-1, 1])).flatten()
    v_inf_wrt_body = DCM@(v_inf_helio.reshape([-1, 1])).flatten()

    h_wrt_body = np.cross(rp_wrt_body, v_inf_wrt_body)

    spec_e = .5*v_inf_mag**2
    v_p_mag = np.sqrt(2*(spec_e + mu/r_p))

    theta_hat = np.cross(h_wrt_body, rp_wrt_body)
    theta_hat /= np.linalg.norm(theta_hat)
    vp_wrt_body = v_p_mag*theta_hat

    return rp_wrt_body, vp_wrt_body, v_inf_wrt_body


def calc_flyby(ls, seq, udp, pop):

    encounters = pop.champion_x[0]
    encounters += np.cumsum(udp._decode_tofs(pop.champion_x))

    r_ps = [0]*(len(ls) - 1)
    deltas = [0]*(len(ls) - 1)
    es = [0]*(len(ls) - 1)
    dvs = [0]*(len(ls) - 1)
    t_fb = [0]*(len(ls) - 1)

    for i in range(len(ls) - 1):
        cur_body = seq[i + 1]
        time = encounters[0]

        v_body = np.array(cur_body.eph(time)[1])

        v_in = np.array(ls[i].get_v2()[0])- v_body
        v_out = np.array(ls[i+1].get_v1()[0]) - v_body

        deltas[i] = np.arccos(v_in.dot(v_out)/np.linalg.norm(v_in)/np.linalg.norm(v_out))

        es[i] = 1/np.sin(deltas[i]/2)
        dvs[i] = pk.fb_vel(v_in, v_out, cur_body)

        a_cur = cur_body.mu_self/np.linalg.norm(v_in)**2
        r_ps[i] = (es[i] - 1)*a_cur/cur_body.radius

        theta_start = -.95*np.arccos(-1/es[i])
        H_start = 2*np.arctanh(np.tanh(theta_start/2)/np.sqrt((1 + es[i])/(es[i] - 1)))
        t_fb[i] = -2*(es[i]*np.sinh(H_start) - H_start)/np.sqrt(cur_body.mu_self/a_cur**3)

    return r_ps, np.rad2deg(deltas)[0], es, dvs, t_fb

def calc_plane(NAIF_ID, when):
    DCM = sppy.tipbod("ECLIPJ2000", NAIF_ID, when).T
    return R.from_matrix(DCM).as_euler('ZXZ')


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

    udp = mga_mod(
        seq=seq,
        vinf=vinf_free,
        t0=[start, end],
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
    uda = pg.mbh(algo=pg.scipy_optimize(), perturb=.9, stop=300)
    algo = pg.algorithm(uda)
    algo.set_verbosity(10)

    # We construct a random population of 20 individuals (the initial guess)
    pop = pg.population(prob, size=1)
    # We solve the problem
    pop = algo.evolve(pop)

    dv_struct = udp._compute_dvs(pop.champion_x)
    seg_list = list((5000*np.ones(len(dv_struct[3]))).astype(int))

    epoch_lists, planet_rs = calc_planets(udp, pop, seg_list, [mer, v, e, m, j, s, u])
    t_list = np.array([epoch.mjd2000 for epoch in epoch_lists])
    eph = udp.get_eph_function(pop.champion_x)

    data = np.vstack([eph(t_list[i])[0] for i in range(len(t_list))])

    all_data = np.hstack([np.reshape(t_list, [len(t_list), 1]), data] + planet_rs)

    small_t = np.linspace(all_data[:, 0][0], all_data[:, 0][-1], 500)
    small_data = np.hstack([small_t.reshape(-1, 1), interp1d(all_data[:, 0], all_data[:, 1:].T)(small_t).T])

    planet_string = ["Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus"]
    col_labels = ['time', 'traj_r_x', 'traj_r_y', 'traj_r_z']
    for planet_str in planet_string:
        col_labels.append(planet_str + '_r_x')
        col_labels.append(planet_str + '_r_y')
        col_labels.append(planet_str + '_r_z')

    all_df = pd.DataFrame(all_data, columns=col_labels)
    all_df.to_csv('data.csv', index=False)

    small_df = pd.DataFrame(small_data, columns=col_labels)
    all_df.to_csv('small_data.csv', index=False)

    plot_all(planet_rs, data, planet_string, ['lightgray', 'cornsilk', 'deepskyblue', 'red', 'darkorange', 'lightgoldenrodyellow', 'lightblue'])

    tot_tof = np.sum(np.array(udp._decode_tofs(pop.champion_x)))
    mjd = pop.champion_x[0] + tot_tof

    v_earth_dep = np.array(list(e.eph(pk.epoch(pop.champion_x[0], "mjd2000"))[1]))
    v_rel_earth = np.array(dv_struct[3][0].get_v1()) - v_earth_dep

    v_depart_earth_wrt_earth = sppy.tipbod("ECLIPJ2000", 399, pop.champion_x[0]) @ v_rel_earth.reshape(-1, 1)

    rp_wrt_earth, vp_wrt_earth, v_inf_wrt_earth_earth = bplane_rp_ecc(v_rel_earth, pk.MU_EARTH, pk.EARTH_RADIUS+200*1000, 399, pop.champion_x[0], "out")
    outbound = pk.ic2par(rp_wrt_earth, vp_wrt_earth, mu=pk.MU_EARTH)


    v_uran_arr = np.array(list(u.eph(pk.epoch(mjd, "mjd2000"))[1]))
    v_sc_helio_arr = dv_struct[3][-1].get_v2()

    #v_sc_uran_arr = change_frame(np.array(v_sc_helio_arr), 799, mjd)

    v_rel_uran = v_sc_helio_arr - v_uran_arr
    rp_wrt_uran, vp_wrt_uran, v_inf_wrt_body = bplane_rp_ecc(v_rel_uran, mu_uranus, r_p_mult2*mu_uranus, 799, mjd, "in")

    fpa = -np.deg2rad(8)

    v_inf_probe = np.sqrt(v_inf_wrt_body.dot(v_inf_wrt_body))
    spec_e = v_inf_probe ** 2 / 2
    probe_sma = mu_uranus / 2 / spec_e

    probe_tof = 6*30*pk.DAY2SEC

    def match_entry(r_p, fpa, r_uranus, sma):
        ecc = r_p*r_uranus/sma + 1
        theta_entry = -np.arccos((r_p*r_uranus*(ecc + 1)/r_uranus - 1)/ecc)
        obj = np.arctan(ecc*np.sin(theta_entry)/(1 + np.cos(theta_entry))) - fpa
        return obj

    r_p_probe = root_scalar(match_entry, args=(fpa, r_uranus, probe_sma), bracket=(0.001, 1)).root * r_uranus

    ecc_probe = r_p_probe / probe_sma + 1
    delta_probe = 2 * np.arcsin(1 / ecc_probe)

    theta_atmo = np.arccos((probe_sma * (ecc_probe**2 - 1)/r_uranus - 1)/ecc_probe)
    H_atmo = 2*np.arctanh(np.tan(theta_atmo/2)/(np.sqrt((ecc_probe + 1)/(ecc_probe - 1))))
    tof_atmo2r_p = (np.sinh(H_atmo) - H_atmo)/np.sqrt(mu_uranus/probe_sma**3)

    r_p_orbiter = r_p_mult2 * r_uranus

    def calc_H_release(H, ecc, mu, sma, tof):
        return ecc*np.sinh(H) - H - np.sqrt(mu/sma**3)*tof


    H_release = root_scalar(calc_H_release, args=(ecc_probe, mu_uranus, probe_sma, -probe_tof), bracket=[-50, 50]).root
    r_release = probe_sma*(ecc_probe*np.cosh(H_release) - 1)
    theta_release = 2*np.arctan(np.sqrt((1 + ecc_probe)/(ecc_probe - 1))*np.tanh(H_release/2))

    v_release_probe = np.sqrt(2*(spec_e + mu_uranus/r_release))
    fpa_release = np.arctan(ecc_probe*np.sin(theta_release)/(1 + ecc_probe*np.cos(theta_release)))

    v_probe_entry = np.sqrt(2*(spec_e + mu_uranus/r_uranus))

    def dv_deflection(ecc_orbiter, r_p_mult2, r_uranus, probe_tof, offset, theta_release, fpa_release, v_release_probe, is_opt):
        if is_opt:
            try:
                orbiter_sma = r_p_mult2*r_uranus/(ecc_orbiter - 1)

                spec_e2 = mu_uranus/(2*orbiter_sma)

                H_release2 = root_scalar(calc_H_release, args=(ecc_orbiter, mu_uranus, orbiter_sma, -probe_tof+offset), bracket=[-50, 50]).root
                r_release2 = orbiter_sma*(ecc_orbiter*np.cosh(H_release2) - 1)
                theta_release2 = 2*np.arctan(np.sqrt((1 + ecc_orbiter)/(ecc_orbiter - 1))*np.tanh(H_release2/2))

                v_release_probe2 = np.sqrt(2*(spec_e2 + mu_uranus/r_release2))
                fpa_release2 = np.arctan(ecc_probe*np.sin(theta_release2)/(1 + ecc_probe*np.cos(theta_release2)))

                delta_arg_p = theta_release2 - theta_release
                ang_between = delta_arg_p + fpa_release2 - fpa_release

                dV_orbiter = np.sqrt(v_release_probe2**2 + v_release_probe**2 - 2*v_release_probe*v_release_probe2*np.cos(ang_between))
            except:
                dV_orbiter=1e9
            return dV_orbiter
        else:
            orbiter_sma = r_p_mult2 * r_uranus / (ecc_orbiter - 1)

            spec_e2 = mu_uranus / (2 * orbiter_sma)

            H_release2 = root_scalar(calc_H_release, args=(ecc_orbiter, mu_uranus, orbiter_sma, -probe_tof + offset),
                                     bracket=[-50, 50]).root
            r_release2 = orbiter_sma * (ecc_orbiter * np.cosh(H_release2) - 1)
            theta_release2 = 2 * np.arctan(np.sqrt((1 + ecc_orbiter) / (ecc_orbiter - 1)) * np.tanh(H_release2 / 2))

            v_release_probe2 = np.sqrt(2 * (spec_e2 + mu_uranus / r_release2))
            fpa_release2 = np.arctan(ecc_probe * np.sin(theta_release2) / (1 + ecc_probe * np.cos(theta_release2)))

            delta_arg_p = theta_release2 - theta_release
            ang_between = delta_arg_p + fpa_release2 - fpa_release

            dV_orbiter = np.sqrt(
                v_release_probe2 ** 2 + v_release_probe ** 2 - 2 * v_release_probe * v_release_probe2 * np.cos(
                    ang_between))

            return dV_orbiter, delta_arg_p

    offset = 3600 + tof_atmo2r_p

    opt_ecc = minimize_scalar(dv_deflection, args=(r_p_mult2, r_uranus, probe_tof, offset, theta_release, fpa_release, v_release_probe, True))
    orbiter_sma = r_p_mult2 * r_uranus / (opt_ecc.x - 1)
    dV_orbiter, delta_arg_p = dv_deflection(opt_ecc.x, r_p_mult2, r_uranus, probe_tof, offset, theta_release, fpa_release, v_release_probe, False)

    print(dV_orbiter)

    orb_elem_wrt_uran = np.array(pk.ic2par(rp_wrt_uran, vp_wrt_uran, mu=mu_uranus))
    orb_elem_wrt_uran[3] -= delta_arg_p

    dV_from_earth = calc_dV_LEO(200, dv_struct[3][0], e, pop)

    deep_space_dv = np.sum(dv_struct[1])

    dV_capture = dv_struct[2]

    point_planet = False
    r_uranus = 25_362_000
    mu_uranus = 5.793939e15

    r_miranda = 471.6*1000/2
    r_ariel = 1157.8*1000/2
    r_umbriel = 1169.4*1000/2
    r_titania = 1576.8*1000/2
    r_oberon = 1522.8*1000/2

    m_miranda = 6400E16
    m_ariel = 125100E16
    m_umbriel = 127500E16
    m_titania = 340000E16
    m_oberon = 307600E16

    G = 6.6743E-11

    mu_miranda = G*m_miranda
    mu_ariel = G*m_ariel
    mu_umbriel = G*m_umbriel
    mu_titania = G*m_titania
    mu_oberon = G*m_oberon

    safe_mult = 1.1

    if point_planet:
        m = spice('MIRANDA', 'URANUS', 'IAU_URANUS', 'NONE', mu_uranus, mu_miranda, 10.0, 10.0)
        a = spice('ARIEL', 'URANUS', 'IAU_URANUS', 'NONE', mu_uranus, mu_ariel, 10.0, 10.0)
        u = spice('UMBRIEL', 'URANUS', 'IAU_URANUS', 'NONE', mu_uranus, mu_umbriel, 10.0, 10.0)
        t = spice('TITANIA', 'URANUS', 'IAU_URANUS', 'NONE', mu_uranus, mu_titania, 10.0, 10.0)
        o = spice('OBERON', 'URANUS', 'IAU_URANUS', 'NONE', mu_uranus, mu_oberon, 10.0, 10.0)

    else:
        m = spice('MIRANDA', 'URANUS', 'IAU_URANUS', 'NONE', mu_uranus, mu_miranda, r_miranda, r_miranda*safe_mult)
        a = spice('ARIEL', 'URANUS', 'IAU_URANUS', 'NONE', mu_uranus, mu_ariel, r_ariel, r_ariel*safe_mult)
        u = spice('UMBRIEL', 'URANUS', 'IAU_URANUS', 'NONE', mu_uranus, mu_umbriel, r_umbriel, r_umbriel*safe_mult)
        t = spice('TITANIA', 'URANUS', 'IAU_URANUS', 'NONE', mu_uranus, mu_titania, r_titania, r_titania*safe_mult)
        o = spice('OBERON', 'URANUS', 'IAU_URANUS', 'NONE', mu_uranus, mu_oberon, r_oberon, r_oberon*safe_mult)

        URANUS = spice('URANUS', 'SUN', 'ECLIPJ2000', 'NONE', pk.MU_SUN, mu_uranus, r_uranus, 4*r_uranus)
    r_uranus = 25_362_000
    mu_uranus = 5.793939e15

    # init_orbit =

    seq2 = [t, a, o, u]

    tof_moon_tour = 3.5

    science_orbit_time = 1
    start_new_mission = .5

    start_end = [mjd + science_orbit_time*365.25, mjd + (start_new_mission + science_orbit_time) * 365.25]
    start2 = pk.epoch(start_end[0])
    end2 = pk.epoch(start_end[1])

    cur_orb = pk.planet.keplerian(start2, [sma, ecc2, orb_elem_wrt_uran[2], orb_elem_wrt_uran[3], orb_elem_wrt_uran[4], orb_elem_wrt_uran[5]],
                                  mu_uranus, 0)

    #fin_orb = pk.planet.keplerian(start2, [sma, ecc2, orb_elem_wrt_uran[2], orb_elem_wrt_uran[3], orb_elem_wrt_uran[4], orb_elem_wrt_uran[5]],
    #                              mu_uranus, 0)

    seq2.insert(0, cur_orb)
    #udp2 = mga_mod_minr(
    #    seq=seq2,
    #    vinf=0,
    #    t0=[start2, end2],
    #    tof=tof_moon_tour * 365.25,
    #    orbit_insertion=False,
    #    tof_encoding='eta',
    #    max_revs=20,
    #    dv_launch_weight=1,
    #    dv_arr_weight=0
    #)

    udp2 = mga_1dsm_mod(
        seq=seq2,
        vinf=[0, .00001],
        t0=[start2, end2],
        tof=tof_moon_tour * 365.25,
        orbit_insertion=False,
        tof_encoding='eta',
        max_revs=20,
        add_vinf_arr=False
    )

    prob2 = pg.problem(udp2)
    # We solve it!!
    # uda = pg.xnes(gen=200, force_bounds=True, xtol=1e-12)
    uda2 = pg.mbh(algo=pg.scipy_optimize(), perturb=.9, stop=300)
    algo2 = pg.algorithm(uda2)
    algo2.set_verbosity(1)

    # We construct a random population of 20 individuals (the initial guess)dv
    pop2 = pg.population(prob2, size=1)
    # We solve the problem
    pop2 = algo2.evolve(pop2)

    dv_struct2 = udp2._compute_dvs(pop2.champion_x)
    seg_list2 = list((5000 * np.ones(len(dv_struct2[1]))).astype(int))

    # pd.DataFrame(data).to_csv('test.csv')

    epoch_lists2, planet_rs2 = calc_planets(udp2, pop2, seg_list2, [cur_orb, m, a, u, t, o])
    t_list2 = np.array([epoch.mjd2000 for epoch in epoch_lists2])
    eph2 = udp2.get_eph_function(pop2.champion_x)

    data2 = np.vstack([eph2(t_list2[i])[0] for i in range(len(t_list2))])

    all_data2 = np.hstack([np.reshape(t_list2, [len(t_list2), 1]), data2] + planet_rs2)

    small_t2 = np.linspace(all_data2[:, 0][0], all_data2[:, 0][-1], 100)
    small_data2 = np.hstack([small_t2.reshape(-1, 1), interp1d(all_data2[:, 0], all_data2[:, 1:].T)(small_t2).T])

    planet_string2 = ['Start', 'Miranda', 'Ariel', 'Umbriel', 'Titania', 'Oberon']
    col_labels2 = ['time', 'traj_r_x', 'traj_r_y', 'traj_r_z']
    for planet_str in planet_string2:
        col_labels2.append(planet_str + '_r_x')
        col_labels2.append(planet_str + '_r_y')
        col_labels2.append(planet_str + '_r_z')

    all_df2 = pd.DataFrame(all_data2, columns=col_labels2)
    all_df2.to_csv('data_moon_tour.csv', index=False)

    plot_all(planet_rs2, data2, planet_string2,
             ['lightblue', 'cornsilk', 'deepskyblue', 'red', 'darkorange', 'green'])

    small_df2 = pd.DataFrame(small_data2, columns=col_labels2)
    small_df2.to_csv('small_data_moon_tour.csv', index=False)

    tot_tof2 = np.sum(np.array(udp2._decode_tofs(pop2.champion_x)))
    mjd2 = pop2.champion_x[0] + tot_tof2

    planet_encounters2 = pop2.champion_x[0] + np.cumsum(udp.eta2direct(pop2.champion_x))

    minpass = min(np.linalg.norm(all_data2[:, 1:4], axis=1)) / r_uranus

    if minpass < 1:
        print("ERROR: ORBITER CRASHES INTO URANUS")

    final_v_pre = dv_struct2[1][-1].get_v2()[0]
    final_v = pk.fb_prop(np.array(final_v_pre),
               np.array(seq2[-1].eph(pop2.champion_x[0] + np.sum(udp2._decode_tofs(pop2.champion_x)))[1]),
               1.1 * r_umbriel, 0, mu_umbriel)
    final_r = dv_struct2[1][-1].get_r2()


    fin_orbit_data = pk.ic2par(final_r, final_v, mu_uranus)
    fin_sma = fin_orbit_data[0]
    fin_ecc = fin_orbit_data[1]
    fin_rest = np.rad2deg(fin_orbit_data[2:])
    r_p_fin = fin_orbit_data[0]*(1 - fin_orbit_data[1])
    r_p_fin_uranus = r_p_fin/r_uranus
    T_fin = 2*np.pi*np.sqrt(fin_orbit_data[0]**3/mu_uranus)/pk.DAY2SEC

    probe_calc = True


#    if probe_calc:


        #plt.plot(x_probe, y_probe)
        #plt.plot(x_orbiter, y_orbiter)
        #plt.plot(r_uranus*np.cos(np.linspace(0, 2*np.pi)), r_uranus*np.sin(np.linspace(0, 2*np.pi)))
        #plt.axis("equal")
        #plt.show()


