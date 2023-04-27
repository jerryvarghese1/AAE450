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

load_spice_kernel("de721_full.bsp")
load_spice_kernel("../current/ura111.bsp")
load_spice_kernel("pck00010.tpc")

e = spice('EARTH', 'SUN', 'ECLIPJ2000', 'NONE', pk.MU_SUN, 3.986e14, 10.0, 10.0)

mer = spice('MERCURY', 'SUN', 'ECLIPJ2000', 'NONE', pk.MU_SUN, 2.2032e13, 10.0, 10.0)
m = spice('MARS', 'SUN', 'ECLIPJ2000', 'NONE', pk.MU_SUN, 4.282e13, 10.0, 10.0)
v = spice('VENUS', 'SUN', 'ECLIPJ2000', 'NONE', pk.MU_SUN, 3.24859e14, 10.0, 10.0)
j = spice('JUPITER BARYCENTER', 'SUN', 'ECLIPJ2000', 'NONE', pk.MU_SUN, 1.266865e17, 10.0, 10.0)
s = spice('SATURN BARYCENTER', 'SUN', 'ECLIPJ2000', 'NONE', pk.MU_SUN, 3.793e16, 10.0, 10.0)

u = spice('URANUS BARYCENTER', 'SUN', 'ECLIPJ2000', 'NONE', pk.MU_SUN, 5.793e15, 10.0, 10.0)

data = pd.read_csv("interplanetary_v5.csv")

time = data['time'].to_numpy()

traj = data[['traj_r_x', 'traj_r_y', 'traj_r_z']].to_numpy()

sun_occlusion = np.zeros(len(time))
mer_occlusion = np.zeros(len(time))
m_occlusion = np.zeros(len(time))
v_occlusion = np.zeros(len(time))
j_occlusion = np.zeros(len(time))
s_occlusion = np.zeros(len(time))
u_occlusion = np.zeros(len(time))

for i in range(len(time)):
    pos_vec = traj[i, :].flatten()
    cur_earth = np.array(e.eph(time[i])[0]).flatten()

    sun_occlusion[i] = np.rad2deg(np.arccos(np.dot(pos_vec/np.linalg.norm(pos_vec), cur_earth/np.linalg.norm(cur_earth))))

    cur_mer = np.array(mer.eph(time[i])[0]).flatten()
    vec = cur_mer - pos_vec
    vec2 = cur_mer - cur_earth

    mer_occlusion[i] = np.rad2deg(np.arccos(np.dot(vec/np.linalg.norm(vec), vec2/np.linalg.norm(vec2))))

    cur_m = np.array(m.eph(time[i])[0]).flatten()
    vec = cur_m - pos_vec
    vec2 = cur_m - cur_earth

    m_occlusion[i] = np.rad2deg(
        np.arccos(np.dot(vec / np.linalg.norm(vec), vec2/np.linalg.norm(vec2))))

    cur_j = np.array(j.eph(time[i])[0]).flatten()
    vec = cur_j - pos_vec
    vec2 = cur_j - cur_earth

    j_occlusion[i] = np.rad2deg(
        np.arccos(np.dot(vec / np.linalg.norm(vec), vec2/np.linalg.norm(vec2))))

    cur_s = np.array(s.eph(time[i])[0]).flatten()
    vec = cur_s - pos_vec
    vec2 = cur_s - cur_earth

    s_occlusion[i] = np.rad2deg(
        np.arccos(np.dot(vec / np.linalg.norm(vec), vec2/np.linalg.norm(vec2))))\

    cur_v = np.array(v.eph(time[i])[0]).flatten()
    vec = cur_v - pos_vec
    vec2 = cur_v - cur_earth

    v_occlusion[i] = np.rad2deg(
        np.arccos(np.dot(vec / np.linalg.norm(vec), vec2 / np.linalg.norm(vec2))))

    cur_u = np.array(u.eph(time[i])[0]).flatten()
    vec = cur_u - pos_vec
    vec2 = cur_u - cur_earth

    u_occlusion[i] = np.rad2deg(
        np.arccos(np.dot(vec / np.linalg.norm(vec), vec2/np.linalg.norm(vec2))))

lab_list = np.array(['Sun', 'Uranus', 'Mercury', 'Mars', 'Venus', 'Jupiter', 'Saturn'])

all_occlusions = np.column_stack([sun_occlusion, u_occlusion, mer_occlusion, m_occlusion, v_occlusion, j_occlusion, s_occlusion])
all_occlusions[np.isnan(all_occlusions)] = 0

interper = interp1d(time, all_occlusions, axis=0)
time2 = np.linspace(time[0], time[-1], len(time)*10)

new_data = interper(time2)
inds = np.where(new_data > 179.8)[0]
planet_inds = np.where(new_data > 179.8)[1]
occ = time2[inds]

tmp1 = np.diff(inds) > 1
seg_inds = np.where(np.diff(inds) > 1)[0] + 1
num_occ = len(tmp1[tmp1])
occ_inds_split = np.split(inds, seg_inds)

occ_planets = lab_list[planet_inds]

for i in range(all_occlusions.shape[1]):
    plt.plot(time2, new_data[:, i], label=lab_list[i])

for i in range(num_occ):
    cur_occ = occ_inds_split[i]
    width = time2[cur_occ[-1]] - time2[cur_occ[0]]
    dates = [pk.epoch(time2[cur_occ[-1]]), pk.epoch(time2[cur_occ[0]])]
    print('Start: %s    End:  %s   Planet: %s'%(dates[1], dates[0], occ_planets[i]))
    rect = matplotlib.patches.Rectangle((time2[cur_occ[0]], 0), width, 180,
                              fill=True,
                              color="purple",
                              linewidth=2)
    plt.gca().add_patch(rect)
plt.xlabel("epoch (mjd2000)")
plt.ylabel("Incidence Angle (deg)")
plt.legend()
plt.show()