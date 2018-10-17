import sys
import numpy as np
import matplotlib.patches as patches
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
'''
usage: python wmap1.py <pmesh> <loc>
<pmesh> is the filename of pmesh file
<loc> is either 'bot' or 'top'
'''


def kart2frac(kart, latmat):
    """
    convert cart coords into frac
    :param kart: a list of cart coords
    :param latmat: [va, vb, vc] in cart
    :return: np array [a, b, c]: frac coords
    """
    p = np.matmul(np.array(kart), np.linalg.inv(np.array(latmat)))
    return p


def xyzarray2frac(x, y, z, latmat):
    """
    convert frac into cart
    :param x:
    :param y:
    :param z:
    :param latmat: [va, vb, vc] in cart
    :return: nx3 mat
    """
    length = min([len(x), len(y), len(z)])
    abc = np.empty((length, 3))
    abc[:] = np.nan
    for i in range(length):
        kart = kart2frac([x[i], y[i], z[i]], latmat)
        abc[i][0] = kart[0]
        abc[i][1] = kart[1]
        abc[i][2] = kart[2]
    return abc


def parse_cif(cif_name='iso.cif'):
    """
    parse lattice vectors from 'iso.cif'
    :return: la 1x3 np array
             lb 1x3 np array
             lc 1x3 np array
             theta_c_rad angle between c axis and z axis
    """
    with open(cif_name) as f_iso:
        content = f_iso.readlines()
    u = np.zeros(6)
    for e in [line.strip().split() for line in content if len(line.strip().split()) == 2]:
        if 'cell_length_a' in e[0]:
            u[0] = float(e[1])
        elif 'cell_length_b' in e[0]:
            u[1] = float(e[1])
        elif 'cell_length_c' in e[0]:
            u[2] = float(e[1])
        elif 'cell_angle_alpha' in e[0]:
            u[3] = float(e[1])
        elif 'cell_angle_beta' in e[0]:
            u[4] = float(e[1])
        elif 'cell_angle_gamma' in e[0]:
            u[5] = float(e[1])
    a, b, c, alpha, beta, gamma = u
    cosdelta_up = np.cos(np.radians(alpha)) - np.cos(np.radians(beta))*np.cos(np.radians(gamma))
    cosdelta_low = np.sin(np.radians(beta))*np.sin(np.radians(gamma))
    cosdelta = cosdelta_up / cosdelta_low
    sindelta = np.sqrt(1-cosdelta**2)
    la = a*np.array([1.0, 0.0, 0.0])
    lb = b*np.array([np.cos(np.radians(gamma)), np.sin(np.radians(gamma)), 0.0])
    lc = c*np.array([np.cos(np.radians(beta)), np.sin(np.radians(beta))*cosdelta,
                     np.sin(np.radians(beta))*sindelta])
    u_lc = lc/np.linalg.norm(lc)
    theta_c_rad = np.arccos(np.clip(np.dot(u_lc, [0, 0, 1]), -1.0, 1.0))
    return la, lb, lc, theta_c_rad


def parse_pmesh(pmesh_name, loc, latmat):
    """

    :param pmesh_name: name of the pmesh file
    :param loc: whether the surface is bot or top
    :param latmat: [la, lb, lc]
    :return: pts nx2 array for x, y coord
             nz normalized z for color
    """
    with open(pmesh_name) as f_pmesh:
        content = f_pmesh.readlines()
    content = [i.strip('\n') for i in content]
    data = []
    for line in content:
        l = list(line.split())
        if len(l) != 1:
            data.append(np.array([float(s) for s in l]))
    sorteddata = sorted(data, key=lambda k: [k[1], k[0]])
    pmesh_data = np.array(sorteddata)
    x = pmesh_data[:, 0]
    y = pmesh_data[:, 1]
    z = pmesh_data[:, 2]
    abc = xyzarray2frac(x, y, z, latmat)
    nc = abc[:, 2]
    if loc == 'bot':
        nz = (nc - min(nc)) / (max(nc) - min(nc))
    elif loc == 'top':
        nz = (max(nc) - nc) / (max(nc) - min(nc))
    else:
        print('wrong loc!')
        sys.exit(1)
    pts = np.column_stack((x, y))
    return pts, nz


def project_pmesh(pts, nz, la, lb):
    x = pts[:, 0]
    y = pts[:, 1]
    boundaries = (min(x), max(x), min(y), max(y))
    grid_x, grid_y = np.mgrid[min(x):max(x):500j, min(y):max(y):500j]
    grid_z = griddata(pts, nz, (grid_x, grid_y), method='cubic')
    fig, ax = plt.subplots()
    ax.imshow(grid_z.T, extent=boundaries, aspect='equal', vmin=0, vmax=1, cmap='bwr', origin='lower')

    # define x y of the left-bot vertex
    anchor_x = (boundaries[0] + boundaries[1]) / 5.0  # parallelogram
    anchor_y = (boundaries[2] + boundaries[3]) / 5.0  # parallelogram

    xp = [anchor_x, anchor_x + la[0], anchor_x + la[0] + lb[0], anchor_x + lb[0]]
    yp = [anchor_y, anchor_y + la[1], anchor_y + la[1] + lb[1], anchor_y + lb[1]]
    ax.add_patch(patches.Polygon(xy=list(zip(xp, yp)), fill=False, linewidth=2))
    ax.set_ylabel(r'Y ($\mathrm{\AA}$)')
    ax.set_xlabel(r'X ($\mathrm{\AA}$)')
    fig.tight_layout()
    fig.savefig('wmap1.png', dpi=200)


Pmesh_File, Location = sys.argv[1], sys.argv[2]
La, Lb, Lc, Theta = parse_cif()
Pts, Nz = parse_pmesh(Pmesh_File, Location, [La, Lb, Lc])
project_pmesh(Pts, Nz, La, Lb)
