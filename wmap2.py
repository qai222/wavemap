import sys
import numpy as np
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
from scipy.interpolate import griddata
'''
usage: python wmap2.py <pmesh_top> <pmesh_bot>
<pmesh_top> is the filename of pmesh file describing the surface of top layer (larger c)
<pmesh_bot> is the filename of pmesh file describing the surface of bot layer (smaller c)
'''

plt.rcParams.update({'font.size': 14})  # use larger font size


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
    # na = abc[:, 0]
    # nb = abc[:, 1]
    nc = abc[:, 2]
    if loc == 'bot':
        nz = (nc - min(nc)) / (max(nc) - min(nc))
    elif loc == 'top':
        nz = (max(nc) - nc) / (max(nc) - min(nc))
    else:
        print('wrong loc!')
        sys.exit(1)
    # pts = np.column_stack((na, nb))
    pts = np.column_stack((x, y))
    return pts, nz


def project_pmesh2(pts_t, pts_b, nz_t, nz_b, la, lb):
    x_t = pts_t[:, 0]
    y_t = pts_t[:, 1]
    bound_t = (min(x_t), max(x_t), min(y_t), max(y_t))

    grid_x, grid_y = np.mgrid[min(x_t):max(x_t):500j, min(y_t):max(y_t):500j]
    grid_z_t = griddata(pts_t, nz_t, (grid_x, grid_y), method='cubic')
    grid_z_b = griddata(pts_b, nz_b, (grid_x, grid_y), method='cubic')
    grid_zsum = (grid_z_t + grid_z_b)

    gs1 = gridspec.GridSpec(ncols=2, nrows=1)
    gs1.update(bottom=0.45, top=0.95,)  # left= 0.00, right = 1.00)
    ax1 = plt.subplot(gs1[0, 0])
    ax2 = plt.subplot(gs1[0, 1])

    gs2 = gridspec.GridSpec(ncols=1, nrows=1)
    gs2.update(top=0.30, bottom=0.10, left=0.15, right=0.6)
    ax4 = plt.subplot(gs2[:])
    ax1.imshow(grid_z_t.T, extent=bound_t, aspect='equal', vmin=0, vmax=1, cmap='bwr')
    ax2.imshow(grid_z_b.T, extent=bound_t, aspect='equal', vmin=0, vmax=1, cmap='bwr')

    anchor_x = (bound_t[0] + bound_t[1]) / 6.0  # parallelogram
    anchor_y = (bound_t[2] + bound_t[3]) / 6.0  # parallelogram
    xp = [anchor_x, anchor_x + la[0], anchor_x + la[0] + lb[0], anchor_x + lb[0]]
    yp = [anchor_y, anchor_y + la[1], anchor_y + la[1] + lb[1], anchor_y + lb[1]]
    for ax in [ax1, ax2]:
        ax.add_patch(patches.Polygon(xy=list(zip(xp, yp)), fill=False, linewidth=2))

    ax1.set_title(r'$Layer$ 2', )
    ax1.set_yticks([0, 5, 10])
    ax1.set_xticks([0, 5, 10])
    ax2.set_yticks([0, 5, 10])
    ax2.set_xticks([0, 5, 10])
    ax2.set_title(r'$Layer$ 1')
    ax1.set_ylabel(r'Y ($\mathrm{\AA}$)')
    ax1.set_xlabel(r'X ($\mathrm{\AA}$)')
    ax2.set_xlabel(r'X ($\mathrm{\AA}$)')

    poly = Polygon(list(zip(xp, yp)))
    grid_zsum_flat = np.ravel(grid_zsum)
    zflat = []
    all_points = np.column_stack((grid_x.ravel(), grid_y.ravel()))
    for i in range(len(all_points)):
        if poly.contains(Point(all_points[i])):
            zflat.append(grid_zsum_flat[i])
    ax4.hist(zflat, 300, alpha=0.75, facecolor='blue')

    ax4.set_ylim([0.0, 1500])
    ax4.set_xlim([0.0, 2.0])
    ax4.set_ylabel(r'$counts$')
    print('sigma', np.std(zflat))
    print('mean', np.mean(zflat))
    plt.savefig('wmap2.png', dpi=200)


Pmesh_File_Top, Pmesh_File_Bot = sys.argv[1], sys.argv[2]
La, Lb, Lc, Theta = parse_cif()
Pts_t, Nz_t = parse_pmesh(Pmesh_File_Top, 'top', [La, Lb, Lc])
Pts_b, Nz_b = parse_pmesh(Pmesh_File_Bot, 'bot', [La, Lb, Lc])
project_pmesh2(Pts_t, Pts_b, Nz_t, Nz_b, La, Lb)
