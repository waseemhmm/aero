import os, re, glob, sys
import numpy as np, pyvista as pv
from scipy.io import savemat
files = sorted(glob.glob('flow_*.vtu'),
               key=lambda f: int(re.search(r'flow_(\d+)\.vtu', f).group(1)))
print(f'found {len(files)} snapshots')
first = pv.read(files[0])
coords = np.asarray(first.points)[:, :2]
N = coords.shape[0]; M = len(files)
rho = np.zeros((N, M)); u = np.zeros((N, M)); v = np.zeros((N, M)); p = np.zeros((N, M))
t = np.zeros(M)
for j, f in enumerate(files):
    m = pv.read(f); pd = m.point_data
    rho[:, j] = np.asarray(pd['Density'])
    if 'Velocity' in pd.keys():
        vel = np.asarray(pd['Velocity'])
        u[:, j] = vel[:, 0]; v[:, j] = vel[:, 1]
    else:
        mom = np.asarray(pd['Momentum'])
        u[:, j] = mom[:, 0] / rho[:, j]
        v[:, j] = mom[:, 1] / rho[:, j]
    p[:, j] = np.asarray(pd['Pressure'])
    t[j] = float(re.search(r'flow_(\d+)\.vtu', f).group(1))
savemat('snaps.mat', dict(rho=rho, u=u, v=v, p=p, t=t.reshape(-1, 1), coords=coords),
        do_compression=True)
print(f'wrote snaps.mat: N={N} M={M}')
