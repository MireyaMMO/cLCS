import numpy as np
import pickle
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cLCS.utils import *

def plot_colourline(x, y, c, cmap, ax=None, transform=None):
  # Plots LCSs using colouredlines to define intensity
    c = cmap((c-.4)/(1.4-.4))
    if ax == None:
        ax = plt.gca()
    for i in np.arange(len(x)-1):
        ax.plot([x[i], x[i+1]], [y[i], y[i+1]],
                c=c[i], transform=transform, lw=.8)
    return

def plot_blacklines(x, y, ax=None, transform=None):
  # Plots LCSs as black lines
    if ax == None:
        ax = plt.gca()
    ax.plot(x, y, 'k', transform=transform, lw=.8)
    return


def cLCSrho_cartopy(dirr, monthvec, fig=None, ax=None, projection=None):
  # Plots cLCS using cartopy
  # check projection line as the one used is for New Zealand due to the +-180
    if not projection:
        projection = ccrs.PlateCarree()
    if not fig:
        fig = plt.figure(figsize=(8, 6), constrained_layout=True)
    if not ax:
        ax = fig.add_subplot(projection=projection)
    f = cfeature.GSHHSFeature(scale='high', levels=[1])
    ax.add_geometries(
        f.geometries(),
        ccrs.PlateCarree(),
        facecolor=cfeature.COLORS['land'],
        edgecolor='black')
    m = '%02d' % monthvec
    month_dirr = dirr+m+'/'
    lon, lat, _, sqrtlda2total, _, _, _, _, _, xspan, yspan, count = pickle.load(
        open(f'{month_dirr}TOT-{m}.p', 'rb'))
    N = count
    # because the landmask was made NaN and the coast is on the left side
    xi = xspan[-1, :]
    yi = yspan[:, -1]
    X0, Y0 = np.meshgrid(xi, yi)
    X0 = X0.ravel()
    Y0 = Y0.ravel()
    pxt, pyt = pickle.load(open(f'{month_dirr}/cLCS_{m}.p', 'rb'))
    lonmin = lon.min()
    latmin = lat.min()
    lonmax = lon.max()
    latmax = lat.max()
    corners = [lonmin, lonmax, latmin, latmax]
    ax.set_extent(corners, crs=ccrs.PlateCarree())
    z = np.log(sqrtlda2total/N)
    nLCS = pxt.shape[0]
    cmap = colormap()
    for kk in range(0, nLCS, 4):  # Change 4 for more or less LCSs
        xs = pxt[kk, :]
        ys = pyt[kk, :]
        xs[np.where(xs < 0)] = 0
        ys[np.where(ys < 0)] = 0
        z[np.where(np.isnan(z))] = 0
        zs = griddata((X0.ravel(), Y0.ravel()),
                      z.ravel(), (xs.ravel(), ys.ravel()))
        [xS, yS] = xy2sph(xs*1e3, lonmin, ys*1e3, latmin)
        zs[np.where(np.isnan(zs))] = 0
        plot_colourline(xS, yS, zs, cmap, transform=ccrs.PlateCarree())
