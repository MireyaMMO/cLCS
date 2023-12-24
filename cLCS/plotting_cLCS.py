import numpy as np
import pickle
import matplotlib.pyplot as plt
import matplotlib
import scipy
import matplotlib.cm as cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from scipy.interpolate import griddata
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cLCS.utils import *


def colormap():
  # Colormap used by Rodrigo Duran cLCSs
    s = np.array([[1.0000,   1.0000,   0.9987, 1],
                  [0.9971,   1.0000,   0.9970, 1],
                  [0.9896,   1.0000,   0.9931, 1],
                  [0.9771,   1.0000,   0.9871, 1],
                  [0.9593,   0.9900,   0.9789, 1],
                  [0.9364,   0.9708,   0.9686, 1],
                  [0.9084,   0.9484,   0.9564, 1],
                  [0.8759,   0.9243,   0.9422, 1],
                  [0.8395,   0.8994,   0.9264, 1],
                  [0.8000,   0.8749,   0.9092, 1],
                  [0.7585,   0.8516,   0.8906, 1],
                  [0.7160,   0.8301,   0.8710, 1],
                  [0.6738,   0.8110,   0.8506, 1],
                  [0.6330,   0.7948,   0.8296, 1],
                  [0.5949,   0.7817,   0.8081, 1],
                  [0.5606,   0.7719,   0.7865, 1],
                  [0.5310,   0.7657,   0.7649, 1],
                  [0.5073,   0.7628,   0.7435, 1],
                  [0.4900,   0.7633,   0.7225, 1],
                  [0.4798,   0.7671,   0.7019, 1],
                  [0.4771,   0.7737,   0.6821, 1],
                  [0.4819,   0.7831,   0.6629, 1],
                  [0.4943,   0.7946,   0.6446, 1],
                  [0.5138,   0.8081,   0.6271, 1],
                  [0.5399,   0.8230,   0.6105, 1],
                  [0.5720,   0.8387,   0.5948, 1],
                  [0.6090,   0.8548,   0.5800, 1],
                  [0.6500,   0.8708,   0.5659, 1],
                  [0.6938,   0.8860,   0.5525, 1],
                  [0.7391,   0.9000,   0.5398, 1],
                  [0.7847,   0.9120,   0.5275, 1],
                  [0.8292,   0.9217,   0.5155, 1],
                  [0.8716,   0.9284,   0.5037, 1],
                  [0.9108,   0.9317,   0.4918, 1],
                  [0.9457,   0.9310,   0.4797, 1],
                  [0.9756,   0.9260,   0.4672, 1],
                  [1.0000,   0.9162,   0.4541, 1],
                  [1.0000,   0.9013,   0.4401, 1],
                  [1.0000,   0.8810,   0.4251, 1],
                  [1.0000,   0.8551,   0.4089, 1],
                  [1.0000,   0.8235,   0.3912, 1],
                  [1.0000,   0.7862,   0.3720, 1],
                  [1.0000,   0.7432,   0.3511, 1],
                  [1.0000,   0.6947,   0.3284, 1],
                  [1.0000,   0.6408,   0.3039, 1],
                  [1.0000,   0.5821,   0.2775, 1],
                  [0.9900,   0.5190,   0.2494, 1],
                  [0.9819,   0.4521,   0.2195, 1],
                  [0.9765,   0.3822,   0.1882, 1],
                  [0.9744,   0.3102,   0.1556, 1],
                  [0.9756,   0.2372,   0.1222, 1],
                  [0.9799,   0.1643,   0.0884, 1],
                  [0.9864,   0.0931,   0.0547, 1],
                  [0.9938,   0.0251,   0.0219, 1],
                  [1.0000,        0,        0, 1],
                  [1.0000,        0,        0, 1],
                  [0.9989,        0,        0, 1],
                  [0.9858,        0,        0, 1],
                  [0.9601,        0,        0, 1],
                  [0.9194,        0,        0, 1],
                  [0.8618,        0,        0, 1],
                  [0.7874,        0,        0, 1],
                  [0.6982,        0,        0, 1],
                  [0.6000,   0.0069,    0.0013, 1]])
    newcmap = LinearSegmentedColormap.from_list('mycmap', s)
    return newcmap


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
