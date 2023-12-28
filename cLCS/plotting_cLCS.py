import numpy as np
import pickle
import matplotlib.pyplot as plt
from scipy.interpolate import griddata, LinearNDInterpolator
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cLCS.utils import *
import os
import logging


logging.basicConfig(level=logging.INFO)
logger = logging


def plot_colourline(x, y, c, cmap, ax=None, lw=0.8, transform=None):
    # Plots LCSs using coloured lines to define strength of attraction
    c = cmap((c - 0.4) / (1.4 - 0.4))
    if ax == None:
        ax = plt.gca()
    for i in np.arange(len(x) - 1):
        ax.plot([x[i], x[i + 1]], [y[i], y[i + 1]], c=c[i], lw=lw, transform=transform)
    return


def plot_lines(x, y, ax=None, color="k", alpha="1", lw=0.8, transform=None):
    # Plots LCSs as black lines
    if ax == None:
        ax = plt.gca()
    ax.plot(x, y, color=color, alpha=alpha, lw=lw, transform=transform)
    return


def cLCSrho_cartopy_colour(
    dirr,
    monthvec,
    colourmap=None,
    lw=0.8,
    fig=None,
    ax=None,
    projection=None,
    line_spacing=4,
    save_fig=False,
    corners=None,
):
    print(f"---- Generating figure with cartopy features")
    if not projection:
        projection = ccrs.PlateCarree()
    if not fig:
        fig = plt.figure(figsize=(8, 6), constrained_layout=True)
    if not ax:
        ax = fig.add_subplot(projection=projection)
    f = cfeature.GSHHSFeature(scale="high", levels=[1])
    ax.add_geometries(
        f.geometries(),
        ccrs.PlateCarree(),
        facecolor=cfeature.COLORS["land"],
        edgecolor="black",
    )
    print(f"---- High-resolution Cartopy features added")
    m = "%02d" % monthvec
    month_dirr = os.path.join(dirr, m)
    TOT_CG_path = os.path.join(month_dirr, f"TOT-{m}.p")
    lon, lat, _, sqrtlda2total, _, _, _, _, _, _, _ , count = pickle.load(
        open(TOT_CG_path, "rb")
    )
    N = count
    pxt, pyt = pickle.load(open(f"{month_dirr}/cLCS_{m}.p", "rb"))
    pxt[np.where(pxt < 0)] = np.nan
    pyt[np.where(pyt < 0)] = np.nan
    z = np.log(sqrtlda2total / N)
    #z[np.where(np.isnan(z))] = 0
    lonmin = lon.min()
    latmin = lat.min()
    lonmax = lon.max()
    latmax = lat.max()
    #LON0, LAT0, _ = projection.transform_points(ccrs.Geodetic(), lon, lat).T
    #z[np.where(np.isnan(z))] = 0
    nonanindex = np.invert(np.isnan(z))
    interpolator = LinearNDInterpolator(
        list(zip(lon[nonanindex].ravel(), lat[nonanindex].ravel())), z[nonanindex].ravel()
    )
    Plon, Plat = xy2sph(pxt*1e3, lonmin, pyt*1e3, latmin)
    #out_xyz = projection.transform_points(ccrs.Geodetic(), Plon, Plat)
    #out_lon = out_xyz[:, :, 0]
    #out_lat = out_xyz[:, :, 1]
    nLCS = pxt.shape[0]
    if not corners:
        corners = [lonmin, lonmax, latmin, latmax]
    ax.set_extent(corners, crs=ccrs.PlateCarree())

    try:
        cmap = get_colourmap(colourmap)
    except:
        logger.warn(
            "colourmap is not defined on utils trying to obtain colourmap from matplotlib"
        )
        cmap = plt.get_cmap(colourmap)
    print(f"---- Squeezeline and associated data loaded")

    for kk in range(0, nLCS, line_spacing): 
        xs = Plon[kk, :]
        ys = Plat[kk, :]
        zs = interpolator(xs, ys)
        #zs[np.where(np.isnan(zs))] = 0
        plot_colourline(xs, ys, zs, cmap, ax=ax, lw=lw, transform=ccrs.PlateCarree())
    if save_fig:
        print(f"---- Saving Figure")
        fig.savefig(f"{month_dirr}/cLCS_{m}.{save_fig}")
    print(f"---- Done")
    return fig, ax


def cLCSrho_cartopy_monochrome(
    dirr,
    monthvec,
    fig=None,
    ax=None,
    color="k",
    alpha="1",
    lw=0.8,
    projection=None,
    line_spacing=1,
    save_fig=False,
    corners=None,
):
    print(f"---- Generating figure with cartopy features")
    if not projection:
        projection = ccrs.PlateCarree()
    if not fig:
        fig = plt.figure(figsize=(8, 6), constrained_layout=True)
    if not ax:
        ax = fig.add_subplot(projection=projection)
    f = cfeature.GSHHSFeature(scale="high", levels=[1])
    ax.add_geometries(
        f.geometries(),
        ccrs.PlateCarree(),
        facecolor=cfeature.COLORS["land"],
        edgecolor="black",
    )
    print(f"---- High-resolution Cartopy features added")
    m = "%02d" % monthvec
    month_dirr = os.path.join(dirr, m)
    TOT_CG_path = os.path.join(month_dirr, f"TOT-{m}.p")
    lon, lat, _, _, _, _, _, _, _, _, _, _ = pickle.load(open(TOT_CG_path, "rb"))
    pxt, pyt = pickle.load(open(f"{month_dirr}/cLCS_{m}.p", "rb"))
    pxt[np.where(pxt < 0)] = np.nan
    pyt[np.where(pyt < 0)] = np.nan
    lonmin = lon.min()
    latmin = lat.min()
    lonmax = lon.max()
    latmax = lat.max()
    Plon, Plat = xy2sph(pxt*1e3, lonmin, pyt*1e3, latmin)
    #out_xyz = projection.transform_points(ccrs.Geodetic(), Plon, Plat)
    #out_lon = out_xyz[:, :, 0]
    #out_lat = out_xyz[:, :, 1]
    if not corners:
        corners = [lonmin, lonmax, latmin, latmax]
    ax.set_extent(corners, crs=ccrs.PlateCarree())
    nLCS = pxt.shape[0]
    for kk in range(0, nLCS, line_spacing):
        plot_lines(
            Plon[kk, :], Plat[kk, :], ax=ax, color=color, alpha=alpha, lw=lw, transform=ccrs.PlateCarree()
        )
    if save_fig:
        print(f"---- Saving Figure")
        fig.savefig(f"{month_dirr}/cLCS_monochrome_{m}.{save_fig}")
    print(f"---- Done")   
    return fig, ax
