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


def plot_colourline(x, y, c, cmap, ax=None, lw=0.8, transform=None, climatology=None):
    # Plots LCSs using coloured lines to define strength of attraction
    if climatology:
        c = cmap((c - 0.4) / (1.4 - 0.4))
    else:
        c = cmap((c - 0.4) / (3.0 - 0.4))
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
    month_or_id,
    colourmap=None,
    lw=0.8,
    fig=None,
    ax=None,
    projection=None,
    line_spacing=4,
    save_fig=False,
    corners=None,
    climatology=True,
    squeezelines_file=None,
    CG_file=None,
    LCS_strength=None,
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
    if isinstance(month_or_id, int):
            month_or_id = "%02d" % month_or_id
            
    month_dirr = os.path.join(dirr, month_or_id)
    if climatology:
        TOT_CG_path = os.path.join(month_dirr, f"TOT-{month_or_id}.p")
        lon, lat, _, sqrtlda2total, _, _, _, _, _, _, _ , count = pickle.load(
            open(TOT_CG_path, "rb")
        )
        N = count
        pxt, pyt = pickle.load(open(f"{month_dirr}/cLCS_{month_or_id}.p", "rb"))
        pxt[np.where(pxt < 0)] = np.nan
        pyt[np.where(pyt < 0)] = np.nan
        z = np.log(sqrtlda2total / N)
        z[np.where(z>3)]=np.nan
        outfile = f"{month_dirr}/cLCS_{month_or_id}"
    else:
        try:
            lon, lat, _, sqrtlda2total, _, _, _, _, _, _, _  = pickle.load(
            open(CG_file, "rb")
        )
            pxt, pyt = pickle.load(open(squeezelines_file, "rb"))
            pxt[np.where(pxt < 0)] = np.nan
            pyt[np.where(pyt < 0)] = np.nan
            z = np.log(sqrtlda2total)
            #z[np.where(z>10)]=np.nan
            outfile = squeezelines_file.split(".")[0]
        except:
            print("Path to CG file and squeezelines file must be provided")
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

    cmap = get_colourmap(colourmap)
    print(f"---- Squeezeline and associated data loaded")
    if LCS_strength=='weak':
        climatology=True
    for kk in range(0, nLCS, line_spacing): 
        xs = Plon[kk, :]
        ys = Plat[kk, :]
        zs = interpolator(xs, ys)
        #zs[np.where(np.isnan(zs))] = 0
        plot_colourline(xs, ys, zs, cmap, ax=ax, lw=lw, transform=ccrs.PlateCarree(), climatology=climatology)
    if save_fig:
        print(f"---- Saving Figure")
        fig.savefig(f"{outfile}.{save_fig}")
    print(f"---- Done")
    return fig, ax


def cLCSrho_cartopy_monochrome(
    dirr,
    month_or_id,
    fig=None,
    ax=None,
    color="k",
    alpha="1",
    lw=0.8,
    projection=None,
    line_spacing=1,
    save_fig=False,
    corners=None,
    climatology=True,
    squeezelines_file=None,
    CG_file=None,
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
    if isinstance(month_or_id, int):
        month_or_id = "%02d" %month_or_id
    month_dirr = os.path.join(dirr, month_or_id)
    if climatology:
        TOT_CG_path = os.path.join(month_dirr, f"TOT-{month_or_id}.p")
        lon, lat, _, _, _, _, _, _, _, _, _ , _ = pickle.load(
            open(TOT_CG_path, "rb")
        )
        pxt, pyt = pickle.load(open(f"{month_dirr}/cLCS_{month_or_id}.p", "rb"))
        outfile = f"{month_dirr}/cLCS_{month_or_id}"
    else:
        try:
            lon, lat, _, _, _, _, _, _, _, _, _  = pickle.load(
            open(CG_file, "rb")
        )
            pxt, pyt = pickle.load(open(squeezelines_file, "rb"))
            outfile = squeezelines_file.split(".")[0]
        except:
            print("Path to CG file and squeezelines file must be provided")       
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
        fig.savefig(f"{outfile}_monochrome.{save_fig}")
    print(f"---- Done")   
    return fig, ax
