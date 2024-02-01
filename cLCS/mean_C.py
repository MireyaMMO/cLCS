import os
import numpy as np
from datetime import datetime, timedelta
import pickle
import numpy.ma as ma
from cLCS.utils import *
import logging
import logging.config
import xarray as xr
import pandas as pd
import calendar

logging.basicConfig(level=logging.INFO)


def create_seed_times(start, end, delta):
    """
    create times at given interval to seed particles
    """
    out = []
    start_t = start
    end_t = datetime.strptime(str(end), "%Y-%m-%d %H:%M:%S")
    while start_t < end:
        out.append(start_t)
        start_t += delta
    return out


class mean_CG(object):
    """
    mean_C
    computes sequential Cauchy-Green tensors, and the average for the full period.
    No vertical motion is considered

    Parameters
    ----------
    dirr: str
        directory where the Climatological files are located
    filename: str
        name of the climatological file. Climatological file used here is a ROMS file. Change reader if the file is not ROMS
    month: int
        month that is going to be analysed
    dx: int
        Spacing between particles on the x-axis in km. Default grid structure (particle per grid point)
    dy: int
        Spacing between particles on the y-axis in km. Default grid structure (particle per grid point)
    T: int
        Timescale to use, use negative timescale if a backward in time simulation is needed
    dt: int
        Default 6 hrs
    frequency_of_deployments: int
        Default 1 day
    time_step_output: int
        Default 1 day
    opendrift_reader: str
        Reader used by opendrift. This depends on the hydrodynamic model being used. For this study the default is 'reader_ROMS_native_MOANA'
    opendrift_model: str
        Module used by opendrift. For the LCS we use OceanDrift (passive particles) however other modules are available (see OpenDrift for more information)
        Default 'OceanDrift'

    Other Parameters:
    ----------------
    log_level: int
        Level of legel wanted for OpenDrift processes, 0 to 100 range.
    vars_dict: dict
        Dictionary containing the names of the variables corresponding to the mask, longitude and latitude (mask_rho, lon_rho, lat_rho for ROMS files)
    max_speed: int
        Maximum speed that can be reached by the particles in m s^{-1}. Default 5.
    horizontal_diffusivity: int
        horizontal_diffusivity in m^{2} s^{-1}. Default 0.1
    advection_scheme: str
        Advection scheme for the particle releases. Default "runge-kutta4" (see OpenDrift for more options)
    save_trajectories: boolean
        Option to save the trajectories obtained from OpenDrift. Default False
    save_daily_CG = boolean
        Option to save the Cauchy Green Tensors for each day. Default False

    Returns
    -------
    Output files:
    TOT-%m.p: pickle file
        Pickle file containing the accumulated values for the particle releases deployed for the climatolgical LCS calculation.
            lon, lat : coordinates (degrees)
            lda2total, sqrtlda2total, ftletotal, C11total, C22total, C12total: Cauchy-Green associated parameters
            T: Used timescale
            xspan, yspan: x,y coordinates (km)
            count: number of deploys, important to obtain averages
    Trajectories_%m%d.nc: NetCDF file (Optional)
        NetCDF file of the particle trajectories. Default: False
    LCS_%m%d_CG.p: pickle File  (Optional)
        Pickle file containing the Cauchy Green Tensors for each day (C11, C22, C12). Default: False

    """

    def __init__(
        self,
        dirr,
        file,
        month,
        domain=False,
        dx0=None,
        dy0=None,
        T=-7,
        dt=6,  # six hour time step
        frequency_of_deployments=1,
        time_step_output=None,  # daily output
        z=0,  # surface
        opendrift_reader="reader_ROMS_native_MOANA",
        opendrift_model="OceanDrift",
        log_level=20,
        vars_dict={"mask": "mask_rho"},
        max_speed=5,
        horizontal_diffusivity=0.1,
        advection_scheme="runge-kutta4",
        save_trajectories=False,
        save_daily_CG=False,
        start_releases=None,
        end_releases=None,
    ):
        # Environment parameters
        self.dirr = dirr
        self.file = file
        self.month = month
        if isinstance(self.month, int):
            self.m = "%02d" % self.month
        self.new_directory = os.path.join(self.dirr, self.m)
        self.domain = domain

        # Particle release parameters
        self.dx0 = dx0
        self.dy0 = dy0
        self.T = T
        self.dt = dt
        self.frequency_of_deployments = frequency_of_deployments
        self.time_step_output = time_step_output
        self.start_releases = start_releases
        self.end_releases = end_releases

        # OpenDrift Configuration Parameters
        self.opendrift_reader = opendrift_reader
        self.opendrift_model = opendrift_model
        self.log_level = log_level
        self.vars_dict = vars_dict
        self.z = z
        self.max_speed = max_speed
        self.horizontal_diffusivity = horizontal_diffusivity
        self.advection_scheme = advection_scheme
        self.save_trajectories = save_trajectories
        self.save_daily_CG = save_daily_CG

        # Cauchy-Green terms
        self.lda2total = 0
        self.sqrtlda2total = 0
        self.ftletotal = 0
        self.C11total = 0
        self.C22total = 0
        self.C12total = 0

        # Log information
        self.logger = logging

    def set_directories(self):
        """Create output directories."""
        self.logger.info("--- Creating output directory")
        if not os.path.isdir(self.new_directory):
            os.makedirs(self.new_directory)

    def get_reader(self, file):
        exec(f"from opendrift.readers import {self.opendrift_reader}")
        reader = eval(self.opendrift_reader).Reader(file)
        return reader
    
    def set_opendrift_configuration(self, file):
        # self.logger.info("--- Setting OpenDrift Configuration")
        exec(
            f"from opendrift.models.{self.opendrift_model.lower()} import {self.opendrift_model}"
        )
        o = eval(self.opendrift_model)(loglevel=self.log_level)
        reader = self.get_reader(file)
        # dynamical landmask if true
        o.set_config("general:use_auto_landmask", False)
        o.max_speed = self.max_speed
        o.add_reader(reader)
        # keep only particles from the "frame" that are on the ocean
        o.set_config("seed:ocean_only", True)
        ###############################
        # PHYSICS of Opendrift
        ###############################
        o.set_config("environment:fallback:x_wind", 0.0)
        o.set_config("environment:fallback:y_wind", 0.0)
        o.set_config("environment:fallback:x_sea_water_velocity", 0.0)
        o.set_config("environment:fallback:y_sea_water_velocity", 0.0)
        o.set_config("environment:fallback:sea_floor_depth_below_sea_level", 10000.0)

        # drift
        o.set_config("environment:fallback:land_binary_mask", 0)
        o.set_config(
            "drift:advection_scheme", self.advection_scheme
        )  # or 'runge-kutta'
        # note current_uncertainty can be used to replicate an horizontal diffusion s
        o.set_config("drift:current_uncertainty", 0.0)

        Kxy = self.horizontal_diffusivity  # m2/s-1
        # using new config rather than current uncertainty
        o.set_config("drift:horizontal_diffusivity", Kxy)

        o.disable_vertical_motion()
        #        o.list_config()
        #        o.list_configspec()
        return o, reader
    
    def get_domain(self, ds, ):
        lonvar = self.vars_dict["lon"]
        latvar = self.vars_dict["lat"]
        maskvar = self.vars_dict["mask"]
        londim = self.vars_dict["lon_dim"]
        latdim = self.vars_dict["lat_dim"]
        xi = np.where((np.unique(self.lon)>=self.domain[0]) & (np.unique(self.lon)<=self.domain[1]))[0]
        yi = np.where((np.unique(self.lat)>=self.domain[2]) & (np.unique(self.lat)<=self.domain[3]))[0]
        ##Assuming a 2D grid as for roms
        try:
            self.lon = ds[lonvar].isel({londim : slice(xi[0], xi[-1]), latdim:slice(yi[0], yi[-1])}).values
            self.lat = ds[latvar].isel({londim : slice(xi[0], xi[-1]), latdim:slice(yi[0], yi[-1])}).values
        except:
            self.lon = ds[lonvar].isel({londim : slice(xi[0], xi[-1])}).values
            self.lat = ds[latvar].isel({latdim:slice(yi[0], yi[-1])}).values
        self.mask = ds[maskvar].isel({londim : slice(xi[0], xi[-1]), latdim:slice(yi[0], yi[-1])}).values
        self.lon_origin = self.lon.min()
        self.lat_origin = self.lat.min()

    def seed_particles(self, ds):
        lonvar = self.vars_dict["lon"]
        latvar = self.vars_dict["lat"]
        # maskvar = self.vars_dict["mask"]
        self.lon = ds[lonvar].values
        self.lat = ds[latvar].values
        # self.mask = ds[maskvar].values
        self.lon_origin = self.lon.min()
        self.lat_origin = self.lat.min()
        if self.domain:
            self.get_domain(ds)
        if self.dx0 and self.dy0:
            lonmax = self.lon.max()
            latmax = self.lat.max()
            xmax, ymax = sph2xy(lonmax, self.lon_origin, latmax, self.lat_origin)
            x = np.arange(0, xmax*1e-3, self.dx0)
            y = np.arange(0, ymax*1e-3, self.dy0)
            self.xspan, self.yspan = np.meshgrid(x,y)
            lon, lat = xy2sph(self.xspan*1e3, self.lon_origin, self.yspan*1e3, self.lat_origin)
            from scipy.interpolate import griddata
            self.mask = griddata((self.lon.ravel(),self.lat.ravel()), self.mask.ravel(), (lon.ravel(), lat.ravel()), method='nearest').reshape(lon.shape)
            self.lon = lon
            self.lat = lat
        else:
            [xspan, yspan] = sph2xy(self.lon, self.lon_origin, self.lat, self.lat_origin)
            self.xspan = xspan * 1e-3  # [km]
            self.yspan = yspan * 1e-3  # [km]
            self.dx0 = np.nanmean(np.diff(self.xspan, axis=1))
            self.dy0 = np.nanmean(np.diff(self.yspan, axis=0))
        self.x = np.copy(self.xspan)
        self.y = np.copy(self.yspan)
        self.Nx0 = self.xspan.shape[0]
        self.Ny0 = self.yspan.shape[1]
        Nxy0 = self.Nx0 * self.Ny0
        X0 = self.xspan.ravel()
        Y0 = self.yspan.ravel()
        # mask the land points from the beginning using mask_rho from roms
        X0[np.where(self.mask.ravel() == 0)] = np.nan
        Y0[np.where(self.mask.ravel() == 0)] = np.nan

        lon0, lat0 = xy2sph(X0 * 1e3, self.lon_origin, Y0 * 1e3, self.lat_origin)

        lon0 = ma.masked_where(np.isnan(lon0), lon0)
        lat0 = ma.masked_where(np.isnan(lat0), lat0)
        nonanindex = np.invert(np.isnan(lon0))  # non-land particles

        return lon0, lat0, nonanindex

    def seed_particles_schism(self, ds):

        self.lon = ds[self.vars_dict["lon"]].values
        self.lat = ds[self.vars_dict["lat"]].values   

        self.lon_origin = self.domain[0]
        self.lat_origin = self.domain[-1]
        lonmax = self.domain[1]
        latmax = self.domain[2]

        xmax, ymax = sph2xy(lonmax, self.lon_origin, latmax, self.lat_origin)
        x = np.arange(0, xmax*1e-3, self.dx0)
        y = np.arange(0, ymax*1e-3, self.dy0)
        self.xspan, self.yspan = np.meshgrid(x,y)
        self.lon, self.lat = xy2sph(self.xspan*1e3, self.lon_origin, self.yspan*1e3, self.lat_origin)

        self.x = np.copy(self.xspan)
        self.y = np.copy(self.yspan)
        self.Nx0 = self.xspan.shape[0]
        self.Ny0 = self.yspan.shape[1]
        Nxy0 = self.Nx0 * self.Ny0
        X0 = self.xspan.ravel()
        Y0 = self.yspan.ravel()

        ## TODO: 
        ## GET LAND INDECES (AS IN "nonnanindex") BASED ON GRID POINTS THAT
        ## FALL OUT OF SCHISM MESH POLYGON. TO DO THIS I WILL NEED TO INSTALL
        ## SCHISM-TOOLS PACKAGE SO WE CAN LEVERAGE SOME TOOLS FOR HANDLING
        ## SCHISM UNSTRUCTURED MESH AND SOME POLYGON OPERATIONS

        ## mask the land points from the beginning using mask_rho from roms
        # X0[np.where(self.mask.ravel() == 0)] = np.nan
        # Y0[np.where(self.mask.ravel() == 0)] = np.nan

        lon0, lat0 = xy2sph(X0 * 1e3, self.lon_origin, Y0 * 1e3, self.lat_origin)

        # lon0 = ma.masked_where(np.isnan(lon0), lon0)
        # lat0 = ma.masked_where(np.isnan(lat0), lat0)
        # nonanindex = np.invert(np.isnan(lon0))  # non-land particles
        
        # return lon0, lat0, nonanindex
        return lon0, lat0

    def obtain_final_particle_position(self, o, lon0, lat0, nonanindex):
        self.logger.info("--- Obtaining the final position of particles")
        print("--- Obtaining the final position of particles")
        lon_OpenDrift = o.history["lon"][:]
        lat_OpenDrift = o.history["lat"][:]
        lon_OpenDrift[np.where(lon_OpenDrift[:] < 0)] = (
            lon_OpenDrift[np.where(o.history["lon"] < 0)][:] + 360
        )
        end_lat = np.zeros(lat0.shape)
        end_lon = np.zeros(lon0.shape)
        end_lat[np.where(end_lat == 0)] = np.nan
        end_lon[np.where(end_lon == 0)] = np.nan
        lon = np.zeros(lon_OpenDrift.shape[0])
        lat = np.zeros(lat_OpenDrift.shape[0])
        # If particle is stranded use last location
        for k in range(lon_OpenDrift.shape[0]):
            notmask = np.where(lon_OpenDrift[k, :].mask == False)[0]
            lon[k] = lon_OpenDrift[k, notmask][-1]
            lat[k] = lat_OpenDrift[k, notmask][-1]
        if self.T > 0:
            end_lon[nonanindex] = lon
            end_lon[np.where(end_lon[:] < 0)] = (
                end_lon[np.where(end_lon[:] < 0)][:] + 360
            )
            end_lat[nonanindex] = lat
        elif (
            self.T < 0
        ):  # Opendrift does something that when run backwards it inverts the order
            end_lon[nonanindex] = lon[::-1]
            end_lon[np.where(end_lon[:] < 0)] = (
                end_lon[np.where(end_lon[:] < 0)][:] + 360
            )
            end_lat[nonanindex] = lat[::-1]
        return end_lon, end_lat

    def obtain_final_particle_position_schism(self, o, lon0, lat0):
        '''
        NOTE: THIS IS A COPY/PASTE FROM "obtain_final_particle_position" ABOVE,
        ONLY DIFFERENCE IS THAT '[nonanindex]' INDEXATION ISN'T USED HERE.
        '''

        self.logger.info("--- Obtaining the final position of particles")
        print("--- Obtaining the final position of particles")
        lon_OpenDrift = o.history["lon"][:]
        lat_OpenDrift = o.history["lat"][:]
        lon_OpenDrift[np.where(lon_OpenDrift[:] < 0)] = (
            lon_OpenDrift[np.where(o.history["lon"] < 0)][:] + 360
        )
        end_lat = np.zeros(lat0.shape)
        end_lon = np.zeros(lon0.shape)
        end_lat[np.where(end_lat == 0)] = np.nan
        end_lon[np.where(end_lon == 0)] = np.nan
        lon = np.zeros(lon_OpenDrift.shape[0])
        lat = np.zeros(lat_OpenDrift.shape[0])
        # If particle is stranded use last location
        for k in range(lon_OpenDrift.shape[0]):
            notmask = np.where(lon_OpenDrift[k, :].mask == False)[0]
            lon[k] = lon_OpenDrift[k, notmask][-1]
            lat[k] = lat_OpenDrift[k, notmask][-1]
        if self.T > 0:
            end_lon = lon
            end_lon[np.where(end_lon[:] < 0)] = (
                end_lon[np.where(end_lon[:] < 0)][:] + 360
            )
            end_lat = lat
        elif (
            self.T < 0
        ):  # Opendrift does something that when run backwards it inverts the order
            end_lon = lon[::-1]
            end_lon[np.where(end_lon[:] < 0)] = (
                end_lon[np.where(end_lon[:] < 0)][:] + 360
            )
            end_lat = lat[::-1]
        return end_lon, end_lat

    def lat_lon_to_x_y(self, end_lon, end_lat):
        [end_x, end_y] = sph2xy(end_lon, self.lon_origin, end_lat, self.lat_origin)
        end_x = end_x * 1e-3
        end_y = end_y * 1e-3
        end_x = end_x.reshape(self.Nx0, self.Ny0)
        end_y = end_y.reshape(self.Nx0, self.Ny0)
        return end_x, end_y

    def calculate_Cauchy_Green(self, end_x, end_y):
        self.logger.info("--- Calculating the Cauchy-Green Tensor")
        print("--- Calculating Cauchy-Green Tensor")
        dxdy, dxdx = np.gradient(end_x, self.dy0, self.dx0)
        dydy, dydx = np.gradient(end_y, self.dy0, self.dx0)
        dxdx0 = np.reshape(dxdx, (self.Nx0, self.Ny0))
        dxdy0 = np.reshape(dxdy, (self.Nx0, self.Ny0))
        dydx0 = np.reshape(dydx, (self.Nx0, self.Ny0))
        dydy0 = np.reshape(dydy, (self.Nx0, self.Ny0))
        C11 = (dxdx0**2) + (dydx0**2)
        C12 = (dxdx0 * dxdy0) + (dydx0 * dydy0)
        C22 = (dxdy0**2) + (dydy0**2)
        detC = (C11 * C22) - (C12**2)
        trC = C11 + C22
        lda2 = np.real(0.5 * trC + np.sqrt(0.25 * trC**2 - detC))
        ftle = np.log(lda2) / (2 * np.abs(self.T))
        return C11, C12, C22, lda2, ftle

    def run(self):
        self.set_directories()  # creates directory for output
        ds = xr.open_dataset(self.file)
        if not self.start_releases:
            start_time = pd.to_datetime(ds[self.vars_dict["time"]][0].values)
        else:
            start_time = self.start_releases
        if not self.end_releases:
            end_time = pd.to_datetime(ds[self.vars_dict["time"]][-1].values)
        else:
            end_time = self.end_releases
            
        if self.T < 0:
            runtime = [
                start_time + timedelta(days=int(np.abs(self.T))),
                end_time + timedelta(days=1),
            ]
        elif self.T > 0:
            runtime = [
                start_time,
                end_time - timedelta(days=int(np.abs(self.T))),
            ]

        time = create_seed_times(
            runtime[0], runtime[1], timedelta(days=self.frequency_of_deployments)
        )
        time_step = timedelta(hours=self.dt)
        duration = timedelta(days=int(np.abs(self.T)))

        reader = self.get_reader(self.file)
        

        ## NOTE: USING SCHISM SPECIFIC METHODS ---------------------------------------------------
        # self.lon0, self.lat0, nonanindex = self.seed_particles(ds) ## NOTE: WORKS ONLY FOR ROMS
        self.lon0, self.lat0 = self.seed_particles_schism(ds)
        ## ---------------------------------------------------------------------------------------

        for count, t in enumerate(time, 1):
            self.logger.info(f"--- {t} Release")
            print(f"--- {t} Release {count}/{len(time)}")
            o, reader = self.set_opendrift_configuration(self.file)
            d = "%02d" % t.day
            if self.save_trajectories:
                namefile = f"{self.new_directory}/Trajectories_{self.m}{d}.nc"
            else:
                namefile = None
            
            ## NOTE: USING SCHISM SPECIFIC METHODS -----------------------------------------------------------------------
            # o.seed_elements(self.lon0[nonanindex], self.lat0[nonanindex], time=t, z=self.z) ## NOTE: WORKS ONLY FOR ROMS
            o.seed_elements(self.lon0, self.lat0, time=t, z=self.z)
            ## -----------------------------------------------------------------------------------------------------------

            # Foward in time
            if self.T > 0:
                o.run(
                    duration=duration,
                    time_step=time_step,
                    time_step_output=self.time_step_output,
                    outfile=namefile,
                )
            # Backward in time
            elif self.T < 0:
                o.run(
                    duration=duration,
                    time_step=-time_step,
                    time_step_output=self.time_step_output,
                    outfile=namefile,
                )

            ## NOTE: USING SCHISM SPECIFIC METHODS -----------------------------------------------------------------------------------
            # end_lon, end_lat = self.obtain_final_particle_position(o, self.lon0, self.lat0, nonanindex) ## NOTE: WORKS ONLY FOR ROMS
            end_lon, end_lat = self.obtain_final_particle_position_schism(o, self.lon0, self.lat0)
            ## -----------------------------------------------------------------------------------------------------------------------
            
            end_x, end_y = self.lat_lon_to_x_y(end_lon, end_lat)

            C11, C12, C22, lda2, ftle = self.calculate_Cauchy_Green(end_x, end_y)

            self.lda2total += lda2
            self.sqrtlda2total += np.sqrt(lda2)
            self.ftletotal += ftle
            self.C11total += C11
            self.C22total += C22
            self.C12total += C12
            if self.save_daily_CG:
                pickle.dump(
                    [self.lon, self.lat, lda2, np.sqrt(lda2), self.T, ftle, C11, C22, C12, self.x, self.y],
                    open(f"{self.new_directory}/LCS_{self.m}{d}-CG.p", "wb"),
                )
            del o, reader
        self.count = count
        pickle.dump(
                [
                    self.lon,
                    self.lat,
                    self.lda2total,
                    self.sqrtlda2total,
                    self.T,
                    self.ftletotal,
                    self.C11total,
                    self.C22total,
                    self.C12total,
                    self.x,
                    self.y,
                    self.count,
                ],
                open(f"{self.new_directory}/TOT-{self.m}.p", "wb"),
            )
        
        print(
            f"Calculation of climatological LCS done for {calendar.month_name[self.month]}"
        )
