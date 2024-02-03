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
    month or id: int or str
        month that is going to be analysed or id for single release
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
        month_or_id,
        climatology=True,
        domain=False,
        dx0=1,
        dy0=1,
        T=-7,
        dt=6*3600,  # six hour time step
        frequency_of_deployments=1,
        time_step_output=86400,  # daily output
        z=0,  # surface
        opendrift_reader="reader_ROMS_native_MOANA",
        opendrift_model="OceanDrift",
        vars_dict=None,
        log_level=20,
        max_speed=5,
        horizontal_diffusivity=0.1,
        advection_scheme="runge-kutta4",
        save_trajectories=False,
        save_daily_CG=False,
        start_release=None,
    ):
        # Environment parameters
        self.dirr = dirr
        self.file = file
        self.month_or_id = month_or_id
        if isinstance(self.month_or_id, int):
            self.month_or_id = "%02d" % self.month_or_id
        self.new_directory = os.path.join(self.dirr, self.month_or_id)
        self.domain = domain
        self.climatology = climatology

        # Particle release parameters
        self.dx0 = dx0
        self.dy0 = dy0
        self.T = T
        self.dt = dt
        self.frequency_of_deployments = frequency_of_deployments
        self.time_step_output = time_step_output
        self.start_release = start_release
        #self.end_releases = end_releases

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
        self.reader = eval(self.opendrift_reader).Reader(file)
    
    def set_opendrift_configuration(self):
        # self.logger.info("--- Setting OpenDrift Configuration")
        exec(
            f"from opendrift.models.{self.opendrift_model.lower()} import {self.opendrift_model}"
        )
        o = eval(self.opendrift_model)(loglevel=self.log_level)
        self.get_reader(self.file)
        
        o.set_config("general:use_auto_landmask", False)
        o.max_speed = self.max_speed
        if 'schism' in self.opendrift_reader:
            from opendrift.readers import reader_global_landmask
            reader_landmask = reader_global_landmask.Reader()
            o.add_reader([self.reader, reader_landmask])
        else:
            o.add_reader(self.reader)
        o.set_config("seed:ocean_only", False) #Particles set on land not moved to ocean
        
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
        )  
        # note current_uncertainty can be used to replicate an horizontal diffusion s
        o.set_config("drift:current_uncertainty", 0.0)

        Kxy = self.horizontal_diffusivity  # m2/s-1
        # using new config rather than current uncertainty
        o.set_config("drift:horizontal_diffusivity", Kxy)

        o.disable_vertical_motion()
        #        o.list_config()
        #        o.list_configspec()
        return o

    # def get_mask(self,ds, type):
    #     maskvar = self.variable_mapping["land_binary_mask"]
    #     self.mask = ds[maskvar].values

    def seed_particles(self, ds):
        if self.domain:
            self.lon_origin = self.domain[0] 
            lonmax = self.domain[1]
            self.lat_origin = self.domain[2]
            latmax = self.domain[3]
        else:
            try:
                lonvar = self.vars_dict["lon"]
                latvar = self.vars_dict["lat"]
                lon = ds[lonvar].values
                lat = ds[latvar].values
                self.lon_origin = lon.min()
                self.lat_origin = lat.min()
                lonmax = lon.max()
                latmax = lat.max()
            except:
                print('Please provide lon and lat variable names for mapping in a dictionary')
        #if self.dx0 and self.dy0:
        xmax, ymax = sph2xy(lonmax, self.lon_origin, latmax, self.lat_origin)
        x = np.arange(0, xmax*1e-3, self.dx0)
        y = np.arange(0, ymax*1e-3, self.dy0)
        self.xspan, self.yspan = np.meshgrid(x,y)
        self.lon, self.lat = xy2sph(self.xspan*1e3, self.lon_origin, self.yspan*1e3, self.lat_origin)
        self.Nx0 = self.xspan.shape[1]
        self.Ny0 = self.yspan.shape[0]
    
    def obtain_final_particle_position(self, o):#, nonanindex):
        self.logger.info("--- Obtaining the final position of particles")
        print("--- Obtaining the final position of particles")
        lon_OpenDrift = o.history["lon"][:]
        lat_OpenDrift = o.history["lat"][:]
        #Particles that were deactivated since step 1 (over-land particles)
        status = o.history["status"][:]
        nonanindex= np.where(status[:,0]==1)[0] 
        lon = np.zeros(lon_OpenDrift.shape[0])
        lat = np.zeros(lat_OpenDrift.shape[0])
        # If particle is stranded use last location
        for k in range(lon_OpenDrift.shape[0]):
            notmask = np.where(lon_OpenDrift[k, :].mask == False)[0]
            if len(notmask)>0:
                lon[k] = lon_OpenDrift[k, notmask][-1]
                lat[k] = lat_OpenDrift[k, notmask][-1]
        lon[np.where(lon[:] < 0)] = (
                lon[np.where(lon[:] < 0)][:] + 360
            )
        if self.T < 0:  
            # Opendrift does something that when run backwards it inverts the order
            lon = lon[::-1]
            lat = lat[::-1]
            status = o.history["status"][::-1]
            nonanindex= np.where(status[:,0]==1)[0] 
        lon[nonanindex]=0
        lat[nonanindex]=0
        return lon, lat

    def lat_lon_to_x_y(self, end_lon, end_lat):
        [end_x, end_y] = sph2xy(end_lon, self.lon_origin, end_lat, self.lat_origin)
        end_x = end_x * 1e-3
        end_y = end_y * 1e-3
        end_x = end_x.reshape(self.Ny0, self.Nx0)
        end_y = end_y.reshape(self.Ny0, self.Nx0)
        return end_x, end_y

    def calculate_Cauchy_Green(self, end_x, end_y):
        self.logger.info("--- Calculating the Cauchy-Green Tensor")
        print("--- Calculating Cauchy-Green Tensor")
        dxdy, dxdx = np.gradient(end_x, self.dy0, self.dx0)
        dydy, dydx = np.gradient(end_y, self.dy0, self.dx0)
        dxdx0 = np.reshape(dxdx, (self.Ny0, self.Nx0))
        dxdy0 = np.reshape(dxdy, (self.Ny0, self.Nx0))
        dydx0 = np.reshape(dydx, (self.Ny0, self.Nx0))
        dydy0 = np.reshape(dydy, (self.Ny0, self.Nx0))

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
        self.get_reader(self.file) 
        if "ROMS" in self.opendrift_reader:
            self.variable_mapping = self.reader.ROMS_variable_mapping
        else:
            try:
                self.variable_mapping = self.reader.variable_mapping
            except:
                print("Can't obtain variable mapping")
        #self.var_dict= reader.variable_mapping
        try:
            ds = xr.open_dataset(self.file)
        except:
            ds = xr.open_mfdataset(self.file)
            
        if self.climatology:
            try:
                start_time = pd.to_datetime(self.reader.times[0])
                end_time = pd.to_datetime(self.reader.times[-1])
            except:
                start_time = self.reader.times[0]
                end_time = self.reader.times[-1]
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
        elif not self.climatology:
            time = self.start_release
            
        time_step = timedelta(seconds=self.dt)
        duration = timedelta(days=int(np.abs(self.T)))

        #reader = self.get_reader(self.file)
        self.seed_particles(ds)
        for count, t in enumerate(time, 1):
            self.logger.info(f"--- {t} Release")
            print(f"--- {t} Release {count}/{len(time)}")
            o  = self.set_opendrift_configuration()
            d = "%02d" % t.day
            if self.save_trajectories:
                namefile = f"{self.new_directory}/Trajectories_{self.month_or_id}{d}.nc"
            else:
                namefile = None
            o.seed_elements(
                self.lon, self.lat, time=t, z=self.z
            )
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

            end_lon, end_lat = self.obtain_final_particle_position(o)
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
                    [self.lon, self.lat, lda2, np.sqrt(lda2), self.T, ftle, C11, C22, C12, self.xspan, self.yspan],
                    open(f"{self.new_directory}/LCS_{self.month_or_id}_{d}-CG.p", "wb"),
                )
            del o, self.reader
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
                    self.xspan,
                    self.yspan,
                    self.count,
                ],
                open(f"{self.new_directory}/TOT-{self.month_or_id}.p", "wb"),
            )
        try:
            print(
                f"Calculation of climatological LCS done for {calendar.month_name[self.month_or_id]}"
            )
        except:
            print(
                f"Calculation of climatological LCS done for {self.month_or_id}"
            )
