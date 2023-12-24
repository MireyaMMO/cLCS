import os
import sys
import numpy as np
from datetime import datetime, timedelta
import pickle
import opendrift
import numpy.ma as ma
from netCDF4 import Dataset
from cLCS.utils import *
import importlib


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

class mean_C(object):
    """
    mean_C
    computes sequential Cauchy-Green tensors, and the average for the full period

    Parameters
    ----------
    dirr: str
        directory where the Climatological files are located 
    filename: str
        name of the climatological file. Climatological file used here is a ROMS file. Change reader if the file is not ROMS
    month: int
        month that is going to be analysed
    T: int
        Timescale to use, use negative timescale if a backward in time simulation is needed 
    dt: int
        Default 6 hrs
    time_step_output: int
        Default 1 day 
    opendrift_reader: str
        Reader used by opendrift. This depends on the hydrodynamic model being used. For this study the default is 'reader_ROMS_native_MOANA'
    opendrift_model: str
        Module used by opendrift. For the LCS we use OceanDrift (passive particles) however other modules are available (see OpenDrift for more information)
        Default 'OceanDrift'
    vars_dict: dict
        Dictionary containing the names of the variables corresponding to the mask, longitude and latitude (mask_rho, lon_rho, lat_rho for ROMS files)


    Returns
    -------
    Output files:
    Trajectories_%m%d.nc: NetCDF file
        NetCDF file of the particle trajectories
    LCS_%m%d_CG.p: pickle File
        Pickle file containing the Cauchy Green Tensors for each day (C11, C22, C12)
    TOT-%m.p: pickle file
        Pickle file containing the accumulated values for the particle releases deployed for the climatolgical LCS calculation. 
            lon, lat : coordinates (degrees)
            lda2total, sqrtlda2total, ftletotal, C11total, C22total, C12total: Cauchy-Green associated parameters
            T: Used timescale
            xspan, yspan: x,y coordinates (km)
            count: number of deploys, important to obtain averages
    """
    def __init__(
      self,
      dirr, 
      climatology_file, 
      month, 
      T, 
      dt=6, #six hour time step
      time_step_output=86400, #daily output
      z=0, #surface
      opendrift_reader='reader_ROMS_native_MOANA', 
      opendrift_model='OceanDrift', 
      vars_dict={'mask': 'mask_rho'}
      ):
      # Import opendrift modules
      self.dirr = dirr
      self.climatology_file = climatology_file
#      self.climatology_file = dirr+filename
      self.month = month 
      self.T = T
      self.dt = dt
      self.time_step_output = time_step_output
      if isinstance(self.month, int):
          self.m = '%02d' % self.month
      self.new_directory = os.path.join(self.dirr, self.m)
      self.opendrift_reader = opendrift_reader
      self.opendrift_model = opendrift_model
      exec(f'from opendrift.readers import {self.opendrift_reader}')
      exec(
          f'from opendrift.models.{self.opendrift_model.lower()} import {self.opendrift_model}')
      self.vars_dict = vars_dict
      self.z = z
      # Cauchy-Green terms
      self.lda2total = 0
      self.sqrtlda2total = 0
      self.ftletotal = 0
      self.C11total = 0
      self.C22total = 0
      self.C12total = 0
      
    def set_directories(self):
        """Create output directories."""
        if not os.path.isdir(self.new_directory):
            os.makedirs(self.new_directory)
            
    def set_opendrift_configuration(self, file):
      o = eval(self.opendrift_model)(loglevel=20)
      reader = eval(self.opendrift_reader)(file)
      # dynamical landmask if true
      o.set_config('general:use_auto_landmask', False)
      o.max_speed = 5.0
      o.add_reader(reader)  # This adds the reader
      # seed
      # keep only particles from the "frame" that are on the ocean
      o.set_config('seed:ocean_only', True)
      ###############################
      # PHYSICS of Opendrift
      ###############################
      o.set_config('environment:fallback:x_wind', 0.0)
      o.set_config('environment:fallback:y_wind', 0.0)
      o.set_config('environment:fallback:x_sea_water_velocity', 0.0)
      o.set_config('environment:fallback:y_sea_water_velocity', 0.0)
      o.set_config(
          'environment:fallback:sea_floor_depth_below_sea_level', 10000.0)

      # drift
      o.set_config('environment:fallback:land_binary_mask', 0)
      o.set_config('drift:advection_scheme', 'runge-kutta4')  # or 'runge-kutta'
      # note current_uncertainty can be used to replicate an horizontal diffusion s
      o.set_config('drift:current_uncertainty', 0.0)

      Kxy = 0.1  # m2/s-1
      # using new config rather than current uncertainty
      o.set_config('drift:horizontal_diffusivity', Kxy)

      o.disable_vertical_motion()
      o.list_config()
      o.list_configspec()
      return o, reader
    
    def seed_particles_full_grid(self, reader, file):
        lon = reader.lon
        lat = reader.lat
        nc = Dataset(file)
        maskvar = self.vars_dict['mask']
        mask = nc.variables[maskvar][:]

        self.lon_origin = np.min(lon)  # [deg]
        self.lat_origin = np.min(lat)  # [deg]

        [xspan, yspan] = sph2xy(lon, self.lon_origin, lat, self.lat_origin)
        self.xspan = xspan*1e-3  # [km]
        self.yspan = yspan*1e-3  # [km]
        self.Nx0 = xspan.shape[0]
        self.Ny0 = yspan.shape[1]
        Nxy0 = self.Nx0*self.Ny0
        self.dx0 = np.nanmean(np.diff(self.xspan, axis=1))
        self.dy0 = np.nanmean(np.diff(self.yspan, axis=0))

        X0 = self.xspan.ravel()
        Y0 = self.yspan.ravel()
        # mask the land points from the beginning using mask_rho from roms
        X0[np.where(mask.ravel() == 0)] = np.nan
        Y0[np.where(mask.ravel() == 0)] = np.nan

        lon0, lat0 = xy2sph(X0*1e3, self.lon_origin, Y0*1e3, self.lat_origin)

        lon0 = ma.masked_where(np.isnan(lon0), lon0)
        lat0 = ma.masked_where(np.isnan(lat0), lat0)
        nonanindex = np.invert(np.isnan(lon0)) #non-land particles
        return lon0, lat0, nonanindex
    
    def obtain_final_particle_position(self, o, lon0, lat0, nonanindex):
        lon_OpenDrift = o.history['lon'][:]
        lat_OpenDrift = o.history['lat'][:]
        lon_OpenDrift[np.where(lon_OpenDrift[:] < 0)] = lon_OpenDrift[np.where(
            o.history['lon'] < 0)][:]+360
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
            end_lon[np.where(end_lon[:] < 0)] = end_lon[np.where(
                end_lon[:] < 0)][:]+360
            end_lat[nonanindex] = lat
        elif self.T < 0:  # Opendrift does something that when run backwards it inverts the order
            end_lon[nonanindex] = lon[::-1]
            end_lon[np.where(end_lon[:] < 0)] = end_lon[np.where(
                end_lon[:] < 0)][:]+360
            end_lat[nonanindex] = lat[::-1]
        return end_lon, end_lat

    def lat_lon_to_x_y(self, end_lon, end_lat):
        [end_x, end_y] = sph2xy(end_lon, self.lon_origin, end_lat, self.lat_origin)
        end_x = end_x*1e-3
        end_y = end_y*1e-3
        end_x = end_x.reshape(self.Nx0, self.Ny0)
        end_y = end_y.reshape(self.Nx0, self.Ny0)
        return end_x, end_y
        
    def calculate_Cauchy_Green(self, end_x, end_y):
          dxdy, dxdx = np.gradient(end_x, self.dy0, self.dx0)
          dydy, dydx = np.gradient(end_y, self.dy0, self.dx0)
          dxdx0 = np.reshape(dxdx, (self.Nx0, self.Ny0))
          dxdy0 = np.reshape(dxdy, (self.Nx0, self.Ny0))
          dydx0 = np.reshape(dydx, (self.Nx0, self.Ny0))
          dydy0 = np.reshape(dydy, (self.Nx0, self.Ny0))
          C11 = (dxdx0**2) + (dydx0**2)
          C12 = (dxdx0*dxdy0) + (dydx0*dydy0)
          C22 = (dxdy0**2) + (dydy0**2)
          detC = (C11*C22) - (C12**2)
          trC = C11 + C22
          lda2 = np.real(.5*trC + np.sqrt(.25*trC**2 - detC))
          ftle = np.log(lda2)/(2*np.abs(self.T))
          return C11, C12, C22, lda2, ftle

    def run(self):
      self.set_directories() # creates directory for output
      o, reader = self.set_opendrift_configuration(self.climatological_file)

      ######################################
      # Defines the times to run (1 per day)
      ######################################
      if self.T < 0:
          runtime = [reader.start_time +
                    timedelta(days=int(np.abs(self.T))), reader.end_time+timedelta(days=1)]
      elif self.T > 0:
          runtime = [reader.start_time, reader.end_time -
                    timedelta(days=int(np.abs(self.T)))]

      time = create_seed_times(runtime[0], runtime[1], timedelta(days=1))
      time_step = timedelta(hours=self.dt)
      duration = timedelta(days=int(np.abs(self.T)))

      self.lon0, self.lat0, nonanindex = self.seed_particles_full_grid(reader, self.climatological_file)
      
      for count, t in enumerate(time):
          # These lines are repeated as this will generate a deploy for each day allowing us to calculate the Cauchy Green tensors needed for the cLCS calculation
          # Set loglevel to 0 for debug information
          o, reader = self.set_opendrift_configuration(self.climatological_file, log_level=self.log_level)
          d = '%02d' % t.day
          namefile = f'{self.new_directory}/Trajectories_{self.m}{d}.nc'
          o.seed_elements(self.lon0[nonanindex], self.lat0[nonanindex],
                          time=t, z=self.z)
          # Foward in time
          if self.T > 0:
              o.run(duration=duration, time_step=time_step,
                    time_step_output=self.time_step_output, outfile=namefile)
          # Backward in time
          elif self.T < 0:
              o.run(duration=duration, time_step=-time_step,
                    time_step_output=self.time_step_output, outfile=namefile)

          end_lon, end_lat = self.obtain_final_particle_position(o, self.lon0, self.lat0, nonanindex)
          
          end_x, end_y = self.lat_lon_to_x_y(end_lon,end_lat)
          
          C11, C12, C22, lda2, ftle = self.calculate_cauchy_green(end_x, end_y)

          self.lda2total += lda2
          self.sqrtlda2total += np.sqrt(lda2)
          self.ftletotal += ftle
          self.C11total += C11
          self.C22total += C22
          self.C12total += C12
          
          pickle.dump([C11, C22, C12], open(f'{self.new_directory}/LCS_{self.m}{d}_CG.p', 'wb'))
      self.count = count
      pickle.dump([self.lon0, self.lat0, self.lda2total, self.sqrtlda2total, self.T, self.ftletotal, self.C11total,
                  self.C22total, self.C12total, self.xspan, self.yspan, self.count], open(f'{self.new_directory}/TOT-{self.m}.p', 'wb'))
