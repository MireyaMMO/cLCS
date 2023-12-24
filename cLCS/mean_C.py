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


def mean_C(dirr,filename,month,T, dt=None,time_step_output=None, opendrift_reader= 'reader_ROMS_native_MOANA', opendrift_model = 'OceanDrift', vars_dict={}):
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
  #Import opendrift modules
  exec(f'from opendrift.readers import {opendrift_reader}')
  exec(f'from opendrift.models.{opendrift_model.lower()} import {opendrift_model}')
  if not dt:
    dt=6

  if not time_step_output:
    time_step_output=86400

  if isinstance(month,int):
    m = '%02d' % month

  os.mkdir(m) #creates directory for output
  o = eval(opendrift_model)(loglevel=20)  # Set loglevel to 0 for all the debug information
  reader = eval(opendrift_reader)(dirr+filename)
  o.set_config('general:use_auto_landmask', False) # dynamical landmask if true
  o.max_speed = 5.0
  o.add_reader(reader) #This adds the reader
  # seed
  o.set_config('seed:ocean_only',True) # keep only particles from the "frame" that are on the ocean
  ###############################
  # PHYSICS of Opendrift
  ###############################
  o.set_config('environment:fallback:x_wind', 0.0)
  o.set_config('environment:fallback:y_wind', 0.0)
  o.set_config('environment:fallback:x_sea_water_velocity', 0.0)
  o.set_config('environment:fallback:y_sea_water_velocity', 0.0)
  o.set_config('environment:fallback:sea_floor_depth_below_sea_level', 10000.0)

  # drift
  o.set_config('environment:fallback:land_binary_mask', 0)
  o.set_config('drift:advection_scheme','runge-kutta4') # or 'runge-kutta'
  o.set_config('drift:current_uncertainty', 0.0 ) # note current_uncertainty can be used to replicate an horizontal diffusion s

  Kxy = 0.1 # m2/s-1
  o.set_config('drift:horizontal_diffusivity',Kxy) # using new config rather than current uncertainty

  o.disable_vertical_motion()
  o.list_config()
  o.list_configspec()
  ######################################
  # Defines the times to run (1 per day)
  ######################################
  if T<0:  
    runtime=[reader.start_time+timedelta(days=int(np.abs(T))), reader.end_time+timedelta(days=1)]
  elif T>0:
    runtime=[reader.start_time, reader.end_time-timedelta(days=int(np.abs(T)))]

  time = create_seed_times(runtime[0], runtime[1], timedelta(days = 1))
  time_step=timedelta(hours=dt)
  duration=timedelta(days=int(np.abs(T)))

  lon=reader.lon
  lat=reader.lat
  climatological_file=dirr+filename
  nc = Dataset(climatological_file)
  maskvar = vars_dict['mask']
  maskrho=nc.variables[maskvar][:]                                                     

  lon0 = np.min(lon) # [deg]
  lat0 = np.min(lat) # [deg]

  [xspan, yspan] = sph2xy(lon, lon0, lat, lat0)
  xspan=xspan*1e-3 #[km]
  yspan=yspan*1e-3 #[km]
  Nx0=xspan.shape[0]
  Ny0=yspan.shape[1]
  Nxy0 = Nx0*Ny0

  dx0 = np.nanmean(np.diff(xspan,axis=1))
  dy0 = np.nanmean(np.diff(yspan,axis=0))

  X0 = xspan.ravel()
  Y0 = yspan.ravel()
  X0[np.where(maskrho.ravel()==0)]=np.nan #mask the land points from the beginning using mask_rho from roms
  Y0[np.where(maskrho.ravel()==0)]=np.nan

  lonNSWE0, latNSWE0 = xy2sph(X0*1e3, lon0, Y0*1e3, lat0)

  lonNSWE0=ma.masked_where(np.isnan(lonNSWE0),lonNSWE0)
  latNSWE0=ma.masked_where(np.isnan(latNSWE0),latNSWE0)
  nonanindex=np.invert(np.isnan(lonNSWE0))
  z=0

  count=0
  lda2total=0
  sqrtlda2total=0
  ftletotal=0
  C11total=0
  C22total=0
  C12total=0

  for i, t in enumerate(time):
    #These lines are repeated as this will generate a deploy for each day allowing us to calculate the Cauchy Green tensors needed for the cLCS calculation
    o = eval(opendrift_model)(loglevel=20)  # Set loglevel to 0 for debug information
    reader_bop = eval(opendrift_reader)(dirr+filename)
    o.set_config('general:use_auto_landmask', False) # dynamical landmask if true
    o.max_speed = 5.0
    o.add_reader(reader_bop) #This adds the reader
    # seed
    o.set_config('seed:ocean_only',True) # keep only particles from the "frame" that are on the ocean
    o.set_config('environment:fallback:x_wind', 0.0)
    o.set_config('environment:fallback:y_wind', 0.0)
    o.set_config('environment:fallback:x_sea_water_velocity', 0.0)
    o.set_config('environment:fallback:y_sea_water_velocity', 0.0)
    o.set_config('environment:fallback:sea_floor_depth_below_sea_level', 10000.0)
    o.set_config('environment:fallback:land_binary_mask', 0)
    o.set_config('drift:advection_scheme','runge-kutta4') # or 'runge-kutta'
    o.set_config('drift:current_uncertainty', 0.0 ) # note current_uncertainty can be used to replicate an horizontal diffusion spd_uncertain = sqrt(Kxy*2/dt)
    # horizontal diffusion
    Kxy = 0.1 # m2/s-1
    o.set_config('drift:horizontal_diffusivity',Kxy) # using new config rather than current uncertainty
    o.disable_vertical_motion()

    d = '%02d' %  time[i].day
    namefile = dirr+m+'/Trajectories_'+m+d+'.nc'
    o.seed_elements(lonNSWE0[nonanindex], latNSWE0[nonanindex],
                      time=time[i], z=z)
    #Foward in time
    if T>0:
      o.run(duration=duration, time_step=time_step, time_step_output = time_step_output, outfile= namefile)
    #Backward in time
    elif T<0:
      o.run(duration=duration, time_step=-time_step, time_step_output = time_step_output, outfile= namefile)

    loni=o.history['lon'][:]
    lati=o.history['lat'][:]
    loni[np.where(loni[:]<0)]=loni[np.where(o.history['lon']<0)][:]+360 
    urlat=np.zeros(latNSWE0.shape)
    urlon=np.zeros(lonNSWE0.shape)
    urlat[np.where(urlat==0)]=np.nan
    urlon[np.where(urlon==0)]=np.nan
    tlon=np.zeros(loni.shape[0])
    tlat=np.zeros(lati.shape[0])
    #If particle is stranded use last location
    for k in range(loni.shape[0]):
        notmask=np.where(loni[k,:].mask == False)[0]
        tlon[k]=loni[k,notmask][-1]
        tlat[k]=lati[k,notmask][-1]
    if T>0:
      urlon[nonanindex]=tlon
      urlon[np.where(urlon[:]<0)]=urlon[np.where(urlon[:]<0)][:]+360
      urlat[nonanindex]=tlat
    elif T<0: #Opendrift does something that when run backwards it inverts the order
      urlon[nonanindex]=tlon[::-1]
      urlon[np.where(urlon[:]<0)]=urlon[np.where(urlon[:]<0)][:]+360
      urlat[nonanindex]=tlat[::-1]

    [xNSWE, yNSWE] = sph2xy(urlon, lon0, urlat, lat0)
    xNSWE=xNSWE*1e-3
    yNSWE=yNSWE*1e-3
    xNSWE=xNSWE.reshape(Nx0,Ny0)
    yNSWE=yNSWE.reshape(Nx0,Ny0)
  ######################################
  # Calculate Cauchy-Green Tensors
  ######################################
    dxdy,dxdx = np.gradient(xNSWE,dy0,dx0)
    dydy,dydx = np.gradient(yNSWE,dy0,dx0)
    dxdx0 = np.reshape(dxdx, (Nx0,Ny0))
    dxdy0 = np.reshape(dxdy, (Nx0,Ny0))
    dydx0 = np.reshape(dydx, (Nx0,Ny0))
    dydy0 = np.reshape(dydy, (Nx0,Ny0))
    C11 = (dxdx0**2) + (dydx0**2)
    C12 = (dxdx0*dxdy0) + (dydx0*dydy0)
    C22 = (dxdy0**2) + (dydy0**2)

    detC = (C11*C22) - (C12**2)
    trC = C11 + C22
    lda2 = np.real(.5*trC + np.sqrt(.25*trC**2 - detC))
    ftle=np.log(lda2)/(2*np.abs(T))

    lda2total=lda2+lda2total
    sqrtlda2total=sqrtlda2total+np.sqrt(lda2)
    ftletotal=ftletotal+ftle
    C11total=C11total+C11
    C22total=C22total+C22
    C12total=C12total+C12
    count=count+1
    pickle.dump([C11,C22,C12], open(dirr+m+'/LCS_'+m+d+'_CG.p', 'wb'))
    del o, xNSWE, yNSWE, loni, lati,tlon,tlat, urlat, urlon, dxdy, dxdx, dydy, dydx, dxdx0, dxdy0, dydx0, dydy0#delete to clear memory 
  pickle.dump([lon,lat,lda2total,sqrtlda2total,T,ftletotal,C11total,C22total,C12total,xspan,yspan,count], open(dirr+m+'/TOT-'+m+'.p', 'wb'))
