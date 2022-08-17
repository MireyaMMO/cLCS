import os
import sys
import numpy as np
from datetime import datetime, timedelta
import pickle
import opendrift
from opendrift.readers import reader_ROMS_native_MOANA
from opendrift.models.oceandrift import OceanDrift
import numpy.ma as ma
from netCDF4 import Dataset

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

def sph2xy(lambda0,lambda1, theta0, theta1):
############# SPH2XY Spherical to curvilinear spherical. ############
##### where X,Y are in meters and LAMBDA0,THETA0 are in degrees#####
   R=6371 * 1e3
   deg2rad = np.pi/180
   x = R * (lambda0 - lambda1) * deg2rad * np.cos(theta1*deg2rad)
   y = R * (theta0 - theta1) * deg2rad
   return x,y

def xy2sph(x, lambda1, y, theta1):
############# XY2SPH Curvilinear spherical to spherical. ############
##### where X,Y are in meters and LAMBDA1,THETA1 are in degrees#####
   R = 6371 * 1e3
   deg2rad = np.pi/180
   lambda0 = lambda1 + x/(R*np.cos(theta1*deg2rad)) / deg2rad
   theta0 = theta1 + y/R / deg2rad
   return lambda0,theta0


def mean_C(dirr,filename,month,T,dt=None,time_step_output=None):
 # mean_C computes sequential Cauchy-Green tensors, and the average for the full period
 # dirr: directory where the Climatological files are located 
 # filename: name of the climatological file. Climatological file used here is a ROMS file. Change reader if the file is not ROMS
 # month: month that is going to be analysed
 # T: Timescale to use, use negative timescale if a backward in time simulation is needed 
 # dt: Default 6 hrs
 # time_step_output: Default 1 day 
   if dt==None:
     dt=6

   if time_step_output==None:
     time_step_output=86400

   if isinstance(month,int):
     m = '%02d' % month

   os.mkdir(m) #creates directory for output
   o = OceanDrift(loglevel=20)  # Set loglevel to 0 for all the debug information
   reader_bop = reader_ROMS_native_MOANA.Reader(dirr+filename)
   o.set_config('general:use_auto_landmask', False) # dynamical landmask if true
   o.max_speed = 5.0
   o.add_reader(reader_bop) #This adds the reader
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
     runtime=[reader_bop.start_time+timedelta(days=int(np.abs(T))), reader_bop.end_time+timedelta(days=1)]
   elif T>0:
     runtime=[reader_bop.start_time, reader_bop.end_time-timedelta(days=int(np.abs(T)))]

   time = create_seed_times(runtime[0], runtime[1], timedelta(days = 1))
   time_step=timedelta(hours=dt)
   duration=timedelta(days=int(np.abs(T)))

   lon=reader_bop.lon
   lat=reader_bop.lat
   ph=dirr+filename
   nc = Dataset(ph)
   maskrho=nc.variables['mask_rho'][:]                                                     

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
     o = OceanDrift(loglevel=20)  # Set loglevel to 0 for debug information
     reader_bop = reader_ROMS_native_MOANA.Reader(dirr+filename)
     o.set_config('general:use_auto_landmask', False) # dynamical landmask if true
     o.max_speed = 5.0
     o.add_reader(reader_bop) #This adds the reader
     # seed
     o.set_config('seed:ocean_only',True) # keep only particles from the "frame" that are on the ocean
     o.set_config('environment:fallback:x_wind', 0.0)
     o.set_config('environment:fallback:y_wind', 0.0)
     o.set_config('environment:fallback:x_sea_water_velocity', 0.0)
     o.set_config('environment:fallback:y_sea_water_velocity', 0.0)
     o.set_config('environment:fallback:sea_floor_depth_below_sea_level', 10000.0
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
       ulrat[nonanindex]=tlat
     elif T<0: #Opendrift does something that when run backwards it inverts the order
       urlon[nonanindex]=tlon[::-1]
       urlon[np.where(urlon[:]<0)]=urlon[np.where(urlon[:]<0)][:]+360
       urlat[nonanindex]=tlat[::-1]

     [xNSWE, yNSWE] = sph2xy(urlon, lon_origin, urlat, lat_origin)
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
