import numpy as np

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