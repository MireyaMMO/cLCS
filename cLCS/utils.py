import numpy as np
import matplotlib.cm as cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

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

def get_colourmap(name):
    if name=='Zissou':
        colors = [(.98, .98, .95),(.23, .60, .69), (.47, .71, .77),(.92, .8, .16), (.88,.68,0),(.95,.10,0),(.79,.08,0)]  # R -> G -> B
        cmap = LinearSegmentedColormap.from_list(name, colors, N=200)
    elif name=='BlueOrange':
        top = cm.get_cmap('Oranges', 128) # r means reversed version
        bottom = cm.get_cmap('Blues_r', 128)# combine it all
        colors = np.vstack((bottom(np.linspace(0, 1, 128)),
                               top(np.linspace(0, 1, 128))))# create a new colormaps with a name of OrangeBlue
        cmap = ListedColormap(colors, name)
    elif name=='Color_blind_1':
        colors = [(.67, .34, .11), (.89, .61, .34),(1, .87, .67), (.67,.72,.86),(.30,.45,.71)]  # R -> G -> B
        cmap = LinearSegmentedColormap.from_list(name, colors, N=200)
    elif name=='Duran_cLCS':
        colors = np.array([[1.0000,   1.0000,   0.9987, 1],
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
        cmap = LinearSegmentedColormap.from_list(name, colors)
    elif name=='RedYellowBlue':
        colors = [(.843, .188, .153), (.988, .553, .349),(.996, .878, .565), (.569,.749,.859),(.271,.459,.706)]  # R -> G -> B
        cmap = LinearSegmentedColormap.from_list(name, colors, N=200)
    elif name=='BlueYellowRed':
        colors = [(.843, .188, .153), (.988, .553, .349),(.996, .878, .565), (.569,.749,.859),(.271,.459,.706)]  # R -> G -> B
        cmap = LinearSegmentedColormap.from_list(name, colors[::-1], N=200)
    elif name=='AlgaeSalmon':
        colors = [(.557, .792, .902), (.165, .616, .561),(.914, .769, .416), (.957,.635,.38),(.906,.435,.318)]  # R -> G -> B
        cmap = LinearSegmentedColormap.from_list(name, colors, N=200)
    elif name=='OceanSun':
        colors = [(.0, .188, .286), (.839, .157, .157),(.969, .498, 0), (.988,.749,.286),(.918,.886,.718)]
        cmap = LinearSegmentedColormap.from_list(name, colors, N=200)
    elif name=='SunOcean':
        colors = [(.0, .188, .286), (.839, .157, .157),(.969, .498, 0), (.988,.749,.286),(.918,.886,.718)]
        cmap = LinearSegmentedColormap.from_list(name, colors[::-1], N=200)
    elif name=='RedBlue':
        colors = [(.792, .0, .125), (.957, .647, .511),(.969, .969, .969), (.573,.773,.871),(.024,.439,.690)]  # R -> G -> B
        cmap = LinearSegmentedColormap.from_list(name, colors, N=200)
    elif name=='BlueRed':
        colors = [(.792, .0, .125), (.957, .647, .511),(.969, .969, .969), (.573,.773,.871),(.024,.439,.690)]  # R -> G -> B
        cmap = LinearSegmentedColormap.from_list(name, colors[::-1], N=200)
    elif name=='PurpleOrange':
        colors = [(.369, .235, .60), (.698, .671, .824),(.969, .969, .969), (.992,.722,.388),(.902,.380,.004)]  # R -> G -> B
        cmap = LinearSegmentedColormap.from_list(name, colors, N=200)
    elif name=='SeaLand':
        colors = [(.004, .522, .443), (.502, .804, .757),(.969, .969, .969), (.875,.761,.490),(.651,.380,.102)]  # R -> G -> B
        cmap = LinearSegmentedColormap.from_list(name, colors, N=200)
    elif name=='Reds':
        colors = [(.996, .941, .851), (.992, .800, .541),(.988, .553, .349), (.843,.188,.122)]  # R -> G -> B
        cmap = LinearSegmentedColormap.from_list(name, colors, N=200)
    return cmap
