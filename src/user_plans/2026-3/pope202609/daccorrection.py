import math

dphi = 2
ctr = -0.1424
pls = -0.1447
mis = -0.1399
dx = (mis-pls)/(2*math.sin(dphi/180.0*math.pi))
print('move z by ',dx)