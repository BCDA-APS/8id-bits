import math

dphi = 2
ctr = -0.224
pls = -0.222
mis = -0.223
dx = (mis-pls)/(2*math.sin(dphi/180.0*math.pi))
print('move z by ',dx)