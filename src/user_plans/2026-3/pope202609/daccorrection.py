import math

dphi = 2
ctr = -1.656
pls = -1.649
mis = -1.655
dx = (mis-pls)/(2*math.sin(dphi/180.0*math.pi))
print('move z by ',dx)