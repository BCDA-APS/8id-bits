import numpy as np
from scipy.optimize import least_squares
rng=np.random.default_rng(1)
E=15.209; LAM=12.398419/E; PX=0.076; HU,HV=1024*PX/2,512*PX/2
A=5.411651
tt=lambda s:2*np.degrees(np.arcsin(2*np.pi/(A/np.sqrt(s))*LAM/(4*np.pi)))
RINGS={"111":tt(3),"200":tt(4),"220":tt(8),"311":tt(11)}

def two_theta(u,v,L,c,v0):
    """detector face normal to the central ray at 2theta=c; returns 2theta per pixel"""
    cr,sr=np.cos(np.radians(c)),np.sin(np.radians(c))
    x=L*sr+u*cr; y=v-v0; z=L*cr-u*sr
    return np.degrees(np.arccos(z/np.sqrt(x*x+y*y+z*z)))

def ring_pixels(t,L,c,v0,n=300):
    """pixels on the detector whose 2theta == t (scan v, solve u)"""
    out=[]
    for v in np.linspace(-HV,HV,n):
        f=lambda u: two_theta(u,v,L,c,v0)-t
        a,b=-HU,HU
        if f(a)*f(b)>0: continue
        for _ in range(60):
            m=0.5*(a+b)
            if f(a)*f(m)<=0: b=m
            else: a=m
        out.append((0.5*(a+b),v))
    return np.array(out)

def run(L,c,ringset,noise_px=0.2,fit_tilt=False,label=""):
    pts=[]
    for nm in ringset:
        p=ring_pixels(RINGS[nm],L,c,0.0)
        if len(p)>=15: pts.append((RINGS[nm],p+rng.normal(0,noise_px*PX,p.shape)))
    if len(pts)<1: return print(f"  {label}: no rings on detector")
    def resid(q):
        Lf,cf,v0=q[:3]
        return np.concatenate([two_theta(p[:,0],p[:,1],Lf,cf,v0)-t for t,p in pts])
    s=least_squares(resid,[L*1.03,c*1.02,0.8])
    J=s.jac; dof=max(len(s.fun)-3,1)
    cov=np.linalg.inv(J.T@J)*(s.fun@s.fun)/dof; e=np.sqrt(np.diag(cov))
    print(f"  {label:26s} rings_on_det={len(pts)}  L={s.x[0]:7.2f}+-{e[0]:.3f}mm "
          f"({100*e[0]/L:.3f}%)  ctr={s.x[1]:6.3f}+-{e[1]:.4f}deg  "
          f"-> dq/q={100*e[0]/L:.3f}% = {3.3*e[0]/L:.5f} A-1")
    return s.x,e

print("CeO2 2theta:",{k:round(v,2) for k,v in RINGS.items()})
print("\n=== L=300 mm, centre 22.0 deg ===")
run(300,22.0,["111","200","220","311"],label="4 rings")
run(300,22.0,["200","220","311"],       label="3 rings (lose 111)")
run(300,22.0,["220","311"],             label="2 rings")
run(300,22.0,["220"],                   label="1 ring (220 only)")
print("\n=== L=350 mm, centre 23.2 deg ===")
run(350,23.2,["200","220","311"],       label="3 rings")
run(350,23.2,["220","311"],             label="2 rings")
print("\n=== degraded ring-finding (1 px scatter) ===")
run(300,22.0,["111","200","220","311"],noise_px=1.0,label="4 rings, 1px noise")
