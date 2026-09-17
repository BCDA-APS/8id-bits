import numpy as np
from scipy.optimize import least_squares
rng=np.random.default_rng(2)
E=15.209; LAM=12.398419/E; PX=0.076; HU,HV=1024*PX/2,512*PX/2
A=5.411651
tt=lambda s:2*np.degrees(np.arcsin(2*np.pi/(A/np.sqrt(s))*LAM/(4*np.pi)))
RINGS=[tt(3),tt(4),tt(8),tt(11)]
qof=lambda t:4*np.pi*np.sin(np.radians(t)/2)/LAM

def tth(u,v,L,c,v0,tu=0.0,tv=0.0,dz=0.0):
    """dz = sample displaced along beam (CeO2 vs DAC). tu,tv = detector tilts (deg)."""
    cr,sr=np.cos(np.radians(c)),np.sin(np.radians(c))
    n =np.array([sr,0,cr]); eu=np.array([cr,0,-sr]); ev=np.array([0,1.0,0])
    for ang,ax in ((tu,ev),(tv,eu)):
        if ang:
            a=np.radians(ang); K=np.array([[0,-ax[2],ax[1]],[ax[2],0,-ax[0]],[-ax[1],ax[0],0]])
            R=np.eye(3)+np.sin(a)*K+(1-np.cos(a))*K@K
            n,eu,ev=R@n,R@eu,R@ev
    P=L*n[None,:]+np.outer(u,eu)+np.outer(v-v0,ev); P[:,2]-=dz
    return np.degrees(np.arccos(P[:,2]/np.linalg.norm(P,axis=1)))

def ringpix(t,L,c,**kw):
    out=[]
    for v in np.linspace(-HV,HV,200):
        f=lambda u: tth(np.array([u]),np.array([v]),L,c,0.0,**kw)[0]-t
        a,b=-HU,HU
        if f(a)*f(b)>0: continue
        for _ in range(50):
            m=.5*(a+b); (b:=m) if f(a)*f(m)<=0 else (a:=m)
        out.append((.5*(a+b),v))
    return np.array(out)

def calibrate(L,c,truth):
    pts=[]
    for t in RINGS:
        p=ringpix(t,L,c,**truth)
        if len(p)>=15: pts.append((t,p+rng.normal(0,0.2*PX,p.shape)))
    def resid(q):
        return np.concatenate([tth(p[:,0],p[:,1],q[0],q[1],q[2])-t for t,p in pts])
    s=least_squares(resid,[L*1.02,c*1.01,0.5])
    # residual q error across the whole detector
    uu,vv=np.meshgrid(np.linspace(-HU,HU,60),np.linspace(-HV,HV,30))
    u,v=uu.ravel(),vv.ravel()
    # TRANSFER the calibration to the real sample: tilts stay (detector is fixed),
    # but the sample sits at dz=0, not where the CeO2 standard was.
    tr_s={k:v2 for k,v2 in truth.items() if k!="dz"}
    true=tth(u,v,L,c,0.0,**tr_s); got=tth(u,v,s.x[0],s.x[1],s.x[2])
    dq=qof(got)-qof(true)
    return s.x,len(pts),np.abs(dq).max(),np.std(dq)

L,C=300.0,22.0
print("Fitting only (L, centre, v-offset) -- tilt assumed perfect, sample assumed co-located\n")
print("  truth                        Lfit(mm)  max |dq| (A-1)   rms dq    verdict")
for lbl,tr in [("perfect",                     {}),
               ("detector tilt 0.5 deg",       {"tu":0.5}),
               ("detector tilt 1.0 deg",       {"tu":1.0}),
               ("detector tilt 2.0 deg",       {"tu":2.0}),
               ("detector tilt 5.0 deg",       {"tu":5.0}),
               ("out-of-plane tilt 1.0 deg",   {"tv":1.0}),
               ("CeO2 0.5 mm off in z",        {"dz":0.5}),
               ("CeO2 1.0 mm off in z",        {"dz":1.0}),
               ("CeO2 2.0 mm off in z",        {"dz":2.0})]:
    x,n,mx,rms=calibrate(L,C,tr)
    v="ok" if mx<0.005 else ("marginal" if mx<0.015 else "NOT OK")
    print(f"  {lbl:28s} {x[0]:7.2f}   {mx:.5f}      {rms:.5f}   {v}")
print(f"\n  (peak FWHM = 0.0495 A-1; 5 GPa gold shift = 0.0155 A-1 at q=3.3)")
