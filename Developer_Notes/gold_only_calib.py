"""Can Au 111 + Au 200 alone calibrate the detector, with no CeO2?
Difference from CeO2: the gold 2theta are NOT known a priori -- a_Au depends on
the unknown pressure. Only the RATIO q200/q111 = 2/sqrt(3) is known."""
import numpy as np
from scipy.optimize import least_squares
rng=np.random.default_rng(3)
E=15.209; LAM=12.398419/E; PX=0.076; HU,HV=1024*PX/2,512*PX/2
CE=5.411651
tof=lambda d:2*np.degrees(np.arcsin(LAM/(2*d)))

def tth(u,v,L,c,v0):
    cr,sr=np.cos(np.radians(c)),np.sin(np.radians(c))
    x=L*sr+u*cr; y=v-v0; z=L*cr-u*sr
    return np.degrees(np.arccos(z/np.sqrt(x*x+y*y+z*z)))

def ringpix(t,L,c,n=200):
    out=[]
    for v in np.linspace(-HV,HV,n):
        f=lambda u: tth(np.array([u]),np.array([v]),L,c,0.0)[0]-t
        a,b=-HU,HU
        if f(a)*f(b)>0: continue
        for _ in range(55):
            m=.5*(a+b)
            if f(a)*f(m)<=0: b=m
            else: a=m
        out.append((.5*(a+b),v))
    return np.array(out)

L,C=350.0,23.0
def mkpts(dspacings,noise=0.2):
    pts=[]
    for d in dspacings:
        t=tof(d); p=ringpix(t,L,C)
        if len(p)>=15: pts.append((d,p+rng.normal(0,noise*PX,p.shape)))
    return pts

print("=== A: CeO2, 3 rings, 2theta KNOWN absolutely ===")
pts=mkpts([CE/np.sqrt(s) for s in (4,8,11)])
r=lambda p: np.concatenate([tth(x[:,0],x[:,1],p[0],p[1],p[2])-tof(d) for d,x in pts])
s=least_squares(r,[L*1.02,C*1.01,0.5]); J=s.jac
e=np.sqrt(np.diag(np.linalg.inv(J.T@J)*(s.fun@s.fun)/(len(s.fun)-3)))
print(f"   L={s.x[0]:7.2f}+-{e[0]:.3f} mm ({100*e[0]/L:.3f}%)  -> dq/q={100*e[0]/L:.3f}%")

print("\n=== B: gold only, Au111+Au200, a_Au FREE (pressure unknown) ===")
for aTrue,Pl in [(4.0782,"0 GPa"),(3.9543,"~20 GPa")]:
    pts=mkpts([aTrue/np.sqrt(3),aTrue/2.0])
    def r2(p):
        L_,c_,v0_,a_=p
        return np.concatenate([tth(x[:,0],x[:,1],L_,c_,v0_)-tof(a_/np.sqrt(3) if abs(d-aTrue/np.sqrt(3))<1e-9 else a_/2.0)
                               for d,x in pts])
    s2=least_squares(r2,[L*1.02,C*1.01,0.5,aTrue*1.01]); J2=s2.jac
    try:
        cov=np.linalg.inv(J2.T@J2)*(s2.fun@s2.fun)/max(len(s2.fun)-4,1)
        e2=np.sqrt(np.diag(cov)); corr=cov[0,3]/np.sqrt(cov[0,0]*cov[3,3])
        print(f"   {Pl:8s} L={s2.x[0]:8.2f}+-{e2[0]:8.2f} mm ({100*e2[0]/L:7.3f}%)"
              f"   a_Au={s2.x[3]:.4f}+-{e2[3]:.4f} A   corr(L,a_Au)={corr:+.4f}")
    except Exception as ex: print("   singular:",ex)

print("\n=== C: with L FIXED by CeO2, do the two gold lines AGREE? ===")
print("   (CeO2 sits at z-offset dz, so the fitted L is wrong by ~dz.")
print("    Each gold line then gives its own a_Au. Disagreement = the alarm.)\n")
def a_from(t,mult):     # mult = sqrt(3) for 111, 2 for 200
    return mult*LAM/(2*np.sin(np.radians(t)/2))
print("   dz(mm)   L_used    a(from 111)  a(from 200)   split      P error")
aT=4.0160          # 8.8 GPa
for dz in [0.0,0.25,0.5,1.0,2.0]:
    Lu=L-dz*np.cos(np.radians(C))          # CeO2 z-offset propagates into fitted L
    out=[]
    for mult,d in ((np.sqrt(3),aT/np.sqrt(3)),(2.0,aT/2.0)):
        rtrue=L*np.tan(np.radians(tof(d)))            # true ring radius
        tapp=np.degrees(np.arctan(rtrue/Lu))          # apparent 2theta with wrong L
        out.append(a_from(tapp,mult))
    split=out[1]-out[0]
    # convert the a-split into an equivalent pressure error (Au, dP/da ~ -3K/a)
    dP=abs(split)*3*167.0/aT
    print(f"   {dz:5.2f}   {Lu:7.2f}   {out[0]:.5f}     {out[1]:.5f}   "
          f"{split:+.5f} A   {dP:5.2f} GPa")
print("\n   -> the 111/200 split is a DIRECT readout of the calibration error")
print("      that CeO2 alone cannot see (sec 4.4).")
