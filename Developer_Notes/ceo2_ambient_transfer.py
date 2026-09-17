"""CeO2 can only be measured at ambient, on a different mount from the DAC.
So the standard is NOT at the sample position. Question: does the two-stage
scheme (CeO2 for map shape + gold for absolute scale) survive that?"""
import numpy as np
from scipy.optimize import least_squares
rng=np.random.default_rng(7)
E=15.209; LAM=12.398419/E; PX=0.076; HU,HV=1024*PX/2,512*PX/2
CE=5.411651; AU=4.0160          # Au at 8.8 GPa
tof=lambda d:2*np.degrees(np.arcsin(LAM/(2*d)))
qof=lambda t:4*np.pi*np.sin(np.radians(t)/2)/LAM

def tth(u,v,L,c,v0,dz=0.0):
    cr,sr=np.cos(np.radians(c)),np.sin(np.radians(c))
    x=L*sr+u*cr; y=v-v0; z=L*cr-u*sr-dz
    return np.degrees(np.arccos(z/np.sqrt(x*x+y*y+z*z)))

def ringpix(t,L,c,dz):
    out=[]
    for v in np.linspace(-HV,HV,160):
        f=lambda u: tth(np.array([u]),np.array([v]),L,c,0.0,dz)[0]-t
        a,b=-HU,HU
        if f(a)*f(b)>0: continue
        for _ in range(55):
            m=.5*(a+b)
            if f(a)*f(m)<=0: b=m
            else: a=m
        out.append((.5*(a+b),v))
    return np.array(out)

L,C=350.0,23.0
uu,vv=np.meshgrid(np.linspace(-HU,HU,80),np.linspace(-HV,HV,40))
ug,vg=uu.ravel(),vv.ravel()
qtrue=qof(tth(ug,vg,L,C,0.0,0.0))        # what the SAMPLE (dz=0) really sees

print("  CeO2   stage1: CeO2-only calib   stage2: + gold transfer from Lambda2M")
print("  dz(mm)   L_fit    max|dq|         L_corr    max|dq|      residual P bias")
for dz in [0.0,0.5,1.0,2.0,5.0]:
    pts=[]
    for s_ in (4,8,11):
        p=ringpix(tof(CE/np.sqrt(s_)),L,C,dz)
        if len(p)>=15: pts.append((CE/np.sqrt(s_),p+rng.normal(0,0.2*PX,p.shape)))
    r=lambda q:np.concatenate([tth(x[:,0],x[:,1],q[0],q[1],q[2])-tof(d) for d,x in pts])
    f1=least_squares(r,[L*1.02,C*1.01,0.5]).x
    e1=np.abs(qof(tth(ug,vg,*f1))-qtrue).max()
    # stage 2: one scalar on L so Au111 lands where Lambda2M says it does
    t_au=tof(AU/np.sqrt(3))
    g=lambda k: (np.interp(0,[0],[0]),)  # placeholder
    def mis(k):
        # radius of Au111 as the (biased) calib would place it vs truth
        rt=L*np.tan(np.radians(t_au))
        return np.degrees(np.arctan(rt/(f1[0]*k[0])))-t_au
    k=least_squares(mis,[1.0]).x[0]
    f2=(f1[0]*k,f1[1],f1[2])
    e2=np.abs(qof(tth(ug,vg,*f2))-qtrue).max()
    dP=abs(e2)/3.3*3*167.0
    print(f"  {dz:5.2f}   {f1[0]:7.2f}   {e1:.5f} A-1   {f2[0]:7.2f}   {e2:.5f} A-1   {dP:5.2f} GPa")
print("\n  (5 GPa gold shift = 0.0155 A-1; sample peak FWHM = 0.0495 A-1)")

print("\n=== stage 2b: transfer using BOTH gold lines (111 AND 200) ===")
print("  a_Au taken from Lambda2M; two lines -> two constraints -> fix L and centre\n")
print("  dz(mm)   L_corr    ctr_corr    max|dq|      residual P bias")
for dz in [0.0,0.5,1.0,2.0,5.0]:
    pts=[]
    for s_ in (4,8,11):
        p=ringpix(tof(CE/np.sqrt(s_)),L,C,dz)
        if len(p)>=15: pts.append((CE/np.sqrt(s_),p+rng.normal(0,0.2*PX,p.shape)))
    r=lambda q:np.concatenate([tth(x[:,0],x[:,1],q[0],q[1],q[2])-tof(d) for d,x in pts])
    f1=least_squares(r,[L*1.02,C*1.01,0.5]).x
    # both gold rings, as they truly appear (sample at dz=0)
    gold=[]
    for m in (np.sqrt(3),2.0):
        t=tof(AU/m); p=ringpix(t,L,C,0.0)
        if len(p)>=15: gold.append((t,p))
    def mis(q):
        return np.concatenate([tth(p[:,0],p[:,1],q[0],q[1],f1[2])-t for t,p in gold])
    f2=least_squares(mis,[f1[0],f1[1]]).x
    e2=np.abs(qof(tth(ug,vg,f2[0],f2[1],f1[2]))-qtrue).max()
    print(f"  {dz:5.2f}   {f2[0]:7.2f}   {f2[1]:7.3f}    {e2:.5f} A-1   {abs(e2)/3.3*3*167:5.2f} GPa")
print("\n  -> the gold PAIR removes the CeO2 position error almost completely.")
print("     Single gold line cannot: one constraint, two biased parameters.")
