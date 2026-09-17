import numpy as np
E=15.209; LAM=12.398419/E
A_CEO2=5.411651                      # NIST SRM 674b
PX=76e-6; NLONG,NSHORT=1024,512
HALF_L=NLONG*PX/2; HALF_S=NSHORT*PX/2
TTH_DAC=30.0                         # DAC full aperture -> 2theta max

q=lambda t: 4*np.pi*np.sin(np.radians(t)/2)/LAM
tth=lambda Q: 2*np.degrees(np.arcsin(np.clip(Q*LAM/(4*np.pi),-1,1)))

print(f"E={E} keV  lambda={LAM:.5f} A")
print(f"DAC 2th_max={TTH_DAC} deg -> q_max={q(TTH_DAC):.3f} A-1")
print(f"detector active: {NLONG*PX*1e3:.1f} x {NSHORT*PX*1e3:.1f} mm\n")

print("=== CeO2 rings within the DAC cone ===")
rings=[]
for hkl in [(1,1,1),(2,0,0),(2,2,0),(3,1,1),(2,2,2),(4,0,0)]:
    s=sum(h*h for h in hkl); d=A_CEO2/np.sqrt(s); Q=2*np.pi/d; t=tth(Q)
    ok = t<=TTH_DAC
    if ok: rings.append((hkl,t,Q))
    print(f"  {hkl}  d={d:.4f} A  q={Q:.4f} A-1  2th={t:7.3f} deg  {'' if ok else '<- outside DAC'}")

print("\n=== coverage vs distance (LONG axis along 2theta) ===")
print("   L(m)   span(deg)  pix d(2th)   dq/pix(A-1) @2th=22")
for L in [0.50,0.45,0.40,0.35,0.30,0.25,0.20]:
    span=2*np.degrees(np.arctan(HALF_L/L))
    dt=np.degrees(PX/L)
    dqdp=2*np.pi*np.cos(np.radians(11))/LAM*np.radians(dt)
    print(f"  {L:5.2f}   {span:7.2f}    {dt:.4f}      {dqdp:.5f}")
print(f"  (short axis along 2th at 0.5 m would give only "
      f"{2*np.degrees(np.arctan(HALF_S/0.5)):.2f} deg -- do not do this)")

print("\n=== how many CeO2 rings land on the detector? ===")
print("  L(m)  center  2th range        q range        CeO2 rings on-detector")
best=None
for L in [0.50,0.40,0.35,0.30,0.25]:
    half=np.degrees(np.arctan(HALF_L/L))
    for c in np.arange(16,26.01,0.5):
        lo,hi=c-half,c+half
        if hi>TTH_DAC+1: continue
        on=[r for r in rings if lo<=r[1]<=hi]
        if best is None or len(on)>best[0]: best=(len(on),L,c,lo,hi,on)
    # report the best centering for this L
    cands=[]
    for c in np.arange(14,26.01,0.25):
        lo,hi=c-half,c+half
        if hi>TTH_DAC+1: continue
        on=[r for r in rings if lo<=r[1]<=hi]
        cands.append((len(on),c,lo,hi,on))
    n,c,lo,hi,on=max(cands,key=lambda x:(x[0],-abs(x[1]-22)))
    names="  ".join("".join(map(str,r[0])) for r in on)
    print(f"  {L:4.2f}  {c:5.1f}   {lo:5.2f}-{hi:5.2f}   {q(max(lo,0)):.2f}-{q(hi):.2f}   [{n}]  {names}")

print("\n=== do the science peaks land in the recommended window? ===")
def bm3(P,K0=167.0,Kp=5.79):          # Au, Anderson 1989
    from scipy.optimize import brentq
    f=lambda x: 1.5*K0*(x**(-7/3)-x**(-5/3))*(1+0.75*(Kp-4)*(x**(-2/3)-1))-P
    return brentq(f,0.5,1.0)
A_AU=4.0782
L,C=0.30,22.0
half=np.degrees(np.arctan(HALF_L/L)); lo,hi=C-half,C+half
print(f"  geometry: L={L} m, center 2th={C} deg  ->  {lo:.2f}-{hi:.2f} deg,"
      f"  q {q(lo):.2f}-{q(hi):.2f} A-1\n")
print("   P(GPa)   Au a(A)   Au111 q   2th     on detector?")
for P in [0,8.8,10.5,20.2,28.5]:
    a=A_AU*bm3(P)**(1/3); Q=2*np.pi/(a/np.sqrt(3)); t=tth(Q)
    print(f"   {P:5.1f}   {a:.4f}   {Q:.4f}   {t:6.2f}   {'YES' if lo<=t<=hi else 'no'}")

print("\n   sample peaks seen this beamtime:")
for label,Q in [("HEA low-q (2th~22 deg)",q(22.0)),
                ("HEA high-q (2th~25 deg)",q(25.0)),
                ("vanished line d=2.006 A",2*np.pi/2.006)]:
    t=tth(Q); print(f"     {label:26s} q={Q:.4f}  2th={t:6.2f}  {'YES' if lo<=t<=hi else 'no'}")

print("\n=== error budget after CeO2 calibration ===")
print("  source                         effect on q at q=3 A-1")
for lbl,rel in [("uncalibrated L off by 1% (their slide)",0.0097),
                ("ring centroid, 0.1 px on r=121 mm",0.076*0.1/121),
                ("L fitted from 4 rings (~0.05%)",5e-4),
                ("CeO2 vs DAC sample z-offset 1.0 mm",1.0/(L*1e3)),
                ("CeO2 vs DAC sample z-offset 0.1 mm",0.1/(L*1e3))]:
    print(f"  {lbl:38s} {rel*3:.4f} A-1   ({rel*100:.3f}%)")
print(f"\n  one pixel at this geometry           {2*np.pi*np.cos(np.radians(11))/LAM*PX/L:.5f} A-1")

print("\n=== azimuthal arc sampled (matters: DAC rings are spotty) ===")
for L in [0.50,0.30,0.25]:
    for c in [22.0]:
        r=L*np.tan(np.radians(c))*1e3
        arc=2*np.degrees(np.arctan((HALF_S*1e3)/r))
        print(f"  L={L:.2f} m  ring radius {r:6.1f} mm  azimuth covered = {arc:5.1f} deg"
              f"  ({arc/360*100:.1f}% of the ring)")
print("  (this beamtime: gold spottiness 5.18, sample <=0.14 -- gold is the")
print("   grainy one, so gold centroids from a short arc are the least reliable)")

print("\n=== physical footprint at the recommended geometry ===")
for L,c in [(0.50,21.8),(0.30,22.0),(0.25,22.0)]:
    r=L*np.tan(np.radians(c))*1e3
    print(f"  L={L:.2f} m: detector centre {r:5.1f} mm off-axis, "
          f"spans r = {r-HALF_L*1e3:5.1f} to {r+HALF_L*1e3:5.1f} mm")
print("\n=== relative counts per pixel (solid angle), 0.5 m = 1.0 ===")
for L in [0.50,0.40,0.30,0.25]:
    print(f"  L={L:.2f} m -> {(0.5/L)**2:.2f}x")
print("\n=== sensitivity to sample z-placement (the cost of moving closer) ===")
for L in [0.50,0.30,0.25]:
    print(f"  L={L:.2f} m: 1 mm z-offset -> {1.0/(L*1e3)*100:.3f}% in q "
          f"= {1.0/(L*1e3)*3:.4f} A-1 at q=3")
