import numpy as np
from scipy.optimize import brentq
E=15.209; LAM=12.398419/E; PX=76e-6; HALF_L=1024*PX/2; HALF_S=512*PX/2
A_CEO2=5.411651; TTH_DAC=30.0
q=lambda t: 4*np.pi*np.sin(np.radians(t)/2)/LAM
tth=lambda Q: 2*np.degrees(np.arcsin(np.clip(Q*LAM/(4*np.pi),-1,1)))
def vv0(P,K0,Kp):
    f=lambda x:1.5*K0*(x**(-7/3)-x**(-5/3))*(1+0.75*(Kp-4)*(x**(-2/3)-1))-P
    return brentq(f,0.3,1.0)

# --- how far does each peak move per 5 GPa? -------------------------------
print("=== peak shift per 5 GPa step (the thing we must resolve) ===")
print("  phase              P range      d(2th) per 5 GPa")
cases=[("Au 111 (K0=167,K'=5.79)", 4.0782/np.sqrt(3),167,5.79),
       ("HEA   (K0=250,K'=4)",     2*np.pi/2.9413, 250,4.0)]
shifts={}
for name,val,K0,Kp in cases:
    d0 = val if name.startswith("Au") else 2*np.pi/val
    if not name.startswith("Au"): d0 = 2*np.pi/(2*np.pi/2.9413)  # d for q=2.9413
    d0 = (4.0782/np.sqrt(3)) if name.startswith("Au") else (2*np.pi/2.9413)
    rows=[]
    for P0 in [0,5,10,15,20,25]:
        d1=d0*vv0(P0,K0,Kp)**(1/3); d2=d0*vv0(P0+5,K0,Kp)**(1/3)
        rows.append(tth(2*np.pi/d2)-tth(2*np.pi/d1))
    shifts[name]=rows
    print(f"  {name:24s} 0->30 GPa   {min(rows):.4f} (high P) .. {max(rows):.4f} (low P) deg")
worst=min(min(v) for v in shifts.values())
print(f"\n  WORST CASE (stiff sample, high pressure): {worst:.4f} deg per 5 GPa")

# --- pixels per 5 GPa vs distance -----------------------------------------
print("\n=== can we resolve it? pixels per 5 GPa step ===")
print("   L(m)  span   px/5GPa(worst)  px/5GPa(Au@0GPa)  CeO2 rings  arc sagitta")
for L in [0.25,0.30,0.35,0.40,0.45,0.50,0.60,0.70,0.80]:
    half=np.degrees(np.arctan(HALF_L/L)); span=2*half
    pxw=np.radians(worst)*L/PX
    pxa=np.radians(max(shifts["Au 111 (K0=167,K'=5.79)"]))*L/PX
    # best centering that still covers Au111(19.9) .. HEA high(25.3)
    best=0
    for c in np.arange(14,27,0.25):
        lo,hi=c-half,c+half
        if lo>19.7 or hi<25.3 or hi>TTH_DAC+1.5: continue
        n=sum(1 for s in [3,4,8,11] if lo<=tth(2*np.pi/(A_CEO2/np.sqrt(s)))<=hi)
        best=max(best,n)
    r=L*np.tan(np.radians(22))*1e3
    sag=r*(1-np.cos(np.arctan(HALF_S*1e3/r)))/ (PX*1e3)
    flag="" if best>=2 and pxw>=5 else ("  <- fails" if best<2 or pxw<5 else "")
    print(f"  {L:5.2f} {span:6.2f}     {pxw:6.1f}          {pxa:6.1f}         {best}"
          f"        {sag:5.1f} px{flag}")

print("\n=== shift detection is DIFFERENTIAL: what actually limits it ===")
print("  Calibration error (L, z-offset, tilt) CANCELS in a pressure difference")
print("  measured on the same fixed detector. What does NOT cancel:")
print("  sample displacement between pressure points (documented this beamtime).\n")
SIG=0.1078   # deg, worst-case 5 GPa shift
print("   L(m)   5GPa signal   err from 50um    err from 100um   err/signal(100um)")
for L in [0.25,0.30,0.35,0.40,0.50,0.60]:
    e50=np.degrees(50e-6/L); e100=np.degrees(100e-6/L)
    print(f"  {L:5.2f}   {SIG:.4f} deg   {e50:.4f} deg      {e100:.4f} deg      {e100/SIG*100:5.1f}%")
print("\n  -> farther is better here: the signal is fixed in degrees,")
print("     the displacement error shrinks as 1/L.")

print("\n=== what peak FWHM would make 5 GPa marginal? ===")
for crit,lbl in [(5,"obvious by eye      (shift >= FWHM/5)"),
                 (10,"resolvable by fit   (shift >= FWHM/10)")]:
    print(f"  {lbl}: need FWHM <= {SIG*crit:.2f} deg")
print(f"  CHECK THIS on scan S00234: FWHM of the 22 deg peak in huber.delta.")

print("\n=== combined score ===")
print("   L(m)  CeO2  px/5GPa  disp.err(100um)  verdict")
for L,n in [(0.25,4),(0.30,4),(0.35,3),(0.40,3),(0.45,2),(0.50,2)]:
    px=np.radians(SIG)*L/PX; e=np.degrees(100e-6/L)/SIG*100
    v=("BEST compromise" if L==0.35 else
       "best calibration" if L==0.30 else
       "best displacement rejection" if L==0.50 else "")
    print(f"  {L:5.2f}   {n}     {px:5.1f}        {e:5.1f}%        {v}")
