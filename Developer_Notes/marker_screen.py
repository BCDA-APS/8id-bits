"""Which internal pressure marker has lines that never collide with the sample?"""
import numpy as np
from scipy.optimize import brentq
E=15.209; LAM=12.398419/E
tof=lambda d:2*np.degrees(np.arcsin(np.clip(LAM/(2*d),-1,1))) if LAM/(2*d)<=1 else 999
WLO,WHI=16.66,29.34          # detector window
def vv0(P,K0,Kp):
    f=lambda x:1.5*K0*(x**(-7/3)-x**(-5/3))*(1+0.75*(Kp-4)*(x**(-2/3)-1))-P
    return brentq(f,0.3,1.0)
def lines(a0,K0,Kp,struct,P):
    a=a0*vv0(P,K0,Kp)**(1/3)
    hkls={"fcc":[(1,1,1),(2,0,0),(2,2,0)],"bcc":[(1,1,0),(2,0,0),(2,1,1)],
          "b1":[(2,0,0),(2,2,0),(2,2,2)]}[struct]
    return {"".join(map(str,h)):tof(a/np.sqrt(sum(x*x for x in h))) for h in hkls}

# sample occupancy across 0-30 GPa: FCC a0=3.591 (paper ambient) + BCC a0=2.871
samp=[]
for P in range(0,31,2):
    samp+=[v for v in lines(3.591,250,4,"fcc",P).values() if WLO<v<WHI]
    samp+=[v for v in lines(2.871,250,4,"bcc",P).values() if WLO<v<WHI]
SLO,SHI=min(samp)-0.4,max(samp)+0.4
print(f"SAMPLE occupies {min(samp):.2f}-{max(samp):.2f} deg over 0-30 GPa")
print(f"  -> forbidden band (with 0.4 deg guard): {SLO:.2f}-{SHI:.2f}")
print(f"  clean zones in window: {WLO:.2f}-{SLO:.2f}  and  {SHI:.2f}-{WHI:.2f}\n")

MK=[("Au  gold",     4.0782,167,5.79,"fcc"),("Pt  platinum",3.9231,277,5.08,"fcc"),
    ("Ag  silver",   4.0853,101,6.12,"fcc"),("Cu  copper",  3.6149,133,5.01,"fcc"),
    ("Ni  nickel",   3.5240,183,4.30,"fcc"),("Pd  palladium",3.8907,195,5.35,"fcc"),
    ("Pb  lead",     4.9508, 46,5.50,"fcc"),("Al  aluminium",4.0495, 73,4.40,"fcc"),
    ("Mo  molybdenum",3.1470,268,4.00,"bcc"),("W   tungsten",3.1652,310,4.00,"bcc"),
    ("Ta  tantalum", 3.3058,194,3.50,"bcc"),("Nb  niobium", 3.3008,170,4.00,"bcc"),
    ("NaCl B1",      5.6402, 24,5.00,"b1"), ("MgO",         4.2117,160,4.15,"b1")]
print("  marker          lines that stay CLEAN and in-window over 0-30 GPa")
for nm,a0,K0,Kp,st in MK:
    keep=[]
    for h in lines(a0,K0,Kp,st,0):
        ts=[lines(a0,K0,Kp,st,P)[h] for P in range(0,31,2)]
        if all(WLO+0.3<t<WHI-0.3 and not (SLO<t<SHI) for t in ts):
            keep.append(f"{h}({min(ts):.1f}-{max(ts):.1f})")
    print(f"  {nm:16s} {len(keep)}  {'  '.join(keep) if keep else '-- none --'}")

print("\n=== ALL in-window lines, and which ones foul the sample ===")
print("  marker         line   2theta over 0-30 GPa   status")
for nm,a0,K0,Kp,st in [("Au  gold",4.0782,167,5.79,"fcc"),("Pt  platinum",3.9231,277,5.08,"fcc"),
                       ("Ta  tantalum",3.3058,194,3.50,"bcc"),("W   tungsten",3.1652,310,4.00,"bcc"),
                       ("Mo  molybdenum",3.1470,268,4.00,"bcc")]:
    any_in=False
    for h in lines(a0,K0,Kp,st,0):
        ts=[lines(a0,K0,Kp,st,P)[h] for P in range(0,31,1)]
        vis=[t for t in ts if WLO<t<WHI]
        if not vis: continue
        any_in=True
        foul=any(SLO<t<SHI for t in vis)
        frac=100*len(vis)/len(ts)
        st2=("** FOULS SAMPLE **" if foul else
             f"clean, visible {frac:.0f}% of range")
        print(f"  {nm:15s} {h:4s}  {min(vis):5.2f}-{max(vis):5.2f}   {st2}")
    if not any_in: print(f"  {nm:15s} --  no lines in window")

print("\n=== pressure sensitivity of the clean line (5 GPa step) ===")
for nm,a0,K0,Kp,st,h in [("Au 111",4.0782,167,5.79,"fcc","111"),
                         ("Ta 110",3.3058,194,3.50,"bcc","110"),
                         ("W  110",3.1652,310,4.00,"bcc","110"),
                         ("Pt 111",3.9231,277,5.08,"fcc","111")]:
    d0=lines(a0,K0,Kp,st,0)[h]-lines(a0,K0,Kp,st,5)[h]
    d25=lines(a0,K0,Kp,st,25)[h]-lines(a0,K0,Kp,st,30)[h]
    print(f"  {nm}: {abs(d0):.4f} deg at low P, {abs(d25):.4f} deg at high P"
          f"   ({abs(d25)/np.degrees(0.076/350):.1f} px worst case)")
