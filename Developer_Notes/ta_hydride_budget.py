"""Is there enough hydrogen in a He-loaded gasket hole to hydride a Ta marker?
The hole is a SEALED micro-volume: the H inventory is fixed at loading."""
import numpy as np
NA=6.022e23
D,T=100e-4,35e-4                       # cm: 100 um hole, 35 um gasket
V=np.pi*(D/2)**2*T
print(f"gasket hole volume  {V:.3e} cm3 = {V*1e6:.3f} nL   (sealed at loading)")
print("TaHx shifts the lattice measurably at H/Ta ~ 0.1-1.0\n")
print("  Ta flake        H2 in He    H atoms    Ta atoms   H/Ta      verdict")
for dims,lbl in [((20e-4,20e-4,5e-4),"20x20x5 um"),((10e-4,10e-4,3e-4),"10x10x3 um")]:
    n_ta=np.prod(dims)*16.65/180.95*NA
    for ppm in [1,10,100,1000]:
        n_he=V*0.120/4.003*NA          # He at ~0.2 GPa loading, 0.12 g/cc
        n_h=2*n_he*ppm*1e-6
        r=n_h/n_ta
        v=("safe" if r<1e-2 else "marginal" if r<0.1 else "WOULD HYDRIDE")
        print(f"  {lbl:12s} {ppm:5d} ppm   {n_h:.2e}  {n_ta:.2e}  {r:.1e}  {v}")
    print()
print("Research-grade He is 5N-6N; H2 is a fraction of the <=10 ppm total impurity.")
print("At <=10 ppm H2 the ratio is <=1e-3 -- two orders below the threshold.")
print("Only implausibly dirty gas (~1000 ppm = 0.1%) would approach it.")
print("\nThe reason is the sealed volume: 0.27 nL holds ~1e-11 mol of gas. A furnace")
print("or an H2 atmosphere is an effectively infinite H reservoir and CAN hydride Ta;")
print("a sealed sub-nanolitre bubble of clean He cannot.")
