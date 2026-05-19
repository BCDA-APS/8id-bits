from hklpy2.user import add_sample
from hklpy2.user import set_diffractometer


def configure_hklpy2(oregistry):
    psic = oregistry["psic"]
    set_diffractometer(psic)

    # undulator = oregistry["undulator"]
    psic.beam.wavelength.put(12.38/15.200)
    try:
        add_sample("STO", a=3.905)
        
    except Exception:
        pass


    # psic.add_reflection((0, 0, 1), {"delta": 10, "chi": 90, "phi": 0, "eta": 5, "mu": 0, "nu": 0}, name="r001")  
    # psic.add_reflection((2, 0, 0), {"delta": 10, "chi": 180, "phi": 0, "eta": 5, "mu": 0, "nu": 0}, name="r110")
    #psic.add_reflection((2, 0, 0), {"delta": 23.9832, "chi": 92.9999, "phi": 0.0322, "eta": 14, "mu": 8.6242, "nu": 17.8778}, name="r112")
    # psic.add_reflection((0, 0, 2), {"delta": 24.0537, "chi": 93, "phi": 0, "eta": 11.2169, "mu": 0, "nu": 0}, name="r002")
    # psic.add_reflection((0.5, 0.5,  1.5 ), {"delta": 18.1616, "chi": 93.1198, "phi": -0.0132, "eta": 10.5397, "mu": 3.2674, "nu": 8.9182}, name="r0.50.51.5_8k")

    psic.add_reflection((0, 0, 2), {"delta": 24.0497, "chi": 93.12, "phi": 0, "eta": 11.08, "mu": 0, "nu": 0}, name="r002_8k")
    psic.add_reflection((1, 1, 2), {"delta": 24.160, "chi": 93.12, "phi": 0, "eta": 13.96, "mu": 8.910, "nu": 17.719}, name="r112")
    # ub_matrix = psic.core.calc_UB("r112", "r002_8k")

    # #9eb8417    0.000    0.000    2.000     0.000    11.080    93.120    -0.016    -0.003    24.140   second
    # #62a748e    1.000    1.000    2.000     8.910    13.960    93.120     0.007    17.719    24.160   first


    # reals: mu=-0.0004, eta=11.2169, chi=92.9999, phi=-0.001, nu=0, delta=24.0537
