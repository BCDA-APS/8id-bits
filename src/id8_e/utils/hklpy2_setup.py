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


    psic.add_reflection((0, 0, 1), {"delta": 10, "chi": 90, "phi": 0, "eta": 5, "mu": 0, "nu": 0}, name="r001")  
    psic.add_reflection((2, 0, 0), {"delta": 10, "chi": 180, "phi": 0, "eta": 5, "mu": 0, "nu": 0}, name="r110")
    psic.add_reflection((2, 0, 0), {"delta": 23.9832, "chi": 92.9999, "phi": 0.0322, "eta": 14, "mu": 8.6242, "nu": 17.8778}, name="r112")
    psic.add_reflection((0, 0, 2), {"delta": 24.0537, "chi": 93, "phi": 0, "eta": 11.2169, "mu": 0, "nu": 0}, name="r002")
    ub_matrix = psic.core.calc_UB("r001", "r200")
    psic.core.calc_UB("r002", "r112")


    reals: mu=-0.0004, eta=11.2169, chi=92.9999, phi=-0.001, nu=0, delta=24.0537
