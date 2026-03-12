from hklpy2.user import add_sample
from hklpy2.user import set_diffractometer


def configure_hklpy2(oregistry):
    psic = oregistry["psic"]
    set_diffractometer(psic)

    undulator = oregistry["undulator"]
    psic.beam.wavelength.put(12.398 / undulator.upstream.value)
    try:
        add_sample("HZO", a=5.22)
        
    except Exception:
        pass


    # psic.add_reflection((1, 1, 1), {"delta": 26.7, "chi": 90.5, "phi": -3.351, "eta": 13.69, "mu": -2.098, "nu": 2.979}, name="r111")  
    # psic.add_reflection((2, 0, 0), {"delta": 17.35, "chi": 90, "phi": 8.25, "eta": 12.9, "mu": -2.648, "nu": 27.45}, name="r200")
    # ub_matrix = psic.core.calc_UB("r111", "r200")
