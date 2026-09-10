from hklpy2.user import add_sample
from hklpy2.user import set_diffractometer


def configure_hklpy2(oregistry):
    psic = oregistry["psic"]
    set_diffractometer(psic)

    # undulator = oregistry["undulator"]
    psic.beam.wavelength.put(12.38/8.800)
    
    add_sample("BTO", a=3.84, b=3.84, c=32.83, replace=True)

    psic.add_reflection((0, 0, 8), {"delta": 19.85, "chi": 90.445, "phi": 0, "eta": 9.925, "mu": 0, "nu": 0}, 
                        name="r008")
    # psic.add_reflection((1, 1, 8), {"delta": 16.63, "chi": 90.445, "phi": 42.32, "eta": 13, "mu": 0, "nu": 26.92}, 
    #                     name="r118")
    psic.add_reflection((1, 1, 9), {"delta": 15.353, "chi": 90.445, "phi": 35.67, "eta": 14.278, "mu": 0, "nu": 31.789}, 
                        name="r119")
    psic.add_reflection((1, 1, 13), {"delta": 30.6085, "chi": 90.355, "phi": 30.7121, "eta": 14.278, "mu": 0, "nu": 34.939}, 
                        name="r1113")
    psic.add_reflection((1.5, 0.5, 6), {"delta": 12.424, "chi": 90.45, "phi": 8.023, "eta": 14.029, "mu": 0, "nu": 35.141}, 
                        name="r1p5p56")
    psic.add_reflection((1.5, 0.5, 7), {"delta": 12.424, "chi": 90.45, "phi": 8.023, "eta": 14.029, "mu": 0, "nu": 35.141}, 
                        name="r1p5p57")
    psic.add_reflection((1.5, 0.5, 9), {"delta": 17.158, "chi": 90.45, "phi": 7.939, "eta": 14.784, "mu": 0, "nu": 35.953}, 
                        name="r1p5p59")
    
    # ub_matrix = psic.core.calc_UB("r008", "r1p5p56")
    # ub_matrix = psic.core.calc_UB("r008", "r1113")

    psic.core.mode = "lifting_detector_omega"

    psic.core.constraints["delta"].limits = (0, 90)
    psic.core.constraints["eta"].limits   = (0, 45)
    psic.core.constraints["nu"].limits    = (0, 40)


    # psic.add_reflection((0, 0, 1), {"delta": 10, "chi": 90, "phi": 0, "eta": 5, "mu": 0, "nu": 0}, name="r001")  
    # psic.add_reflection((2, 0, 0), {"delta": 10, "chi": 180, "phi": 0, "eta": 5, "mu": 0, "nu": 0}, name="r110")
    #psic.add_reflection((2, 0, 0), {"delta": 23.9832, "chi": 92.9999, "phi": 0.0322, "eta": 14, "mu": 8.6242, "nu": 17.8778}, name="r112")
    # psic.add_reflection((0, 0, 2), {"delta": 24.0537, "chi": 93, "phi": 0, "eta": 11.2169, "mu": 0, "nu": 0}, name="r002")
    # psic.add_reflection((0.5, 0.5,  1.5 ), {"delta": 18.1616, "chi": 93.1198, "phi": -0.0132, "eta": 10.5397, "mu": 3.2674, "nu": 8.9182}, name="r0.50.51.5_8k")
    # #9eb8417    0.000    0.000    2.000     0.000    11.080    93.120    -0.016    -0.003    24.140   second
    # #62a748e    1.000    1.000    2.000     8.910    13.960    93.120     0.007    17.719    24.160   first


    # reals: mu=-0.0004, eta=11.2169, chi=92.9999, phi=-0.001, nu=0, delta=24.0537
