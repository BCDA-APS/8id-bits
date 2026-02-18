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
