"""
Auxilary HKL functions.

.. autosummary::
    ~newsample
    ~sampleChange
    ~sampleRemove
    ~_sampleList
    ~list_reflections
    ~or_swap
    ~setor0
    ~setor1
    ~set_orienting
    ~del_reflection
    ~list_orienting
    ~or0
    ~or1
    ~compute_UB
    ~calc_UB
    ~setmode
    ~ca
    ~ubr
    ~br
    ~uan
    ~an
    ~_wh
    ~setlat
    ~setaz
    ~freeze
    ~freeze_psi
    ~wh
"""
"""
Provide a simplified UI for hklpy2 diffractometer users.

The user must define a diffractometer instance, then
register that instance here calling `set_diffractometer(instance)`.

FUNCTIONS
"""

__all__ = """
    newsample
    sampleChange
    sampleRemove
    _sampleList
    list_reflections
    or_swap
    setor0
    setor1
    set_orienting
    del_reflection
    list_orienting
    or0
    or1
    compute_UB
    calc_UB
    setmode
    ca
    ubr
    br
    uan
    an
    _wh
    setlat
    setaz
    freeze
    freeze_psi
""".split()

try:
    from hklpy2.user import set_diffractometer, get_diffractometer, cahkl, add_sample
    from bluesky import RunEngine, RunEngineInterrupted
    from bluesky.utils import ProgressBarManager
    import asyncio
    from bluesky.plan_stubs import mv
    from apsbits.core.instrument_init import oregistry
    import numpy as np
except ModuleNotFoundError:
    print("gi module is not installed, the hkl_utils functions will not work!")

RE = RunEngine({}, loop=asyncio.new_event_loop())
pbar_manager = ProgressBarManager()

class Geometries:
    """Register the diffractometer geometries."""

    _psic = None
    _psic_sim = None
    _psic_q2 = None
    _psic_psi = None

    def __init__(self, psic_nm):
        self.psic = oregistry.find(psic_nm)
#        self.sim = oregistry.find(sim_nm)
#        self.q2 = oregistry.find(q2_nm)
#        self.psi = oregistry.find(psi_nm)

    # @property
    # def sim(self):
    #     return self._psic_sim

    # @sim.setter
    # def sim(self, value):
    #     self._psic_sim = value

    @property
    def psic(self):
        return self._psic

    @psic.setter
    def psic(self, value):
        self._psic = value

    # @property
    # def q2(self):
    #     return self._psic_q2

    # @q2.setter
    # def q2(self, value):
    #     self._psic_q2 = value

    # @property
    # def psi(self):
    #     return self._psic_psi

    # @psi.setter
    # def psi(self, value):
    #     self._psic_psi = value


geometries = Geometries("psic")
# geometries.psic = oregistry.find("psic")
# geometries.sim = oregistry.find("psic_sim")
# geometries.q2 = oregistry.find("psic_q")
# geometries.psi = oregistry.find("psic_psi")
    
set_diffractometer(geometries.psic)

def newsample():
    """
    Interactively add a new sample with lattice constants.
    
    Prompts the user for:
      - sample name
      - lattice constants (a, b, c, alpha, beta, gamma)
    
    Uses defaults if the user just presses Enter.
    """
    # Ask for sample name
    name = input("Sample name? ").strip()
    if not name:
        print("Sample name cannot be empty.")
        return

    # Default values
    defaults = {
        "a": 5.0, "b": 5.0, "c": 5.0,
        "alpha": 90.0, "beta": 90.0, "gamma": 90.0,
    }

    def ask(name, default):
        val = input(f"{name} ({default})? ") or default
        try:
            return float(val)
        except ValueError:
            print(f"Invalid {name}, keeping {default}")
            return default

    # Collect lattice parameters
    a = ask("a", defaults["a"])
    b = ask("b", defaults["b"])
    c = ask("c", defaults["c"])
    alpha = ask("alpha", defaults["alpha"])
    beta  = ask("beta", defaults["beta"])
    gamma = ask("gamma", defaults["gamma"])

    # Call your existing function
    add_sample(name, a=a, b=b, c=c, alpha=alpha, beta=beta, gamma=gamma)

    print(f"Added sample '{name}' with lattice constants: "
          f"a={a}, b={b}, c={c}, alpha={alpha}, beta={beta}, gamma={gamma}")

def _sampleList():
    """List all samples currently defined in hklpy2; specify  current one."""
    _geom_ = get_diffractometer()
    samples = _geom_.samples
    print("")
    for nm, sample in samples.items():
#        orienting_refl = samples[x].reflections
        print(f"Sample = {nm}")
        print("Lattice:", end=" ")
#        print(*samples[x].lattice._fields, sep=", ", end=" = ")
        print(sample.lattice)
        print(
            "======================================================================"
        )
    print("\nCurrent sample: " + _geom_.sample.name)


def list_reflections(all_samples=False):
    """
    Lists all reflections defined in hklpy for six-circle geometry.

    Parameters
    ----------
    all_samples : bool, optional
        If True, list reflections for all samples.
        If False, only the current sample. Defaults to False.
    """
    _geom_ = get_diffractometer()
    samples = _geom_.samples.values() if all_samples else [_geom_.sample]

    for sample in samples:
        print(f"Sample: {sample.name}")

        orienting_refl = sample.reflections.order  # orientation reflection keys

        # column widths
        refl_width = 8
        pseudo_width = 9   # for h, k, l
        real_width = 10    # for motor positions

        # Header
        header = (
            f"\n{'#':>{refl_width}}"
            + "".join(f"{m:>{pseudo_width}}" for m in _geom_.pseudo_positioners._fields).upper()
            + "".join(f"{k:>{real_width}}" for k in _geom_.real_positioners._fields)
            + "   orienting"
        )
        print(header)

        # Rows
        for key, ref in sample.reflections.items():
            h, k, l = list(ref.pseudos.values())
            pos = list(ref.reals.values())

            tag = ""
            if orienting_refl and key == orienting_refl[0]:
                tag = "first"
            elif len(orienting_refl) > 1 and key == orienting_refl[1]:
                tag = "second"

            row = (
                f"{key:>{refl_width}}"
                f"{h:{pseudo_width}.3f}{k:{pseudo_width}.3f}{l:{pseudo_width}.3f}"
                + "".join(f"{v:{real_width}.3f}" for v in pos)
                + (f"   {tag}" if tag else "")
            )
            print(row)

        if len(samples) > 1 and all_samples:
            print("=" * 107)



def sampleChange(sample_key=None):
    """
    Change selected sample in hklpy2.

    Parameters
    ----------
    sample_key : string, optional
        Name of the sample as set in hklpy. If None it will ask for which
        sample.
    """
    _geom_ = get_diffractometer()
    if sample_key is None:
        d = _geom_.samples.keys()
        print("Sample keys:", list(d))
        sample_key = (
            input("\nEnter sample key [{}]: ".format(_geom_.sample.name))
            or _geom_.calc.sample.name
        )
    try:
        _geom_.sample = sample_key  # define the current sample
        print("\nCurrent sample: " + _geom_.sample.name)
        # to be done: check if orienting reflections exist
#        compute_UB()

    except KeyError:
        print("Not a valid sample key")

def sampleRemove(sample_key=None):
    """
    Remove selected sample in hklpy.

    Parameters
    ----------
    sample_key : string, optional
        Name of the sample as set in hklpy. If None it will ask for which
        sample.
    """
    _geom_ = get_diffractometer()
    if sample_key is None:
        d = _geom_.samples.keys()
        print("Sample keys:", list(d))
        sample_key = (
            input("\nEnter sample key to remove [{}]: ".format(_geom_.sample.name))
        )
        if sample_key == _geom_.sample.name:
            print("The current sample cannot be removed.")
            sample_key = " "
    try:
        _geom_.samples.pop(sample_key)  # remove the sample
        print("\n[{}] sample is removed.".format(sample_key))

    except KeyError:
        print("Not a valid sample key")

def or_swap():
    """Swaps the two orientation reflections in hklpy."""
    _geom_ = get_diffractometer()
    sample = _geom_.sample
    sample.reflections.swap()
    list_reflections()
    print("Computing UB!")
    sample.core.calc_UB(
        sample.reflections.order[0], sample.reflections.order[1]
    )
    _geom_.forward(1, 0, 0)


def calc_UB(r1, r2, wavelength=None, output=False):
    """
    Compute the UB matrix with two reflections.

    Parameters
    ----------
    r1, r2 : hklpy reflections
        Orienting reflections from hklpy.
    wavelength : float, optional
        This is not used...
    output : boolean
        Toggle to decide whether to print the UB matrix.
    """
#    _check_geom_selected()
    _geom_ = get_diffractometer()    
    _geom_.sample.core.calc_UB(r1, r2)
    if output:
        print(_geom_.sample.UB)

def compute_UB():
    """
    Calculates the UB matrix.

    This fixes one issue with the hklpy calc_UB in that using wh() right after
    will not work, it needs to run one calculation first.

    Parameters
    ----------
    h, k, l : float, optional
        Values of H, K, L positions for current reflection. If None, it will ask
        for it.
    """
    _geom_ = get_diffractometer()
    sample = _geom_.sample
    print("Computing UB!")
    calc_UB(
        sample.reflections.order[0], sample.reflections.order[1]
    )
    _geom_.forward(1, 0, 0)

def setor0():
    """
    Sets the primary orientation in hklpy2.

    Parameters
    ----------
    diffractometer real motors : float, optional
        Values of motor positions for current reflection. It will ask
        for it.
    h, k, l : float, optional
        Values of H, K, L positions for current reflection. It will ask
        for it.
    """
    _geom_ = get_diffractometer()
    sample = _geom_.sample
    orienting_refl = sample.reflections.order
    motors = _geom_.real_positioners._fields
    if len(orienting_refl) > 0:
        for key, ref in sample.reflections.items():
            if key == orienting_refl[0]:
                pos = list(ref.reals.values())
                old_h, old_k, old_l = list(ref.pseudos.values())
    else:
        pos = [None]* len(motors)
        for i in range(0,len(_geom_.real_positioners._fields)):
            pos[i]=0
        old_h = 0
        old_k = 0
        old_l = 0

    print("Enter primary-reflection angles:")
    or0pos = [None] * len(_geom_.real_positioners._fields)
    for i in range(0,len(_geom_.real_positioners._fields)):
        temppos = input("{} = [{:6.2f}]: ".format(motors[i], pos[i])) or pos[i]
        or0pos[i] = float(temppos)
    h = input("H = [{}]: ".format(old_h)) or old_h
    k = input("K = [{}]: ".format(old_k)) or old_k
    l = input("L = [{}]: ".format(old_l)) or old_l

    _geom_.add_reflection(
        (float(h),
        float(k),
        float(l)),
        or0pos
    )

    if len(orienting_refl) > 1:
        sample.reflections.order.pop(0)
        sample.reflections.order.pop(-1)
        sample.reflections.order.insert(
            0, list(sample.reflections.keys())[-1]
        )

    if len(orienting_refl) > 2:
        print("Computing UB!")
        sample.core.calc_UB(
            sample.reflections.order[0],
            sample.reflections.order[1],
        )
        _geom_.forward(1, 0, 0)



def setor1():
    """
    Sets the primary orientation in hklpy2.

    Parameters
    ----------
    diffractometer real motors : float, optional
        Values of motor positions for current reflection. It will ask
        for it.
    h, k, l : float, optional
        Values of H, K, L positions for current reflection. It will ask
        for it.
    """
    _geom_ = get_diffractometer()
    sample = _geom_.sample
    orienting_refl = sample.reflections.order
    motors = _geom_.real_positioners._fields
    if len(orienting_refl) > 1:
        for key, ref in sample.reflections.items():
            if key == orienting_refl[1]:
                pos = list(ref.reals.values())
                old_h, old_k, old_l = list(ref.pseudos.values())
    else:
        pos = [None]* len(motors)
        for i in range(0,len(_geom_.real_positioners._fields)):
            pos[i]=0
        old_h = 0
        old_k = 0
        old_l = 0

    print("Enter primary-reflection angles:")
    or0pos = [None] * len(_geom_.real_positioners._fields)
    for i in range(0,len(_geom_.real_positioners._fields)):
        temppos = input("{} = [{:6.2f}]: ".format(motors[i], pos[i])) or pos[i]
        or0pos[i] = float(temppos)
    h = input("H = [{}]: ".format(old_h)) or old_h
    k = input("K = [{}]: ".format(old_k)) or old_k
    l = input("L = [{}]: ".format(old_l)) or old_l

    _geom_.add_reflection(
        (float(h),
        float(k),
        float(l)),
        or0pos
    )

    if len(orienting_refl) > 2:
        sample.reflections.order.pop(1)
        sample.reflections.order.pop(-1)
        sample.reflections.order.insert(
            1, list(sample.reflections.keys())[-1]
        )

    if len(orienting_refl) > 1:
        print("Computing UB!")
        sample.core.calc_UB(
            sample.reflections.order[0],
            sample.reflections.order[1],
        )
        _geom_.forward(1, 0, 0)


def or0(h=None, k=None, l=None):
    """
    Sets the primary orientation in hklpy2 using the current motor positions.

    Parameters
    ----------
    h, k, l : float, optional
        Values of H, K, L positions for current reflection. If None, it will ask
        for it.
    """
    _geom_ = get_diffractometer()
    sample = _geom_.sample
    orienting_refl = sample.reflections.order
    if not h and not k and not l:
        if len(orienting_refl) > 0:
            for key, ref in sample.reflections.items():
                if key == orienting_refl[0]:
                    hr, kr, lr = list(ref.pseudos.values())
        else:
            hr = 2
            kr = 0
            lr = 0
        h = (input("H ({})? ".format(hr)) if not h else h) or hr
        k = (input("K ({})? ".format(kr)) if not k else k) or kr
        l = (input("L ({})? ".format(lr)) if not l else l) or lr
    _geom_.add_reflection(
        (float(h),
        float(k),
        float(l)),
        _geom_.real_position
    )

    if len(orienting_refl) > 1:
        sample.reflections.order.pop(0)
        sample.reflections.order.pop(-1)
        sample.reflections.order.insert(
            0, list(sample.reflections.keys())[-1]
        )

    if len(orienting_refl) > 2:
        print("Computing UB!")
        sample.core.calc_UB(
            sample.reflections.order[0],
            sample.reflections.order[1],
        )
        _geom_.forward(1, 0, 0)

def or1(h=None, k=None, l=None):
    """
    Sets the secondary orientation in hklpy2 using the current motor positions.

    Parameters
    ----------
    h, k, l : float, optional
        Values of H, K, L positions for current reflection. If None, it will ask
        for it.
    """
    _geom_ = get_diffractometer()
    sample = _geom_.sample
    orienting_refl = sample.reflections.order
    if not h and not k and not l:
        if len(orienting_refl) > 1:
            for key, ref in sample.reflections.items():
                if key == orienting_refl[1]:
                    hr, kr, lr = list(ref.pseudos.values())
        else:
            hr = 0
            kr = 2
            lr = 0
        h = (input("H ({})? ".format(hr)) if not h else h) or hr
        k = (input("K ({})? ".format(kr)) if not k else k) or kr
        l = (input("L ({})? ".format(lr)) if not l else l) or lr
    _geom_.add_reflection(
        (float(h),
        float(k),
        float(l)),
        _geom_.real_position
    )

    if len(orienting_refl) > 2:
        sample.reflections.order.pop(1)
        sample.reflections.order.pop(-1)
        sample.reflections.order.insert(
            1, list(sample.reflections.keys())[-1]
        )

    if len(orienting_refl) > 1:
        print("Computing UB!")
        sample.core.calc_UB(
            sample.reflections.order[0],
            sample.reflections.order[1],
        )
        _geom_.forward(1, 0, 0)

def set_orienting():
    """
    Change the primary and secondary orienting reflections
    to existing reflections in the reflection list (hklpy2).

    WARNING: Works only with six-circle geometry (fix in progress).
    """
    _geom_ = get_diffractometer()
    sample = _geom_.sample
    orienting_refl = sample.reflections.order

    # Column widths for alignment
    idx_width = 3
    pseudo_width = 9
    real_width = 10

    # Print header
    header = (
        f"\n{'#':>{idx_width}}"
        + "".join(f"{m:>{pseudo_width}}" for m in _geom_.pseudo_positioners._fields).upper()
        + "".join(f"{k:>{real_width}}" for k in _geom_.real_positioners._fields)
        + "   orienting"
    )
    print(header)

    # Find current orienting reflection indices (defaults: 0,1)
    or0_old, or1_old = 0, 1
    if len(orienting_refl) > 0:
        or0_old = list(sample.reflections.keys()).index(orienting_refl[0])
    if len(orienting_refl) > 1:
        or1_old = list(sample.reflections.keys()).index(orienting_refl[1])

    # Print reflection list
    for i, (key, ref) in enumerate(sample.reflections.items()):
        h, k, l = list(ref.pseudos.values())
        pos = list(ref.reals.values())

        tag = ""
        if len(orienting_refl) > 0 and key == orienting_refl[0]:
            tag = "first"
        elif len(orienting_refl) > 1 and key == orienting_refl[1]:
            tag = "second"

        row = (
            f"{i:>{idx_width}}"
            f"{h:{pseudo_width}.3f}{k:{pseudo_width}.3f}{l:{pseudo_width}.3f}"
            + "".join(f"{v:{real_width}.3f}" for v in pos)
            + (f"   {tag}" if tag else "")
        )
        print(row)

    # Prompt user for orienting reflections
    try:
        or0 = int(input(f"\nFirst orienting ({or0_old})? ") or or0_old)
        or1 = int(input(f"Second orienting ({or1_old})? ") or or1_old)
    except ValueError:
        print("Invalid input, keeping old orienting reflections.")
        or0, or1 = or0_old, or1_old

    keys = list(sample.reflections.keys())

    # Safely update orienting reflections
    if 0 <= or0 < len(keys):
        sample.reflections.order[0:1] = [keys[or0]]
    if 0 <= or1 < len(keys):
        if len(sample.reflections.order) > 1:
            sample.reflections.order[1:2] = [keys[or1]]
        else:
            sample.reflections.order.append(keys[or1])

    # Recompute UB matrix
    print("Computing UB!")
    sample.core.calc_UB(sample.reflections.order[0], sample.reflections.order[1])
    _geom_.forward(1, 0, 0)

def del_reflection():
    """
    Delete a reflection from the current sample's reflection list.

    - Displays reflections like set_orienting()
    - User inputs index number of reflection to delete
    - Orienting reflections (first/second) cannot be deleted
    """
    _geom_ = get_diffractometer()
    sample = _geom_.sample
    orienting_refl = sample.reflections.order

    # Column widths for alignment
    idx_width = 3
    pseudo_width = 9
    real_width = 10

    # Print header
    header = (
        f"\n{'#':>{idx_width}}"
        + "".join(f"{m:>{pseudo_width}}" for m in _geom_.pseudo_positioners._fields).upper()
        + "".join(f"{k:>{real_width}}" for k in _geom_.real_positioners._fields)
        + "   orienting"
    )
    print(header)

    # Print reflection list
    keys = list(sample.reflections.keys())
    for i, (key, ref) in enumerate(sample.reflections.items()):
        h, k, l = list(ref.pseudos.values())
        pos = list(ref.reals.values())

        tag = ""
        if len(orienting_refl) > 0 and key == orienting_refl[0]:
            tag = "first"
        elif len(orienting_refl) > 1 and key == orienting_refl[1]:
            tag = "second"

        row = (
            f"{i:>{idx_width}}"
            f"{h:{pseudo_width}.3f}{k:{pseudo_width}.3f}{l:{pseudo_width}.3f}"
            + "".join(f"{v:{real_width}.3f}" for v in pos)
            + (f"   {tag}" if tag else "")
        )
        print(row)

    # Prompt for index to delete
    try:
        idx = int(input("\nEnter reflection number to delete: "))
    except ValueError:
        print("Invalid input. Aborting delete.")
        return

    # Safety checks
    if idx < 0 or idx >= len(keys):
        print("Invalid index. Aborting delete.")
        return
    if keys[idx] in orienting_refl[:2]:
        print("Cannot delete orienting reflections (first/second). Aborting delete.")
        return

    # Delete reflection
    ref_key = keys[idx]
    del sample.reflections[ref_key]
    print(f"Deleted reflection '{ref_key}'.")

def list_orienting():
    """
    List the first and second orienting reflections only
    from the current sample's reflection list.
    """
    _geom_ = get_diffractometer()
    sample = _geom_.sample
    orienting_refl = sample.reflections.order

    if len(orienting_refl) < 2:
        print("Not enough orienting reflections defined.")
        return

    # Column widths for alignment
    idx_width = 3
    pseudo_width = 9
    real_width = 10

    # Print header
    header = (
        f"\n{'#':>{idx_width}}"
        + "".join(f"{m:>{pseudo_width}}" for m in _geom_.pseudo_positioners._fields).upper()
        + "".join(f"{k:>{real_width}}" for k in _geom_.real_positioners._fields)
        + "   orienting"
    )
    print(header)

    # Print only orienting reflections
    keys = list(sample.reflections.keys())
    for tag, key in zip(["first", "second"], orienting_refl[:2]):
        if key not in sample.reflections:
            continue
        idx = keys.index(key)
        ref = sample.reflections[key]
        h, k, l = list(ref.pseudos.values())
        pos = list(ref.reals.values())

        row = (
            f"{idx:>{idx_width}}"
            f"{h:{pseudo_width}.3f}{k:{pseudo_width}.3f}{l:{pseudo_width}.3f}"
            + "".join(f"{v:{real_width}.3f}" for v in pos)
            + f"   {tag}"
        )
        print(row)

def setmode(mode=None):
    """
    Set the mode of the currently selected diffractometer.

    Parameters
    ----------
    mode : string, optional
        Mode to be selected. If None, it will ask.
    """
    _geom_=get_diffractometer()
    current_mode = _geom_.core.mode
    for index, item in enumerate(_geom_.core.modes):
        print("{:2d}. {}".format(index + 1, item))
        if current_mode == item:
            current_index = index
    if mode:
        _geom_.core.mode = _geom_.core.modes[int(mode) - 1]
        print("\nSet mode to {}".format(mode))
    else:
        mode = input("\nMode ({})? ".format(current_index + 1)) or (
            current_index + 1
        )
        _geom_.core.mode = _geom_.core.modes[int(mode) - 1]

    # Freeze the appropriate detector angle at 0 based on mode geometry.
    # Local axis names are derived from the solver→local axis name mapping
    # (E6C solver: 'gamma' = horizontal detector, 'delta' = vertical detector).
    selected_mode = _geom_.core.mode
    if "vertical" in selected_mode:
        solver_det = "gamma"   # horizontal detector angle frozen at 0
    elif "horizontal" in selected_mode:
        solver_det = "delta"   # vertical detector angle frozen at 0
    else:
        solver_det = None

    if solver_det is not None:
        solver_reals = _geom_.core.solver_real_axis_names
        local_reals = list(_geom_.real_positioners._fields)
        solver_to_local = dict(zip(solver_reals, local_reals))
        det_angle = solver_to_local.get(solver_det)
        if det_angle is not None:
            _geom_.core.presets = {det_angle: 0}
            print(f"Preset {det_angle}=0 for mode '{selected_mode}'")
        else:
            _geom_.core.presets = {}
    else:
        _geom_.core.presets = {}

def freeze():
    """
    Interactively set preset values for the constant axes of the current mode.

    For each angle held constant in the current mode, prompts for a value
    (showing the current preset as default). Updates ``_geom_.core.presets``.
    """
    _geom_ = get_diffractometer()
    current_mode = _geom_.core.mode
    constant_axes = _geom_.core.constant_axis_names

    if not constant_axes:
        print(f"Mode '{current_mode}' has no constant (frozen) axes.")
        return

    print(f"Mode: {current_mode}")
    current_presets = _geom_.core.presets
    new_presets = {}
    for axis in constant_axes:
        current_val = current_presets.get(axis, 0.0)
        raw = input(f"  {axis} ({current_val})? ") or current_val
        try:
            new_presets[axis] = float(raw)
        except ValueError:
            print(f"  Invalid value for '{axis}', keeping {current_val}")
            new_presets[axis] = float(current_val)

    _geom_.core.presets = new_presets
    print("Frozen angles:")
    for axis, val in new_presets.items():
        print(f"  {axis} = {val}")


def ca(h, k, l):
    """
    Calculate the motors position of a reflection.

    Parameters
    ----------
    h, k, l : float
        H, K, and L values.
    """
    _geom_ = get_diffractometer()
    pos = cahkl(h, k, l)
    print("\n   Calculated Positions:")
    print(
        "\n   H K L = {:5f} {:5f} {:5f}".format(
            h,
            k,
            l,
        )
    )
    print(
        f"\n   Lambda (Energy) = {_geom_.beam.wavelength.get():6.4f} \u212b"
        f" ({_geom_.beam.energy.get():6.4f}) keV"
    )
    print(
        f"\n{''.join(f'{k:>10}' for k in _geom_.real_positioners._fields)}"
        f"\n{''.join(f'{v:>10.3f}' for v in pos)}"
    )


def _wh():
    import numpy as np
    """
    Retrieve information on the current reciprocal space position.

    WARNING: This function will only work with six circles. This will be fixed
    in future releases.
    """
    _geom_ = get_diffractometer()
    # _geom_for_psi_ = geometries.psi
    # _geom_for_psi_.sample.UB = _geom_.sample.UB
    # _geom_for_q_ = geometries.q2
    print(
        f"\n   {' '.join(_geom_.pseudo_positioners._fields).upper()}"
        f" = {', '.join([f'{v.position:5f}' for v in _geom_.pseudo_positioners])}"
    )
    print(
        f"\n   Lambda (Energy) = {_geom_.beam.wavelength.get():6.4f} \u212b"
        f" ({_geom_.beam.energy.get():6.4f}) keV"
    )
    print(
        f"\n{''.join(f'{k:>10}' for k in _geom_.real_positioners._fields)}"
        f"\n{''.join(f'{v.position:>10.3f}' for v in _geom_.real_positioners)}"
    )
    # print(
    #     "\n   PSI = {:5.4f} ".format(
    #         _geom_for_psi_.inverse(0).psi,
    #     )
    # )
    # _h2, _k2, _l2 = _geom_for_psi_.core.extras.values()
    # print(
    #     "   PSI reference vector = {:3.3f} {:3.3f} {:3.3f}".format(
    #         _h2,
    #         _k2,
    #         _l2,
    #     )
    # )
    # tth_from_q = 2*np.emath.arcsin(_geom_for_q_.inverse(0).q/4/np.pi*_geom_.beam.wavelength.get())*180/np.pi
    # print(
    #     "\n   Q = {:5f}  tth = {:5f}".format(
    #         _geom_for_q_.inverse(0).q, tth_from_q
    #     )
    # )
 
def _ensure_idle():
    if  RE.state != 'idle':
        print('The RunEngine invoked by magics cannot be resumed.')
        print('Aborting...')
        RE.abort()

def ubr(h, k, l):
    """
    Move the motors to a reciprocal space point.

    Parameters
    ----------
    h, k, l : float
        H, K, and L values.

    Returns
    -------
    """
    _geom_ = get_diffractometer()
    plan = mv(
        _geom_.h, float(h), _geom_.k, float(k), _geom_.l, float(l)
    )
    RE.waiting_hook = pbar_manager
    try:
        RE(plan)
    except RunEngineInterrupted:
        ...
    RE.waiting_hook = None
    _ensure_idle()
    return None


def br(h, k, l):
    """
    Move the motors to a reciprocal space point.

    Parameters
    ----------
    h, k, l : float
        H, K, and L values.

    Returns
    -------
    Generator for the bluesky Run Engine.

    Example
    -------
    RE(br(2, 2, 0))
    """
    _geom_ = get_diffractometer()
    yield from bps.mv(
        _geom_.h, float(h), _geom_.k, float(k), _geom_.l, float(l)
    )

def uan(*args):
    """
    Moves the delta and eta motors.

    WARNING: This function will only work with six circles. This will be fixed
    in future releases.

    Parameters
    ----------
    delta, eta: float, optional??
        Delta and eta motor angles to be moved to.

    Returns
    -------
    """
    _geom_ = get_diffractometer()
    if len(args) != 2:
        delta, eta = args
        raise ValueError("Usage: uan(delta/tth,eta/th)")
    else:
        delta, eta = args
        if len(_geom_.real_position) == 6:
            print("Moving to (delta,eta)=({},{})".format(delta, eta))
            plan = mv(_geom_.delta, delta, _geom_.eta, eta)
        elif len(_geom_.real_position) == 4:
            print("Moving to (tth,th)=({},{})".format(delta, eta))
            plan = mv(_geom_.tth, delta, _geom_.th, eta)
    RE.waiting_hook = pbar_manager
    try:
        RE(plan)
    except RunEngineInterrupted:
        ...
    RE.waiting_hook = None
    _ensure_idle()
    return None



def an(delta=None, eta=None):
    """
    Moves the delta and theta motors.

    WARNING: This function will only work with six circles. This will be fixed
    in future releases.

    Parameters
    ----------
    delta, th: float, optional??
        Delta and th motor angles to be moved to.

    Returns
    -------
    Generator for the bluesky Run Engine.
    """
    _geom_ = get_diffractometer()
    if len(args) != 2:
        delta, eta = args
        raise ValueError("Usage: uan(delta/tth,eta/th)")
    else:
        delta, eta = args
        if len(_geom_.real_position) == 6:
            print("Moving to (delta,eta)=({},{})".format(delta, eta))
            yield from bps.mv(_geom_.delta, delta, _geom_.eta, eta)
        elif len(_geom_.real_position) == 4:
            print("Moving to (tth,th)=({},{})".format(delta, eta))
            yield from bps.mv(_geom_.tth, delta, _geom_.th, eta)

def setlat(*args):
    """
    Set the lattice constants for the current sample.

    Parameters
    ----------
    a, b, c, alpha, beta, gamma : float, optional
        Lattice constants. 
        - If all 6 are provided as arguments, they are used directly.
        - If none are provided, the user will be prompted interactively.
    
    Notes
    -----
    If two orienting reflections are defined, UB will be recomputed automatically.
    """
    _geom_ = get_diffractometer()
    sample = _geom_.sample

    # Current lattice values
    param_names = sample.lattice.system_parameter_names(0)
    current_values = [getattr(sample.lattice, p) for p in param_names]

    # Case 1: direct arguments
    if len(args) == 6:
        new_values = args

    # Case 2: interactive input
    elif len(args) == 0:
        new_values = []
        for name, current in zip(param_names, current_values):
            try:
                val = input(f"Lattice {name} ({current})? ") or current
                new_values.append(float(val))
            except ValueError:
                print(f"Invalid input for {name}, keeping old value {current}.")
                new_values.append(current)

    else:
        raise ValueError("Provide either 0 or 6 arguments: a, b, c, alpha, beta, gamma.")

    # Apply new lattice parameters
    for name, val in zip(param_names, new_values):
        setattr(sample.lattice, name, float(val))

    # Recompute UB if orienting reflections exist
    if len(sample.reflections.order) > 1:
        print("Computing UB...")
        sample.core.calc_UB(sample.reflections.order[0], sample.reflections.order[1])
        _geom_.forward(1, 0, 0)

    # Final confirmation
    print("\nUpdated lattice parameters:")
    for name in param_names:
        print(f"  {name} = {getattr(sample.lattice, name):.4f}")

def setaz(*args):
    """
    Set the azimuthal reference reflection (h2, k2, l2).
    
    Works differently depending on geometry:
      - 4-circle: uses 'psi_constant'
      - 6-circle: uses 'psi_constant_vertical' and 'psi_constant_horizontal'

    Parameters
    ----------
    h2, k2, l2 : int or float, optional
        Miller indices of the azimuthal reference reflection.
        - If provided, these are used directly.
        - If not provided, the user will be prompted interactively.
    """
    _geom_ = get_diffractometer()
    _geom_for_psi_ = geometries.psi

    # Determine geometry type
    npos = len(_geom_.real_position)
    if npos not in (4, 6):
        raise ValueError(f"Function not available in mode '{_geom_.core.mode}'")

    # Save current mode
    mode_temp = _geom_.core.mode

    # Set mode based on geometry
    if npos == 4:
        _geom_.core.mode = "psi_constant"
    else:  # 6-circle
        _geom_.core.mode = "psi_constant_vertical"

    # Current values
    _h2, _k2, _l2, psi = list(_geom_.core.extras.values())

    # Parse input
    if len(args) == 3:
        h2, k2, l2 = args
    elif len(args) == 0:
        def ask(name, current):
            val = input(f"{name} = ({current})? ") or current
            try:
                return float(val)
            except ValueError:
                print(f"Invalid {name}, keeping {current}")
                return current

        h2 = ask("H", _h2)
        k2 = ask("K", _k2)
        l2 = ask("L", _l2)
    else:
        raise ValueError("Either no arguments or exactly 3 arguments (h, k, l) must be provided.")

    # Apply new extras to both geometries
    extras = {"h2": h2, "k2": k2, "l2": l2}
    _geom_.core.extras = extras
    _geom_for_psi_.core.extras = extras

    # For 6-circle, also set horizontal psi mode
    if npos == 6:
        _geom_.core.mode = "psi_constant_horizontal"
        _geom_.core.extras = extras

    # Print result
    print(f"Azimuth = {h2} {k2} {l2} with Psi fixed at {psi}")

    # Restore original mode if desired
    _geom_.core.mode = mode_temp

def freeze_psi(*args):
    _geom_ = get_diffractometer()
    _geom_for_psi_ = geometries.psi
    if (_geom_.core.mode == "psi_constant" or 
        _geom_.core.mode == "psi_constant_vertical" or
        _geom_.core.mode == "psi_constant_horizontal"
    ):
        _h2, _k2, _l2, psi = list(_geom_.core.extras.values())
        if len(args) == 0:
            psi =_geom_for_psi_.inverse(0).psi
        elif len(args) == 1:
            psi = args[0]
        else:
            raise ValueError(
                "either no argument or azimuth needs to be provided."
            )
        _geom_.core.extras = {'psi': psi}
        print("Psi = {}".format(psi))
    else:
        raise ValueError(
            "Function not available in mode '{}'".format(
                _geom_.core.mode
            )
        )
 
