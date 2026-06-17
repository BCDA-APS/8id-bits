import bluesky.plan_stubs as bps
from apsbits.core.instrument_init import oregistry

ur5 = oregistry["ur5"]

# Coordinates: [x, y, z, roll, pitch, yaw]
ur5_home = [454.5, -314.9, 500, -179.74, 3.23, -0.12]

vial_up = [454.5, -584.9, 350, -179.74, 3.23, -0.12]
vial_down = [454.5, -584.9, 240, -179.74, 3.23, -0.12]

pipette_lock_up = [397.86, 98.4, 100, -179.74, 3.23, -0.12]
pipette_lock_down = [397.86, 98.4, 6.41, -179.74, 3.23, -0.12]

pipette_unlock_up = [397.86, 188.43, 200, -179.74, 3.23, -0.12]
pipette_unlock_down = [397.86, 188.43, 6.41, -179.74, 3.23, -0.12]

capillary_up = [473.85, -518.94, 300, -179.74, 3.23, -0.12]
capillary_down = [473.85, -518.94, 248.5, -179.74, 3.23, -0.12]


# ---------------------------------------------------------------------------
# Low-level helpers
# ---------------------------------------------------------------------------

def move_to(coords):
    """Linear move to [x, y, z, roll, pitch, yaw]; waits for AsyncMoveDone 1→0→1."""
    yield from bps.mv(ur5.pose, coords)


def home_pipette():
    yield from bps.abs_set(ur5.pipette.home, 1, wait=True)


def aspirate(volume=20):
    yield from bps.abs_set(ur5.pipette.set_volume, volume, wait=True)
    yield from bps.abs_set(ur5.pipette.aspirate, 1, wait=True)


def dispense(volume=20):
    yield from bps.abs_set(ur5.pipette.set_volume, volume, wait=True)
    yield from bps.abs_set(ur5.pipette.dispense, 1, wait=True)


# ---------------------------------------------------------------------------
# Workflow steps
# ---------------------------------------------------------------------------

def pickup_pipette():
    """
    ur5_home → pipette_lock_up → pipette_lock_down → pipette_unlock_down
             → pipette_unlock_up → home pipette → ur5_home
    """
    yield from move_to(ur5_home)
    yield from move_to(pipette_lock_up)
    yield from move_to(pipette_lock_down)
    yield from move_to(pipette_unlock_down)
    yield from move_to(pipette_unlock_up)
    yield from home_pipette()
    yield from move_to(ur5_home)


def pipette_liquid(volume=20):
    """
    ur5_home → vial_up → vial_down → aspirate → vial_up → ur5_home
    """
    yield from move_to(ur5_home)
    yield from move_to(vial_up)
    yield from move_to(vial_down)
    yield from aspirate(volume)
    yield from move_to(vial_up)
    yield from move_to(ur5_home)


def dispense_liquid(volume=20):
    """
    ur5_home → capillary_up → capillary_down → dispense → capillary_up → ur5_home
    """
    yield from move_to(ur5_home)
    yield from move_to(capillary_up)
    yield from move_to(capillary_down)
    yield from dispense(volume)
    yield from move_to(capillary_up)
    yield from move_to(ur5_home)


def dock_pipette():
    """
    ur5_home → pipette_unlock_up → pipette_unlock_down → pipette_lock_down
             → pipette_lock_up → ur5_home
    """
    yield from move_to(ur5_home)
    yield from move_to(pipette_unlock_up)
    yield from move_to(pipette_unlock_down)
    yield from move_to(pipette_lock_down)
    yield from move_to(pipette_lock_up)
    yield from move_to(ur5_home)


# ---------------------------------------------------------------------------
# Full workflow
# ---------------------------------------------------------------------------

def run_pipetting_workflow(volume=20):
    """
    Full pipetting workflow submitted to the Run Engine:
    1. Move to home
    2. Pick up pipette
    3. Aspirate from vial
    4. Dispense into capillary
    5. Dock pipette
    """
    yield from move_to(ur5_home)
    yield from pickup_pipette()
    yield from pipette_liquid(volume)
    yield from dispense_liquid(volume)
    yield from dock_pipette()
