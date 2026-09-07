import time
from apsbits.core.instrument_init import oregistry

ur5 = oregistry["ur5"]

# Coordinates: [x, y, z, roll, pitch, yaw]
ur5_home = [454.5, -314.9, 500, -179.74, 3.23, -0.12]

vial_up = [448.48, -594.82, 350, -179.74, 3.23, -0.12]
vial_down = [448.48, -594.82, 240, -179.74, 3.23, -0.12]

pipette_lock_up = [397.86, 98.4, 100, -179.74, 3.23, -0.12]
pipette_lock_down = [397.86, 98.4, 6.41, -179.74, 3.23, -0.12]

pipette_unlock_up = [397.86, 188.43, 200, -179.74, 3.23, -0.12]
pipette_unlock_down = [397.86, 188.43, 6.41, -179.74, 3.23, -0.12]

capillary_up = [473.95, -519.38, 300, -179.74, 3.23, -0.12]
capillary_down = [473.95, -519.38, 248.5, -179.74, 3.23, -0.12]


# ---------------------------------------------------------------------------
# Workflow steps
# ---------------------------------------------------------------------------

def home():
    ur5.auto_move_l.put(0)
    ur5.pose.set(ur5_home)


def pickup_pipette():
    """
    ur5_home → pipette_lock_up → pipette_lock_down → pipette_unlock_down
             → pipette_unlock_up → home pipette → ur5_home
    """
    print()
    print("[pickup_pipette] Starting")
    home()
    ur5.pose.set(pipette_lock_up)
    ur5.pose.set(pipette_lock_down)
    ur5.pose.set(pipette_unlock_down)
    ur5.pose.set(pipette_unlock_up)
    home()

    ur5.pipette.home.put(1)
    time.sleep(1)
    print("[pickup_pipette] Done")
    print()


def pipette_liquid(volume=20):
    """
    ur5_home → vial_up → vial_down → aspirate → vial_up → ur5_home
    """
    print()
    print(f"[pipette_liquid] Starting  (volume={volume}uL)")
    home()
    ur5.pose.set(vial_up)
    ur5.pose.set(vial_down)

    ur5.pipette.set_volume.put(volume)
    time.sleep(1)
    ur5.pipette.aspirate.put(1)
    time.sleep(3)

    ur5.pose.set(vial_up)
    home()
    print("[pipette_liquid] Done")
    print()


def dispense_liquid(volume=20):
    """
    ur5_home → capillary_up → capillary_down → dispense → capillary_up → ur5_home
    """
    print()
    print(f"[dispense_liquid] Starting  (volume={volume}uL)")
    home()
    ur5.pose.set(capillary_up)
    ur5.pose.set(capillary_down)

    ur5.pipette.set_volume.put(volume)
    time.sleep(1)
    ur5.pipette.dispense.put(1)
    time.sleep(5)

    ur5.pose.set(capillary_up)
    home()
    print("[dispense_liquid] Done")
    print()


def dock_pipette():
    """
    ur5_home → pipette_unlock_up → pipette_unlock_down → pipette_lock_down
             → pipette_lock_up → ur5_home
    """
    print()
    print("[dock_pipette] Starting")
    home()
    ur5.pose.set(pipette_unlock_up)
    ur5.pose.set(pipette_unlock_down)
    ur5.pose.set(pipette_lock_down)
    ur5.pose.set(pipette_lock_up)
    home()
    print("[dock_pipette] Done")
    print()


# ---------------------------------------------------------------------------
# Full workflow
# ---------------------------------------------------------------------------

def run_pipetting_workflow(volume=30):
    """
    1. Move to home
    2. Pick up pipette
    3. Aspirate from vial
    4. Dispense into capillary
    5. Dock pipette
    """
    print()
    print("=== Pipetting workflow start ===")
    home()
    pickup_pipette()
    pipette_liquid(volume)
    dispense_liquid(volume)
    dock_pipette()
    print("=== Pipetting workflow complete ===")
    print()
