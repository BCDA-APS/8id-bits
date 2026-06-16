"""UR5 robotic arm control device.

Ophyd wrapper around the Universal Robots UR5 PVs served by the
``urRobot`` EPICS support module. Exposes the six joint readbacks and
command setpoints, the TCP pose (Cartesian + Euler) readbacks and command
setpoints, motion parameters, and the move/stop triggers needed to drive
the arm from Bluesky.
"""

from ophyd import Component
from ophyd import Device
from ophyd import EpicsSignal
from ophyd import EpicsSignalRO
from ophyd.status import SubscriptionStatus


def _motion_complete_callback():
    """Build a SubscriptionStatus callback that requires a 1->0->1 cycle.

    ``Control:AsyncMoveDone`` reads 1 when idle. A fresh ``.set()`` must
    only resolve on the *new* move's 0->1 transition, not on the tail end
    of a prior move that was still completing when we subscribed. Tracking
    whether we have seen the 1->0 falling edge first closes that race.
    """
    saw_motion_start = [False]

    def _done(*, old_value, value, **kw):
        if old_value == 1 and value == 0:
            saw_motion_start[0] = True
            return False
        return saw_motion_start[0] and old_value == 0 and value == 1

    return _done


class _UR5JointGroup(Device):
    """Multi-axis joint positioner.

    ``set((j1, j2, j3, j4, j5, j6))`` writes all six J*Cmd setpoints,
    triggers ``Control:moveJ.PROC``, and returns a Status that completes
    when ``Control:AsyncMoveDone`` cycles 1->0->1 (motion started then
    finished). The 1->0 requirement ensures the Status does not resolve
    on the tail of a prior move that happened to still be in flight.
    """

    j1 = Component(EpicsSignal, "Control:J1Cmd", kind="omitted")
    j2 = Component(EpicsSignal, "Control:J2Cmd", kind="omitted")
    j3 = Component(EpicsSignal, "Control:J3Cmd", kind="omitted")
    j4 = Component(EpicsSignal, "Control:J4Cmd", kind="omitted")
    j5 = Component(EpicsSignal, "Control:J5Cmd", kind="omitted")
    j6 = Component(EpicsSignal, "Control:J6Cmd", kind="omitted")
    actuate = Component(EpicsSignal, "Control:moveJ.PROC", kind="omitted")
    done = Component(EpicsSignalRO, "Control:AsyncMoveDone", kind="omitted")

    set_timeout = 60.0

    def set(self, target):
        targets = tuple(target)
        if len(targets) != 6:
            raise ValueError(
                f"UR5 joint group needs 6 targets (j1..j6), got {len(targets)}"
            )
        for cpt, val in zip(
            (self.j1, self.j2, self.j3, self.j4, self.j5, self.j6), targets
        ):
            cpt.put(val)
        status = SubscriptionStatus(
            self.done,
            _motion_complete_callback(),
            timeout=self.set_timeout,
        )
        self.actuate.put(1)
        return status


class _UR5PoseGroup(Device):
    """Cartesian TCP positioner.

    ``set((x, y, z, roll, pitch, yaw))`` writes all six Pose*Cmd setpoints,
    triggers ``Control:moveL.PROC``, and returns a Status that completes
    when ``Control:AsyncMoveDone`` cycles 1->0->1 (motion started then
    finished). The 1->0 requirement ensures the Status does not resolve
    on the tail of a prior move that happened to still be in flight.
    """

    x = Component(EpicsSignal, "Control:PoseXCmd", kind="omitted")
    y = Component(EpicsSignal, "Control:PoseYCmd", kind="omitted")
    z = Component(EpicsSignal, "Control:PoseZCmd", kind="omitted")
    roll = Component(EpicsSignal, "Control:PoseRollCmd", kind="omitted")
    pitch = Component(EpicsSignal, "Control:PosePitchCmd", kind="omitted")
    yaw = Component(EpicsSignal, "Control:PoseYawCmd", kind="omitted")
    actuate = Component(EpicsSignal, "Control:moveL.PROC", kind="omitted")
    done = Component(EpicsSignalRO, "Control:AsyncMoveDone", kind="omitted")

    set_timeout = 60.0

    def set(self, target):
        targets = tuple(target)
        if len(targets) != 6:
            raise ValueError(
                f"UR5 pose group needs 6 targets (x,y,z,roll,pitch,yaw), "
                f"got {len(targets)}"
            )
        for cpt, val in zip(
            (self.x, self.y, self.z, self.roll, self.pitch, self.yaw), targets
        ):
            cpt.put(val)
        status = SubscriptionStatus(
            self.done,
            _motion_complete_callback(),
            timeout=self.set_timeout,
        )
        self.actuate.put(1)
        return status


class _PipetteGroup(Device):
    """Tricontinent OEM pipette mounted on the UR5.

    Wraps the ``Pipette:*`` records from the IOC's
    ``tricontinent_pipette.db``. ``home`` and ``aspirate`` are ``bo`` write
    triggers; the ``set_*`` signals are scalar setpoints in their record's
    engineering units; ``calc_plunger_position`` / ``actual_plunger_position``
    are the calculated and encoder-derived plunger positions in increments.
    """

    home = Component(EpicsSignal, "Pipette:Home")
    set_volume = Component(EpicsSignal, "Pipette:SetVolume")
    dispense = Component(EpicsSignal, "Pipette:Dispense")
    aspirate = Component(EpicsSignal, "Pipette:Aspirate")
    set_speed = Component(EpicsSignal, "Pipette:SetSpeed")
    set_dead_volume = Component(EpicsSignal, "Pipette:SetDeadVolume")
    calc_plunger_position = Component(EpicsSignalRO, "Pipette:CPPI")
    actual_plunger_position = Component(EpicsSignalRO, "Pipette:APPS")
    set_plunger = Component(EpicsSignal, "Pipette:SetPlunger")


class UR5(Device):
    """Universal Robots UR5 arm exposed via the urRobot RTDE support module."""

    joint1_readback = Component(EpicsSignalRO, "Receive:Joint1")
    joint1_setpoint = Component(EpicsSignal, "Control:J1Cmd")
    joint2_readback = Component(EpicsSignalRO, "Receive:Joint2")
    joint2_setpoint = Component(EpicsSignal, "Control:J2Cmd")
    joint3_readback = Component(EpicsSignalRO, "Receive:Joint3")
    joint3_setpoint = Component(EpicsSignal, "Control:J3Cmd")
    joint4_readback = Component(EpicsSignalRO, "Receive:Joint4")
    joint4_setpoint = Component(EpicsSignal, "Control:J4Cmd")
    joint5_readback = Component(EpicsSignalRO, "Receive:Joint5")
    joint5_setpoint = Component(EpicsSignal, "Control:J5Cmd")
    joint6_readback = Component(EpicsSignalRO, "Receive:Joint6")
    joint6_setpoint = Component(EpicsSignal, "Control:J6Cmd")

    pose_x_readback = Component(EpicsSignalRO, "Receive:PoseX")
    pose_x_setpoint = Component(EpicsSignal, "Control:PoseXCmd")
    pose_y_readback = Component(EpicsSignalRO, "Receive:PoseY")
    pose_y_setpoint = Component(EpicsSignal, "Control:PoseYCmd")
    pose_z_readback = Component(EpicsSignalRO, "Receive:PoseZ")
    pose_z_setpoint = Component(EpicsSignal, "Control:PoseZCmd")
    pose_roll_readback = Component(EpicsSignalRO, "Receive:PoseRoll")
    pose_roll_setpoint = Component(EpicsSignal, "Control:PoseRollCmd")
    pose_pitch_readback = Component(EpicsSignalRO, "Receive:PosePitch")
    pose_pitch_setpoint = Component(EpicsSignal, "Control:PosePitchCmd")
    pose_yaw_readback = Component(EpicsSignalRO, "Receive:PoseYaw")
    pose_yaw_setpoint = Component(EpicsSignal, "Control:PoseYawCmd")

    linear_speed = Component(EpicsSignal, "Control:LinearSpeed")
    joint_acceleration = Component(EpicsSignal, "Control:JointAcceleration")

    move_j = Component(EpicsSignal, "Control:moveJ.PROC", kind="omitted")
    move_l = Component(EpicsSignal, "Control:moveL.PROC", kind="omitted")
    stop_j = Component(EpicsSignal, "Control:stopJ.PROC", kind="omitted")
    stop_l = Component(EpicsSignal, "Control:stopL.PROC", kind="omitted")
    auto_move_j = Component(EpicsSignal, "Control:AutoMoveJ")
    auto_move_l = Component(EpicsSignal, "Control:AutoMoveL")

    joints = Component(_UR5JointGroup, "")
    pose = Component(_UR5PoseGroup, "")
    pipette = Component(_PipetteGroup, "")
