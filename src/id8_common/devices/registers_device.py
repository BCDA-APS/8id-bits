"""
EPICS PVs as Storage registers.
"""

from ophyd import Component
from ophyd import Device
from ophyd import EpicsSignal


class EpicsPvStorageRegisters(Device):
    """A device class for managing EPICS PV storage registers.

    This class provides functionality for storing and retrieving values
    in EPICS PV registers. It is used for temporary storage of values
    and parameters during beamline operation.
    """

    cycle_name = Component(EpicsSignal, "StrReg1", string=True)
    geometry = Component(EpicsSignal, "StrReg2", string=True)
    mount_point = Component(EpicsSignal, "StrReg3", string=True)
    experiment_name = Component(EpicsSignal, "StrReg4", string=True)
    analysis_machine = Component(EpicsSignal, "StrReg5", string=True)

    workflow_name = Component(EpicsSignal, "StrReg6", string=True)
    use_subfolder = Component(EpicsSignal, "StrReg7", string=True)

    header = Component(EpicsSignal, "StrReg11", string=True)
    sample_name = Component(EpicsSignal, "StrReg12", string=True)
    sample_move = Component(EpicsSignal, "StrReg13", string=True)
    inner_motor = Component(EpicsSignal, "StrReg14", string=True)
    outer_motor = Component(EpicsSignal, "StrReg15", string=True)

    det_name = Component(EpicsSignal, "StrReg16", string=True)
    det_mode = Component(EpicsSignal, "StrReg17", string=True)
    qmap_file = Component(EpicsSignal, "StrReg18", string=True)
    analysis_type = Component(EpicsSignal, "StrReg19", string=True)

    measurement_num = Component(EpicsSignal, "Reg1")
    acq_time = Component(EpicsSignal, "Reg2")
    acq_period = Component(EpicsSignal, "Reg3")
    num_frames = Component(EpicsSignal, "Reg4")
    num_repeats = Component(EpicsSignal, "Reg5")

    sample_index = Component(EpicsSignal, "Reg6")
    inner_center = Component(EpicsSignal, "Reg7")
    outer_center = Component(EpicsSignal, "Reg8")
    inner_range = Component(EpicsSignal, "Reg9")
    outer_range = Component(EpicsSignal, "Reg10")

    inner_pts = Component(EpicsSignal, "Reg11")
    outer_pts = Component(EpicsSignal, "Reg12")

    sample1_pos = Component(EpicsSignal, "Reg16")
    sample2_pos = Component(EpicsSignal, "Reg17")
    sample3_pos = Component(EpicsSignal, "Reg18")
    sample4_pos = Component(EpicsSignal, "Reg19")
    sample5_pos = Component(EpicsSignal, "Reg20")

    sample6_pos = Component(EpicsSignal, "Reg21")
    sample7_pos = Component(EpicsSignal, "Reg22")
    sample8_pos = Component(EpicsSignal, "Reg23")
    sample9_pos = Component(EpicsSignal, "Reg24")
    sample10_pos = Component(EpicsSignal, "Reg25")

    sample11_pos = Component(EpicsSignal, "Reg26")
    sample12_pos = Component(EpicsSignal, "Reg27")
    sample13_pos = Component(EpicsSignal, "Reg28")
    sample14_pos = Component(EpicsSignal, "Reg29")
    sample15_pos = Component(EpicsSignal, "Reg30")

    sample16_pos = Component(EpicsSignal, "Reg31")
    sample17_pos = Component(EpicsSignal, "Reg32")
    sample18_pos = Component(EpicsSignal, "Reg33")
    sample19_pos = Component(EpicsSignal, "Reg34")
    sample20_pos = Component(EpicsSignal, "Reg35")

    sample21_pos = Component(EpicsSignal, "Reg36")
    sample22_pos = Component(EpicsSignal, "Reg37")
    sample23_pos = Component(EpicsSignal, "Reg38")
    sample24_pos = Component(EpicsSignal, "Reg39")
    sample25_pos = Component(EpicsSignal, "Reg40")

    sample26_pos = Component(EpicsSignal, "Reg41")

    det_pixel_size = Component(EpicsSignal, "Reg48")
    eiger_det_x0 = Component(EpicsSignal, "Reg49")
    eiger_det_y0 = Component(EpicsSignal, "Reg50")

    eiger_db_x0 = Component(EpicsSignal, "Reg51")
    eiger_db_y0 = Component(EpicsSignal, "Reg52")
    rigaku_det_x0 = Component(EpicsSignal, "Reg53")
    rigaku_det_y0 = Component(EpicsSignal, "Reg54")
    rigaku_db_x0 = Component(EpicsSignal, "Reg55")

    rigaku_db_y0 = Component(EpicsSignal, "Reg56")
    current_det_x0 = Component(EpicsSignal, "Reg57")
    current_det_y0 = Component(EpicsSignal, "Reg58")
    current_db_x0 = Component(EpicsSignal, "Reg59")
    current_db_y0 = Component(EpicsSignal, "Reg60")

