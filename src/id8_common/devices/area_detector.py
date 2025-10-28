"""
EPICS area_detector definitions for ID8.
"""

import logging

from apstools.devices import AD_EpicsFileNameHDF5Plugin
from apstools.devices import CamMixin_V34
from ophyd import ADComponent
from ophyd import EpicsSignal
from ophyd import EpicsSignalRO
from ophyd import EpicsSignalWithRBV
from ophyd.areadetector import AreaDetector
from ophyd.areadetector import CamBase
from ophyd.areadetector import EigerDetectorCam
from ophyd.areadetector import SimDetectorCam
from ophyd.areadetector.plugins import CodecPlugin_V34
from ophyd.areadetector.plugins import FileBase
from ophyd.areadetector.plugins import ImagePlugin_V34
from ophyd.areadetector.plugins import OverlayPlugin_V34
from ophyd.areadetector.plugins import PluginBase_V34
from ophyd.areadetector.plugins import ProcessPlugin_V34
from ophyd.areadetector.plugins import PvaPlugin_V34
from ophyd.areadetector.plugins import ROIPlugin_V34
from ophyd.areadetector.plugins import StatsPlugin_V34
from ophyd.areadetector.plugins import TransformPlugin_V34
from ophyd.ophydobj import Kind
from ophyd.status import Status

logger = logging.getLogger(__name__)
logger.info(__file__)

PLUGINS__CLEAR_STAGE_SIGS = "image process1 transform1".split()


def ad_setup(det: AreaDetector, iconfig: dict) -> None:
    """not a plan: Steps to prepare an area detector object."""
    # just in case these things are not defined in the class source code
    det.cam.stage_sigs["wait_for_plugins"] = "Yes"
    for nm in det.component_names:
        obj = getattr(det, nm)
        if "blocking_callbacks" in dir(obj):  # is it a plugin?
            obj.stage_sigs["blocking_callbacks"] = "No"
    for nm in PLUGINS__CLEAR_STAGE_SIGS:
        getattr(det, nm).stage_sigs = {}

    plugin = det.hdf1  # for convenience below
    plugin.kind = Kind.config | Kind.normal  # Ensure plugin's read is called.
    plugin.stage_sigs.move_to_end("capture", last=True)

    # if iconfig.get("ALLOW_AREA_DETECTOR_WARMUP", False):
    #     if det.connected:
    #         if not AD_plugin_primed(plugin):
    #             AD_prime_plugin2(plugin)


class CamBase_V34(CamBase):
    """
    Updates to CamBase since v22.

    PVs removed from AD now.
    """

    pool_max_buffers = None


class FileBase_V34(FileBase):
    """
    Updates to FileBase since v22.

    PVs removed from AD now.
    """

    file_number_sync = None
    file_number_write = None


class EigerDetectorCam_V34(CamMixin_V34, EigerDetectorCam):
    """Revise EigerDetectorCam for ADCore revisions."""

    initialize = ADComponent(EpicsSignal, "Initialize", kind="config")

    # These components not found on Eiger 4M at 8-ID-I
    file_number_sync = None
    file_number_write = None
    fw_clear = None
    link_0 = None
    link_1 = None
    link_2 = None
    link_3 = None
    dcu_buff_free = None
    offset = None


class BasicCam_V34(CamMixin_V34, CamBase_V34):
    """Basic camera support, such as the HHL Mirrors in 8-ID-A."""

    offset = None


class Lambda2MCam(CamBase_V34):
    """
    Support for the Lambda 2M cam controls.

    https://x-spectrum.de/products/lambda/
    """

    _html_docs = ["Lambda2MCam.html"]

    firmware_version = ADComponent(EpicsSignalRO, "FirmwareVersion_RBV", kind="omitted")
    operating_mode = ADComponent(EpicsSignalWithRBV, "OperatingMode", kind="config")
    serial_number = ADComponent(EpicsSignalRO, "SerialNumber_RBV", kind="omitted")
    temperature = ADComponent(EpicsSignalWithRBV, "Temperature", kind="config")
    wait_for_plugins = ADComponent(
        EpicsSignal,
        "WaitForPlugins",
        string=True,
        kind="config",
    )

    energy_threshold = ADComponent(EpicsSignalWithRBV, "EnergyThreshold", kind="config")
    dual_mode = ADComponent(EpicsSignalWithRBV, "DualMode", string=True, kind="config")
    dual_threshold = ADComponent(EpicsSignalWithRBV, "DualThreshold", kind="config")

    EXT_TRIGGER = 0


class Rigaku3MCam(CamBase_V34):
    """Support for the RigakuSi3M camera controls."""

    _html_docs = ["Rigaku3MCam.html"]
    wait_for_plugins = ADComponent(
        EpicsSignal, "WaitForPlugins", string=True, kind="config"
    )

    # sparse_enable = ADComponent(EpicsSignal, "SparseEnable", string=True)
    fast_file_name = ADComponent(EpicsSignalWithRBV, "FileName", string=True)
    fast_file_path = ADComponent(EpicsSignalWithRBV, "FilePath", string=True)
    image_mode = ADComponent(EpicsSignalWithRBV, "ImageMode", string=True)
    output_resolution = ADComponent(EpicsSignalWithRBV, "OutputResolution", string=True)
    dual_threshold = ADComponent(EpicsSignal, "DualThreshold")
    upper_threshold = ADComponent(EpicsSignalWithRBV, "UpperThreshold")
    lower_threshold = ADComponent(EpicsSignalWithRBV, "LowerThreshold")

    output_control = ADComponent(EpicsSignalWithRBV, "OutputControl", string=True)
    trigger_edge = ADComponent(EpicsSignalWithRBV, "TriggerEdge", string=True)
    # trigger_mode replaces same (SignalWithRBV) in CamBase
    trigger_mode = ADComponent(EpicsSignalWithRBV, "TriggerMode", string=True)

    def __init__(self, *args, **kwargs):
        """Initialize the Rigaku3MCam.

        Args:
            *args: Variable length argument list.
            **kwargs: Arbitrary keyword arguments.

        After initialization, sets the kind attribute for various components
        to 'config'.
        """
        super().__init__(*args, **kwargs)

        # kind attribute must be set after Components are initialized
        attrs = """
            dual_threshold
            fast_file_name
            fast_file_path
            image_mode
            lower_threshold
            output_control
            output_resolution
            trigger_edge
            trigger_mode
            upper_threshold
        """.split()
        for attr in attrs:
            getattr(self, attr).kind = "config"


class SimDetectorCam_V34(CamMixin_V34, SimDetectorCam):
    """Revise SimDetectorCam for ADCore revisions."""


class ID8_PluginMixin(PluginBase_V34):
    """Remove property attribute found in AD IOCs now."""

    _asyn_pipeline_configuration_names = None


class ID8_CodecPlugin(ID8_PluginMixin, CodecPlugin_V34):
    """Remove property attribute found in AD IOCs now."""


class ID8_EpicsFileNameHDF5Plugin(ID8_PluginMixin, AD_EpicsFileNameHDF5Plugin):
    """Remove property attribute not found in AD IOCs now."""

    @property
    def _plugin_enabled(self):
        """Check if the plugin is enabled.

        Returns:
            bool: True if the plugin is enabled, False otherwise.
        """
        return self.stage_sigs.get("enable") in (1, "Enable")

    def generate_datum(self, *args, **kwargs):
        """Generate datum if plugin is enabled.

        Args:
            *args: Variable length argument list.
            **kwargs: Arbitrary keyword arguments.
        """
        if self._plugin_enabled:
            super().generate_datum(*args, **kwargs)

    def read(self):
        """Read data from the plugin if enabled.

        Returns:
            dict: Dictionary containing readings if plugin is enabled,
                empty dictionary otherwise.
        """
        if self._plugin_enabled:
            readings = super().read()
        else:
            readings = {}
        return readings

    def stage(self):
        """Stage the plugin if enabled.

        Returns:
            list: List of staged objects if plugin is enabled,
                empty list otherwise.
        """
        if self._plugin_enabled:
            staged_objects = super().stage()
        else:
            staged_objects = []
        return staged_objects

    def trigger(self):
        """Trigger the plugin if enabled.

        Returns:
            Status: Status object indicating completion.
        """
        if self._plugin_enabled:
            trigger_status = super().trigger()
        else:
            trigger_status = Status(self)
            trigger_status.set_finished()
        return trigger_status


class ID8_ImagePlugin(ID8_PluginMixin, ImagePlugin_V34):
    """Remove property attribute found in AD IOCs now."""


class ID8_OverlayPlugin(ID8_PluginMixin, OverlayPlugin_V34):
    """Remove property attribute found in AD IOCs now."""


class ID8_ProcessPlugin(ID8_PluginMixin, ProcessPlugin_V34):
    """Remove property attribute found in AD IOCs now."""


class ID8_PvaPlugin(ID8_PluginMixin, PvaPlugin_V34):
    """Remove property attribute found in AD IOCs now."""


class ID8_ROIPlugin(ID8_PluginMixin, ROIPlugin_V34):
    """Remove property attribute found in AD IOCs now."""


class ID8_StatsPlugin(ID8_PluginMixin, StatsPlugin_V34):
    """Remove property attribute found in AD IOCs now."""


class ID8_TransformPlugin(ID8_PluginMixin, TransformPlugin_V34):
    """Remove property attribute found in AD IOCs now."""
