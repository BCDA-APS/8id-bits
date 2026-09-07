
from apsbits.core.instrument_init import oregistry
from bluesky import plan_stubs as bps
from legacy.id8_i.plans.master_plan import run_measurement_info
from id8_common.devices import *
from legacy.id8_i.plans import *
import subprocess


def test_aeon():

    yield from bps.sleep(0.1)

    remote_cmd = f"LD_PRELOAD=libvma.so /home/xspadmin/tpx_env/daq_aps/acquire_event_based.py -s 10 -ds 10000 -th1 7476 -th3 7585 -fp /mnt/ramdisk/test/test2.bin"

    subprocess.run(["ssh", "root@164.54.116.167", remote_cmd], check=True)

    local_1 = f"/ramdisk/test/test21-0.bin"
    local_2 = f"/ramdisk/test/test22-0.bin"
    local_3 = f"/ramdisk/test/test23-0.bin"


    mount_point = "/gdata/dm/8ID/8IDI/2025-3/timepix202512/data/"
    folder_name = "/A001_test/A001_test_001/"
    subprocess.run(["mkdir", "-p", f"{mount_point}/{folder_name}"], check=True)

    dest_1 = f"{mount_point}/{folder_name}/A001_test_001.tpx.000"
    dest_2 = f"{mount_point}/{folder_name}/A001_test_001.tpx.001"
    dest_3 = f"{mount_point}/{folder_name}/A001_test_001.tpx.002"

    # subprocess.run(["mv", "-f", "-t", dest, local_1, local_2, local_3], check=True)
    subprocess.run(["mv", "-f", local_1, dest_1], check=True)
    subprocess.run(["mv", "-f", local_2, dest_2], check=True)
    subprocess.run(["mv", "-f", local_3, dest_3], check=True)
    

