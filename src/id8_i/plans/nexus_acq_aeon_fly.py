
from apsbits.core.instrument_init import oregistry
from bluesky import plan_stubs as bps
import subprocess

from id8_common.devices import *
from id8_i.utils.dm_util import dm_run_job
from id8_i.utils.dm_util import dm_setup
from id8_i.utils.nexus_utils import create_nexus_format_metadata

from id8_i.plans.sample_info_unpack import gen_folder_prefix
from id8_i.plans.sample_info_unpack import mesh_grid_move
from id8_i.plans.shutter_logic import blockbeam
from id8_i.plans.shutter_logic import post_align
from id8_i.plans.shutter_logic import showbeam
from id8_i.plans.shutter_logic import shutteroff


rigaku3M = oregistry["rigaku3M"]
pv_registers = oregistry["pv_registers"]
sample = oregistry["sample"]


def aeon_acquire_fly(file_header, file_name, fly_speed = 0.1):

    ## Setup variables based on info from Reg
    cycle_name = pv_registers.cycle_name.get()
    exp_name = pv_registers.experiment_name.get()
    use_subfolder = pv_registers.use_subfolder.get()

    if use_subfolder == 'Yes':
        file_path = f"{exp_name}/data/{file_header}/{file_name}"
    elif use_subfolder == 'No':   
        file_path = f"{exp_name}/data/{file_name}"
    else: 
        print("Sub folder options can only be either Yes or No")

    folder_path = f"/gdata/dm/8ID/8IDI/{cycle_name}/{file_path}"
    yield from bps.mv(pv_registers.file_name, file_name)
    yield from bps.mv(
        pv_registers.metadata_full_path,
        f"{folder_path}/{file_name}_metadata.hdf",
    )

    extra_acq_time = 6
    sample_y_orig = sample.y.user_setpoint.get()

    total_travel = fly_speed * (1.0 + extra_acq_time) # Acquisition time is fixed at 1 second

    # Acquisition. ssh into Aeon server and save data to Beryl ramdisk
    # Acquisition time is set to 1 second by default
    yield from showbeam()

    yield from bps.mv(sample.y.velocity, fly_speed)
    yield from bps.mvr(sample.y.user_setpoint, total_travel)

    remote_cmd = f"LD_PRELOAD=libvma.so /home/xspadmin/tpx_env/daq_aps/acquire_event_based-shuttercontrol.py -s 1 -ds 100000 -th1 7476 -th3 7585 -fp /mnt/ramdisk/{file_name}.bin"
    subprocess.run(["ssh", "root@164.54.116.167", remote_cmd], check=True)

    yield from blockbeam()

    yield from bps.mv(sample.y.velocity, 10)
    yield from bps.mv(sample.y, sample_y_orig)

    # Data transfer from Beryl ramdisk to /gdata 
    local_1 = f"/ramdisk/{file_name}1-0.bin"
    local_2 = f"/ramdisk/{file_name}2-0.bin"
    local_3 = f"/ramdisk/{file_name}3-0.bin"
    
    subprocess.run(["mkdir", "-p", f"{folder_path}"], check=True)

    dest_1 = f"{folder_path}/{file_name}.tpx.000"
    dest_2 = f"{folder_path}/{file_name}.tpx.001"
    dest_3 = f"{folder_path}/{file_name}.tpx.002"

    subprocess.run(["mv", "-f", local_1, dest_1], check=True)
    subprocess.run(["mv", "-f", local_2, dest_2], check=True)
    subprocess.run(["mv", "-f", local_3, dest_3], check=True)
    

def aeon_acq_series_fly(
    num_reps=2,
    wait_time=0,
    fly_speed=0.1,
    sample_move=False,
):

    # try:
    yield from post_align()
    yield from shutteroff()
    workflowProcApi, dmuser = dm_setup()
    folder_prefix = gen_folder_prefix()

    for ii in range(num_reps):
        if sample_move:
            yield from mesh_grid_move()

        file_header = f"{folder_prefix}"
        file_name = f"{folder_prefix}_r{ii+1:05d}"

        print(file_name)

        print(f"\nStarting Measurement {file_name}")
        yield from aeon_acquire_fly(file_header, file_name, fly_speed)
        print(f"Measurement {file_name} Complete")

        metadata_fname = pv_registers.metadata_full_path.get()
        create_nexus_format_metadata(metadata_fname, det=rigaku3M)

        dm_run_job(workflowProcApi, dmuser)

        yield from bps.sleep(wait_time)

    # except Exception as e:
    #     print(f"Error occurred during measurement: {e}")
    # finally:
    #     pass
