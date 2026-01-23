
import yaml 
import re

from id8_i.plans.nexus_acq_eiger_int import eiger_acq_int_series
from id8_i.plans.sample_info_unpack import select_sample
from id8_i.plans.select_detector import select_detector
from id8_i.plans.scan_8idi import att
from id8_i.plans.nexus_acq_eiger_ext import eiger_acq_ext_trig

def run_measurement_info(file_name="measurement_info.yaml"):
    """Run measurement based on measurement_info.yaml file."""

    file_path = "/home/beams10/8IDIUSER/bluesky/src/user_plans/"

    try:
        with open(file_path + file_name, "r") as f:
            measurement_info = yaml.safe_load(f)

        for block_key, block_value in measurement_info.items():
                        
            match = re.search(r"_(\d+)", block_key)
            sam_index = int(match.group(1))
            det_name = block_value.get("detector")
            att_list = block_value.get("att_list")
            acq_time_list = block_value.get("acq_time_list")
            acq_period_list = block_value.get("acq_period_list")
            num_frames_list = block_value.get("num_frames_list")
            wait_time_list = block_value.get("wait_time_list")
            num_reps_list = block_value.get("num_reps_list")
            sample_move_yes_list = block_value.get("sample_move_yes_list")

            print(f"\n --- Measurement Block {block_key} ---")

            print(f"Sample index: {sam_index}")
            select_sample(sam_index)

            print(f"Detector name: {det_name}")
            select_detector(det_name)

            for ii in range(len(att_list)):
                print(f'\n At Attenuation Ratio {att_list[ii]}:\n')
                att(att_list[ii])

                for jj in range(len(acq_time_list[ii])):
                    acq_time = acq_time_list[ii][jj]
                    acq_period = acq_period_list[ii][jj]
                    num_frames = num_frames_list[ii][jj]
                    num_reps = num_reps_list[ii][jj]
                    wait_time = wait_time_list[ii][jj]
                    sample_move_yes = sample_move_yes_list[ii][jj]

                    print(f"    Acquisition Time: {acq_time}")
                    print(f"    Acquisition Period: {acq_period}")
                    print(f"    Number of Frames: {num_frames}")
                    print(f"    Number of Repeats: {num_reps}")

                    if det_name == "eiger4M":
                        if acq_time == acq_period:
                            print(f"Using eiger internal series")
                            eiger_acq_int_series(
                                acq_time=acq_time,
                                num_frames=num_frames,
                                num_reps=num_reps,
                                wait_time=wait_time,
                                sample_move=sample_move_yes,
                            )
                        elif acq_time != acq_period and acq_period >= 0.1:
                            print(f"Using eiger external trigger series")
                            eiger_acq_ext_trig(
                                acq_time=acq_time,
                                acq_period=acq_period,
                                num_frames=num_frames,
                                num_reps=num_reps,
                                wait_time=wait_time,
                                sample_move=sample_move_yes,
                            )
                        else:
                            print("Error: use acquition period larger than 0.1 s for Eiger Ext Trig mode")
                    
                    else:
                        print("Detector name must be eiger4M, rigaku3M, or aeon750k")

    except KeyboardInterrupt as err:
        raise RuntimeError("\n Bluesky plan stopped by user (Ctrl+C).") from err
    except Exception as e:
        print(f"Error occurred during measurement: {e}")
        raise Exception from e
    finally:
        pass

def run_round_robin(num_loops=1, filename="measurement_info.yaml"):

    for _ in range(num_loops):
        run_measurement_info(filename)