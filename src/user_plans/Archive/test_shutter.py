    
import subprocess

file_name = "A001_001"
remote_cmd = f"LD_PRELOAD=libvma.so /home/xspadmin/tpx_env/daq_aps/acquire_event_based-shuttercontrol.py -s 1 -ds 100000 -th1 7476 -th3 7585 -fp /mnt/ramdisk/{file_name}.bin"
subprocess.run(["ssh", "xspadmin@164.54.116.167", remote_cmd], check=True)