import time
from apsbits.core.instrument_init import oregistry

pv_registers = oregistry["pv_registers"]

try:
    # Loop through a range (increased to 100 to give you time to press Ctrl+C)
    for i in range(100):
        pv_registers.sample0_pos.put(i)
        print(f"Set sample0_pos to {i}")
        
        # Sleep to make the loop visible and give time to interrupt
        time.sleep(0.5)

except KeyboardInterrupt:
    print("\nLoop interrupted by user (Ctrl+C).")
    print("\nRunning clean-up macros")
    
    # Execute the cleanup command
    pv_registers.measurement_num.put(100)
    print("\nDone.")