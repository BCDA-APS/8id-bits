
import epics
import time

from apsbits.core.instrument_init import oregistry

tetramm1 = oregistry["tetramm1"]
tetramm2 = oregistry["tetramm2"]
tetramm3 = oregistry["tetramm3"]
tetramm4 = oregistry["tetramm4"]
pv_registers = oregistry["pv_registers"]
huber = oregistry["huber"]


print("\n 1. monitor the pv: S09-FOFB:ShortHistControlModeM. ")
print('\n start the below when it goes from 1 to 2 (i will get more details \n')
while True:
    #time.sleep(0.05)
    beam_status = epics.caget("S09-FOFB:ShortHistControlModeM")
    
    if beam_status != 2:
        time.sleep(0.05)
    if beam_status == 2:
        break

# horizontal edge (x,y): (5.94, 54)
print('\n 2. huber move xy to a value (knife x) \n')
huber.sample_x.move(5.939)
huber.sample_y.move(54)

print('\n 3. tetraam acq \n')
tetramm_acq_series(det=tetramm3, filename='fofboff_run_hori_h10hz', num_capture=10)

# vertical edge (x,y):(6.0, 54.182)
print('\n 4. huber move xy to a value (knife y) \n')
#huber.sample_x.move(6)
#no crl
huber.sample_x.move(6.0)
huber.sample_y.move(54.185)


print('\n 5. tetraam acq \n')
tetramm_acq_series(det=tetramm3, filename='fofboff_run_vert_h10hz', num_capture=10)

# knife out (x,y): (6, 54)
print('\n 6. huber move xy to a value (no knife) \n')

huber.sample_x.move(6)
#huber.sample_x.move(6.5)
huber.sample_y.move(54)
tetramm_acq_series(det=tetramm3, filename='fofboff_run_noKE_h10hz', num_capture=10)

print('\n 7. move crl-1 lenses out')
# epics.caput("8iddCRL:1xD:1:sortedIndex", 0)


print('\n 8. tetraam acq (with no knife - need to verify with Xianbo)')
# tetramm_acq_series(det=tetramm3, filename='fofbon_run_full', num_capture=10)


