import time

#RE(bps.mv(huber.mu,8.084))
#time.sleep(90)
#keysight.func.set("SQU").wait()
#time.sleep(5)

lambda2M.hdf1.file_template.put('%s%s.h5')

time.sleep(5)
keysight.frequency.set(15000)
keysight.amplitude.set(5.0)
keysight.burst_count.set(24000)
keysight.output.set(1).wait
lambda_acq_ext_series(acq_time=1, acq_period=8, num_frames=3030,num_reps=1,frames_before_voltage=30)
keysight.output.set(0).wait