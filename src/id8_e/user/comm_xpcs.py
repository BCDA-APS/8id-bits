import time

#RE(bps.mv(huber.mu,8.084))
#time.sleep(90)
#keysight.func.set("SQU").wait()
#time.sleep(5)

lambda2M.hdf1.file_template.put('%s%s.h5')

lambda_acq_int_series(acq_time=1.00, num_frames=100,num_reps=1)

