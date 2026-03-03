
import time

####
# huber.sample_x.move(0.211)
# huber.sample_y.move(42.185)
# att(20)
select_sample(1)
pv_registers.qmap_file.put('lambda2m-qmap-A0015_CsPbBr3_PeakB_0p2GPA.hdf')
huber.delta.move(9.95)
lambda_acq_ext_series(acq_time=1.00, acq_period=2.00, num_frames=600, num_reps=1)

select_sample(2)
pv_registers.qmap_file.put('lambda2m-qmap-C0013_CsPbBr3_PeakC_0p2GPA.hdf')
huber.delta.move(14.70)
lambda_acq_ext_series(acq_time=1.00, acq_period=2.00, num_frames=600, num_reps=1)

select_sample(3)
pv_registers.qmap_file.put('lambda2m-qmap-C0014_CsPbBr3_PeakA_0p2GPA.hdf')
huber.delta.move(20.20)
lambda_acq_ext_series(acq_time=1.00, acq_period=2.00, num_frames=600, num_reps=1)

###
# huber.sample_x.move(0.201)
# huber.sample_y.move(42.20)

# select_sample(1)
# pv_registers.qmap_file.put('lambda2m-qmap-A0015_CsPbBr3_PeakB_0p2GPA.hdf')
# huber.delta.move(9.85)
# lambda_acq_ext_series(acq_time=1.00, acq_period=2.00, num_frames=1500, num_reps=1)

# select_sample(2)
# pv_registers.qmap_file.put('lambda2m-qmap-C0013_CsPbBr3_PeakC_0p2GPA.hdf')
# huber.delta.move(14.55)
# lambda_acq_ext_series(acq_time=1.00, acq_period=2.00, num_frames=600, num_reps=1)

# select_sample(3)
# pv_registers.qmap_file.put('lambda2m-qmap-C0014_CsPbBr3_PeakA_0p2GPA.hdf')
# huber.delta.move(19.98)
# lambda_acq_ext_series(acq_time=1.00, acq_period=2.00, num_frames=1500, num_reps=1)

# ###
# huber.sample_x.move(0.1272)
# huber.sample_y.move(42.25)

# select_sample(1)
# pv_registers.qmap_file.put('lambda2m-qmap-A0015_CsPbBr3_PeakB_0p2GPA.hdf')
# huber.delta.move(9.85)
# lambda_acq_ext_series(acq_time=1.00, acq_period=2.00, num_frames=1500, num_reps=1)

# select_sample(2)
# pv_registers.qmap_file.put('lambda2m-qmap-C0013_CsPbBr3_PeakC_0p2GPA.hdf')
# huber.delta.move(14.21)
# lambda_acq_ext_series(acq_time=1.00, acq_period=2.00, num_frames=1500, num_reps=1)

# select_sample(3)
# pv_registers.qmap_file.put('lambda2m-qmap-C0014_CsPbBr3_PeakA_0p2GPA.hdf')
# huber.delta.move(19.98)
# lambda_acq_ext_series(acq_time=1.00, acq_period=2.00, num_frames=1500, num_reps=1)
