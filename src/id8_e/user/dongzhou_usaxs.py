
import time

huber.eta.move(0)

att(50)

select_sample(1)
huber.delta.move(-1.3)
pv_registers.analysis_type.put('Multitau')
pv_registers.det_name.put('lambda2M')
pv_registers.qmap_file.put('lambda_m1p3Delta_Sq360_Dq36_lin.hdf')
lambda_acq_int_series(acq_time=1.0, num_frames=5, num_reps=1, sample_move=False)


select_sample(2)
huber.delta.move(4)
pv_registers.analysis_type.put('Multitau')
pv_registers.det_name.put('lambda2M')
pv_registers.qmap_file.put('lambda_4p0Delta_Sq360_Dq36_lin.hdf')
lambda_acq_int_series(acq_time=1.0, num_frames=5, num_reps=1, sample_move=False)


att(50)
select_sample(3)
huber.delta.move(20.0)
pv_registers.analysis_type.put('Multitau')
pv_registers.det_name.put('eiger4M')
pv_registers.qmap_file.put('eiger4m_qmap_Sq360_Sphi16_Dq36_Dphi1_log.hdf')
eiger_acq_int_series(acq_time=0.3, num_frames=5000, num_reps=1, sample_move=False)


