import time
lambda2M.hdf1.file_template.put('%s%s.h5')


#RE(bps.mv(huber.mu,8.084))
#time.sleep(90)

##pv_registers.qmap_file.put('lambda2m_qmap_S30x1_D3x1.hdf')
pv_registers.qmap_file.put('eiger4m_qmap_S30x1_D3x1_new.hdf')

#RE(bps.mv(huber.sample_x,-3.8))
pv_registers.sample_name.put('delete2')
att(1)

###WAXS
#lambda_acq_int_series(acq_time=1,num_frames=100,num_reps=1)
#eiger_acq_int_series(acq_time=1,num_frames=10,num_reps=1)
###XPCS


eiger_acq_ext_trig(acq_time=0.000995, acq_period=0.001, num_frames=20000,num_reps=1,segment_length=500)


