import time
#RE(bps.mv(huber.mu,8.084))
#time.sleep(90)
##pv_registers.qmap_file.put('lambda2m_qmap_S30x1_D3x1.hdf')
#RE(bps.mv(huber.sample_x,-3.8))
###WAXS
#lambda_acq_int_series(acq_time=1,num_frames=100,num_reps=1)
# eiger_acq_int_series(acq_time=0.01,num_frames=10000,num_reps=1)
###XPCS
pv_registers.sample_name.put(f'code_test_030926')
eiger_acq_ext_trig(acq_time=0.000995, acq_period=0.001, num_frames=10000,num_reps=1,segment_length=10000)
#eiger_acq_ext_trig(acq_time=0.009995, acq_period=0.01, num_frames=100000,num_reps=1,segment_length=1000)
# eiger_acq_ext_trig(acq_time=0.00995, acq_period=0.01, num_frames=100000,num_reps=1,segment_length=100)

# for j in range(10):
#     # slot time set
#     lambda2M.hdf1.file_template.put('%s%s.h5')
#     pv_registers.qmap_file.put('eiger4m_qmap_S30x1_D3x1_new.hdf')
#     #set name attenuation and run
#     pv_registers.sample_name.put(f'dry_3_5mil_m50C_att10_{j:d}')
#     att(10)
#     eiger_acq_ext_trig(acq_time=0.00995, acq_period=0.01, num_frames=100000,num_reps=1,segment_length=100)
#     #
#     # fast time set 
#     #
#     lambda2M.hdf1.file_template.put('%s%s.h5')
#     pv_registers.qmap_file.put('eiger4m_qmap_S30x1_D3x1_new.hdf')
#     #set name attenuation and run
#     pv_registers.sample_name.put(f'dry_3_5mil_m50C_att100_{j:d}')
#     att(100)
#     eiger_acq_ext_trig(acq_time=0.0995, acq_period=0.1, num_frames=10000,num_reps=1,segment_length=100)


