'''
temperature-dependent XPCS macro for STO50 
'''

# set_temp_lakeshore(40, 1)
# lakeshore2.setpoint_out1.put(40)
# time.sleep(5)

# '''rock sample in eta and leave at peak'''
# rock_and_move(RE, cat, huber, lambda2M, start=-0.5, stop=0.5, num=80, exp_time=1, att_ratio=20)

# '''translate in y and leave at peak'''
# center_y(RE, cat, huber, lambda2M, start=-0.05, stop=0.05, num=80, exp_time=1, att_ratio=20)

# lambda_acq_int_series(acq_time=1, num_frames=10, num_reps=1)



# moving to (0.5, 0.5, 2.5)
# mu.move(3.2674)
# nu.move(8.9182)
# delta.move(18.1616)
# rock_and_move(RE, cat, huber, lambda2M, start=-0.5, stop=0.5, num=40, exp_time=1, att_ratio=10, method="com")
# att(20)
# lambda_acq_int_series(acq_time=1, num_frames=3600, num_reps=1)


# # set_temp_lakeshore(80, 3600)  # this command doesn't work somehow after restart bluesky
# lakeshore2.setpoint_out1.put(80)
# time.sleep(3600)

# # moving to (1, 1, 2)
# delta.move(24.538)
# nu.move(17.5754)
# mu.move(13.92)
# rock_and_move(RE, cat, huber, lambda2M, start=-0.5, stop=0.5, num=40, exp_time=1, att_ratio=10, method="com")
# att(10)
# lambda_acq_int_series(acq_time=1, num_frames=1800, num_reps=1)

# # moving to (0.5, 0.5, 2.5)
# mu.move(3.2674)
# nu.move(8.9182)
# delta.move(18.1616)
# rock_and_move(RE, cat, huber, lambda2M, start=-0.5, stop=0.5, num=40, exp_time=1, att_ratio=10, method="com")
att(20)
lambda_acq_int_series(acq_time=1, num_frames=3600, num_reps=1)


# set_temp_lakeshore(120, 3600)
lakeshore2.setpoint_out1.put(120)
time.sleep(3600)

# moving to (1, 1, 2)
delta.move(24.538)
nu.move(17.5754)
mu.move(13.92)
rock_and_move(RE, cat, huber, lambda2M, start=-0.5, stop=0.5, num=40, exp_time=1, att_ratio=10, method="com")
att(10)
lambda_acq_int_series(acq_time=1, num_frames=1800, num_reps=1)

# moving to (0.5, 0.5, 2.5)
mu.move(3.2674)
nu.move(8.9182)
delta.move(18.1616)
rock_and_move(RE, cat, huber, lambda2M, start=-0.5, stop=0.5, num=40, exp_time=3, att_ratio=1, method="max")
att(20)
lambda_acq_int_series(acq_time=1, num_frames=3600, num_reps=1)



# set_temp_lakeshore(140, 3600)
lakeshore2.setpoint_out1.put(140)
time.sleep(3600)

# moving to (1, 1, 2)
delta.move(24.538)
nu.move(17.5754)
mu.move(13.92)
rock_and_move(RE, cat, huber, lambda2M, start=-0.5, stop=0.5, num=40, exp_time=1, att_ratio=10, method="com")
att(10)
lambda_acq_int_series(acq_time=1, num_frames=1800, num_reps=1)

# moving to (0.5, 0.5, 2.5)
mu.move(3.2674)
nu.move(8.9182)
delta.move(18.1616)
rock_and_move(RE, cat, huber, lambda2M, start=-0.5, stop=0.5, num=40, exp_time=5, att_ratio=1, method="max")
att(20)
lambda_acq_int_series(acq_time=1, num_frames=3600, num_reps=1)



# set_temp_lakeshore(200, 3600)
lakeshore2.setpoint_out1.put(200)
time.sleep(3600)

# moving to (1, 1, 2)
delta.move(24.538)
nu.move(17.5754)
mu.move(13.92)
rock_and_move(RE, cat, huber, lambda2M, start=-0.5, stop=0.5, num=40, exp_time=1, att_ratio=10, method="com")
att(10)
lambda_acq_int_series(acq_time=1, num_frames=3600, num_reps=1)

