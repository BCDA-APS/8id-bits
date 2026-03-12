
##to align a new sample
att(1e6)

RE(bps.mv(huber.delta,0.0))
RE(bps.mv(huber.nu,0.0))

RE(dscan(motor=huber.sample_x, rel_begin=-1, rel_end=1, num_pts=20, count_time=1.0, det=lambda2M, att_ratio=1e6))
RE(dscan(motor=huber.sample_x, rel_begin=-1, rel_end=1, num_pts=20, count_time=1.0, det=lambda2M, att_ratio=1e6))


RE(bps.mv(huber.sample_x,0.0))
RE(bps.mv(huber.sample_y,38.4))

##to measure on lambda2m at a new angle (away from the direct beam)
att(1)
RE(bps.mv(huber.nu,13.7))
RE(bps.mv(huber.delta,10.0))
###collect lambda acq


go back to eiger 

RE(bps.mv(huber.nu,13.7))
RE(bps.mv(huber.delta,10