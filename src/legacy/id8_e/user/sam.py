
total_runs = 30

runs = range(0, total_runs)

for run in runs: 
    print("run ", run)
    RE(bps.mvr(huber.sample_x, 0.008))
    eiger_acq_int_series(acq_time=0.5, num_frames=3600, num_reps=1)
