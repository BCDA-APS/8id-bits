def auto_attenuate(
    det,
    pilot_exptime: float = 0.05,
    rate_limit: float = 4e5,
    filter_factor: float = 5.0,
    retry_max: int = 10,
    grace_factor: float = 0.25,
):
    """Find the optimal attenuation using short pilot exposures.

    Takes pilot frames and reads det.stats1.max_value as the feedback signal.
    Adjusts filter_8ide.transmission iteratively until the count rate
    (max_pixel / pilot_exptime) falls within the acceptable band:

        rate > rate_limit                       -> reduce transmission by filter_factor, retry
        rate_limit*grace_factor < rate < limit  -> acceptable, stop
        rate <= rate_limit * grace_factor       -> raise transmission toward 0.75*rate_limit, retry

    Starts from maximum attenuation (transmission=1e-10) for safety.
    Leaves the filter at the converged setting and restores the detector
    acquire_time. Does not save data (HDF capture is not started).

    Args:
        det:            eiger4M or lambda2M
        pilot_exptime:  duration of each test frame (s)
        rate_limit:     max acceptable count rate (max pixel cts/s)
        filter_factor:  transmission multiplier per step, must be > 1
        retry_max:      max iterations before giving up
        grace_factor:   lower rate bound = rate_limit * grace_factor

    Example::

        auto_attenuate(eiger4M, pilot_exptime=0.05, rate_limit=4e5)
        RE(dscan(sample.x, -1, 1, 100, count_time=1.0))
    """
    is_eiger = ("eiger" in det.name.lower()) or ("eiger" in det.prefix.lower())
    low_rate = rate_limit * grace_factor

    orig_acq_time = det.cam.acquire_time.get()
    orig_acq_period = det.cam.acquire_period.get()

    det.cam.acquire_time.put(pilot_exptime)
    det.cam.acquire_period.put(pilot_exptime)

    if is_eiger:
        det.cam.trigger_mode.put("Internal Series")
        det.cam.num_images.put(1)
        det.cam.num_triggers.put(1)
        det.cam.manual_trigger.put("Disable")
    else:
        # lambda2M
        det.cam.trigger_mode.put("Internal")
        det.cam.num_images.put(1)

    det.stats1.enable.put(1)
    det.stats1.compute_statistics.put(1)

    # Start from maximum attenuation (minimum transmission) for safety
    filter_beam.transmission.move(1e-6)
    time.sleep(0.5)

    showbeam()
    try:
        for attempt in range(retry_max):
            det.cam.acquire.put(1)
            t0 = time.time()
            timeout = pilot_exptime * 5 + 2
            while det.cam.acquire.get() == 1:
                time.sleep(0.02)
                if time.time() - t0 > timeout:
                    print("  WARNING: pilot frame timed out")
                    break

            max_cts = det.stats1.max_value.get()
            rate = max_cts / pilot_exptime
            current_trans = filter_beam.transmission.readback.get()

            print(
                f"Attempt {attempt + 1}: trans={current_trans:.4f}"
                f"max_cts={max_cts:.0f}  rate={rate:.0f} cts/s"
            )

            if rate > rate_limit:
                new_trans = current_trans / filter_factor
                print(f"    Rate too high -> reducing transmission to {new_trans:.6f}")
                filter_beam.transmission.move(new_trans)
            elif rate < low_rate:
                if rate > 0:
                    new_trans = current_trans * (0.75 * rate_limit / rate)
                else:
                    new_trans = current_trans * filter_factor
                print(f"    Rate too low -> raising transmission to {new_trans:.4f}")
                filter_beam.transmission.move(new_trans)
            else:
                print(f"    Rate in [{low_rate:.0f}, {rate_limit:.0f}] cts/s -- converged.")
                break
        else:
            print(f"WARNING: auto_attenuate did not converge in {retry_max} attempts")
    finally:
        blockbeam()
        det.cam.acquire_time.put(orig_acq_time)
        det.cam.acquire_period.put(orig_acq_period)

    trans = filter_beam.transmission.readback.get()
    atten = filter_beam.attenuation.readback.get()
    print(f"  Final: transmission={trans:.4f}  attenuation={atten}")