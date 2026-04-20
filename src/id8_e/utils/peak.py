import numpy as np
from bluesky.plan_stubs import mv, mvr
from id8_e.plans.scan_8ide import dscan


def rock_and_move(RE, cat, huber, lambda2M, start=-0.2, stop=0.2, num=80, exp_time=1, att_ratio=20, method="max"):
    RE(dscan(huber.mu, start, stop, num, exp_time, det=lambda2M, att_ratio=att_ratio))

    last_scan_id = cat.values().last().metadata["start"]["scan_id"]
    run = cat[last_scan_id]
    ds = run.primary.read()

    motor_vals = ds["huber_mu"].values
    det_vals = ds["lambda2M_stats1_total"].values
    # det_vals = ds["eiger4M_stats1_total"].values

    if method == "max":
        idx = det_vals.argmax()
        peak_pos = motor_vals[idx]
    elif method == "com":
        peak_pos = np.sum(motor_vals * det_vals) / np.sum(det_vals)

    print(f"Peak at huber_mu = {peak_pos:.4f}  ({method})")
    RE(mv(huber.mu, peak_pos))

def center_x(RE, cat, huber,lambda2M, start=-0.2, stop=0.2, num=80, exp_time=1, att_ratio=20, method="max"):
    RE(dscan(huber.sample_x, start, stop, num, exp_time, det=lambda2M, att_ratio=att_ratio))

    last_scan_id = cat.values().last().metadata["start"]["scan_id"]
    run = cat[last_scan_id]
    ds = run.primary.read()

    motor_vals = ds["huber_sample_x"].values
    det_vals = ds["lambda2M_stats1_total"].values

    if method == "max":
        idx = det_vals.argmax()
        peak_pos = motor_vals[idx]
    elif method == "com":
        peak_pos = np.sum(motor_vals * det_vals) / np.sum(det_vals)

    print(f"Peak at sample_x = {peak_pos:.4f}  ({method})")
    RE(mv(huber.sample_x, peak_pos))

def center_y(RE, cat, huber, lambda2M, start=-0.2, stop=0.2, num=80, exp_time=1, att_ratio=20, method="max"):
    RE(dscan(huber.sample_y, start, stop, num, exp_time, det=lambda2M, att_ratio=att_ratio))

    last_scan_id = cat.values().last().metadata["start"]["scan_id"]
    run = cat[last_scan_id]
    ds = run.primary.read()

    motor_vals = ds["huber_sample_y"].values
    det_vals = ds["lambda2M_stats1_total"].values

    if method == "max":
        idx = det_vals.argmax()
        peak_pos = motor_vals[idx]
    elif method == "com":
        peak_pos = np.sum(motor_vals * det_vals) / np.sum(det_vals)

    print(f"Peak at sample_y = {peak_pos:.4f}  ({method})")
    RE(mv(huber.sample_y, peak_pos))

def center_delta(RE, cat, huber, lambda2M, start=-0.2, stop=0.2, num=80, exp_time=1, att_ratio=20, method="max"):
    RE(dscan(huber.delta, start, stop, num, exp_time, det=lambda2M, att_ratio=att_ratio))

    last_scan_id = cat.values().last().metadata["start"]["scan_id"]
    run = cat[last_scan_id]
    ds = run.primary.read()

    motor_vals = ds["huber_delta"].values
    det_vals = ds["lambda2M_stats1_total"].values

    if method == "max":
        idx = det_vals.argmax()
        peak_pos = motor_vals[idx]
    elif method == "com":
        peak_pos = np.sum(motor_vals * det_vals) / np.sum(det_vals)

    print(f"Peak at sample_y = {peak_pos:.4f}  ({method})")
    RE(mv(huber.delta, peak_pos))

def center_nu(RE, cat, huber, lambda2M, start=-0.2, stop=0.2, num=80, exp_time=1, att_ratio=20, method="max"):
    RE(dscan(huber.nu, start, stop, num, exp_time, det=lambda2M, att_ratio=att_ratio))

    last_scan_id = cat.values().last().metadata["start"]["scan_id"]
    run = cat[last_scan_id]
    ds = run.primary.read()

    motor_vals = ds["huber_nu"].values
    det_vals = ds["lambda2M_stats1_total"].values

    if method == "max":
        idx = det_vals.argmax()
        peak_pos = motor_vals[idx]
    elif method == "com":
        peak_pos = np.sum(motor_vals * det_vals) / np.sum(det_vals)

    print(f"Peak at sample_y = {peak_pos:.4f}  ({method})")
    RE(mv(huber.nu, peak_pos))


