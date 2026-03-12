def mv_samx(x):
    return RE(bps.mv(huber.sample_x,x))

def mv_samy(y):
    return RE(bps.mv(huber.sample_y,y))

def mv_nu(nu):
    return RE(bps.mv(huber.nu,nu))



def mv_delta(delta):
    return RE(bps.mv(huber.delta,delta))

def go_align():
    att(1e6)
    mv_delta(0)
    mv_nu(0)

def go_eiger():
    mv_delta(10)
    mv_nu(13.7)
    att(1)


def dscan_samx(rel_begin, rel_end, num_pts, count_time):
    RE(dscan(
        motor=huber.sample_x,
        rel_begin=rel_begin,
        rel_end=rel_end,
        num_pts=num_pts,
        count_time=count_time,
        det=lambda2M,
        att_ratio=1e6
    ))

def dscan_samy(rel_begin, rel_end, num_pts, count_time):
    RE(dscan(
        motor=huber.sample_y,
        rel_begin=rel_begin,
        rel_end=rel_end,
        num_pts=num_pts,
        count_time=count_time,
        det=lambda2M,
        att_ratio=1e6
    ))

def wm_samx():
    print(f"{huber.sample_x.position:4.2f}")

def wm_samy():
    print(f"{huber.sample_y.position:4.2f}")