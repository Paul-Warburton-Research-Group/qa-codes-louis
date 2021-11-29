import numpy as np

# Transverse Field Schedules
def A(t,p):
    if p["stype"] == "linear":
        if t < 0:
            return 1.0
        elif t >= 0 and t <= p["tan"]:
            return (1.0-t/p["tan"])
        else:
            return 0.0
    elif p["stype"] == "hsine":
        if t < 0:
            return 0.0
        elif t >= 0 and t <= p["tan"]:
            return np.sin(t/p["tan"]*np.pi)**p["dist"]
        else:
            return 0.0
    elif p["stype"] == "interf":
        if t < 0:
            return 1.0
        elif t >= 0 and t/p["tan"] < p["smg"]+p["sw"]/2:
            return 1.0-t/p["tan"]
        elif t/p["tan"] >= p["smg"]+p["sw"]/2 and t/p["tan"] < p["smg"]+1.5*p["sw"]:
            return t/p["tan"]-(p["smg"]+p["sw"]/2)*2+1.0
        elif t/p["tan"] >= p["smg"]+1.5*p["sw"] and t/p["tan"] <= p["smg"]+1.5*p["sw"]+1-p["smg"]+p["sw"]/2:
            return -t/p["tan"] + (p["smg"]+1.5*p["sw"]+1-p["smg"]+p["sw"]/2)
        elif t/p["tan"] > p["smg"]+1.5*p["sw"]+1-p["smg"]+p["sw"]/2:
            return 0.0

# Problem Field Schedules
def B(t,p):
    if p["stype"] == "linear":
        if t < 0:
            return 0.0
        elif t >= 0 and t <= p["tan"]:
            return (t/p["tan"])
        else:
            return 1.0
    elif p["stype"] == "hsine":
        if t < 0:
            return 1.0
        elif t >= 0 and t <= p["tan"]:
            return 1-np.sin(t/p["tan"]*np.pi)
        else:
            return 1.0
    elif p["stype"] == "interf":
        if t < 0:
            return 0.0
        elif t >= 0 and t/p["tan"] < p["smg"]+p["sw"]/2:
            return t/p["tan"]
        elif t/p["tan"] >= p["smg"]+p["sw"]/2 and t/p["tan"] < p["smg"]+1.5*p["sw"]:
            return -t/p["tan"]+(p["smg"]+p["sw"]/2)*2
        elif t/p["tan"] >= p["smg"]+1.5*p["sw"] and t/p["tan"] <= p["smg"]+1.5*p["sw"]+1-p["smg"]+p["sw"]/2:
            return 1+t/p["tan"] - (p["smg"]+1.5*p["sw"]+1-p["smg"]+p["sw"]/2)
        elif t/p["tan"] > p["smg"]+1.5*p["sw"]+1-p["smg"]+p["sw"]/2:
            return 1.0
