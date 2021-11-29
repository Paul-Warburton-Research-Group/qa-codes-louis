import qutip as qt
import numpy as np
def solve():
    a = np.zeros((512,512))
    a[501] = a[32] = a[332] = a[154] = a[451] = 1
    aq = qt.Qobj(a,isherm=True)
    return aq.eigenstates()
