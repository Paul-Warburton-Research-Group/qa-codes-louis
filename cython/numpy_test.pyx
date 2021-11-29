import numpy as np
cimport numpy as np
local_dtype = np.float64
ctypedef np.float64_t local_dtype_t

def solve(np.int n):
    cdef np.ndarray a = np.zeros([n, n], dtype=local_dtype)
    a[501] = a[32] = a[332] = a[154] = a[451] = 1
    return np.linalg.eigvalsh(a)

