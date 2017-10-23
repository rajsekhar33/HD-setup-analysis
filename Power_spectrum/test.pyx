import numpy as np
cimport numpy as np
DTYPEi = np.int
DTYPE = np.double
ctypedef np.double_t DTYPE_t
ctypedef np.int_t DTYPEi_t

def test_2(double[:] u):
    cdef int i
    cdef int N = u.shape[0]
    cdef double[:] u1= np.zeros(N)
    for i in range(N):
        u1[i] = 2*u[i]
    return [u1[i] for i in range(N)]
