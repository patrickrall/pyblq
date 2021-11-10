

from pyblq import *
from pyblq.polys import rect

n = 5

def a(k): # Jordan Wigner transformation
    out = 1
    for i in range(k): out = out @ Blq([[1,0],[0,1]])
    out = out @ Blq([[0,1],[0,0]])
    for i in range(n-k-1): out = out @ Blq([[1,0],[0,-1]])
    return out

# Assemble Hamiltonian
t = 10
h = 5
H = t * (a(n-1) * ~a(0) + a(0) * ~a(n-1))
for k in range(0,n-1): H += h * (a(k) * ~a(k+1) + a(k+1) * ~a(k))
H += h * (a(0) * ~a(0))
for k in range(1,n): H += h * (a(k) * ~a(k))

# Normalize
norm = abs(H)
H /= norm

# Make positive
H = (H + Blq([[1,0],[0,1]]).tensorpower(n))/2

eps = 0.1

def thresh_poly(E):
    assert E - eps > 0
    assert E + eps < 1
    return rect((E-eps,E+eps))

n_discard = discard(2)
for i in range(n-1): n_discard = n_discard @ discard(2)

out = H * (~n_discard)
out = out.svt(thresh_poly(0.5))


#test = n_discard * H * (~n_discard)
#print(test.kraus())



