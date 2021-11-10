
from pyblq import Blq, ket, discard
from pyblq.polys import rect

I = Blq([[1,0],[0,1]])
Z = Blq([[1,0],[0,-1]])
a0 = (ket(0) * ~ket(1)) @ I
a1 = Z @ (ket(0) * ~ket(1))

H = 0.1 * a0 * ~a0 + 0.3 * (a0 * ~a1 + a1 * ~a0)
norm = abs(H)
H = (H/norm + I @ I)/2

state = H * ~(discard() @ discard())/2
def threshold(E,eps):
    p = rect((E-eps,E+eps))
    energy_threshold = (discard() @ discard()) * state.svt(p)
    probability_threshold = energy_threshold.svt(rect((0.0,4**(-1/2))))
    print(probability_threshold)
    return probability_threshold.sample() is not None

Emin, Emax, eps = 0, 1, 0.2
while Emax - Emin > 2*eps:
    if threshold( (Emax+Emin)/2 , eps): Emax = (Emax+Emin)/2
    else: Emax = (Emax+Emin)/2
print((Emax+Emin)/2)

