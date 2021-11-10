
from pyblq import Blq, ket, discard, maxmixed
from pyblq.polys import rect

I = Blq([[1,0],[0,1]])
Z = Blq([[1,0],[0,-1]])
a0 = (ket(0) * ~ket(1)) @ I
a1 = Z @ (ket(0) * ~ket(1))

H = 0.1 * a0 * ~a0 + 0.3 * (a0 * ~a1 + a1 * ~a0)

def threshold(E,eps):
    p = rect({(-H.scale, E-eps):0, (E+eps, H.scale):1})
    energy_threshold = (discard() @ discard()) * H.svt(p) * (maxmixed() @ maxmixed())
    probability_threshold = energy_threshold.svt(rect({0:0, (1/2,1):1},parity='odd'))
    return probability_threshold.sample() is not None

Emin, Emax, eps = -H.scale, H.scale, 0.2
while Emax - Emin > 2*eps:
    if threshold((Emax+Emin)/2, eps): Emax = (Emax+Emin)/2
    else: Emax = (Emax+Emin)/2
print((Emax+Emin)/2)

