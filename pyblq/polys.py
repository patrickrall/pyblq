
import numpy as np
from numpy.polynomial import Chebyshev as Cheb
from scipy.special import jv, iv

################## Hamiltonian Simulation

# equation 52 of 1806.01838
# R = floor( r( e*|t|/2 , (5/4)*eps )/2 )
# where r(t',eps') is the smallest r > t such that eps' = (t'/r)**r
def ja_R(t,eps):
    eps_p = (5/4)*eps
    t_p = np.e*abs(t)/2

    # R is the largest integer below this threshold
    def thresh(R):
        r = 2*R
        if r == 0: return False
        return eps_p >= (t_p/r)**r

    # binary search an interval Rmin,Rmax containing R
    Rmin = int(np.ceil(t_p))
    Rmax = 2*Rmin
    while True:
        if thresh(Rmax): break
        Rmin *= 2
        Rmax *= 2

    # binary search within that interval
    while (Rmax - Rmin) > 1:
        Rmed = int((Rmax + Rmin)/2)
        if thresh(Rmed): Rmax = Rmed
        else: Rmin = Rmed

    return Rmin

# convenience function for making chebyshev 'monomials'
def cheb(n):
    assert int(n) == n
    n = int(n)
    assert n >= 0
    coefs = np.zeros(n+1)
    coefs[-1] = 1
    return Cheb(coefs)


# Lemma 57 of 1806.01838
# Can also just specify the degree 2*R+1
#   to do so, set eps=None
def ja_sin(t,eps=0.01,R=None):
    assert eps is None or R is None
    assert not (eps is None and R is None)
    if R is None: R = ja_R(t,eps)

    out = Cheb([0])
    for k in range(R+1):
        out += (-1)**k * jv(2*k+1,t) * cheb(2*k+1)

    return 2*out

# Lemma 57 of 1806.01838
# Can also just specify the degree 2*R
#   to do so, set eps=None
def ja_cos(t,eps=0.01,R=None):
    assert eps is None or R is None
    assert not (eps is None and R is None)
    if R is None: R = ja_R(t,eps)

    # note that Lemma 57 actually has a typo in it
    # it should be | cos(tx) - J_0(t) - 2 sum_k^R (-1)^k J_2k(t) T_2k(x) | < eps
    # but second minus sign is written as a plus sign
    out = Cheb([jv(0,t)])
    for k in range(1,R+1):
        out += 2*(-1)**k * jv(2*k,t) * cheb(2*k)

    return out

################# Thresholding


# Corollary 4 of 1707.05391
def erf(k,eps=0.01):
    assert eps < 1/np.sqrt(np.pi) # probably should be even less than this

    # bound on eps_exp from Lemma 13, equation 66
    def exp_degree(eps_p,beta_p):
        t = np.ceil( max( beta_p * np.exp(2) , np.log( 4/eps_p)  )  )
        return np.ceil( np.sqrt(  2 * t * np.log(4/eps_p)  ) )

    # Corollary 3, equation 69
    def gauss_degree(eps_p, gamma_p):
        return exp_degree(eps_p  , gamma_p**2 / 2 ) * 2

    # Corollary 4, equation 70
    # need to find n such that 4k/(sqrt(pi) n ) * gauss_error(k, n-1) < eps
    # if we have converged then n = gauss_degree(k, eps * sqrt(pi) * n / (4*k)) - 1
    n = 1 # initial guess for n
    while True:
        if np.isnan(n): assert False
        new_n = gauss_degree(eps * np.sqrt(np.pi) * n / (4*k), k) - 1
        if n > 1 and new_n >= n: break
        n = new_n

    out = Cheb([0,iv(0,k**2/2)])

    for j in range(1,int((n-1)/2)+1):
        out += (-1)**j * iv(j, k**2/2) * ( cheb(2*j+1)/(2*j+1) - cheb(2*j-1)/(2*j-1) )

    return out*2*k*np.exp(-k**2/2)/np.sqrt(np.pi)



# Lemma 10 of 1707.05391
# approximates sign(x) for |x| > kappa/2
def sign(kappa,eps=0.01):
    eps_half = eps/2
    k = (np.sqrt(2)/kappa) * np.log(2/(np.pi * eps_half**2))**(1/2)
    return erf(k, eps_half)

# General rectangle function for the region 0 < x < 1
# If odd=True then P(0)=0 and P is odd
#   otherwise P(1)=1 and P is even
# Input a list of tuples (0.2,0.3) indicating regions
#   where the polynomial switches between 0 and 1
def rect(*switch_regions,odd=True,eps=0.01):

    if odd:
        out = Cheb([0])
        low = True
    else:
        out = Cheb([1])
        low = False

    if len(switch_regions) == 0:
        return out

    eps /= len(switch_regions)

    pos = 0
    for (l,r) in switch_regions:
        assert pos <= l and l < r and r <= 1
        m = (l+r)/2
        kappa = r-l
        # can do a little bit better than doubling the degree
        # as done in Corollary 5
        rescale = max(abs(1-m), abs(m+1))
        kappa /= rescale
        Sign = sign(kappa, eps=eps)(Cheb([-m/rescale,1/rescale]))

        if odd: Poly = (Sign - Sign(Cheb([0,-1])))/2
        else: Poly = (Sign + Sign(Cheb([0,-1])))/2 + 1

        if low: out += Poly
        else: out -= Poly
        low = not low

    return out


# Theorem 10 of 1707.05391
# approximates k*x for |x| < 1/k
def lin(k,eps=0.01):
    assert False # Finish me
    #mid = (1 + 1/k)/2
    mid = 1/k + 0.1
    kappa = mid - 1/k
    stretch = 5
    Sign = sign(kappa/stretch, eps=eps)(Cheb([-mid/stretch,1/stretch]))
    Rect = -(Sign(Cheb([0,-1])) + Sign)/2
    return Rect * Cheb([0,k])

    return rect((1/k,1),odd=False,eps=eps) * Cheb([0,k])


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    xs = np.linspace(-1,1,100)
    Rect = rect((0.4,0.6),eps=0.4)
    #Lin = lin(3, eps=0.4)
    plt.plot(xs,Rect(xs))
    #plt.plot(xs,3*xs)
    plt.grid()
    plt.show()



# Method one, for isometries
# nudge scale down to nearest solution
# apply Chebyshev

# Method two, for arbitary matrices
# use lin polynomial

# also implement hermitian svt


################ Other

# TODO: I will implement these in a future version.

# approximates x**d using degree O(sqrt(d))
# Theorem 3.3 of "Faster Algorithms via Approximation Theory" Vishnoi Sachadeva
# (also stated as Theorem 63 of 1806.01838)
def fast_forwarding(d,eps=0.01):
    pass

# approximates e**(-beta(1-x))
# Lemma 4.2 of "Faster Algorithms via Approximation Theory" Vishnoi Sachadeva
# (also stated as Corollary 64 of 1806.01838)
# Doesn't this use a taylor series?
def exponential(eps=0.01):
    pass


# Corollary 66 of 1806.01838
# Read the corollary statement for a full description
# coef_lambda computes the coefficients a_l = coef_lambda(l)
# B is the bound satisfying sum_l |a_l| (r + delta)^l <= B
def taylor_series(coef_lambda,B,eps=0.01):
    pass

# Odd polynomial that approximates f(x) = (3/4) (eps / x) for |x| > eps
# Corollary 69 of 1806.01838
def reciprocal(eps=0.01):
    pass

# approximates (2/pi) * arcsin(x)
# Lemma 70 of 1806.01838
def arcsin(eps=0.01):
    pass







