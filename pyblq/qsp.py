

import numpy as np
from numpy.polynomial import Chebyshev as Cheb
from numpy.polynomial import Polynomial as Poly
from numpy.polynomial.polynomial import polyvalfromroots

# this implementation is based on: 1806.01838, 1806.10236, 2003.02831, 2105.02859
# as well as github.com/ichuang/pyqsp and github.com/alibaba-edu/angle-sequence

# new method is in preparation at this time:
# https://simons.berkeley.edu/talks/optimization-based-approach-quantum-signal-processing-and-its-energy-landscape
# new implementations should consider using that new method once available
# https://github.com/qsppack/QSPPACK

# for debugging
def plot_poly(*Ps):
    import matplotlib.pyplot as plt
    xs = np.linspace(-2,2,500)
    sigmas = np.linspace(-2,2,500)
    ws = np.exp(1j*np.arccos(sigmas))
    delta = 0
    for (P,c) in zip(Ps,plt.cm.rainbow(np.linspace(0,1,len(Ps)))):
        if isinstance(P,list):
            P = P[0]
            d = int((len(P.coef)-1)/2)
            plt.plot(sigmas,[delta+np.real(P(w)*w**(-d)) for w in ws],c=c)
            plt.plot(sigmas,[delta+np.imag(P(w)*w**(-d)) for w in ws],c=c,ls="--")
        else:
            plt.plot(xs,np.real(delta+P(xs)),c=c)
            plt.plot(xs,np.imag(delta+P(xs)),c=c,ls="--")
        delta += 0.01

    plt.ylim(-3,3)
    plt.grid()
    plt.show()

# set any tiny entries to 0
# prune high degree terms with coefficient close to 0
def sanitize_polynomial(P):
    d = len(P.coef)-1
    while np.allclose(P.coef[d],0): d -= 1
    new_P_coefs = []
    for i in range(d+1):
        if np.allclose(P.coef[i],0): new_P_coefs.append(0)
        else:
            new_P_coefs.append(P.coef[i])
    return Cheb(new_P_coefs)

# Given numpy polynomial with fixed parity, do quantum signal processing.
# Inputs:
#   P - A numpy polynomial. Must be either completly even or odd.
#   tol - An error parameter for several np.allclose checks in the algorithm.
#           Don't think of this as a guarantee that the output poly will be within tol of P.
#           Better check that yourself.
#   completion_mode - how to pick the polynomial that determines the other matrix elements
#           None - Try to do PQ completion, if that fails do FG completion.
#           "PQ" - Only possible if |P(x)| <= 1 for |x| <= 1
#                                   |P(x)| >= for |x| >= 1
#                                   if P is even P(ix)P*(ix) > 1 for all real x
#                  Output polynomial is P(x) exactly.
#           "FG" - Only possible if P has real coefficients and |P(x)| <= 1 for |x| <= 1.
#                  Output polynomial has an extra imaginary part P(x) + i P'(x).
#   allow_scale
#           True - Attempt to rescale P so that PQ or FG completion is possible.
#           False - Force implementation of P as is.
#   verify
#           True - Enable various assertions that test if things are going smoothly.
#           False - Return *anything* even if it is a mess.

# Outputs:
# scale - the scale of the SVT block-encoding to actually achieve the polynomial
# global_phi - a global phase correction on the final unitary
# phis - a list of deg(P)-1 many angles to implement the SVT circuit
#      - put these in between the applications of U
# completion_mode - either "PQ" or "FG" depending on which mode was used.

# Flowchart of possiblities:
# if PQ completion is possible:
#   returns angle sequence implementing P
# else:
#   if P has nonvanishing imaginary part:
#       Raises ValueError: completion is impossible.
#   else:
#       returns angle sequence implementing P plus some imaginary part
def quantum_signal_processing(P,tol=1e-6,completion_mode=None,allow_scale=True,verify=False):
    assert completion_mode in [None, "PQ", "FG"]

    if not isinstance(P, Cheb):
        P = P.convert(kind=Cheb)


    # check parity
    d = len(P.coef)-1
    for i in range(d+1):
       if i % 2 != d % 2 and not np.allclose(P.coef[i],0):
            raise ValueError("Polynomial for Quantum Signal Processing has mixed parity!")

    # degree-zero case
    if d == 0:
        return P.coef[0], 1, [], "PQ"


    # Note: in general we store a Laurent polynomial as a Poly object
    # Storing Fl(w) = F(w) * w^d. That way f_n = F.coef[d+n].
    # The little l stands for Laurent.

    ####################################
    ############## Attempt PQ completion


    # follow Theorem 4 of 1806.01838
    # compute |P(+-1)|. if not equal, abort PQ completion
    # scale = |P(1)|
    if allow_scale:
        scale = abs(P(1))
        pq_possible = np.allclose(scale, abs(P(-1)), atol=tol)
    else:
        scale = 1
        pq_possible = np.allclose(abs(P(1)), abs(P(-1)), atol=tol)

    if not pq_possible and completion_mode == "PQ":
        raise ValueError("PQ completion failed: |P(1)| = " + str(scale) + " =/= " + str(abs(P(-1))) + " = |P(-1)|")

    if completion_mode == "FG":
        pq_possible = False

    if pq_possible:
        # Pbar = P / scale
        # A(x) = 1 - Pbar(x) Pbar*(x)
        # assert A is even and real
        # At(x^2) = A(x)
        Pbar = P / scale
        A = 1 - Pbar * Cheb(Pbar.coef.conj())
        A = A.convert(kind=Poly)
        atcoef = []
        for i in range(len(A.coef)):
            if verify: assert np.allclose(A.coef[i], np.real(A.coef[i]), atol=tol)
            if i % 2 == 1:
                if verify: assert np.allclose(A.coef[i], 0, atol=tol)
            else:
                atcoef.append(A.coef[i])
        At = Poly(atcoef)

    # check:
    #    At(1) = 0,
    #    if d is even then At(0) = 0
    #   so, except for those roots:
    #    At's real roots have even multiplicity
    #    At's complex roots come in cc pairs
    # if not true, abort PQ completion

    if pq_possible:
        pq_possible = np.allclose(At(1),0,atol=tol)
        if not pq_possible and completion_mode == "PQ":
            raise ValueError("PQ completion failed: At(1) = " + str(At(1)) + " =/= 0")

    if pq_possible and d % 2 == 0:
        pq_possible = np.allclose(At(0),0,atol=tol)
        if not pq_possible and completion_mode == "PQ":
            raise ValueError("PQ completion failed: At(0) = " + str(At(0)) + " =/= 0")

    if pq_possible:
        raw_roots = At.roots()
        paired_roots = []
        unpaired_roots = []
        for r in raw_roots:
            if np.allclose(r,1,atol=tol): continue

            # if d is odd then r=0 is also an extra root
            if d % 2 == 0:
                if np.allclose(r,0,atol=tol): continue

            found = None
            for s in unpaired_roots:
                if np.allclose(r, s.conj(),atol=tol):
                    found = s
                    break
            if found:
                unpaired_roots.remove(found)
                paired_roots.append([found,r])
            else:
                unpaired_roots.append(r)

        pq_possible = (len(unpaired_roots) == 0)
        if not pq_possible and completion_mode == "PQ":
            raise ValueError("PQ completion failed: At has unpaired roots:", unpaired_roots)

    if pq_possible:
        completion_mode = "PQ" # definitely can do this now

        # odd P case:
        #   At(y) = (1 - y) K^2 prod_r (y - r)(y - r*)
        #   W(y) = K prod_r (y - r)
        #   Q(x) = W(x^2)
        # even P case:
        #   At(y) = y (1 - y) K^2 prod_r (y - r)(y - r*)
        #   W(y) = K prod_r (y - r)
        #   Q(x) = x * W(x^2)

        while True:
            test_x = np.random.random()
            if test_x == 0: continue
            if test_x not in raw_roots: break

        if d % 2 == 1:
           K_squared = At(test_x) / polyvalfromroots(test_x, raw_roots)
        else:
           K_squared = At(test_x) / polyvalfromroots(test_x, raw_roots + [0])

        if verify: assert np.allclose(K_squared, np.real(K_squared),atol=tol)
        if verify: assert np.real(K_squared) <= 0
        K = np.sqrt(-K_squared)

        if len(paired_roots) == 0:
            Q = Cheb([K])
        else:
            W = K * Cheb.fromroots([root_pair[0] for root_pair in paired_roots])
            Q = W(Cheb([0.5,0,0.5]))

        if d % 2 == 0:
            Q = Q * Cheb([0,1])


        # Convert PQ to FG via 2105.02859 appendix B 3.
        # say Pbar(sigma) = sum_{n=0}^d pbar_n T_n(sigma)
        # say Q(sigma) = sum_{n=1}^d q_n U_{n-1}(sigma)

        QU = 0 # QU(sigma) = sum_{n=1}^d q_n sigma^n
        for n in range(len(Q.coef)):
            if n == 0:
                # T_0 = U_0 = U_{1-1}
                QU += Poly([0,Q.coef[0]])
            elif n == 1:
                # T_1 = 0.5 * U_1 = 0.5 * U_{2-1}
                QU += 0.5 * Q.coef[n] * Poly([0,0,1])
            else:
                # T_n(x) = 0.5*( U_n(x) - U_{n-2}(x) )
                # = 0.5*( U_{(n+1)-1}(x) - U_{(n-1)-1}(x) )
                QU += 0.5 * Q.coef[n] * Poly([0,1])**(n+1)
                QU -= 0.5 * Q.coef[n] * Poly([0,1])**(n-1)

        pbarcoef = Pbar.coef
        qcoef = QU.coef

        # f_0 = Re(p_0) and g_0 = Im(p_0)
        # for n > 1
        # f_n = Re(pbar_n + q_n)/2
        # f_{-n} = Re(pbar_n - q_n)/2
        # g_n = Im(pbar_n - q_n)/2
        # g_{-n} = Im(pbar_n + q_n)/2

        fcoef = np.zeros(2*d+1)
        gcoef = np.zeros(2*d+1)
        for n in range(-d,d+1):
            if n == 0:
                fcoef[d] = np.real(pbarcoef[0])
                gcoef[d] = np.imag(pbarcoef[0])
            else:
                fcoef[d+n] = np.real(pbarcoef[n] + qcoef[n])/2
                fcoef[d-n] = np.real(pbarcoef[n] - qcoef[n])/2
                gcoef[d+n] = np.imag(pbarcoef[n] - qcoef[n])/2
                gcoef[d-n] = np.imag(pbarcoef[n] + qcoef[n])/2

        Fl = Poly(fcoef)
        Gl = Poly(gcoef)



    ##########################################
    ############# Do FG completion if PQ completion failed

    if not pq_possible:
        completion_mode = "FG"

        # "Capitalization" in 2003.02831 capitalizes F directly
        # via F(w) += tol/4 * w**d  + tol/4 * w**(-d)
        # Capitalization seems to make PQ completion worse, so we don't do it.
        l = np.zeros(d+1)
        l[-1] = 1
        P += tol/2 * Cheb(l)

        # check that P is real. otherwise completion is impossible.
        for n in range(d+1):
            if not np.allclose(P.coef[n], np.real(P.coef[n])):
                raise ValueError("Cannot perform quantum signal processing for complex polynomial P, since there is no Q such that |P|^2 + (1-x^2)|Q|^2 = 1.\nP =", P)
            P.coef[n] = np.real(P.coef[n])

        extrema = (P * Cheb(P.coef.conj())).deriv(1).roots()
        extrema = [x for x in extrema if np.allclose(np.imag(x),0)] # only want real roots
        if allow_scale:
            # scale = maximum |P(x)| for |x| <= 1
            # Pbar = P / scale

            scale = max([ abs(P(x)) for x in list(extrema)+[1,-1] ])
        else:
            for x in list(extrema)+[1,-1]:
                if abs(P(x)) > 1:
                    raise ValueError("Cannot perform quantum signal processing for polynomial P, since it violates |P(x)| < 1 for |x| < 1. Specifically P("+str(x)+") = "+str(P(x)))
            scale = 1

        Pbar = P  / scale
        pbarcoef = Pbar.coef

        # Pbar(x) = sum_{n=0}^d p_n T_n(x)
        # make f_n from pbar_n. The following choice makes |A(w)| = |Pbar(sigma)|
        #       so that |P(sigma)| < 1 implies |A(w)| < 1.
        # f_0 = pbar_0
        # for n > 0 select f_n = f_{-n} = pbar_n/2
        # F(w) = sum_{n=-d}^d f_n w^n
        fcoef = np.zeros(2*d+1)
        for n in range(-d,d+1):
            if n == 0:
                fcoef[d] = np.real(pbarcoef[0])
            else:
                fcoef[d+n] = np.real(pbarcoef[n])/2
                fcoef[d-n] = np.real(pbarcoef[n])/2

        Fl = Poly(fcoef)
        # follow "Completion via root finding" in 2003.02831
        # and "3.2 Complementing Polynomials" in 1806.10236
        # or Lemma 6 in 1806.01838.
        # Honestly, I got none of these methods to work exactly as written.
        # But I found a method that's inspired by Lemma 6 that seems to work.


        # F(w) = F(w) w^(-d)
        # F(w^-1) = Poly(Fl.coef[::-1]) w^(-d)
        # A(w) = 1 - F(w) F(w^-1) = 1 - Fl w^(-d) Poly(Fl.coef[-1]) w^(-d) = 1 - Fl Poly(Fl.coef[::-1]) w^(-2d)
        # AP = A * w^(2d) = w^(2d) - Fl Poly(Fl.coef[::-1])
        l = np.zeros(2*d+1)
        l[-1] = 1
        Al = Poly(l) - Fl * Poly(Fl.coef[::-1])

        # Roots of Al(w) come in quads [r, -r, 1/r, -1/r]

        # Want A(w) = G(w) G(1/w)
        # Select Gl(w) = G(w) w^d = \Pi_r (w^2 - r)
        # Gl(1/w) = G(1/w) w^d = w^2d \Pi_r (w^-2 - r) = \Pi_r (1 - rw^2)
        # G(w)G(1/w) = Gl(w)Gl(1/w) =  \Pi_r (w^2 - r)(1 - rw^2)
        raw_roots = Al.roots()
        root_quads = []
        for r in raw_roots:
            found = False
            for quad in root_quads:
                if len(quad) == 4: continue
                if np.allclose(r, quad[0], atol=tol) \
                        or np.allclose(r, 1/quad[0], atol=tol) \
                        or np.allclose(r, -quad[0], atol=tol) \
                        or np.allclose(r, -1/quad[0], atol=tol):
                    quad.append(r)
                    found = True
                    break
            if not found:
                root_quads.append([r])

        G0l = 1
        for quad in root_quads:
            if verify: assert len(quad) == 4
            G0l *= Poly([-quad[0]**2,0,1])

        # Find a random test point to determine K_squared times all the other constants
        while True:
            test_x = np.random.random()
            if test_x == 0: continue
            if test_x not in raw_roots: break

        K_squared = Al(test_x) / (G0l * Poly(G0l.coef[::-1]))(test_x)
        if verify: assert np.allclose(K_squared, np.real(K_squared),atol=tol)
        if verify: assert np.real(K_squared) >= 0
        K = np.real(np.sqrt(K_squared))

        # prune the terms with incorrect parity, and take the real part
        Gl_coefs = np.zeros(len(G0l.coef))
        for i in range(len(G0l.coef)):
            if i % 2 == 1: continue
            Gl_coefs[i] = np.real(G0l.coef[i])

        Gl = K * Poly(np.real(Gl_coefs))

    #####################################
    ########### Decomposition via halving

    # easier to refer to polynomials in pairs F and Fp (prime)
    # such that together F(w) + i X Fp(w) is unitary
    Fpl = Gl

    # iterate over all indices with the same parity as dim
    def iter(dim):
        n = -dim
        while n <= dim:
            yield n
            n += 2

    # position of -dim <= n <= dim in a 2(dim+1) vector
    def idx(dim,n):
        assert dim % 2 == n % 2
        return (dim+n)//2

    # same but primed
    def idxp(dim,n):
        assert dim % 2 == n % 2
        return (dim+1)+(dim+n)//2

    def laurent_poly_to_array(Al,Apl,d):
        A_arr = np.zeros(2*(d+1))
        for n in iter(d):
            A_arr[idx(d,n)] = np.real(Al.coef[d+n])
            A_arr[idxp(d,n)] = np.real(Apl.coef[d+n])
        return A_arr

    def array_to_laurent_poly(A_arr,d):
        Al = Poly([A_arr[idx(d,dn-d)] if d % 2 == (dn-d) % 2 else 0 for dn in range(2*d+1)])
        Apl = Poly([A_arr[idxp(d,dn-d)] if d % 2 == (dn-d) % 2 else 0  for dn in range(2*d+1)])
        return Al, Apl

    def norm(Al,Apl):
        return Al*Poly(Al.coef[::-1]) + Apl*Poly(Apl.coef[::-1])

    def halve(d,A_arr,l):
        assert len(A_arr) == 2*(d+1)

        # U = A(w) + iX Ap(w)
        # V^dagger = B(w) + iX Bp(w)
        # V^dagger U = C(w) + iX Cp(w)

        # so therefore:
        #  A(w) + iX Ap(w) = U = V (V^dagger U)
        #                 = (B(w) + iX Bp(w))^dagger (C(w) + iX Cp(w))

        # M is 2(d+l+1) by 2(l+1) matrix
        # so that C = M B. Defined by:
        # c_k = sum_(n+m=k) b_m a_n - bp_-m ap_n
        # cp_k = sum_(n+m=k) bp_m a_n + b_-m ap_n
        # obtain these by multiplying out (B(w) + iX Bp(w))*(A(w) + iX Ap(w))
        # notice that 2003.02831 does B(w) + Bp(w) iX whereas we do B(w) + iX Bp(w)
        # so their equations are a bit different.
        M = np.zeros((2*(d+l+1), 2*(l+1)))

        for k in iter(d+l):
            for m in iter(l):
                n = k - m
                if not abs(n) <= d: continue
                M[idx(d+l,k),idx(l,m)] = A_arr[idx(d,n)]
                M[idx(d+l,k),idxp(l,-m)] = -A_arr[idxp(d,n)]
                M[idxp(d+l,k),idxp(l,m)] = A_arr[idx(d,n)]
                M[idxp(d+l,k),idx(l,-m)] = A_arr[idxp(d,n)]

        # Pi is a 2+4*l   by 2(d+l+1) matrix. We will demand Pi C = [1,0,0,0,...].
        # that only contains the c_k and cp_k for d-l < |k| <= d+l (these should vanish)
        # Pi[0,:] = all 1's for c_k for -(d-l) <= k <= d-l   (sum of these must be 1)
        # Pi[1,:] = all 1's for cp_k for -(d-l) <= k <= d-l  (sum of these must be 0)
        Pi = [np.zeros(2*(d+l+1)), np.zeros(2*(d+l+1))]
        for k in iter(d+l):
            if not abs(k) <= (d-l): continue
            Pi[0][idx(d+l,k)] = 1
            Pi[1][idxp(d+l,k)] = 1

        for k in iter(d+l):
            if abs(k) <= d-l: continue
            row = np.zeros(2*(d+l+1))
            row[idx(d+l,k)] = 1
            Pi.append(row)
            row = np.zeros(2*(d+l+1))
            row[idxp(d+l,k)] = 1
            Pi.append(row)
        Pi = np.array(Pi)

        # Require C,Cp are of degree l, and that C(1) = 1 and Cp(1) = 0
        # 2003.02831 actually requires B(1) = 1 and Bp(1) = 0 instead,
        # but doing this to C is more convenient and also makes the solution unique
        # C(1) = 1  -> sum_k C_k = 1
        # Cp(1) = 0 -> sum_k Cp_k = 0

        # solve Pi M B = [1,0,0,0,0,0....] for B
        constraint = np.zeros(2+4*l)
        constraint[0] = 1

        B_arr = np.linalg.lstsq(Pi @ M, constraint, rcond=None)[0]


        # we have that (B(w) + iX Bp(w))^dagger = B(w^-1) - i Bp(w^-1) X = B(w^-1) - i X Bp(w^)
        # so to convert B to the inverse Bd, reverse the order of Bl and negate Bpl
        Bd_arr = np.zeros(2*(l+1))
        for m in iter(l):
            Bd_arr[idx(l,m)] = B_arr[idx(l,-m)]
            Bd_arr[idxp(l,m)] = -B_arr[idxp(l,m)]

        # this C has dimension  2*(d+l+1), but actually
        # k such that d-l <= |k| < d-l vanish.
        C_arr_tmp = M @ B_arr
        C_arr = np.zeros(2*(d-l+1))
        for k in iter(d-l):
            C_arr[idx(d-l,k)] = C_arr_tmp[idx(d+l,k)]
            C_arr[idxp(d-l,k)] = C_arr_tmp[idxp(d+l,k)]

        # Test if we succeeded.
        if verify:
            # Make the polynomials
            Al, Apl = array_to_laurent_poly(A_arr,d)
            Bl, Bpl = array_to_laurent_poly(B_arr,l)
            Bdl, Bdpl = array_to_laurent_poly(Bd_arr,l)
            Cl, Cpl = array_to_laurent_poly(C_arr,d-l)

            # check if Bl,Bpl is unitary:
            Il = Bdl*Bl - Poly(Bdpl.coef[::-1])*Bpl
            Ipl = Bdpl*Bl + Poly(Bdl.coef[::-1])*Bpl
            Il.coef[2*l] -= 1
            assert np.allclose(Il.coef,0,atol=tol)
            assert np.allclose(Ipl.coef,0,atol=tol)

            # check if Cl,Cpl is unitary:
            Cd_arr = np.zeros(2*((d-l)+1))
            for m in iter(d-l):
                Cd_arr[idx(d-l,m)] = C_arr[idx(d-l,-m)]
                Cd_arr[idxp(d-l,m)] = -C_arr[idxp(d-l,m)]
            Cdl, Cdpl = array_to_laurent_poly(Cd_arr,d-l)
            Il = Cdl*Cl - Poly(Cdpl.coef[::-1])*Cpl
            Ipl = Cdpl*Cl + Poly(Cdl.coef[::-1])*Cpl
            Il.coef[2*(d-l)] -= 1
            assert np.allclose(Il.coef,0,atol=tol)
            assert np.allclose(Ipl.coef,0,atol=tol)

            # check Bd,Bdpl, and Cl,Cpl multiply to give Al,Apl
            Dl = Bdl*Cl - Poly(Bdpl.coef[::-1])*Cpl
            Dpl = Bdpl*Cl + Poly(Bdl.coef[::-1])*Cpl
            D_arr = laurent_poly_to_array(Dl,Dpl,d)
            assert np.allclose(A_arr, D_arr, atol=tol)

        return Bd_arr, C_arr

    # recursively decompose A(w) + i X Ap(w) into a product of
    # degree one polynomials
    def decompose(d,A_arr):
        assert len(A_arr) == 2*(d+1)

        if d == 1:
            return [A_arr]
        else:
            l = d // 2
            l = int(np.ceil(d/2))
            B_arr,C_arr = halve(d,A_arr,l)

            assert len(B_arr) == 2*(l+1)
            assert len(C_arr) == 2*((d-l)+1)
            return decompose(l,B_arr) + decompose(d-l,C_arr)

    ###################o

    # check that F, Fp have the correct parity
    for n in range(-d,d+1):
        if n % 2 != d % 2:
            assert np.allclose(Fl.coef[d+n],0)
            assert np.allclose(Fpl.coef[d+n],0)

    # pack FP and GP into single F array, keeping only terms of the right parity
    F_arr = laurent_poly_to_array(Fl,Fpl,d)

    deg1_polys = decompose(d,F_arr)


    # for each A in decompose(d, F):
    # A(w) + iX Ap(w) = e^(i phi X) tildew e^(i varphi X)
    # such that:
    # cos(phi + varphi) = a_1 + a_-1
    # cos(phi - varphi) = a_1 - a_-1
    # sin(phi + varphi) = ap_1 + ap_-1
    # sin(phi - varphi) = ap_1 - ap_-1


    # to solve, set phi+varphi = +- arccos(a_1 + a_-1)
    # and use the sign of ap_1 + ap_-1 to figure out the sign
    # assert that everything is consistent. Similarly for phi-varphi
    phi_pairs = []
    for A_arr in deg1_polys:
        assert len(A_arr) == 4
        a, am, ap, apm = A_arr

        if verify: # check if unitary
            Al, Apl = array_to_laurent_poly(A_arr,1)
            Adl, Adpl = array_to_laurent_poly([am,a,-ap,-apm],1)
            Il = Adl*Al - Poly(Adpl.coef[::-1])*Apl
            Ipl = Adpl*Al + Poly(Adl.coef[::-1])*Apl
            Il.coef[2] -= 1
            assert np.allclose(Il.coef,0,atol=tol)
            assert np.allclose(Ipl.coef,0,atol=tol)

        # s = "sum" = phi+varphi
        # d = "diff" = phi-varphi
        a_sum = a+am
        if a_sum < -1: a_sum = -1
        if a_sum > 1: a_sum = 1
        a_diff = a-am
        if a_diff < -1: a_diff = -1
        if a_diff > 1: a_diff = 1

        su = np.arccos(a_sum)
        di = np.arccos(a_diff)
        if ap+apm < 0: su = -su
        if ap-apm < 0: di = -di

        # These asserts seem to not succeed reliably, even if A is unitary.
        # assert np.allclose(np.sin(su), ap+apm, atol=tol)
        # assert np.allclose(np.sin(di), ap-apm, atol=tol)

        phi, varphi = (su+di)/2, (su-di)/2
        phi_pairs.append((phi,varphi))


    # H (A(w) + iX Ap(w)) H = e^(i phi Z) e^(i theta X) e^(i varphi Z)
    # convert to reflection convention
    # e^(i X theta) = diag(1, i) refl(sigma) diag(1, i)
    # diag(1, i) = e^(i pi/4) * e^(- i Z pi/4)
    refl_phi_pairs = []
    refl_phase = 1
    for (phi,varphi) in phi_pairs:
        phi -= np.pi/4
        if phi < 0: phi += 2*np.pi
        varphi -= np.pi/4
        if varphi < 0: varphi += 2*np.pi
        refl_phi_pairs.append((phi,varphi))
        refl_phase *= 1j

    # the outer phis contribute to a global phase
    phase = np.exp(1j*(refl_phi_pairs[0][0] + refl_phi_pairs[-1][1])) * refl_phase

    # add adjacent varphi,phi
    phis = []
    for (phi,varphi) in refl_phi_pairs:
        if len(phis) > 0: phis[-1] += phi
        if len(phis) < len(phi_pairs)-1: phis.append(varphi)

    # reverse to match matrix multiplication order
    phis = list(reversed(phis))

    return scale, phase, phis, completion_mode


def cheby_angles(d):
    global_phi = 1j**(1-d)
    phis = [np.pi/2]*(d-1)
    return global_phi, phis


# computes P(sigma) where P(sigma) =  <+| (F(w) + i X Fp(w)) |+>
def angles_to_poly(phis):
    zero = Poly([0,0,1])/np.sqrt(2)
    one = Poly([1,0,0])/np.sqrt(2)

    for i in range(len(phis)):
        # diag(1, -i) refl(sigma) diag(1, -i) =  e^(i X theta)

        zero, one = (zero+one)/np.sqrt(2), (zero-one)/np.sqrt(2)
        one *= -1j
        zero *= np.exp(1j*(phis[i]))
        one *= np.exp(-1j*(phis[i]))
        one *= -1j
        zero, one = (zero+one)/np.sqrt(2), (zero-one)/np.sqrt(2)

        # implements e^(i theta Z)
        zero *= Poly([0,0,1])
        one *= Poly([1,0,0])

    F_plus_iFp = (zero + one)/np.sqrt(2)

    d = len(phis)+1
    assert len(F_plus_iFp.coef) % 2 == 1
    dprime = int((len(F_plus_iFp.coef) - 1)/2)
    assert dprime <= d

    P_coefs = np.zeros(d+1).astype(complex)
    for i in range(-dprime,dprime+1,1):
        P_coefs[abs(i)] += F_plus_iFp.coef[dprime+i]
    return Cheb(P_coefs)


# evaluates the matrix H (F(w) + i X Fp(w)) H which takes the form:
# [   P(sigma)                   i Q*(sigma) sqrt(1-sigma^2)  ]
# [ i Q(sigma) sqrt(1-sigma^2)              P*(sigma)         ]
def angles_to_matrix(phis,sigma):
    signal = np.array([[sigma,np.sqrt(1-sigma**2)],[np.sqrt(1-sigma**2), -sigma]])

    out = signal

    for i in range(len(phis)):
        out = signal @ np.diag([np.exp(1j*phis[i]),np.exp(-1j*phis[i])]) @  out

    return out


########################################33


def test_qsp(d,PQ=False):
    assert d == int(d)
    d = int(d)
    assert d >= 0

    if PQ: # pq completion mode

        if d % 2 == 1: # p is odd
            # random complex odd pbar such that:
            # |pbar(x)| <= 1 for real |x| <= 1
            # |pbar(x)| >= 1 for real |x| >= 1

            # keep trying until we find something good
            while True:
                coef = np.random.random(d+1) + 1j*np.random.random(d+1)
                for i in range(d+1):
                    if i % 2 != 1: coef[i] = 0
                Pbar = Poly(coef)
                Pbar = Pbar / abs(Pbar(1))

                extrema = Pbar.deriv(1).roots()
                all_good = True
                for x in extrema:
                    if not np.allclose(np.imag(x),0): continue
                    if abs(x) > 1:
                        # |pbar(x)| >= 1 for real |x| >= 1
                        if abs(Pbar(x)) < 1:
                            all_good = False
                            break
                    else:
                        # |pbar(x)| <= 1 for real |x| <= 1
                        if abs(Pbar(x)) > 1:
                            all_good = False
                            break

                if all_good: break

            P = Pbar * np.random.random()*5


        else: # p is even
            # random complex odd q such that
            # (1 - x^2) |q(x)|^2 <= 1 for real |x| <= 1
            # this is basically going through Theorem 4 of 1806.01838 in reverse

            coef = np.random.random(d+1) + 1j*np.random.random(d+1)
            for i in range(d+1):
                if i % 2 != 1: coef[i] = 0
            Q = Poly(coef)
            Test = Poly([1,0,-1]) * Q * Poly(coef.conj())

            extrema = Test.deriv(1).roots()
            scale = 1
            for x in extrema:
                val = abs(Test(x))
                if Test(x) > 1:
                    if val > scale: scale = val
            Qbar = Q / np.sqrt(scale)

            # A(x) = 1 - (1 - x^2) Qbar(x) Qbar*(x)
            # A is real and even
            # At(x^2) = A(x)
            A = 1 - Poly([1,0,-1]) * Qbar * Poly(Qbar.coef.conj())
            atcoef = []
            for i in range(len(A.coef)):
                assert np.allclose(A.coef[i], np.real(A.coef[i]))
                if i % 2 == 1:
                    assert np.allclose(A.coef[i], 0)
                else:
                    atcoef.append(A.coef[i])
            At = Poly(atcoef)

            # At(y) >= 0 for all real y. Therefore at's real roots have even multiplicity
            # At coeffs real. Therefore At's complex roots come in conjugate pairs.
            # At(y) = K^2 prod_r (y - r) (y - r*)
            raw_roots = At.roots()

            paired_roots = []
            unpaired_roots = []
            for r in raw_roots:
                found = None
                for s in unpaired_roots:
                    if np.allclose(r, s.conj()):
                        found = s
                        break
                if found:
                    unpaired_roots.remove(found)
                    paired_roots.append([found,r])
                else:
                    unpaired_roots.append(r)

            assert len(unpaired_roots) == 0

            while True:
                test_x = np.random.random()
                if test_x not in raw_roots: break

            K_squared = At(test_x) / polyvalfromroots(test_x, raw_roots)
            assert np.allclose(K_squared, np.real(K_squared))
            assert np.real(K_squared) > 0
            K = np.sqrt(K_squared)

            # let W(y) = K prod_r (y -r)
            # let Pbar(x) = W(x^2)
            # P = Pbar * random number
            W = K * Poly.fromroots([root_pair[0] for root_pair in paired_roots])
            Pbar = W(Poly([0,0,1]))

            P = Pbar * np.random.random()*5

    else: # fg completion mode

        coef = np.random.random(d+1)
        for i in range(d+1):
            if i % 2 != d % 2: coef[i] = 0

        P = Poly(coef)

    if PQ: mode = "PQ"
    else: mode = "FG"
    # mode = None

    scale, phase, phis, PQ = quantum_signal_processing(P,completion_mode=mode,tol=1e-3,verify=False)

    import matplotlib.pyplot as plt
    xs = np.linspace(-1,1,100)
    out_real = []
    out_imag = []
    for x in xs:
        y = angles_to_matrix(phis,x)[0,0] * scale * phase
        out_real.append(np.real(y))
        out_imag.append(np.imag(y))
    plt.plot(xs,out_real,c="b")
    plt.plot(xs,out_imag,c="b",ls="--")

    plt.plot(xs,0.01+np.real(P(xs)),c="r")
    plt.plot(xs,0.01+np.imag(P(xs)),c="r",ls="--")
    plt.grid()

    plt.show()


if __name__ == "__main__":
    print(quantum_signal_processing(Poly([0,0,1]),completion_mode="FG",tol=1e-3,verify=False))
    # test_qsp(13,PQ=False)


