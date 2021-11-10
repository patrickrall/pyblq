from pyblq import *
from pyblq.qsp import quantum_signal_processing, angles_to_poly, angles_to_matrix

# scalar X
# row/column vector X
# unitary X
# isometry X
# square matrix X
# non-square matrix X

#th = 0.16
#s, c = np.sin(th), np.cos(th)
#print(np.array([[c,-s],[s,c],[0,0]]))
#plus = Blq(np.array([[c,-s],[s,c],[0,0]]))
#print(plus)
#print(plus.kraus())



from numpy.polynomial import Chebyshev as Cheb
from numpy.polynomial import Polynomial as Poly
import matplotlib.pyplot as plt

#P = Cheb([0,1,0,1,0,0.3])
#P = Cheb([0.5,0,0.5])

b = Blq([[1],[0.2]])
b /= abs(b)

b = b.svt(Cheb([0,0,1]))
b = b.svt(Cheb([0,0,1]))

print(b)
print(b.sample())

assert False


def repeated_squaring(base_dim, exp_bits):

    for reps in range(exp_bits):
        rep_term = Blq("""[base:n @ power:n <- base:n]
            power += base
            repeat reps:
                power = power*power
        """,n=base_dim,reps=reps)




b = repeated_squaring(5,2)
#print(b)
#print(b.kraus())
print(b)
print(b.kraus())


xs = np.linspace(-1,1,100)

if False:
    scale, phase, phis, kind = quantum_signal_processing(P,completion_mode="PQ")
    print(kind)

    P = angles_to_poly(phis) * phase
    plt.plot(xs, [np.real(P(x)) for x in xs], label="poly-re")
    plt.plot(xs, [np.imag(P(x))+0.01 for x in xs], label="poly-im")
    plt.plot(xs, [np.real(angles_to_matrix(phis,x)[0,0] * phase)+0.02 for x in xs], label="mat-re")
    plt.plot(xs, [np.imag(angles_to_matrix(phis,x)[0,0] * phase)+0.03 for x in xs], label="mat-im")

    plt.legend()
    plt.show()


if False:
    out = []
    for x in xs:
        if True:
            b = Blq(np.array([[x,np.sqrt(1-x**2)],[np.sqrt(1-x**2),-x]]))
            b = Blq([[1,0]])  * b * Blq([[1],[0]])
            b = b.svt(P)
            out.append(b.kraus()[0][0,0])
        if False:
            b = Blq(np.array([[x,np.sqrt(1-x**2)],[np.sqrt(1-x**2),-x]]))
            b = Blq("""[x:2 <- x:2]
            declare y : 2
            if (x == 1):
                y <- b * |y
            scalar (~ket(0)) * |y
            """,b=b)
            b = b.svt(P)
            out.append(b.kraus()[0][1,1])
        if False:
            b = Blq(np.diag([1,x])).svt(P)
            out.append(b.kraus()[0][1,1])

    plt.plot(xs, np.real(out), label="svt-re")
    plt.plot(xs, np.imag(out), label="svt-im")
    plt.plot(xs, [np.real(P(x))+0.01 for x in xs], label="poly-re")
    plt.plot(xs, [np.imag(P(x))+0.01 for x in xs], label="poly-im")
    #plt.plot(xs, [np.real(angles_to_matrix(phis,x)[0,0] * phase*scale)+0.02 for x in xs], label="mat-re")
    #plt.plot(xs, [np.imag(angles_to_matrix(phis,x)[0,0] * phase*scale)+0.03 for x in xs], label="mat-im")

    plt.legend()
    plt.grid()
    plt.show()

if False:
    x = 0.5
    if True:
        b = Blq(np.array([[x,np.sqrt(1-x**2)],[np.sqrt(1-x**2),-x]]))
        b = Blq([[1,0]])  * b * Blq([[1],[0]])
        b = b.svt(P)
    if False:
        b = Blq(np.array([[x,np.sqrt(1-x**2)],[np.sqrt(1-x**2),-x]]))
        b = Blq("""[x:2 <- x:2]
        declare y : 2
        if (x == 1):
            y <- b * |y
        scalar (~ket(0)) * |y
        """,b=b)

        b = b.svt(P)
    if False:
        b = Blq(np.diag([1,x]))
        # b = Blq([[0,1]])  * b * Blq([[0],[1]])
        b = b.svt(P)

    print(b)
    print(b.kraus())
    print(x, P(x))




#swap1 = Blq("""[x:2 @ y:2 <- x:2 @ y:2  ] x @ y <- |y @ |x """)
#I = Blq("[z:2 <- z:2] pass")
#print((swap1 @ I) * (I @ swap1) * (swap1 @ I))


# test =  Blq([[1,1],[0,1j]])
#test = my_prog
#print(test)
#print(test.kraus())


#swap2 = Blq("""[x:2 @ y:2 @ z:2 <- x:2 @ y:2 @ z:2 ] pass
#y @ z <- |z @ |y
#""")

# todo: multiply
# test array feature of init

#############
# Z = Blq([[1,0],[0,-1]])
# X = Blq([[0,1],[1,0]])
# Z = ket(0) @ ~ket(0) - ket(1) @ ~ket(1)
# X = ket(1) @ ~ket(0) - ket(0) @ ~ket(1)
# I = Blq("[x:2 -> x:2]")

# H = 0.2 * Z @ I + 0.3 * I @ Z + 0.5 * X @ X

