
from .tokens import tokenize
from .make_ast import make_ast
from .scoping import scope_analysis
from .runtime import build_block
from .qac import *
from .qac_string import qac_to_string
from .qac_sim import *
from .expose import expose
from .svt import svt


from .debug import recur_print

from numpy.polynomial import Chebyshev as Cheb

import numpy as np

class Blq():
    def __init__(self,*args,**kwargs):
        self.qac = []
        self.sim_cache = None

        # empty block
        if len(args) == 0:
            if kwargs:
                raise ValueError("Can't take any kwargs for an empty block.")
            self.qac.append(QAC())
            return

        if len(args) == 1:
            if isinstance(args[0],Blq):
                if kwargs:
                    raise ValueError("Can't take any kwargs when initializing via a Blq object.")
                for qac in args[0].qac:
                    self.qac.append(block_copy(qac))
                return

            if isinstance(args[0],list)\
                    or isinstance(args[0],np.ndarray)\
                    or isinstance(args[0],np.matrix):

                if kwargs:
                    raise ValueError("Can't take any kwargs when initializing via a matrix.")

                mat = np.array(args[0]).astype(complex)

                if len(mat.shape) != 2:
                    raise ValueError("Initializing via ndarray requires a matrix, but input has shape "+str(mat.shape)+".")

                self.qac = [matrix_to_qac(mat)]

                return

        # a list of strings, blocks, and complex-able types
        # pass a single block for the copy constructor.
        args = list(args)
        for i in range(len(args)):
            if not isinstance(args[i], str) and not isinstance(args[i], Blq):
                try: args[i] = complex(args[i])
                except: raise ValueError("Argument "+str(i+1)+" must either be a string, a block, or cast to a complex.")

        for key in kwargs.keys():
            if not isinstance(kwargs[key], Blq):
                try: kwargs[key] = complex(kwargs[key])
                except: raise ValueError("Key word argument "+key+" must either be a block or cast to a complex.")

        tokens = tokenize(args)
        in_decls,out_decls,instrs = make_ast(tokens, args)

        scope_analysis(in_decls,out_decls,instrs,args,kwargs)

        qac = build_block(in_decls,out_decls,instrs,args,kwargs)

        self.qac.append(qac)

    def __repr__(self):
        if len(self.qac) > 1: qac = block_add(*[block_copy(qac) for qac in self.qac])
        else: qac = self.qac[0]

        s = "|" + ",".join([str(x.dim) for x in qac.unnamed_outputs])
        s += "> "+str(qac.scale)+ " <"
        s += ",".join([str(x.dim) for x in qac.unnamed_inputs]) + "|"
        return s

    def __str__(self):
        if len(self.qac) > 1: qac = block_add(*[block_copy(qac) for qac in self.qac])
        else: qac = self.qac[0]
        return qac_to_string(qac)

    def __abs__(self):
        return sum([q.scale for q in self.qac])

    ####### Simulation

    def sample(self):
        if self.sim_cache is None:
            if len(self.qac) > 1: qac = block_add(*[block_copy(qac) for qac in self.qac])
            else: qac = block_copy(self.qac[0])

            if len(qac.unnamed_inputs) > 0:
                raise ValueError("Can only sample the output of an isometry, but this block has shape "+qac.shape())

            if not np.allclose(qac.scale, 1):
                raise ValueError("Can only sample from an unscaled block, but this block has scale "+str(qac.scale))

            self.sim_cache = get_qac_state(qac)

        state, registers, discards = self.sim_cache
        return sample_state(state, registers, discards)

    def kraus(self):
        if len(self.qac) > 1: qac = block_add(*[block_copy(qac) for qac in self.qac])
        else: qac = block_copy(self.qac[0])

        if len(qac.unnamed_inputs) not in [0,1] or len(qac.unnamed_outputs) not in [0,1]:
            raise ValueError("Can only compute the Kraus form of a block with at most one input register and at most one output register, but this block has shape "+qac.shape())

        return get_kraus(qac)


    ####### Expose and SVT

    def expose_declare(self):
        if len(self.qac) > 1: qac = block_add(*[block_copy(qac) for qac in self.qac])
        else: qac = block_copy(self.qac[0])

        qac, regs = expose(qac, which=["declare"])
        out = Blq()
        out.qac = [qac]
        return out, [reg.dim for reg in regs["declare"]]

    def expose_discard(self):
        if len(self.qac) > 1: qac = block_add(*[block_copy(qac) for qac in self.qac])
        else: qac = block_copy(self.qac[0])

        qac, regs = expose(qac, which=["discard"])
        out = Blq()
        out.qac = [qac]
        return out, [reg.dim for reg in regs["discard"]]

    def expose_maxmixed(self):
        if len(self.qac) > 1: qac = block_add(*[block_copy(qac) for qac in self.qac])
        else: qac = block_copy(self.qac[0])

        qac, regs = expose(qac, which=["maxmixed"])
        out = Blq()
        out.qac = [qac]
        return out, [reg.dim for reg in regs["maxmixed"]]

    def expose_zero(self):
        if len(self.qac) > 1: qac = block_add(*[block_copy(qac) for qac in self.qac])
        else: qac = block_copy(self.qac[0])

        qac, regs = expose(qac, which=["zero"])
        out = Blq()
        out.qac = [qac]
        return out, [reg.dim for reg in regs["zero"]]

    def svt(self, P,allow_scale=True):

        # convert to chebyshev polynomial
        if not isinstance(P, Cheb):
            try:
                P = P.convert(kind=Cheb)
            except:
                raise ValueError("Input "+repr(P)+" is not a numpy polynomial object.")

        # check fixed parity
        first_nonzero = None
        for i in range(len(P.coef)):
            if not np.allclose(P.coef[i],0):
                first_nonzero = i
                break

        for i in range(len(P.coef)):
            if i % 2 != first_nonzero % 2 and not np.allclose(P.coef[i],0):
                raise ValueError("Input polynomial has mixed parity: has terms of degree "+str(i)+" and of degree "+str(first_nonzero)+".")

        if len(self.qac) > 1: qac = block_add(*[block_copy(qac) for qac in self.qac])
        else: qac = block_copy(self.qac[0])

        qac = svt(qac,P,allow_scale=allow_scale)
        out = Blq()
        out.qac = [qac]
        return out

    ####### addition -> linear combination

    def __add__(self,other):
        if not isinstance(other,Blq):
            raise ValueError("Can only add blocks, but "+repr(other)+" is not a block.")

        shape1 = self.qac[0].shape()
        shape2 = other.qac[0].shape()
        if shape1 != shape2:
            raise ValueError("Block shapes '"+shape1+"' and '"+shape2+"' do not match.")

        def has_discard(instrs):
            for instr in instrs:
                if instr["kind"] == "qac_discard": return True
                if instr["kind"] == "qac_maxmixed": return True
                if instr["kind"] == "qac_if":
                    if has_discard(instr["instructions"]): return True
            return False

        for block in [self,other]:
            for qac in self.qac:
                if has_discard(qac.instrs):
                    raise ValueError("Block "+repr(block)+" has incoherent instructions (discard/maxmixed) and cannot be used in an LCU circuit. To bypass this, use block.expose_discard and block.expose_maxmixed.")

        out = Blq()
        out.qac = [block_copy(qac) for qac in self.qac]
        out.qac += [block_copy(qac) for qac in other.qac]
        return out

    def __sub__(self,other):
        return self - other

    def __rsub__(self,other):
        return other - self

    def __neg__(self):
        return self * -1

    ####### multiplication

    def __mul__(self,other):
        if not isinstance(other,Blq):
            try: other = complex(other)
            except: raise ValueError("Block product must be either with another block or a scalar.")
            out = Blq(self)
            out.qac = [block_copy(out.qac[i]) for i in range(len(out.qac))]
            out.qac = [block_tensor(out.qac[i], out.qac[i].block_scalar(other))
                                for i in range(len(out.qac))]
            return out

        rshape = self.qac[0].rshape()
        lshape = other.qac[0].lshape()
        if rshape != lshape:
            raise ValueError("Right shape '"+rshape+"' does not match left shape '"+lshape+"'.")

        if len(self.qac) > 1: qac1 = block_add(*[block_copy(qac) for qac in self.qac])
        else: qac1 = block_copy(self.qac[0])
        if len(other.qac) > 1: qac2 = block_add(*[block_copy(qac) for qac in other.qac])
        else: qac2 = block_copy(other.qac[0])

        out = Blq()
        out.qac = [block_mul(qac1,qac2)]
        return out

    def __rmul__(self,other):
        # we know that other is not a block, otherwise it would have called __mul__
        try: other = complex(other)
        except: raise ValueError("Block product must be either with another block or a scalar.")

        out = Blq(self)
        out.qac = [block_copy(out.qac[i]) for i in range(len(out.qac))]
        out.qac = [block_tensor(out.qac[i], out.qac[i].block_scalar(other))
                            for i in range(len(out.qac))]
        return out

    ####### division

    def __truediv__(self,other):
        try: other = complex(other)
        except: raise ValueError("Can only divide blocks by a scalar.")

        out = Blq(self)
        out.qac = [block_copy(out.qac[i]) for i in range(len(out.qac))]
        out.qac = [block_tensor(out.qac[i], out.qac[i].block_scalar(1/other))
                            for i in range(len(out.qac))]
        return out

    ####### matmul -> tensor products and scalar multiplication

    def __matmul__(self,other):
        if isinstance(other,Blq):
            if len(self.qac) > 1: qac1 = block_add(*[block_copy(qac) for qac in self.qac])
            else: qac1 = block_copy(self.qac[0])
            if len(other.qac) > 1: qac2 = block_add(*[block_copy(qac) for qac in other.qac])
            else: qac2 = block_copy(other.qac[0])

            out = Blq()
            out.qac = [block_tensor(qac1,qac2)]
            return out

        try: other = complex(other)
        except: raise ValueError("Block tensor product must be either with another block or a scalar.")
        out = Blq(self)
        out.qac = [block_copy(out.qac[i]) for i in range(len(out.qac))]
        out.qac = [block_tensor(out.qac[i], out.qac[i].block_scalar(other))
                            for i in range(len(out.qac))]
        return out

    def __rmatmul__(self,other):
        # we know that other is not a block, otherwise it would have called __matmul__
        try: other = complex(other)
        except: raise ValueError("Block tensor product must be either with another block or a scalar.")

        out = Blq(self)
        out.qac = [block_copy(out.qac[i]) for i in range(len(out.qac))]
        out.qac = [block_tensor(out.qac[i], out.qac[i].block_scalar(other))
                            for i in range(len(out.qac))]
        return out

    def tensorpower(self,n):
        assert int(n) == n
        n = int(n)
        assert n >= 0

        if n == 0: return Blq("[ <- ] pass")

        if len(self.qac) == 1: qac = block_copy(self.qac[0])
        else: qac = block_add(*[block_copy(qac) for qac in self.qac])

        out = block_copy(qac)
        for i in range(n-1):
            out = block_tensor(out, block_copy(qac))

        out_blq = Blq(self)
        out_blq.qac = [out]
        return out_blq

    ######### conjugate transpose
    def __invert__(self):
        out = Blq()
        out.qac = [block_adjoint(block_copy(self.qac[i])) for i in range(len(self.qac))]
        return out

####################

def ket(value,d=2):
    assert int(d) == d
    assert d > 1
    assert value == int(value)
    assert value >= 0
    assert value < d
    return Blq("[x:",d," <- ] x += ",value)


def discard(d=2):
    return Blq("[ <- x:",d,"] discard x ")

