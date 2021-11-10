

global register_count
register_count = 0
global register_refs
register_refs = {}
global register_dims
register_dims = {}

def Register():
    def __init__(self, dim):
        assert int(dim) == dim
        assert dim > 1

        register_count += 1
        self.id = register_count
        register_refs[self.id] = None
        register_dims[self.id] = dim

    def trace(self):
        out = self.id
        while register_refs[out] != None:
            out = register_refs[out]
        return out

    @property
    def dim(self):
        global register_dims
        return register_dims[self.trace()]

    # TODO: Do I ever need this?
    # @dim.setter
    # def dim(self,dim):
    #     register_dims[self.trace()] = dim

    def __eq__(self,other):
        if not isinstance(other,Register): return False
        return self.trace() == other.trace()

    def substitute(self,other):
        assert isinstance(other,Register)
        assert self.dim == other.dim
        target = other.trace()
        if target == self.trace(): return

        global register_refs
        key = self.id
        while register_refs[key] != None:
            key = register_refs[key]
        register_refs[key] = target


########################################


class QAASM():

    def __init__(self):
        self.scale = 1

        self.inputs = {} # "x":<reg>
        self.outputs = {}  # "y":<reg>
        self.rshape = [] # ["x"]
        self.lshape = [] # ["y"]

    def declare(self,):




##############  creation

def block_from_reg(b, reg):
    pass

def block_from_scalar(b, scalar):
    pass

def block_from_consume(reg):
    pass

def block_from_create(reg):
    pass


############ sequencing

def qaasm_declare(b, dim):
    reg = Register(dim)
    b.instrs.append({"kind":"qaasm_declare", "reg":reg, "dim":dim})
    b.outputs.append(b)
    return reg

def qaasm_increment(b, reg, expn, deps):
    assert reg in b.outputs
    for dep in deps:
        if dep in b.outputs: continue
        if dep not in b.references:
            b.references.append(dep)

    b.instrs.append({"kind":"qaasm_increment", "reg":reg, "expn": expn})

def qaasm_discard(b,reg):
    assert reg in b.outputs
    b.instrs.append({"kind":"qaasm_discard", "reg":reg})
    b.outputs.remove(reg)

def qaasm_zero(b,reg):
    assert reg in b.outputs
    b.instrs.append({"kind":"qaasm_zero", "reg":reg})
    b.outputs.remove(reg)

def qaasm_unitary(b,reg, matrix):
    assert reg in b.outputs
    b.instrs.append({"kind":"qaasm_unitary", "reg":reg, "mat":matrix})

def qaasm_phase(b,ph):
    assert abs(ph) == 1
    b.instrs.append({"kind":"qaasm_unitary", "phase":ph})


# __enter__ copies the correct register into the target, __exit__ uncomputes
# {"target":<reg>, "key":<reg>, "regs":[<reg>,<reg>]}
class qaasm_array_idx():
    def __init__(self,b,dat):
        self.dat = dat
        assert dat["key"] in b.outputs
        for dep in dat["regs"]:
            if dep in b.outputs: continue
            if dep not in b.references:
                b.references.append(dep)

    def flip_flag(self,i):
        b.instrs.append({
            "kind":"qaasm_increment", "reg":self.flag,
            "expn": {"kind": "boolean_expn",
                "terms":[{"kind":"value_expn", "value":i},
                    "==",
                    {"kind":"register_expn", "reg":self.dat["key"]}]
                }
        })

    def copy_target(self,i,flip=False):
        def opt_negate(x):
            if not flip: return x
            return {"kind": "negate_expn", "expn": x}

        b.instrs.append({"kind":"qaasm_if", "cond":flag, "instructions":[
            {
                "kind":"qaasm_increment", "ref": self.dat["target"], "expn":opt_negate({
                    "kind":"register_expn": "reg":self.dat["regs"][i]
                })
            }
            ]})

    def __enter__(self):
        b.instrs.append({"kind":"qaasm_declare", "reg":self.dat["target"], "dim":self.dat["target"].dim})
        b.outputs.append(b)
        self.flag = qaasm_declare(b,2)

        for i in range(len(self.dat["regs"])):
            flipEq(i)
            copy_target(i)
            flipEq(i)


    def __exit__(self):
        for i in range(len(self.dat["regs"])):
            flipEq(i)
            copy_target(i,flip=True)
            flipEq(i)

        qaasm_zero(b,self.dat["target"])
        qaasm_zero(b,self.flag)


def qaasm_postselect(b,p):
    flag = qaasm_declare(b,2)
    alpha = p**(1/2)
    beta = (1-p)**(1/2)
    qaasm_unitary(b, [[a,-b],[b, a]])
    qaasm_zero(b,flag)


class qaasm_if():

    def __init__(b, reg):
        self.b = b
        self.reg = reg
        pass

    def __enter__(self):
        self.tmpblq = Blq()

        for reg in self.b.outputs:
            self.tmpblq.outputs.append(reg)

        # don't have enough info on how to repair:.

        return self.tmpblq

    def __exit__(self):
        pass

def uncompute(b,reg):
    raise NotImplementedError
    pass


##############  arithmetic
def block_add(b1,b2):
    pass

def block_mul(b1,b2):
    pass

def block_tensor(b1,b2):
    pass

def block_adjoint(b1,b2):
    pass



