from .qac_string import qac_to_string
import numpy as np


global register_count
register_count = 0
global register_refs
register_refs = {}
global register_dims
register_dims = {}

# a data structure such that:
#    identity unique upon initialization
#    can be merged with other registers
#    stores dimension
#    name hint: a list of names this register has
#       to make auto-generated qac more readable

class Register():
    def __init__(self, dim):
        assert int(dim) == dim
        assert dim > 1

        global register_count
        global register_refs
        global register_dims
        register_count += 1
        self.id = register_count
        register_refs[self.id] = None
        register_dims[self.id] = dim

        self.name_hints = []

    def trace(self):
        out = self.id
        global register_refs
        while register_refs[out] != None:
            out = register_refs[out]
        return out

    def __hash__(self):
        return self.trace()

    @property
    def dim(self):
        global register_dims
        return register_dims[self.trace()]

    def __eq__(self,other):
        if not isinstance(other,Register): return False
        return self.trace() == other.trace()

    def substitute(self,other):
        assert isinstance(other,Register)
        assert self.dim == other.dim
        target = other.trace()
        if target == self.trace(): return

        key = self.id
        global register_refs
        while register_refs[key] != None:
            key = register_refs[key]
        register_refs[key] = target

###################################

# Idea: registers in instructions are inherently named. registers in algebra are inherently unnamed
# consume, create, cast are sort-of the boundary between named and unnamed registers

# two types of QAC registers: scoped / unnamed
# types of registers is determined by their presence in the bookkeeping dictionaries,
# self.scope_inputs/outputs and self.unnamed_inputs/outputs
# not by the register objects themselves, or the Qaasm instructions.

# blocks as inputs to if statements can't have any unnamed registers, must encode 1x1 matrices.
# algebraic blocks can't have any scoped registers as output (but can as input, that's what consumption is)
# user-level blocks can't have any scoped registers and referenced registers.

# every QAC object needs to have their own registers for things.
# this includes registers imported from parents.

class QAC():

    def __init__(self,parent=None):
        self.scale = 1
        self.instrs = []

        # scope and unnamed registers are guaranteed unique to this qac object

        # unordered
        self.scope_inputs  = {}  # "x":<reg>, or "x":[<reg>,<reg>] if array
        self.scope_outputs  = {}

        # ordered
        self.unnamed_inputs = [] # <reg>, <reg> (no arrays allowed)
        self.unnamed_outputs = []

        assert isinstance(parent,QAC) or (parent is None)
        self.parent = parent


    def __str__(self):
        return qac_to_string(self)

    ### shape strings
    def scope_lshape(self):
        d = self.scope_outputs
        out = []
        for k in sorted(d.keys()):
            s = k+":"
            if isinstance(d[k],list): s += str(d[k][0].dim)+"@@"+str(len(d[k]))
            else: s += str(d[k].dim)
            out.append(s)
        if len(out) == 0: return "Unit"
        return " @ ".join(out)

    def scope_rshape(self):
        d = self.scope_inputs
        out = []
        for k in sorted(d.keys()):
            s = k+":"
            if isinstance(d[k],list): s += str(d[k][0].dim)+"@@"+str(len(d[k]))
            else: s += str(d[k].dim)
            out.append(s)
        if len(out) == 0: return "Unit"
        return " @ ".join(out)

    def scope_shape(self):
        return self.scope_lshape() + " <- " + self.scope_rshape()

    def unnamed_lshape(self):
        if len(self.unnamed_outputs) == 0: return "Unit"
        return " @ ".join([str(reg.dim) for reg in self.unnamed_outputs])

    def unnamed_rshape(self):
        if len(self.unnamed_inputs) == 0: return "Unit"
        return " @ ".join([str(reg.dim) for reg in self.unnamed_inputs])

    def unnamed_shape(self):
        return self.unnamed_lshape() + " <- " + self.unnamed_rshape()

    def lshape(self):
        scope = self.scope_lshape()
        unnamed = self.unnamed_lshape()
        if scope == "Unit": return unnamed
        if unnamed == "Unit": return scope
        return unnamed +" @ "+scope

    def rshape(self):
        scope = self.scope_rshape()
        unnamed = self.unnamed_rshape()
        if scope == "Unit": return unnamed
        if unnamed == "Unit": return scope
        return unnamed +" @ "+scope

    def shape(self):
        return self.lshape() + " <- " + self.rshape()

    ##############################
    # Getting keys from parents

    # get key from parent scope
    def parentGet(self,key):
        if self.parent is None:
            raise KeyError()
        if key in self.parent.scope_outputs:
            return self.parent.scope_outputs[key]
        return self.parent.parentGet(key)

    # check if parent has key
    def parentHas(self,key):
        if self.parent is None: return False
        if key in self.parent.scope_outputs: return True
        return self.parent.parentHas(key)


    # takes a register from the parent scope
    # and puts it into scope_inputs so it becomes mutable in this block.
    # necessarily makes a copy, because every block needs to have unique registers.
    def promote(self, name):

        assert self.parentHas(name)

        # Check that 'name' is not in scope.
        # I can't check that it wasn't in scope temporarily.
        # That's probably fine, right? The scope analyzer will have caught that.
        assert name not in self.scope_inputs
        assert name not in self.scope_outputs

        prvreg = self.parentGet(name)

        if isinstance(prvreg, list):
            regs = []
            for i in range(len(prvreg)):
                reg = Register(prvreg[i].dim)
                reg.name_hints += prvreg[i].name_hints
                regs.append(reg)

            self.scope_inputs[name] = regs
            self.scope_outputs[name] = regs
        else:
            reg = Register(prvreg.dim)
            reg.name_hints += prvreg.name_hints
            self.scope_inputs[name] = reg
            self.scope_outputs[name] = reg

    # temporary name for a register
    # that is not in scope_outputs
    def tmpname(self, name):
        i = 0
        while True:
            i += 1
            tmpname = name+str(i)
            if tmpname in self.scope_outputs: continue
            if self.parentHas(tmpname): continue
            return tmpname

    ############################################## BASIC INSTRUCTIONS
    # DECLARE, DISCARD, ZERO
    # UNITARY, PHASE
    # INCREMENT

    def declare(self, name, dim, slots=None):
        assert name not in self.scope_outputs
        if slots is None:
            reg = Register(dim)
            reg.name_hints.append(name)
            self.instrs.append({"kind":"qac_declare", "reg":reg, "dim":dim})
            self.scope_outputs[name] = reg
            return

        assert int(slots) == slots
        assert slots > 0
        regs = []
        for i in range(slots):
            reg = Register(dim)
            reg.name_hints.append(name)
            reg.name_hints.append(str(i))
            self.instrs.append({"kind":"qac_declare", "reg":reg, "dim":dim})
            regs.append(reg)
        self.scope_outputs[name] = regs

    def discard(self, name):
        if name not in self.scope_outputs: self.promote(name)
        assert name in self.scope_outputs
        reg = self.scope_outputs[name]
        if isinstance(reg,list):
            for i in range(len(reg)):
                self.instrs.append({"kind":"qac_discard", "reg":reg[i]})
        else:
            self.instrs.append({"kind":"qac_discard", "reg":reg})
        del self.scope_outputs[name]


    def zero(self, name):
        if name not in self.scope_outputs: self.promote(name)
        assert name in self.scope_outputs
        reg = self.scope_outputs[name]
        if isinstance(reg,list):
            for i in range(len(reg)):
                self.instrs.append({"kind":"qac_zero", "reg":reg[i]})
        else:
            self.instrs.append({"kind":"qac_zero", "reg":reg})
        del self.scope_outputs[name]

    def unitary(self, name, unitary):
        if name not in self.scope_outputs: self.promote(name)
        assert name in self.scope_outputs

        # assert is a unitary matrix in list of lists form
        assert isinstance(unitary,list)
        for col in unitary:
            assert isinstance(col,list)
            assert len(col) == len(unitary)
            for i in range(len(col)):
                col[i] = complex(col[i])

        for i in range(len(unitary)):
            for j in range(len(unitary)):
                dot = 0
                for k in range(len(unitary)):
                    dot += unitary[i][k] * unitary[k][j].conj()
                if i == k: assert dot == 1
                else: assert dot == 0

        self.instrs.append({"kind":"qac_unitary", "reg":name, "mat":unitary})

    def phase(self, value):
        self.instrs.append({"kind":"qac_phase", "value":phase})


    ##################################################### INCREMENT


    def increment(self, regname, expn):
        # optimization: don't add zero.
        if expn["kind"] == "value_expn" and expn["value"] == 0: return

        if regname not in self.scope_outputs: self.promote(regname)
        register = self.scope_outputs[regname]

        # Do array reference decompiling.
        # {"target":<reg>, "key":<reg>, "sources":[<reg>,<reg>,<reg>]}
        array_refs = []

        # this recursive function does two things:
        # turn {"kind":"named_register_expn"} into {"kind":"register_expn"}
        # populate array_refs
        def process_expn(expn):
            assert expn["kind"] != "register_expn"

            # {"kind": "named_register_expn", "name":"x", "key": None, <int>, or "k" }
            # transform to: {"kind": "register_expn", "register":<reg> }
            if expn["kind"] == "named_register_expn":
                name, key = expn["name"], expn["key"]

                out = {"kind": "register_expn"}

                if name not in self.scope_outputs: self.promote(name)
                out["register"] = self.scope_outputs[name]

                if key is None:
                    assert not isinstance(out["register"], list)
                    return out

                if isinstance(key,int):
                    assert isinstance(out["register"], list)
                    assert key >= 0
                    assert key < len(out["register"])
                    out["register"] = out["register"][key]
                    return out

                if isinstance(key,str):
                    assert isinstance(out["register"], list)

                    if key not in self.scope_outputs: self.promote(key)
                    key = self.scope_outputs[key]

                    assert not isinstance(key, list)
                    assert key.dim <= len(out["register"])

                    target = Register(out["register"][0].dim)
                    target.name_hints.append(name)

                    array_ref = {"target":target, "key":key, "sources":[]}
                    for i in range(key.dim):
                        array_ref["sources"].append(out["register"][i])
                    array_refs.append(array_ref)

                    # TODO: test me!
                    out["register"] = target

                    return out

                assert False # unreachable

            # {"kind": "value_expn", "value":5j}
            if expn["kind"] == "value_expn":
                return {"kind": "value_expn", "value": expn["value"]}

            # {"kind": "sum_expn", "terms":[<linexp>] }
            if expn["kind"] == "sum_expn":
                return {
                    "kind": "sum_expn",
                    "terms":[process_expn(exp) for exp in expn["terms"]]
                }

            # {"kind": "negate_expn", "expn":<linexp> }
            if expn["kind"] == "negate_expn":
                return {"kind": "negate_expn", "expn":process_expn(expn["expn"]) }

            # {"kind": "adjoint_expn", "expn":<linexp> }
            if expn["kind"] == "adjoint_expn":
                return {"kind": "adjoint_expn", "expn":process_expn(expn["expn"]) }

            # {"kind": "product_expn", "terms":[<linexp>] }
            if expn["kind"] == "product_expn":
                return {
                    "kind": "product_expn",
                    "terms":[process_expn(exp) for exp in expn["terms"]]
                }

            # {"kind": "division_expn", "dividend":<linexp>, "divisor":5j }
            if expn["kind"] == "division_expn":
                return {"kind": "division_expn", "dividend":process_expn(expn["dividend"]), "divisor":expn["divisor"] }

            # {"kind": "modulo_expn", "dividend":<linexp>, "divisor":5 }
            if expn["kind"] == "modulo_expn":
                return {"kind": "modulo_expn", "dividend":process_expn(expn["dividend"]), "divisor":expn["divisor"] }

            # {"kind": "boolean_expn", "terms":[<linexp>, <string>, <linexp>, <string>, ...] }
            if expn["kind"] == "boolean_expn":
                out = {"kind": "boolean_expn", "terms":[] }
                for i in range(len(expn["terms"])):
                    if i % 2 == 0: out["terms"].append(process_expn(expn["terms"][i]))
                    else: out["terms"].append(expn["terms"][i])
                return out

            assert False # unreachable

        expn = process_expn(expn)

        # make a flag for doing if statements
        if len(array_refs) > 0:
            flag = Register(2)
            flag.name_hints.append("increment_flag")
            self.instrs.append({"kind":"qac_declare", "reg":flag, "dim":2})

        # helper function for doing array references
        def flip_flag(ref,i):
            self.instrs.append({"kind":"qac_increment", "reg":flag, "expn":{
                                "kind":"boolean_expn", "terms":[
                                    {"kind":"value_expn", "value": i}, "==",
                                    {"kind":"register_expn", "register":ref["key"]}]}})

        # declare array_ref targets, populate them.
        for ref in array_refs:
            self.instrs.append({"kind":"qac_declare", "reg":ref["target"], "dim":ref["target"].dim})
            for i in range(ref["key"].dim):
                flip_flag(ref,i)
                self.instrs.append({"kind":"qac_if", "cond":flag, "instructions":[
                    {"kind":"qac_increment", "reg":ref["target"], "expn":{"kind":"register_expn", "register":ref["sources"][i]}}
                ]})
                flip_flag(ref,i)

        # do the increment.
        self.instrs.append({"kind":"qac_increment", "reg":register, "expn":expn})

        # uncompute target registers of array_refs
        for ref in array_refs:
            for i in range(ref["key"].dim):
                flip_flag(ref,i)
                self.instrs.append({"kind":"qac_if", "cond":flag, "instructions":[
                    {"kind":"qac_increment", "reg":ref["target"], "expn":{"kind":"negate_expn", "expn":{"kind":"register_expn", "register":ref["sources"][i]}}}
                ]})
                flip_flag(ref,i)
            self.instrs.append({"kind":"qac_zero", "reg":ref["target"]})

        # remove the flag
        if len(array_refs) > 0:
            self.instrs.append({"kind":"qac_zero", "reg":flag})


    def increment_array_fixed(self, regname, idx, expn):
        if regname not in self.scope_outputs: self.promote(regname)
        assert isinstance(self.scope_outputs[regname], list)
        assert isinstance(idx, int)

        tmpname = self.tmpname(regname+"_increment")
        self.scope_outputs[tmpname] = self.scope_outputs[regname][idx]

        self.increment(tmpname, expn)

        del self.scope_outputs[tmpname]

    def increment_array_variable(self, regname, keyname, expn):
        if regname not in self.scope_outputs: self.promote(regname)
        reg = self.scope_outputs[regname]
        assert isinstance(reg, list)

        if keyname not in self.scope_outputs: self.promote(keyname)
        key = self.scope_outputs[keyname]
        assert not isinstance(key, list)
        assert key.dim <= len(self.scope_outputs[regname])

        tmpname = self.tmpname(regname+"_increment")
        tmp = Register(reg[0].dim)
        tmp.name_hints.append(regname+"_increment")
        self.scope_outputs[tmpname] = tmp
        self.instrs.append({"kind":"qac_declare", "reg":tmp, "dim":tmp.dim})

        flag = Register(2)
        flag.name_hints.append("increment_flag")
        self.instrs.append({"kind":"qac_declare", "reg":flag, "dim":2})

        def flip_flag(i):
            self.instrs.append({"kind":"qac_increment", "reg":flag, "expn":{
                                "kind":"boolean_expn", "terms":[
                                    {"kind":"value_expn", "value": i}, "==",
                                    {"kind":"register_expn", "register":key}]}})

        for i in range(key.dim):
            flip_flag(i)
            self.instrs.append({"kind":"qac_if", "cond":flag, "instructions":[
                {"kind":"qac_increment", "reg":tmp, "expn":{"kind":"register_expn", "register":reg[i]}}
            ]})
            flip_flag(i)

        self.increment(tmpname, expn)

        for i in range(key.dim):
            flip_flag(i)
            self.instrs.append({"kind":"qac_if", "cond":flag, "instructions":[
                {"kind":"qac_increment", "reg":tmp, "expn":{"kind":"negate_expn", "expn":{"kind":"register_expn", "register":reg[i]}}}
            ]})
            flip_flag(i)

        self.instrs.append({"kind":"qac_zero", "reg":flag})

        self.instrs.append({"kind":"qac_zero", "reg":tmp})
        del self.scope_outputs[tmpname]


    ##################################################### SYMMETRIZATION

    # If this is a block of code with no unnamed registers,
    # then make all the scope_inputs and scope_outputs registers match one another.
    # Required for if and repeat statements, but also some other operations
    def symmetrize(self):
        # assert qac.scope_input and qac.scope_outputs to have the same set of keys
        for key in self.scope_inputs.keys():  assert key in self.scope_outputs
        for key in self.scope_outputs.keys(): assert key in self.scope_inputs
        # and that the dimensions match
        for key in self.scope_inputs.keys():
            assert self.scope_inputs[key].dim == self.scope_outputs[key].dim

        to_sub = {}
        for key in self.scope_inputs.keys():
            if self.scope_inputs[key] == self.scope_outputs[key]:
                # already match
                continue

            if self.scope_inputs[key] not in self.scope_outputs.values():
                # relabel output to input at end
                self.instrs.insert(0,{"kind":"qac_rename",
                                      "source": self.scope_outputs[key],
                                      "target": self.scope_inputs[key]})
                self.scope_outputs[key] = self.scope_inputs[key]
                continue

            if self.scope_outputs[key] not in self.scope_inputs.values():
                # relabel output to input at beginning
                self.instrs.append({"kind":"qac_rename",
                                    "source": self.scope_outputs[key],
                                    "target": self.scope_inputs[key]})
                self.scope_inputs[key] = self.scope_outputs[key]
                continue

            # make a new register.
            reg = Register(self.scope_inputs[key].dim)
            reg.name_hints += self.scope_inputs[key].name_hints

            self.instrs.insert(0,{"kind":"qac_rename",
                                  "source": reg,
                                  "target": self.scope_inputs[key]})

            self.instrs.append({"kind":"qac_rename",
                                "source": self.scope_outputs[key],
                                "target": reg})

            to_sub[key] = reg

        for key in to_sub.keys():
            self.scope_inputs[key] = to_sub[key]
            self.scope_outputs[key] = to_sub[key]

    ##################################################### IF STATEMENTS

    def if_statement(self, ctrlname, qac):

        assert len(qac.unnamed_inputs) == 0
        assert len(qac.unnamed_outputs) == 0
        qac.symmetrize()

        assert qac.parent == self

        # promote any qac.inputs if needed, and wire them up
        for key in qac.scope_inputs.keys():
            if key not in self.scope_outputs: self.promote(key)
            qac.scope_inputs[key].substitute(self.scope_outputs[key])

        # obtain the condition register
        if ctrlname not in self.scope_outputs: self.promote(ctrlname)
        cond = self.scope_outputs[ctrlname]

        # can't use an entire array as a control
        assert not isinstance(cond,list)

        if_instrs = []

        # take into account scale via a postselect
        if qac.scale != 1:
            if qac.scale < 1:
                p = qac.scale
            else:
                self.scale *= qac.scale
                p = 1/qac.scale

            post_flag = Register(2)
            post_flag.name_hints.append("postselection_flag")
            if_instrs.append({"kind":"qac_declare", "reg":post_flag, "dim":2})
            a1, a2 = p**(1/2), (1-p)**(1/2)
            mat = [[a1, -a2], [a2, a1]]
            if_instrs.append({"kind":"qac_unitary", "reg":post_flag, "mat":mat})
            if_instrs.append({"kind":"qac_zero", "reg":post_flag})

        if_instrs += qac.instrs

        self.instrs.append({"kind":"qac_if", "cond":cond, "instructions": if_instrs})


    # controls [{"expn":<expn>, "qac":<qac>}]
    def compound_if_statement(self, controls, else_qac=None):
        # simple case: there is just one condition
        if len(controls) == 1:
            flag = Register(2)
            flag.name_hints.append("if_flag")
            self.instrs.append({"kind":"qac_declare", "reg":flag, "dim":2})
            flagname = self.tmpname("if_flag")
            self.scope_outputs[flagname] = flag

            self.increment(flagname, controls[0]["expn"])
            self.if_statement(flagname, controls[0]["qac"])

            if else_qac is None:
                # just uncompute the flag
                self.increment(flagname, controls[0]["expn"])
            else:
                # flip the flag

                self.increment(flagname, {"kind":"value_expn", "value":1})

                self.if_statement(flagname, else_qac)

                # uncompute the flag accounting for the extra flip

                self.increment(flagname, {"kind":"sum_expn", "terms":[{"kind":"value_expn", "value":1},
                    controls[0]["expn"]]})


            self.instrs.append({"kind":"qac_zero", "reg":flag})
            del self.scope_outputs[flagname]
            return

        # complex case: multiple instruction blocks

        # a flag counting how many expressions have been true so far
        prv_count = Register(len(controls)+1)
        prv_count.name_hints.append("if_prv_count")
        self.instrs.append({"kind":"qac_declare", "reg":prv_count})

        # allocate a flag for all the if statements, give it a flagname
        flag = Register(2)
        post_flag.name_hints.append("if_flag")
        self.instrs.append({"kind":"qac_declare", "reg":flag})
        flagname = self.tmpname("if_flag")
        self.scope_outputs[flagname] = flag

        for i in range(len(controls)):
            # if controls[i]["expn"] and prv_count == 0, perform the instructions
            self.increment(flagname, {"kind":"product_expn", "terms":[
                                    controls[i]["expn"],
                                    {"kind":"boolean_expn", "terms":[
                                        {"kind":"value_expn", "value": 0}, "==",
                                        {"kind":"register_expn", "register":prv_count}]}
                                ]})

            self.if_statement(flagname, controls[i]["qac"])
            self.increment(flagname, {"kind":"product_expn", "terms":[
                                    controls[i]["expn"],
                                    {"kind":"boolean_expn", "terms":[
                                        {"kind":"value_expn", "value": 0}, "==",
                                        {"kind":"register_expn", "register":prv_count}]}
                                ]})

            # if controls[i]["expn"], increment prv_count
            self.increment(flagname,controls[i]["expn"])
            self.instrs.append({"kind":"qac_if", "cond":flag, "instructions": [
                {"kind":"qac_increment", "reg":prv_count, "expn":{"kind":"value_expn", "value": 1}}
            ]})
            self.increment(flagname,controls[i]["expn"])

        if else_expn is not None:
            # if prv_count is still zero, perform the else expn
            self.increment(flagname, {"kind":"boolean_expn", "terms":[
                                        {"kind":"value_expn", "value": 0}, "==",
                                        {"kind":"register_expn", "register":prv_count}]})
            self.if_statement(flagname, else_expn)
            self.increment(flagname, {"kind":"boolean_expn", "terms":[
                                        {"kind":"value_expn", "value": 0}, "==",
                                        {"kind":"register_expn", "register":prv_count}]})

        # uncompute prv_count
        for i in range(len(controls)):
            self.increment(flagname,controls[i]["expn"])
            self.instrs.append({"kind":"qac_if", "cond":flag, "instructions": [
                {"kind":"qac_increment", "reg":prv_count, "expn":{
                    "kind":"negate_expn", "expn":{"kind":"value_expn", "value": 1}}}
            ]})
            self.increment(flagname,controls[i]["expn"])

        # remove the ancillae
        del self.scope_outputs[flagname]
        self.instrs.append({"kind":"qac_zero", "reg":flag})
        self.instrs.append({"kind":"qac_zero", "reg":prv_count})


    ##################################################### REPEAT STATEMENTS

    def repeat_fixed(self, qac, count):
        assert qac.parent == self
        assert int(count) == count
        assert count >= 0

        if count == 0: return

        assert len(qac.unnamed_inputs) == 0
        assert len(qac.unnamed_outputs) == 0

        qac.symmetrize()

        for i in range(count):
            self.scale *= qac.scale

            if i == 0: tmp = qac
            else: tmp = block_copy(qac)

            # promote any qac.inputs if needed, and wire them up
            for key in tmp.scope_inputs.keys():
                if key not in self.scope_outputs: self.promote(key)
                tmp.scope_inputs[key].substitute(self.scope_outputs[key])

            for instr in tmp.instrs:
                self.instrs.append(instr)


    def repeat_variable(self, qac, regname):
        assert len(qac.unnamed_inputs) == 0
        assert len(qac.unnamed_outputs) == 0
        assert qac.parent == self

        qac.symmetrize()

        # obtain the count register
        if regname not in self.scope_outputs: self.promote(regname)
        count = self.scope_outputs[regname]
        assert not isinstance(count, list)

        flag = Register(2)
        flag.name_hints.append("repeat_flag")
        flagname = self.tmpname("repeat_flag")
        self.scope_outputs[flagname] = flag
        self.instrs.append({"kind":"qac_declare", "reg":flag, "dim":2})

        for i in range(count.dim-1):
            if i == 0: tmp = qac
            else: tmp = block_copy(qac)

            # promote any qac.inputs if needed, and wire them up
            for key in tmp.scope_inputs.keys():
                if key not in self.scope_outputs: self.promote(key)
                tmp.scope_inputs[key].substitute(self.scope_outputs[key])

            self.instrs.append({"kind":"qac_increment", "reg":flag, "expn":{
                                "kind":"boolean_expn", "terms":[
                                    {"kind":"value_expn", "value": i}, "<",
                                    {"kind":"register_expn", "register":count}]}})

            self.if_statement(flagname, tmp)

            self.instrs.append({"kind":"qac_increment", "reg":flag, "expn":{
                                "kind":"boolean_expn", "terms":[
                                    {"kind":"value_expn", "value": i}, "<",
                                    {"kind":"register_expn", "register":count}]}})

        self.instrs.append({"kind":"qac_zero", "reg":flag})
        del self.scope_outputs[flagname]



    ############################################################ SCALAR AND INIT

    def scalar_instr(self, qac):
        assert len(qac.unnamed_outputs) == 0
        assert len(qac.unnamed_inputs) == 0
        assert qac.parent == self

        # promote any qac.inputs if needed, and wire them up
        for key in qac.scope_inputs.keys():
            if key not in self.scope_outputs: self.promote(key)
            qac.scope_inputs[key].substitute(self.scope_outputs[key])

            # block consumes the variables in my scope
            del self.scope_outputs[key]

        # and outputs new variables that were actually only referenced
        for key in qac.scope_outputs.keys():
            assert key not in self.scope_outputs
            self.scope_outputs[key] = qac.scope_outputs[key]

        self.scale *= qac.scale

        for instr in qac.instrs:
            self.instrs.append(instr)


    def init_names(self,target_names, qac):
        assert qac.parent == self
        assert len(qac.unnamed_inputs) == 0
        assert len(qac.unnamed_outputs) == len(target_names)

        # promote targets and check dimensions
        for i in range(len(target_names)):
            targ = target_names[i]
            if targ not in self.scope_outputs: self.promote(targ)
            assert not isinstance(self.scope_outputs[targ], list)
            assert qac.unnamed_outputs[i].dim == self.scope_outputs[targ].dim

        #########################
        # promote any qac.inputs if needed, and wire them up
        for key in qac.scope_inputs.keys():
            if key not in self.scope_outputs: self.promote(key)
            qac.scope_inputs[key].substitute(self.scope_outputs[key])

        # block consumes the variables in my scope
        for key in qac.scope_inputs.keys():
            del self.scope_outputs[key]

        # and outputs new variables that were actually only referenced
        for key in qac.scope_outputs.keys():
            assert key not in self.scope_outputs
            self.scope_outputs[key] = qac.scope_outputs[key]

        self.scale *= qac.scale

        for instr in qac.instrs:
            self.instrs.append(instr)

        # Align the scope
        # 1. discard all the targets that havent been consumed
        # 2. make qac.unnamed_outputs[i] be called target_names[i] in scope
        for i in range(len(target_names)):
            if target_names[i] in self.scope_outputs:
                self.instrs.append({"kind":"qac_discard", "reg":self.scope_outputs[target_names[i]]})
            self.scope_outputs[target_names[i]] = qac.unnamed_outputs[i]
            if target_names[i] not in qac.unnamed_outputs[i].name_hints:
                qac.unnamed_outputs[i].name_hints.append(target_names[i])

    # targets = [{"name":<name>, "key":<int>, <name> or None}]
    def init_instr(self, targets, qac):

        # this register only gets declared if we actually need it.
        flag = None
        def flip_flag(key_reg,i):
            self.instrs.append({"kind":"qac_increment", "reg":flag, "expn":{
                                "kind":"boolean_expn", "terms":[
                                    {"kind":"value_expn", "value": i}, "==",
                                    {"kind":"register_expn", "register":key_reg}]}})

        target_names = []
        for target in targets:
            if target["key"] is None:
                target_names.append(target["name"])
                continue

            if isinstance(target["key"], int):
                regname = target["name"]
                if regname not in self.scope_outputs: self.promote(regname)
                reg = self.scope_outputs[regname]
                assert isinstance(reg,list)
                assert target["key"] >= 0
                assert target["key"] < len(reg)

                tmpname = self.tmpname(target["name"])
                target["tmp"] = tmpname
                self.unnamed_outputs[tmpname] = reg[target["key"]]
                target_names.append(tmpname)
                continue

            if isinstance(target["key"], str):
                regname = target["name"]
                if regname not in self.scope_outputs: self.promote(regname)
                reg = self.scope_outputs[regname]
                assert isinstance(reg,list)

                keyname = target["key"]
                if keyname not in self.scope_outputs: self.promote(keyname)
                key = self.scope_outputs[keyname]
                assert not isinstance(key,list)
                assert key.dim <= len(reg)

                # declare the flag if needed
                if flag is None:
                    flag = Register(2)
                    flag.name_hints.append("init_array_flag")
                    self.instrs.append({"kind":"qac_declare", "reg":flag, "dim":2})

                tmpname = self.tmpname(target["name"])
                target["tmp"] = tmpname
                tmpreg = Register(reg[0].dim)
                self.instrs.append({"kind":"qac_declare", "reg":tmpreg, "dim":tmpreg.dim})
                self.unnamed_outputs[tmpname] = tmpreg
                target_names.append(tmpname)

                # conditionally swap the target register onto tmpreg
                for i in range(key.dim):
                    flip_flag(key,i)
                    self.instrs.append({"kind":"qac_if", "cond":flag, "instructions":[
                        {"kind":"qac_swap", "reg1":tmpreg, "reg2":reg[i]}
                    ]})
                    flip_flag(key,i)

            assert False # unreachable

        self.init_names(target_names, qac)

        for target in targets:
            if isinstance(target["key"], int):
                del self.scope_outputs[target["tmp"]]

            if isinstance(target["key"], str):
                reg = self.scope_outputs[target["name"]]
                key = self.scope_outputs[target["key"]]
                tmpreg = self.scope_outputs[target["tmp"]]

                for i in range(key.dim):
                    flip_flag(key,i)
                    self.instrs.append({"kind":"qac_if", "cond":flag, "instructions":[
                        {"kind":"qac_swap", "reg1":tmpreg, "reg2":reg[i]}
                    ]})
                    flip_flag(key,i)

                del self.scope_outputs[target["tmp"]]

        # get rid of the flag
        if flag is not None:
            self.instrs.append({"kind":"qac_zero", "reg":flag})


    ############################################### ASSIGN INSTRUCTIONS


    def assign(self, regname, expn, undo=None):

        if regname not in self.scope_outputs: out.promote(regname)
        reg = self.scope_outputs[regname]

        # find a temporary name
        i = 0
        while True:
            i += 1
            tmpname = "tmp"+str(i)
            if tmpname in self.scope_outputs: continue
            if self.parentHas(tmpname): continue
            break

        # initialize the temporary register
        # and compute the new value into it
        tmpreg = Register(reg.dim)
        tmpreg.name_hints += reg.name_hints
        self.instrs.append({"kind":"qac_declare", "reg":tmpreg, "dim": reg.dim})
        self.scope_outputs[tmpname] = tmpreg
        self.increment(tmpname, expn)

        if undo is None:

            # trace out the old register
            self.instrs.append({"kind":"qac_discard", "reg":reg})

            # make the temporary register occupy the previous registister's name in scope
            self.scope_outputs[regname] = tmpreg
            del self.scope_outputs[tmpname]

        else:
            # declare tmp
            # tmp += expn(x)
            # x -= undo(tmp)
            # zero x
            # rename x to tmp

            # equivalent to:

            # declare tmp
            # tmp += expn(x)
            # swap x and tmp
            # tmp -= undo(x)
            # zero tmp

            # initialize the temporary register
            # and compute the new value into it
            tmpreg = Register(reg.dim)
            tmpreg.name_hints += reg.name_hints
            self.instrs.append({"kind":"qac_declare", "reg":tmpreg, "dim": reg.dim})
            self.scope_outputs[tmpname] = tmpreg
            self.increment(tmpname, expn)

            # swap the register's names in scope.
            self.scope_outputs[tmpname], self.scope_outputs[regname] = self.scope_outputs[regname], self.scope_outputs[tmpname]

            # decrement tmp value
            self.increment(tmpname, {"kind":"negate_expn", "expn":undo})

            # zero out register with temporary name
            self.instrs.append({"kind":"qac_zero", "reg":self.scope_outputs[tmpname]})
            del self.scope_outputs[tmpname]


    def assign_array_fixed(self, arrayname, idx, expn):
        if arrayname not in self.scope_outputs: out.promote(arrayname)
        arrayreg = self.scope_outputs[arrayname]
        assert isinstance(array,list)

        idx = complex(idx)
        assert idx == int(idx.real)
        idx = int(idx.real)
        assert idx <= len(array)

        # find a temporary name
        i = 0
        while True:
            i += 1
            tmpname = "tmp"+str(i)
            if tmpname in self.scope_outputs: continue
            if self.parentHas(tmpname): continue
            break

        # also give array[idx] that temporary name
        self.scope_outputs[tmpname] = array[idx]

        self.assign_instr(self,tmpname,expn)

        # now replace array[idx] with the new register
        array[idx] = scope.scope_outputs[tmpname]
        del scope.scope_outputs[tmpname]


    def assign_array_variable(self, arrayname, keyname, expn):

        # there are two ways of implementing this.
        # one:
        #    conditionally swap the correct register into a temporary register
        #    perform the assignment on the temporary register
        #    undo the conditional swap
        # two:
        #    for each register in the array, conditionally perform the assignment.

        # I believe option one is better, because that way we avoid lots of conditional discards.
        # But this remains to be tested in practice.

        # I don't really have an undo version of this function, because that
        # insnt properly supported lexially.


        if arrayname not in self.scope_outputs: self.promote(arrayname)
        arrayreg = self.scope_outputs[arrayname]
        assert isinstance(array,list)

        if keyname not in self.scope_outputs: self.promote(keyname)
        keyreg = self.scope_outputs[keyname]

        assert keyreg.dim <= len(array)

        # find a temporary name
        i = 0
        while True:
            i += 1
            tmpname = "tmp"+str(i)
            if tmpname in self.scope_outputs: continue
            if self.parentHas(tmpname): continue
            break

        tmpreg = Register(array[0].dim)
        tmpreg.name_hints += array[0].name_hints
        self.instrs.append({"kind":"qac_declare", "reg":tmpreg, "dim": tmpreg.dim})
        self.scope_outputs[tmpname] = tmpreg

        flag = Register(2)
        flag.name_hints.append("array_assign_flag")
        self.instrs.append({"kind":"qac_declare", "reg":flag, "dim": 2})

        for i in range(keyreg.dim):
            self.instrs.append({"kind":"qac_increment", "reg":flag, "expn":{
                                "kind":"boolean_expn", "terms":[
                                    {"kind":"value_expn", "value": i}, "==",
                                    {"kind":"register_expn", "register":keyreg}]}})

            self.instrs.append({"kind":"qac_if", "cond":flag, "instructions":[
                {"kind":"qac_swap", "reg1":tmpreg, "reg2": array[i]}
            ]})

            self.instrs.append({"kind":"qac_increment", "reg":flag, "expn":{
                                "kind":"boolean_expn", "terms":[
                                    {"kind":"value_expn", "value": i}, "==",
                                    {"kind":"register_expn", "register":keyreg}]}})

        self.assign_instr(self,tmpname,expn)
        tmpreg = self.scope_outputs[tmpname] # assign_instr exchanges the register

        for i in range(keyreg.dim):
            self.instrs.append({"kind":"qac_increment", "reg":flag, "expn":{
                                "kind":"boolean_expn", "terms":[
                                    {"kind":"value_expn", "value": i}, "==",
                                    {"kind":"register_expn", "register":keyreg}]}})

            self.instrs.append({"kind":"qac_if", "cond":flag, "instructions":[
                {"kind":"qac_swap", "reg1":tmpreg, "reg2": array[i]}
            ]})

            self.instrs.append({"kind":"qac_increment", "reg":flag, "expn":{
                                "kind":"boolean_expn", "terms":[
                                    {"kind":"value_expn", "value": i}, "==",
                                    {"kind":"register_expn", "register":keyreg}]}})



        self.instrs.append({"kind":"qac_zero", "reg":flag})
        self.instrs.append({"kind":"qac_zero", "reg":tmpreg})
        del self.scope_outputs[tmpname]



    ############################################## MAKING NEW BLOCKS WITH THIS ONE AS A PARENT


    def block_scalar(parent,z):
        z = complex(z)
        out = QAC()
        out.parent = parent
        out.scale *= abs(z)
        if (z/abs(z) != 1):
            out.instrs.append({
                "kind":"qac_phase",
                "value": z/abs(z)
            })
        return out

    def block_create(parent, expn, dim):
        out = QAC()
        out.parent = parent
        assert int(dim) == dim
        assert dim > 1

        # find a temporary name
        i = 0
        while True:
            tmpname = "tmp"+str(i)
            if not out.parentHas(tmpname): break
            i += 1

        # Temporarily make a named register
        reg = Register(dim)
        reg.name_hints.append("create")
        out.instrs.append({"kind":"qac_declare", "reg":reg, "dim":dim})
        out.scope_outputs[tmpname] = reg

        # that way I can just invoke increment
        out.increment(tmpname, expn)

        # and then un-name it.
        del out.scope_outputs[tmpname]
        out.unnamed_outputs.append(reg)

        # This is a little bit of a cludge - ideally I'd have a version of increment
        # that supports the target being unnamed.

        return out

    def block_consume(parent, name):
        out = QAC()
        out.parent = parent
        out.promote(name)

        out.unnamed_outputs.append(out.scope_outputs[name])
        del out.scope_outputs[name]
        return out

    def block_cast(parent, name, key):
        out = QAC()
        out.parent = parent
        out.promote(name)

        if key is None:
            reg = out.scope_outputs[name]
        elif isinstance(key,int):

            regarray = out.scope_outputs[name]
            assert key >= 0 and key < len(regarray)

            reg = regarray[key]

        else:
            assert isinstance(key,str)

            regarray = out.scope_outputs[name]

            out.promote(key)
            keyreg = out.scope_outputs[key]

            assert keyreg.dim <= len(regarray)

            reg = Register(regarray[0].dim)
            reg.name_hints.append("cast")
            out.instrs.append({"kind":"qac_declare", "reg":reg, "dim":reg.dim})

            array_flag = Register(2)
            array_flag.name_hints.append("cast_array_flag")
            out.instrs.append({"kind":"qac_declare", "reg":array_flag, "dim":2})

            for i in range(keyreg.dim):
                self.instrs.append({"kind":"qac_increment", "reg":array_flag, "expn":{
                                    "kind":"boolean_expn", "terms":[
                                    {"kind":"value_expn", "value": i}, "==",
                                    {"kind":"register_expn", "register":keyreg}]}})

                self.instrs.append({"kind":"qac_increment", "reg":reg, "expn":{
                                    "kind":"register_expn", "register":regarray[i]}})

                self.instrs.append({"kind":"qac_increment", "reg":array_flag, "expn":{
                                    "kind":"boolean_expn", "terms":[
                                    {"kind":"value_expn", "value": i}, "==",
                                    {"kind":"register_expn", "register":keyreg}]}})

        out.scale = (reg.dim-1)

        # to use as a control
        cond_flag = Register(2)
        cond_flag.name_hints.append("cast_control_flag")
        out.instrs.append({"kind":"qac_declare", "reg":cond_flag, "dim":2})

        # to postselect
        post_flag = Register(2)
        post_flag.name_hints.append("cast_postselect_flag")
        out.instrs.append({"kind":"qac_declare", "reg":post_flag, "dim":2})

        for i in range(reg.dim):
            p = i/out.scale
            a1, a2 = p**(1/2), (1-p)**(1/2)
            mat = [[a1, -a2], [a2, a1]]

            out.instrs.append({"kind":"qac_increment", "reg":cond_flag, "expn":{
                               "kind":"boolean_expn", "terms":[
                                    {"kind":"value_expn", "value": i}, "==",
                                    {"kind":"register_expn", "register":reg}]}})

            out.instrs.append({"kind":"qac_if", "cond":cond_flag, "instructions":[
                {"kind":"qac_unitary", "reg":post_flag, "mat":mat}
            ]})

            out.instrs.append({"kind":"qac_increment", "reg":cond_flag, "expn":{
                               "kind":"boolean_expn", "terms":[
                                    {"kind":"value_expn", "value": i}, "==",
                                    {"kind":"register_expn", "register":reg}]}})

        out.instrs.append({"kind":"qac_zero", "reg":post_flag})
        out.instrs.append({"kind":"qac_zero", "reg":cond_flag})

        if key is not None and isinstance(key,str):
            for i in range(keyreg.dim):
                self.instrs.append({"kind":"qac_increment", "reg":array_flag, "expn":{
                                    "kind":"boolean_expn", "terms":[
                                    {"kind":"value_expn", "value": i}, "==",
                                    {"kind":"register_expn", "register":keyreg}]}})

                self.instrs.append({"kind":"qac_increment", "reg":reg, "expn":{
                                    "kind":"negate_expn", "expn":{
                                    "kind":"register_expn", "register":regarray[i]}}})

                self.instrs.append({"kind":"qac_increment", "reg":array_flag, "expn":{
                                    "kind":"boolean_expn", "terms":[
                                    {"kind":"value_expn", "value": i}, "==",
                                    {"kind":"register_expn", "register":keyreg}]}})

            out.instrs.append({"kind":"qac_zero", "reg":reg, "dim":reg.dim})
            out.instrs.append({"kind":"qac_zero", "reg":array_flag, "dim":2})

        return out

    def block_boolean_cast(parent, boolean_expn):
        assert boolean_expn["kind"] == "boolean_expn"

        out = QAC()
        out.parent = parent

        # find a temporary name
        i = 0
        while True:
            i += 1
            tmpname = "tmp"+str(i)
            if out.parentHas(tmpname): continue
            break

        tmpreg = Register(2)
        tmpreg.name_hints.append("bool_cast_flag")
        out.instrs.append({"kind":"qac_declare", "reg":tmpreg, "dim": 2})
        out.scope_outputs[tmpname] = tmpreg
        out.increment(tmpname, {"kind":"boolean_expn", "terms":[
                                    boolean_expn, "==",
                                    {"kind":"value_expn", "value": 0}
                                ]})

        out.zero(tmpname)

        return out

    def block_userspace(parent, blq):
        if len(blq.qac) > 1: qac = block_add(*[block_copy(qac) for qac in blq.qac])
        else: qac = block_copy(blq.qac[0])
        assert qac.parent == None
        qac.parent = parent
        return qac

###############################

# In general, the inputs to these functions can't be used again for anything else.
# That is because they make the basically indispensable assumption that all the registers
# in the input blocks are unique. But register manipulation merges the register identities,
# so that assumption is no longer true after it is used.

def block_add(*blocks):
    assert len(blocks) > 1

    b0 = blocks[0]

    out = QAC()

    for b in blocks:
        while b.parent in blocks:
            b.parent = b.parent.parent
        if b.parent is not None:
            out.parent = b.parent
            break

    # check some things
    for b in blocks:
        if b.parent is None: b.parent = out.parent
        assert b.parent == out.parent # all blocks must have the same parent.

        for i in range(len(b.unnamed_inputs)):
            assert b.unnamed_inputs[i].dim == b0.unnamed_inputs[i].dim

        for i in range(len(b.unnamed_outputs)):
            assert b.unnamed_outputs[i].dim == b0.unnamed_outputs[i].dim

    # wire up named inputs
    for b in blocks:
        b.symmetrize()
        for key in b.scope_inputs.keys():
            if key not in out.scope_outputs: out.promote(key)
            out.scope_inputs[key].substitute(b.scope_outputs[key])

    ##################################################
    # We will be creating lots of qac_if statements for the LCU circuit.
    # Each if_statement can't change the scope: it must have the same registers before and after.
    # The given blocks can have more inputs than outputs, or permute their registers, etc.
    # So we can't use the given blocks' registers as a starting point. We need
    # an official set of registers that are present at the start and end of each if statement,
    # and a template to map the given blocks onto them.

    for i in range(len(b0.unnamed_inputs)):
        r = Register(b0.unnamed_inputs[i].dim)
        for b in blocks:
            for name in b.unnamed_inputs[i].name_hints:
                if name not in r.name_hints:
                    r.name_hints.append(name)
        out.unnamed_inputs.append(r)

    for i in range(len(b0.unnamed_outputs)):
        r = Register(b0.unnamed_outputs[i].dim)
        for b in blocks:
            for name in b.unnamed_outputs[i].name_hints:
                if name not in r.name_hints:
                    r.name_hints.append(name)
        out.unnamed_outputs.append(r)

    #### make control register
    ctrl = Register(len(blocks))
    ctrl.name_hints.append("sum_ctrl")
    out.instrs.append({"kind":"qac_declare", "reg":ctrl, "dim":len(blocks)})

    out.scale = 0
    for b in blocks:
        out.scale += b.scale

    # set the first column
    mat = [[0 for j in range(len(blocks))] for j in range(len(blocks))]
    for i in range(len(blocks)):
        mat[i][0] = np.sqrt(blocks[i].scale/out.scale)

    # complete the matrix using gram-schmidt
    for i in range(1,len(blocks)):
        col = [0 for j in range(len(blocks))]
        col[i] = 1
        for j in range(i):
            for k in range(len(blocks)):
                col[k] -= mat[k][j] * mat[i][j]
        norm = 0
        for k in range(len(blocks)):
            norm += col[k]*col[k]
        for k in range(len(blocks)):
            col[k] /= np.sqrt(norm)
            mat[k][i] = col[k]

    out.instrs.append({"kind":"qac_unitary", "reg":ctrl, "mat":mat})

    # make a flag for the if statements
    flag = Register(2)
    flag.name_hints.append("sum_flag")
    out.instrs.append({"kind":"qac_declare", "reg":flag, "dim":2})

    # declare the outputs, so we can swap to them
    for j in range(len(out.unnamed_outputs)):
        out.instrs.append({"kind":"qac_declare",
            "reg":out.unnamed_outputs[j], "dim":out.unnamed_outputs[j].dim})

    # apply select
    for i in range(len(blocks)):


        out.instrs.append({"kind":"qac_increment", "reg":flag, "expn":{
                        "kind":"boolean_expn", "terms":[
                        {"kind":"value_expn", "value": i}, "==",
                        {"kind":"register_expn", "register":ctrl}]}})

        instrs = []

        for j in range(len(blocks[i].unnamed_inputs)):
            instrs.append({"kind":"qac_declare",
                "reg":blocks[i].unnamed_inputs[j], "dim":blocks[i].unnamed_inputs[j].dim})


        for j in range(len(out.unnamed_inputs)):
            instrs.append({"kind":"qac_swap", "reg1":out.unnamed_inputs[j],
                                              "reg2":blocks[i].unnamed_inputs[j]})

        instrs += blocks[i].instrs

        for j in range(len(out.unnamed_outputs)):
            instrs.append({"kind":"qac_swap", "reg1":blocks[i].unnamed_outputs[j],
                                              "reg2":out.unnamed_outputs[j]})

        for j in range(len(blocks[i].unnamed_outputs)):
            instrs.append({"kind":"qac_zero",
                "reg":blocks[i].unnamed_outputs[j]})

        out.instrs.append({"kind":"qac_if", "cond":flag, "instructions":instrs})

        out.instrs.append({"kind":"qac_increment", "reg":flag, "expn":{
                        "kind":"boolean_expn", "terms":[
                        {"kind":"value_expn", "value": i}, "==",
                        {"kind":"register_expn", "register":ctrl}]}})

    # remove the input registers
    for j in range(len(out.unnamed_inputs)):
        out.instrs.append({"kind":"qac_zero", "reg":out.unnamed_inputs[j]})

    # get rid of the flag
    out.instrs.append({"kind":"qac_zero", "reg":flag})

    # uncompute the ctrl register, by inverting the unitary
    inv = []
    for i in range(len(blocks)):
        row = []
        for j in range(len(blocks)):
            row.append(mat[j][i].conjugate())
        inv.append(row)


    out.instrs.append({"kind":"qac_unitary", "reg":ctrl, "mat":inv})
    out.instrs.append({"kind":"qac_zero", "reg":ctrl})

    return out


def block_mul(b2,b1):
    if b1.parent == b2: b1.parent = b2.parent
    if b2.parent == b1: b2.parent = b1.parent
    if b1.parent is None: b1.parent = b2.parent
    if b2.parent is None: b2.parent = b1.parent
    assert b1.parent == b2.parent

    # check that keys consumed by b1 are not consumed again by b2
    for key in b1.scope_inputs.keys():
        if key in b1.scope_outputs: # b1 might have changed the dimension
            if key in b2.scope_inputs:
                assert b2.scope_inputs[key].dim == b1.scope_outputs[key].dim
            continue
        assert key not in b2.scope_inputs

    # check that b2.unnamed_inputs matches b1.unnamed_outputs in dimension
    assert len(b2.unnamed_inputs) == len(b1.unnamed_outputs)
    for i in range(len(b2.unnamed_inputs)):
        assert b2.unnamed_inputs[i].dim == b1.unnamed_outputs[i].dim

    out = QAC()
    out.parent = b1.parent

    # set up inputs
    for key in b1.scope_inputs.keys():
        out.scope_inputs[key] = b1.scope_inputs[key]
    for key in b2.scope_inputs.keys():
        if key not in out.scope_inputs:
            out.scope_inputs[key] = b2.scope_inputs[key]
    out.unnamed_inputs += b1.unnamed_inputs

    # wire scope_io up by substitution
    for key in b2.scope_inputs.keys():
        b2.scope_inputs[key].substitute(b1.scope_outputs[key])

    # wire unnamed_io up by substitution
    for i in range(len(b2.unnamed_inputs)):
        b2.unnamed_inputs[i].substitute(b1.unnamed_outputs[i])

    # set up outputs
    for key in b1.scope_outputs.keys():
        if key not in b2.scope_inputs:
            out.scope_outputs[key] = b1.scope_outputs[key]
    for key in b2.scope_outputs.keys():
        out.scope_outputs[key] = b2.scope_outputs[key]

    out.unnamed_outputs += b2.unnamed_outputs

    out.instrs += b1.instrs
    out.instrs += b2.instrs

    out.scale = b1.scale * b2.scale

    return out


def block_exponent(b,n):
    assert n == int(n)
    assert n >= 0

    out = QAC()
    out.parent = b.parent
    b.parent = out

    # populate with the identity.
    for key in b.scope_inputs:
        reg = Register(b.scope_inputs[key].dim)
        reg.name_hints += b.scope_inputs[key].name_hints
        out.scope_inputs[key] = reg
        out.scope_outputs[key] = reg

    for i in range(len(b.unnamed_inputs)):
        reg = Register(b.unnamed_inputs[i].dim)
        reg.name_hints += b.unnamed_inputs[i].name_hints
        out.unnamed_inputs[i] = reg
        out.unnamed_outputs[i] = reg

    out.repeat_fixed(b,n)

    return out


def block_exponent_index(b,parent,idxname):

    out = QAC()
    out.parent = parent
    b.parent = out

    # populate with the identity.
    for key in b.scope_inputs:
        reg = Register(b.scope_inputs[key].dim)
        reg.name_hints += b.scope_inputs[key].name_hints
        out.scope_inputs[key] = reg
        out.scope_outputs[key] = reg

    for i in range(len(b.unnamed_inputs)):
        reg = Register(b.unnamed_inputs[i].dim)
        reg.name_hints += b.unnamed_inputs[i].name_hints
        out.unnamed_inputs[i] = reg
        out.unnamed_outputs[i] = reg

    out.repeat_variable(b,idxname)

    return out


def block_tensor(b2,b1):
    if b1.parent == b2: b1.parent = b2.parent
    if b2.parent == b1: b2.parent = b1.parent
    if b1.parent is None: b1.parent = b2.parent
    if b2.parent is None: b2.parent = b1.parent
    assert b1.parent == b2.parent

     # Require that the consumed or created scoped variables in b1,b2 are disjoint.
    for key in b1.scope_inputs.keys():
        if key in b1.scope_outputs: continue
        assert key not in b2.scope_inputs
        assert key not in b2.scope_outputs
    for key in b2.scope_inputs.keys():
        if key in b2.scope_outputs: continue
        assert key not in b1.scope_inputs
        assert key not in b1.scope_outputs
    for key in b1.scope_outputs.keys():
        if key in b1.scope_inputs: continue
        assert key not in b2.scope_inputs
        assert key not in b2.scope_outputs
    for key in b2.scope_outputs.keys():
        if key in b2.scope_inputs: continue
        assert key not in b1.scope_inputs
        assert key not in b1.scope_outputs

    out = QAC()
    out.parent = b1.parent

    for key in b1.scope_inputs:
        out.scope_inputs[key] = b1.scope_inputs[key]
    for key in b1.scope_outputs:
        out.scope_outputs[key] = b1.scope_outputs[key]
    for key in b2.scope_inputs:
        if key not in out.scope_inputs:
            out.scope_inputs[key] = b2.scope_inputs[key]
        else:
            out.scope_outputs[key].substitute(b2.scope_inputs[key])
    for key in b2.scope_outputs:
        if key not in out.scope_outputs:
            out.scope_outputs[key] = b2.scope_outputs[key]

    out.unnamed_inputs = b2.unnamed_inputs + b1.unnamed_inputs
    out.unnamed_outputs = b2.unnamed_outputs + b1.unnamed_outputs

    out.instrs += b1.instrs
    out.instrs += b2.instrs

    out.scale = b1.scale * b2.scale

    return out


def block_tensorexponent(b,n):
    assert n == int(n)
    assert n >= 0

    if n == 0:
        out = QAC()
        out.parent = b.parent
        return out

    if n == 1:
        return b


    out = block_copy(b)
    for i in range(1,n):
        out = block_tensor(out,block_copy(b))

    return out


def block_adjoint(b):
    b.scope_inputs, b.scope_outputs = b.scope_outputs, b.scope_inputs
    b.unnamed_inputs, b.unnamed_outputs = b.unnamed_outputs, b.unnamed_inputs

    def reverse_instructions(old_instrs):
        new_instrs = []

        for instr in old_instrs[::-1]:
            if instr["kind"] == "qac_declare":
                new_instrs.append({"kind":"qac_zero", "reg":instr["reg"]})

            if instr["kind"] == "qac_zero":
                new_instrs.append({"kind":"qac_declare", "reg":instr["reg"], "dim":instr["reg"].dim})

            if instr["kind"] == "qac_discard":
                new_instrs.append({"kind":"qac_maxmixed", "reg":instr["reg"], "dim":instr["reg"].dim})

                # account for scale
                d = instr["reg"].dim
                b.scale *= np.sqrt(d)

            if instr["kind"] == "qac_maxmixed":
                new_instrs.append({"kind":"qac_discard", "reg":instr["reg"]})

                # account for scale
                d = instr["reg"].dim
                b.scale /= np.sqrt(d)

            if instr["kind"] == "qac_increment":
                if instr["expn"]["kind"] == "negate_expn":
                    new_expn = instr["expn"]["expn"]
                else:
                    new_expn = {"kind":"negate_expn", "expn":instr["expn"]}

                new_instrs.append({"kind":"qac_increment", "reg":instr["reg"], "expn":new_expn})

            if instr["kind"] == "qac_unitary":
                inv = []
                for i in range(len(instr["mat"])):
                    row = []
                    for j in range(len(instr["mat"])):
                        row.append(instr["mat"][j][i].conjugate())
                    inv.append(row)
                new_instrs.append({"kind":"qac_unitary", "reg":instr["reg"], "mat":inv})

            if instr["kind"] == "qac_phase":
                new_instrs.append({"kind":"qac_phase", "value":instr["value"].conjugate()})

            if instr["kind"] == "qac_swap":
                new_instrs.append(instr)

            if instr["kind"] == "qac_rename":
                new_instrs.append({"kind":"qac_rename",
                                   "source":instr["target"],
                                   "target":instr["source"]})

            if instr["kind"] == "qac_if":
                new_instrs.append({
                    "kind":"qac_if",
                    "cond":instr["cond"],
                    "instructions": reverse_instructions(instr["instructions"])
                })

        return new_instrs

    b.instrs = reverse_instructions(b.instrs)
    return b


def block_copy(b):

    register_lookup = {}
    def copy_reg(reg):
        tr = reg.trace()
        if tr not in register_lookup:
            register_lookup[tr] = Register(reg.dim)
            register_lookup[tr].name_hints += reg.name_hints
        return register_lookup[tr]

    out = QAC()
    out.parent = b.parent

    out.scale = b.scale

    out.unnamed_inputs = [copy_reg(reg) for reg in b.unnamed_inputs]
    out.unnamed_outputs = [copy_reg(reg) for reg in b.unnamed_outputs]
    out.scope_inputs = {key:copy_reg(reg) for (key,reg) in b.scope_inputs.items()}
    out.scope_outputs = {key:copy_reg(reg) for (key,reg) in b.scope_outputs.items()}

    def copy_expn(expn):
        assert expn["kind"] != "named_register_expn"

        # {"kind": "register_expn", "register":<reg> }
        if expn["kind"] == "register_expn":
            return {"kind": "register_expn", "register":copy_reg(expn["register"]) }

        # {"kind": "value_expn", "value":5j}
        if expn["kind"] == "value_expn":
            return {"kind": "value_expn", "value": expn["value"] }

        # {"kind": "sum_expn", "terms":[<linexp>] }
        if expn["kind"] == "sum_expn":
            return {"kind": "sum_expn",
                    "terms":[copy_expn(exp) for exp in expn["terms"]]
                    }

        # {"kind": "negate_expn", "expn":<linexp> }
        if expn["kind"] == "negate_expn":
            return {
                "kind": "negate_expn",
                "expn": copy_expn(expn["expn"])
            }

        # {"kind": "product_expn", "terms":[<linexp>] }
        if expn["kind"] == "product_expn":
            return {"kind": "product_expn",
                    "terms":[copy_expn(exp) for exp in expn["terms"]]
                    }

        # {"kind": "division_expn", "dividend":<linexp>, "divisor":5j }
        if expn["kind"] == "division_expn":
            return {
                    "kind": "division_expn",
                    "dividend": copy_expn(expn["dividend"]),
                    "divisor": expn["divisor"],
                    }

        # {"kind": "modulo_expn", "dividend":<linexp>, "divisor":5 }
        if expn["kind"] == "modulo_expn":
            return {
                    "kind": "modulo_expn",
                    "dividend": copy_expn(expn["dividend"]),
                    "divisor": expn["divisor"],
                    }

        # {"kind": "boolean_expn", "terms":[<linexp>, <string>, <linexp>, <string>, ...] }
        if expn["kind"] == "boolean_expn":
            terms = []
            for i in range(len(expn["terms"])):
                if i % 2 == 0: terms.append(copy_expn(expn["terms"][i]))
                else: terms.append(expn["terms"][i])
            return { "kind": "boolean_expn", "terms":terms }

        assert False # unreachable


    def copy_instructions(old_instrs):
        new_instrs = []

        for instr in old_instrs:
            if instr["kind"] == "qac_declare":
                new_instrs.append({"kind":"qac_declare",
                                   "reg":copy_reg(instr["reg"]),
                                   "dim": instr["reg"].dim
                                   })
                continue

            if instr["kind"] == "qac_zero":
                new_instrs.append({"kind":"qac_zero", "reg":copy_reg(instr["reg"])})
                continue

            if instr["kind"] == "qac_discard":
                new_instrs.append({"kind":"qac_discard", "reg":copy_reg(instr["reg"])})
                continue

            if instr["kind"] == "qac_maxmixed":
                new_instrs.append({"kind":"qac_maxmixed",
                                   "reg":copy_reg(instr["reg"]),
                                   "dim": instr["reg"].dim
                                   })
                continue

            if instr["kind"] == "qac_increment":
                new_instrs.append({"kind":"qac_increment",
                                   "reg":copy_reg(instr["reg"]),
                                   "expn":copy_expn(instr["expn"])
                                   })
                continue

            if instr["kind"] == "qac_unitary":
                mat_copy = []
                for i in range(instr["reg"].dim):
                    row = []
                    for j in range(instr["reg"].dim):
                        row.append(instr["mat"][i][j])
                    mat_copy.append(row)

                new_instrs.append({"kind":"qac_unitary",
                                   "reg":copy_reg(instr["reg"]),
                                   "mat":mat_copy
                                   })
                continue

            if instr["kind"] == "qac_phase":
                new_instrs.append({"kind":"qac_phase", "value":instr["value"].conjugate()})
                continue

            if instr["kind"] == "qac_swap":
                new_instrs.append({"kind":"qac_swap",
                                    "reg1":copy_reg(instr["reg1"]),
                                    "reg2":copy_reg(instr["reg2"])
                                    })
                continue

            if instr["kind"] == "qac_rename":
                new_instrs.append({"kind":"qac_rename",
                                    "source":copy_reg(instr["source"]),
                                    "target":copy_reg(instr["target"])
                                    })
                continue

            if instr["kind"] == "qac_if":
                new_instrs.append({"kind":"qac_if",
                                   "cond":copy_reg(instr["cond"]),
                                   "instructions":copy_instructions(instr["instructions"])
                                   })
                continue

            assert False # unreachable

        return new_instrs

    out.instrs = copy_instructions(b.instrs)

    return out


def matrix_to_qac(mat):
    assert isinstance(mat, np.ndarray)
    assert mat.dtype == np.complex128
    assert len(mat.shape) == 2

    u,s,vh = np.linalg.svd(mat)

    qac = QAC()

    l,r = mat.shape
    if l == 1 and r == 1:
        z = complex(mat[0,0])
        qac.scale = abs(z)
        if (z/abs(z) != 1):
            qac.instrs.append({
                "kind":"qac_phase",
                "value": z/abs(z)
            })
        return qac

    if r == 1: # column vector
        qac.scale = s[0]
        u *= vh[0,0]
        out_reg = Register(l)
        qac.instrs.append({"kind":"qac_declare", "reg":out_reg, "dim":l})
        if not np.allclose(u,np.eye(u.shape[0])):
            qac.instrs.append({"kind":"qac_unitary", "reg":out_reg, "mat":u.tolist()})
        qac.unnamed_outputs.append(out_reg)
        return qac

    if l == 1: # row vector
        qac.scale = s[0]
        vh *= u[0,0]
        in_reg = Register(r)
        qac.unnamed_inputs.append(in_reg)

        if not np.allclose(vh,np.eye(vh.shape[0])):
            qac.instrs.append({"kind":"qac_unitary", "reg":in_reg, "mat":vh.tolist()})
        qac.instrs.append({"kind":"qac_zero", "reg":in_reg})
        return qac

    in_reg = Register(r)
    qac.unnamed_inputs.append(in_reg)

    if np.allclose(s,s[0]):
        if l == r: # unitary

            qac.scale = s[0]
            mat /= s[0]

            assert np.allclose((mat @ mat.conj().T), np.eye(l))

            if not np.allclose(mat,np.eye(mat.shape[0])):
                qac.instrs.append({"kind":"qac_unitary", "reg":in_reg, "mat":mat.tolist()})
            qac.unnamed_outputs.append(in_reg)

            return qac

        else:
            # isometry
            qac.scale = s[0]

            # vh to r-dim register
            if not np.allclose(vh,np.eye(vh.shape[0])):
                qac.instrs.append({"kind":"qac_unitary", "reg":in_reg, "mat":vh.tolist()})

             # move from r-dim to l-dim register
            out_reg = Register(l)
            qac.instrs.append({"kind":"qac_declare", "reg":out_reg, "dim":l})
            qac.instrs.append({"kind":"qac_increment", "reg":out_reg, "expn":{"kind":"register_expn", "register":in_reg}})
            qac.instrs.append({"kind":"qac_increment", "reg":in_reg, "expn":{"kind":"negate_expn", "expn":{"kind":"register_expn", "register":out_reg}}})
            qac.instrs.append({"kind":"qac_zero", "reg":in_reg})

            # u on l-dim register
            if not np.allclose(u,np.eye(u.shape[0])):
                qac.instrs.append({"kind":"qac_unitary", "reg":out_reg, "mat":u.tolist()})

            qac.unnamed_outputs.append(out_reg)

            return qac

    qac.scale = max(s)
    s /= qac.scale

    flag = Register(2)
    flag.name_hints.append("block_flag")
    qac.instrs.append({"kind":"qac_declare", "reg":flag, "dim":2})

    cond = Register(2)
    cond.name_hints.append("block_cond")
    qac.instrs.append({"kind":"qac_declare", "reg":cond, "dim":2})

    if not np.allclose(vh,np.eye(vh.shape[0])):
        qac.instrs.append({"kind":"qac_unitary", "reg":in_reg, "mat":vh.tolist()})

    for i in range(len(s)):
        if not np.allclose(s[i],1):
            qac.instrs.append({"kind":"qac_increment", "reg":cond, "expn":{
                            "kind":"boolean_expn", "terms":[
                                {"kind":"value_expn", "value": i}, "==",
                                {"kind":"register_expn", "register":in_reg}]}})

            # reflection
            rot = [[s[i],  np.sqrt(1-s[i]**2)],[ np.sqrt(1-s[i]**2), -s[i]]]
            qac.instrs.append({"kind":"qac_if", "cond":cond, "instructions":[
                {"kind":"qac_unitary", "reg":flag, "mat":rot}
            ]})

            qac.instrs.append({"kind":"qac_increment", "reg":cond, "expn":{
                            "kind":"boolean_expn", "terms":[
                                {"kind":"value_expn", "value": i}, "==",
                                {"kind":"register_expn", "register":in_reg}]}})

    qac.instrs.append({"kind":"qac_zero", "reg":flag})
    qac.instrs.append({"kind":"qac_zero", "reg":cond})

    if l == r:
        # square non-unitary matrix
        if not np.allclose(u,np.eye(u.shape[0])):
            qac.instrs.append({"kind":"qac_unitary", "reg":in_reg, "mat":u.tolist()})

        qac.unnamed_outputs.append(in_reg)
        return qac

    # most general case: non-square non-isometric matrix

    # move from r-dim to l-dim register
    out_reg = Register(l)
    qac.instrs.append({"kind":"qac_declare", "reg":out_reg, "dim":l})
    qac.instrs.append({"kind":"qac_increment", "reg":out_reg, "expn":{"kind":"register_expn", "register":in_reg}})
    qac.instrs.append({"kind":"qac_increment", "reg":in_reg, "expn":{"kind":"negate_expn", "expn":{"kind":"register_expn", "register":out_reg}}})
    qac.instrs.append({"kind":"qac_zero", "reg":in_reg})

    # u on l-dim register
    if not np.allclose(u,np.eye(u.shape[0])):
        qac.instrs.append({"kind":"qac_unitary", "reg":out_reg, "mat":u.tolist()})

    qac.unnamed_outputs.append(out_reg)

    return qac

