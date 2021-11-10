

def qac_to_string(qac):

    out = ""

    if qac.scale != 1:
        out += "# Scale: "+str(qac.scale)+"\n"

    scope = {} # dict of registers currently in scope
    # values are None or a string for that register


    if len(qac.scope_inputs.keys()) != 0:
        for key in qac.scope_inputs.keys():
            reg = qac.scope_inputs[key]
            if isinstance(reg,list):
                for i in range(len(reg)):
                    assert reg[i].trace() not in scope
                    scope[reg[i].trace()] = None
                    out += "scope_input "+reg_to_string(reg[i],scope)+" : "+str(reg[i].dim)
                    out += " # "+key+"["+str(i)+"]\n"
            else:
                assert reg.trace() not in scope
                scope[reg.trace()] = None
                out += "scope_input "+reg_to_string(reg,scope)+" : "+str(reg.dim)
                out += " # "+key+"\n"
        out += "\n"

    if len(qac.unnamed_inputs) != 0:
        for reg in qac.unnamed_inputs:
            assert reg.trace() not in scope
            scope[reg.trace()] = None
            out += "input "+reg_to_string(reg,scope)+" : "+str(reg.dim)+"\n"
        out += "\n"

    for instr in qac.instrs:
        out += instr_to_string(instr,scope)

    out += "\n"

    if len(qac.scope_outputs.keys()) != 0:
        for key in qac.scope_outputs.keys():
            reg = qac.scope_outputs[key]
            if isinstance(reg,list):
                for i in range(len(reg)):
                    out += "output "+reg_to_string(reg[i],scope)+" : "+str(reg[i].dim)
                    out += " # "+key+"["+str(i)+"]\n"
                    del scope[reg[i].trace()]
            else:
                out += "scope_output "+reg_to_string(reg,scope)+" : "+str(reg.dim)
                out += " # "+key+"\n"
                del scope[reg.trace()]
        out += "\n"

    if len(qac.unnamed_outputs) != 0:
        for reg in qac.unnamed_outputs:
            out += "output "+reg_to_string(reg,scope)+" : "+str(reg.dim)+"\n"
            del scope[reg.trace()]

    return out

def instr_to_string(instr,scope):
    def test_print(d):
        if hasattr(d, "name_hints"):
            return (d.trace(), d.name_hints)
        if isinstance(d,list):
            return [test_print(val) for val in d]
        if isinstance(d,dict):
            return {key:test_print(val) for key,val in d.items()}
        return d
    # print(scope)
    # print(test_print(instr))

    if instr["kind"] == "qac_declare":
        #declare <identifier> : <positive integer literal>
        # {"kind":"qac_declare", "reg":<register>, "dim":<int>}
        assert instr["reg"].trace() not in scope
        scope[instr["reg"].trace()] = None
        return "declare "+reg_to_string(instr["reg"],scope)+" : "+str(instr["dim"])+"\n"

    if instr["kind"] == "qac_discard":
        # discard <identifier>
        # {"kind":"qac_discard", "reg":<register>}
        out = "discard "+reg_to_string(instr["reg"],scope)+"\n"
        del scope[instr["reg"].trace()]
        return out

    if instr["kind"] == "qac_maxmixed":
        # maxmixed <identifier> : <positive integer literal>
        # {"kind":"qac_maxmixed", "reg":<register>, "dim":<int>}
        assert instr["reg"].trace() not in scope
        scope[instr["reg"].trace()] = None
        return "maxmixed "+reg_to_string(instr["reg"],scope)+"\n"

    if instr["kind"] == "qac_zero":
        # zero <identifier>
        # {"kind":"qac_zero", "reg":<register>}
        out = "zero "+reg_to_string(instr["reg"],scope)+"\n"
        del scope[instr["reg"].trace()]
        return out

    if instr["kind"] == "qac_increment":
        # <identifier> += <expn>
        # {"kind":"qac_increment", "reg":<register>, "expn":<expn>}
        return reg_to_string(instr["reg"],scope)+" += "+expn_to_string(instr["expn"],scope)+"\n"

    if instr["kind"] == "qac_unitary":
        # unitary [[0,1],[1,0]] <identifier>
        # {"kind":"qac_unitary", "reg":<register>, "mat":<matrix>}
        return "unitary "+str(instr["mat"])+" "+reg_to_string(instr["reg"],scope)+"\n"

    if instr["kind"] == "qac_phase":
        # phase <complex lit with mag 1>
        # {"kind":"qac_phase", "value":<complexnr>}
        return "phase "+str(instr["value"])+"\n"

    if instr["kind"] == "qac_swap":
        # swap <reg1> <reg2>
        # {"kind":"qac_swap", "reg1":<register>, "reg2": <register>}
        return "swap "+reg_to_string(instr["reg1"],scope)+" "+reg_to_string(instr["reg2"],scope)+"\n"

    if instr["kind"] == "qac_rename":
        # rename <source> <target>
        # {"kind":"qac_rename", "source":<register>, "target": <register>}
        assert instr["source"].trace() in scope
        assert instr["target"].trace() not in scope
        scope[instr["target"].trace()] = None
        out = "rename "+reg_to_string(instr["source"],scope)+" "+reg_to_string(instr["target"],scope)+"\n"
        del scope[instr["source"].trace()]
        return out

    if instr["kind"] == "qac_if":
        # if <identifier>: <instr \n> or pass
        # {"kind":"qac_if", "cond":<register>, "instructions":[<instrs>] }
        out = "if "+reg_to_string(instr["cond"],scope)+":\n"
        for sub_instr in instr["instructions"]:
            ss = instr_to_string(sub_instr,scope)
            for s in ss.split("\n"):
                if s == "": continue
                out += "    "+s+"\n"
        if len(instr["instructions"]) == 0:
            out += "    pass\n"
        return out

    assert False # unreachable


def expn_to_string(expn,scope):

    if expn["kind"] == "register_expn":
        # <identifier>        # {"kind": "register_expn", "register":<reg> }
        return reg_to_string(expn["register"],scope)

    if expn["kind"] == "value_expn":
        # <complex literal>   # {"kind": "value_expn", "value":5j}
        if expn["value"].imag == 0: return str(expn["value"].real)
        if expn["value"].real == 0: return str(expn["value"])
        return "("+str(expn["value"])+")"

    if expn["kind"] == "sum_expn":
        # <expn> + <expn>     # {"kind": "sum_expn", "terms":[<linexp>] }
        return "(" + " + ".join([expn_to_string(x,scope) for x in expn["terms"]])+")"

    if expn["kind"] == "negate_expn":
        # -<expn>             # {"kind": "negate_expn", "expn":<linexp> }
        return "-"+expn_to_string(expn["expn"],scope)

    if expn["kind"] == "adjoint_expn":
        # ~<expn>             # {"kind": "adjoint_expn", "expn":<linexp> }
        return "~"+expn_to_string(expn["expn"],scope)

    if expn["kind"] == "product_expn":
        # <expn> * <expn>     # {"kind": "product_expn", "terms":[<linexp>] }
        return "(" + " * ".join([expn_to_string(x,scope) for x in expn["terms"]])+")"

    if expn["kind"] == "division_expn":
        # <expn> / <expn>     # {"kind": "division_expn", "dividend":<linexp>, "divisor":5j }
        return "("+expn_to_string(expn["dividend"],scope) +"/("+str(expn["divisor"])+"))"

    if expn["kind"] == "modulo_expn":
        # <expn> % <expn>     # {"kind": "modulo_expn", "dividend":<linexp>, "divisor":5 }
        return "("+expn_to_string(expn["dividend"],scope) +"%("+str(expn["divisor"])+"))"

    if expn["kind"] == "boolean_expn":
        # {"kind": "boolean_expn", "terms":[<linexp>, <string>, <linexp>, <string>, ...] }
        # <string> is one of ==, !=, >, <, >=, <=
        out = "("
        for i in range(len(expn["terms"])):
            if i % 2 == 0:
                out += " " + expn_to_string(expn["terms"][i],scope) + " "
            else:
                out += expn["terms"][i]
        return out + ")"

    if expn["kind"] == "named_register_expn":
        return str(expn)

    assert False # unreachable


def reg_to_string(reg,scope):
    if reg.trace() not in scope:
        print("Warning: register with trace " + str(reg.trace()) + " is undeclared.")
        return "?"
    if scope[reg.trace()] is not None: return scope[reg.trace()]

    # for debugging
    if True and len(reg.name_hints) > 1:
        combined_hint = "_".join(reg.name_hints)
        if combined_hint not in scope.values():
            scope[reg.trace()] = combined_hint
            return combined_hint

    for hint in reg.name_hints:
        if hint in scope.values(): continue
        scope[reg.trace()] = hint
        return hint

    if len(reg.name_hints) == 0: base = "reg"
    else: base = reg.name_hints[0]

    i = 1
    while base+str(i) in scope.values(): i+=1
    scope[reg.trace()] = base+str(i)
    return base+str(i)
