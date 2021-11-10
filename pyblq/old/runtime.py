
from .qaasm import *


# returns qaasm_expn. Evaluates as much as possible.
# perhaps need to wrap the thing to keep track of dependent registers
#   and also deal with arrays

# alternative: realize that a qaasm expn can only ever end up in a qaasm increment instruction.
# the wrapper could be a list of instructions itself.
#  {"expn":<qaasm_expn>, "depends":[<reg>,<reg>,<reg>],
#                        "arrays":[{"target":<reg>, "key":<reg>, "regs":[<reg>,<reg>]}] }
# returns expn, depends, arrays
def process_qaasm_expn(ast):

    if ast["kind"] == "symbol_expression":
        name = ast["identifier"]["value"]

        if name in kwargs:
            assert complex(kwargs[name]) == kwargs[name]
            return {"kind": "value_expn", "value": kwargs[name]}, [], []

        assert name in scope

        if ast["key"] is None:
            assert isinstance(scope[name], Register)
            return {"kind": "register_expn", "register":scope[name]}, [scope[name]], []

        key, deps, arrs = process_qaasm_expn(ast["key"])
        assert len(arrs) == 0

        if key["kind"] == "value_expn":
            assert len(deps) == 0
            v = int(key["value"].real)
            if v < 0 or v >= len(scope[name]):
                raise IndexError("Array index '"+str(v)+"' out of range"+at(ast["loc"]))

            return {"kind": "register_expn", "register":scope[name][v]}, [scope[name][v]], []

        assert key["kind"] == "register_expn"

        assert len(deps) == 1
        assert deps[0] == key["register"]

        keydim = key["register"].dim
        if keydim >= len(scope[name]):
            keyname = ast["key"]["identifier"]["value"]]
            raise IndexError("Array index register '"+keyname+"' dimension "+str(keydim)+" out of range"+\
                    " for array '"+key+"' of length "+str(len(scope[name]))+error_at(ast["loc"],args))

        dim = scope[name][0].dim
        reg = Register(dim)
        out = {"kind":"register_expn", "register":reg}
        array = {"target": reg, "key": key["register"], "regs":[scope[name][i] for i in range(keydim)]}

        return out, deps+array["regs"], [array]

    if ast["kind"] == "scalar_expression":
        return {"kind": "value_expn", "value": ast["value"]}, [], []

    assert ast["kind"] != "block_expression"
    assert ast["kind"] != "consume_expression"
    assert ast["kind"] != "create_expression"
    assert ast["kind"] != "adjoint_expression"

    if ast["kind"] == "parenthetical_expression":
        return process_qaasm_expn(ast["expn"])

    if ast["kind"] == "negate_expression":
        child, deps, arrs = process_qaasm_expn(ast["expn"])
        if child["kind"] == "value_expn":
            return {"kind": "value_expn", "value": -child["value"]}, [], []

        return {"kind": "negate_expn", "expn": child }, deps, arrs

    if ast["kind"] == "boolean_expression":
        # {"kind": "boolean_expn", "terms":[<linexp>, <string>, <linexp>, <string>, ...] }

        terms = []
        out_deps = []
        out_arrs = []
        for i in range(len(ast["terms"])):
            if i % 2 == 1:
                assert ast["terms"][i] in ["==", "!=", "<", ">", ">=", "<="]
                terms.append(ast["terms"][i])
            else:
                expn, deps, arrs = process_qaasm_expn(ast["terms"][i])
                out_arrs += arrs
                for dep in deps:
                    if dep not in out_deps:
                        out_deps.append(dep)
                terms.append(expn)

        # try to pre-evaluate entire expression
        all_values_known = True
        for i in range(len(ast["terms"])):
            if i % 2 == 1:
                if terms[i-1]["kind"] == "value_expn" and terms[i+1]["kind"] == "value_expn":
                    if terms[i] == "==": value = (terms[i-1]["value"] == terms[i+1]["value"])
                    elif terms[i] == "!=": value = (terms[i-1]["value"] != terms[i+1]["value"])
                    elif terms[i] == "<": value = (terms[i-1]["value"] < terms[i+1]["value"])
                    elif terms[i] == ">": value = (terms[i-1]["value"] > terms[i+1]["value"])
                    elif terms[i] == ">=": value = (terms[i-1]["value"] >= terms[i+1]["value"])
                    else:
                        assert terms[i] == "<="
                        value = (terms[i-1]["value"] <= terms[i+1]["value"])

                    if not value:
                        return {"kind": "value_expn", "value": complex(0)}, [], []
                else:
                    all_values_known = False

        if all_values_known:
            return {"kind": "value_expn", "value": complex(1)}, [], []

        return {"kind": "boolean_expn", "terms":terms }, out_deps, out_arrs

    if ast["kind"] == "sum_expression":

        scalar = None
        terms = []
        out_deps = []
        out_arrs = []
        for term in ast["terms"]:
            value, deps, arrs = process_qaasm_expn(term)
            out_arrs += arrs
            for dep in deps:
                if dep not in out_deps:
                    out_deps.append(dep)
            if value["kind"] == "value_expn":
                if scalar is None:
                    scalar = value["value"]
                else:
                    scalar += value["value"]
            else:
                terms.append(value)

        if len(terms) == 0:
            return {"kind": "value_expn", "value": scalar}, [], []

        if scalar is not None:
            terms.append({"kind": "value_expn", "value": scalar})

        return {"kind": "sum_expn", "terms":terms }, out_deps, out_arrs

    assert ast["kind"] != "tensorproduct_expression"

    if ast["kind"] == "product_expression":
        scalar = None
        terms = []
        out_deps = []
        out_arrs = []
        for term in ast["terms"]:
            value,deps,arrs = process_qaasm_expn(term)
            out_arrs += arrs
            for dep in deps:
                if dep not in out_deps:
                    out_deps.append(dep)
            if value["kind"] == "value_expn":
                if scalar is None:
                    scalar = value["value"]
                else:
                    scalar *= value["value"]
            else:
                terms.append(value)

        if len(terms) == 0:
            return {"kind": "value_expn", "value": scalar}, [], []

        if scalar is not None:
            terms.append({"kind": "value_expn", "value": scalar})

        return {"kind": "product_expn", "terms":terms }, out_deps, out_arrs

    if ast["kind"] == "division_expression":
        dividend, deps, arrs = process_qaasm_expn(ast["dividend"])
        divisor, _, _ = process_qaasm_expn(ast["divisor"])

        assert divisor["kind"] == "value_expn"

        if dividend["kind"] == "value_expn":
            return {"kind": "value_expn", "value": dividend["value"]/divisor["value"]}, [], []

        return {"kind":"division_expn", "dividend": dividend, "divisor":divisor["value"]}, deps, arrs

    if ast["kind"] == "modulo_expression":
        dividend, deps, arrs = process_qaasm_expn(ast["dividend"])
        divisor, _, _ = process_qaasm_expn(ast["divisor"])

        assert divisor["kind"] == "value_expn"
        v = divisor["value"]
        bad = (v.real != v)
        if not bad: bad = int(v.real) != v
        if not bad: bad = int(v.real) < 1

        if v.real != v:
            raise IndexError("Modulo divisor dimension "+str(v)+" must be a positive integer"+error_at(ast["loc"],args))

        v = int(v.real)

        if dividend["kind"] == "value_expn":
            # Honestly, I don't know what generalization of modulo to complex numbers I should pick.
            # There are some sensible candidates but none of them are obvious or standard.
            # But I need a modulo operation, so I'm going with this simple thing for now - should be replaced later:
            # Insist the divisor is a positive integer, and shift the real part of the dividend into the range [0,divisor).
            # According to the wikipedia article on modulo, languages vary wildly in their implementation of this operation.
            # Sadly the python implementation won't work for me because I need it to be well defined for any complex dividend.
            out = dividend["value"]
            while out.real < 0: out += v
            while out.real >= v.real: out -= v
            return {"kind": "value_expn", "value": out}, [], []

        return {"kind":"modulo_expn", "dividend": dividend, "divisor":v}

    if ast["kind"] == "exponent_expression":
        base, out_deps, out_arrs = process_qaasm_expn(ast["base"])
        exponent, deps, arrs = process_qaasm_expn(ast["exponent"])

        out_arrs += arrs
        for dep in deps:
            if dep not in out_deps:
                out_deps.append(dep)

        if base["kind"] == "value_expn" and exponent["kind"] == "value_expn":
            return {"kind": "value_expn", "value": base["value"] ** exponent["value"]}, [], []

        return {"kind":"exponent_expn", "base": base, "divisor": exponent}, out_deps, out_arrs

    assert ast["kind"] != "tensorexponent_expression"

    assert False # should be unreachable



########################################################################


def build_block(out,in_decls,out_decls,instrs,args,kwargs):

    scope = {}

    for decl in in_decls:
        # TODO: correct?
        scope[name] = Register(in_decls["dim"])

    #################
    # returns Blq
    def process_block_expn(ast,scope):
        if ast["kind"] == "symbol_expression":

            if name in kwargs:
                assert isinstance(kwargs[name], blq)
                # don't need to copy, since all other evaluations copy
                return kwargs[name]

            assert name in scope
            expn, deps, arrs = process_qaasm_expn(ast)

            assert expn["kind"] == "register_expn"
            assert len(arrs) in [0,1]

            out = Blq()
            out.scale = expn["register"].dim

            if len(arrs) == 1:
                arr = qaasm_array_idx(b,arrs[0])
                arr.__enter__()

            for i in range(expn["register"].dim):
                tmp = Blq()

            qaasm_postselect

            if len(arrs) == 1:
                arr.__exit__()

            return

        if ast["kind"] == "scalar_expression":
            return

        if ast["kind"] == "block_expression":
            return

        if ast["kind"] == "consume_expression":
            return

        if ast["kind"] == "create_expression":
            return

        if ast["kind"] == "adjoint_expression":
            return

        if ast["kind"] == "parenthetical_expression":
            return

        if ast["kind"] == "negate_expression":
            return

        if ast["kind"] == "boolean_expression":
            return

        if ast["kind"] == "sum_expression":
            return

        if ast["kind"] == "tensorproduct_expression":
            return

        if ast["kind"] == "product_expression":
            return

        if ast["kind"] == "division_expression":
            return

        if ast["kind"] == "modulo_expression":
            return

        if ast["kind"] == "exponent_expression":
            return

        if ast["kind"] == "tensorexponent_expression":
            return

        assert False # should be unreachable



    def process_instruction(ast):
        if ast["kind"] == "declare_instruction":
            return

        if ast["kind"] == "discard_instruction":
            return

        if ast["kind"] == "uncompute_instruction":
            return

        if ast["kind"] == "scalar_instruction":
            return

        if ast["kind"] == "pass_instruction":
            return

        if ast["kind"] == "repeat_instruction":
            return

        if ast["kind"] == "if_instruction":
            return

        if ast["kind"] == "init_instruction":
            return

        if ast["kind"] == "assign_instruction":
            return

        if ast["kind"] == "increment_instruction":
            return

        if ast["kind"] == "decrement_instruction":
            return

        assert False # should be unreachable

    #################

    for instr in instrs:
        process_instruction(instr)

    # package up outQaasm somehow?

    return out
