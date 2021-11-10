
from .debug import recur_print
from .qac import *
from .tokens import error_at

# Turn an ast_expression into a qac_expression
#    Pre-evaluate as much of the expression as possible.
#    Also, doesn't use {"register_expn"}, instead uses:
# {"kind": "named_register_expn", "name":"x", "key": None, <int>, or "k" }
def evaluate_qac_expn(ast,parent,kind,args,kwargs):
    # assert kind in ["comp", "index", "scalar", "block"]
    if kind == "block":
        # Either returns {"named_register_expn" (maybe with key)}, {"value_expn"} or a QAC object.
        mode = "block"
    elif kind == "scalar":
        # Returns anything that isnt QAC.
        mode = "scalar"
    elif kind == "comp":
        out = evaluate_qac_expn(ast,parent,"scalar",args,kwargs)
        if isinstance(out, QAC) or out["kind"] != "value_expn":
            raise SyntaxError("Expected value to be known at compile time"+error_at(ast["loc"],args))
        return out
    else:
        assert kind == "index"
        out = evaluate_qac_expn(ast,parent,"scalar",args,kwargs)

        if isinstance(out, QAC) or out["kind"] not in ["value_expn", "named_register_expn"]:
            raise SyntaxError("Expected integer or register"+error_at(ast["loc"],args))

        if out["kind"] == "value_expn" and (int(out["value"].real) != out["value"] or out["value"].real < 0):
            raise SyntaxError("Expected value '"+str(out["value"])+"' to be a nonnegative integer"+error_at(ast["loc"],args))

        if out["kind"] == "named_register_expn" and out["key"] is not None:
            raise SyntaxError("Expected register to not be an array"+error_at(ast["loc"],args))

        return out


    # {"kind":"symbol_expression", "identifier": <identifier>, "key": <index-expn>/None, "loc":<loc> }
    if ast["kind"] == "symbol_expression":
        name = ast["identifier"]["value"]
        if name in kwargs:
            if isinstance(kwargs[name],complex):
                return {"kind": "value_expn", "value": kwargs[name]}
            else:
                if mode == "scalar":
                    raise SyntaxError("Expected number for keyword argument '"+name+"' but got a block"+error_at(ast["loc"],args))

                return parent.block_userspace(kwargs[name])

        if name not in parent.scope_outputs:
            parent.promote(name)

        if ast["key"] is None:
            return {"kind": "named_register_expn", "name":name, "key":None}

        key = evaluate_qac_expn(ast["key"],parent,"index",args,kwargs)

        # assert key["kind"] in ["value_expn", "named_register_expn"]

        if key["kind"] == "value_expn":
            v = int(key["value"].real)
            if v < 0 or v >= len(parent.scope_outputs[name]):
                raise IndexError("Array index '"+str(v)+"' out of range"+error_at(ast["loc"],args))

            return {"kind": "named_register_expn", "name":name, "key":v}

        if key["kind"] == "named_register_expn":
            assert key["key"] is None
            return {"kind": "named_register_expn", "name":name, "key":key["name"]}

        assert False # should be unreachable

    # {"kind":"scalar_expression", "value": <value>, "loc":<loc> }
    if ast["kind"] == "scalar_expression":
        return {"kind": "value_expn", "value": ast["value"]}

    # {"kind":"block_expression", "value": <value>, "loc":<loc> }
    if ast["kind"] == "block_expression":
        assert mode == "block"
        return parent.block_userspace(ast["value"])

    # {"kind":"consume_expression", "identifier": <identifier>, "loc":<loc> }
    if ast["kind"] == "consume_expression":
        assert mode == "block"
        name = ast["identifier"]["value"]
        return parent.block_consume(name)

    # {"kind":"create_expression", "expn": <scalar-expn>, "dim":<comp-expn> or None, "loc":<loc> }
    if ast["kind"] == "create_expression":
        assert mode == "block"

        expn = evaluate_qac_expn(ast["expn"],parent,"scalar",args,kwargs)

        if ast["dim"] is None:
            dim = {"kind":"value_expn", "value":complex(2)}
            dim = dim["value"]
        else:
            dim = obtain_value(ast["dim"], "int", args, kwargs)

        # assert dim["kind"] == "value_expn"
        if int(dim.real) != dim:
            raise IndexError("Register size '"+str(dim)+"' is not an integer"+error_at(ast["loc"],args))
        dim = int(dim.real)

        if dim <= 1:
            raise IndexError("Register size '"+str(dim)+"' must be at least 2"+error_at(ast["loc"],args))

        return parent.block_create(expn, dim)

    # {"kind":"adjoint_expression", "expn": <block-expn>, "loc":<loc> }
    if ast["kind"] == "adjoint_expression":
        if mode == "block":
            expn = evaluate_qac_expn(ast["expn"],parent,"block",args,kwargs)

            if isinstance(expn, QAC):
                assert kind == "block"
                return block_adjoint(expn)

            assert isinstance(expn, dict)
            if expn["kind"] == "named_register_expn":
                # it's real, so taking the adjoint does nothing.
                return expn

            if expn["kind"] == "value_expn":
                expn["value"] = expn["value"].conjugate()
                return expn

            assert False # unreachable

        if mode == "scalar":
            expn = evaluate_qac_expn(ast["expn"],parent,"scalar",args,kwargs)

            if expn["kind"] == "value_expn":
                expn["value"] = expn["value"].conjugate()
                return expn

            return {"kind": "adjoint_expn", "expn": expn }

        assert False # unreachable.


    # {"kind":"parenthetical_expression", "expn": <expn>, "loc":<loc> }
    if ast["kind"] == "parenthetical_expression":
        return evaluate_qac_expn(ast["expn"],parent,kind,args,kwargs)


    # {"kind":"negate_expression", "expn": <expn>, "loc":<loc> }
    if ast["kind"] == "negate_expression":
        if mode == "block":
            expn = evaluate_qac_expn(ast["expn"],parent,"block",args,kwargs)

            if isinstance(expn, QAC):
                assert kind == "block"
                ph = parent.block_scalar(-1)
                return block_tensor(expn,ph)

            assert isinstance(expn, dict)
            if expn["kind"] == "named_register_expn":
                return parent.block_cast(expn["name"], expn["key"])

            if expn["kind"] == "value_expn":
                expn["value"] = -expn["value"]
                return expn

            assert False # unreachable

        if mode == "scalar":
            expn = evaluate_qac_expn(ast["expn"],parent,"scalar",args,kwargs)

            if expn["kind"] == "value_expn":
                expn["value"] = -expn["value"]
                return expn

            return {"kind": "negate_expn", "expn": expn }

        assert False # unreachable.


    # {"kind":"boolean_expression", "terms": [<expn>, <symbol_token>, <expn>, ... <expn> ], "loc":<loc>  }
    if ast["kind"] == "boolean_expression":
        terms = []
        for i in range(len(ast["terms"])):
            if i % 2 == 1:
                assert ast["terms"][i]["kind"] in ["==", "!=", "<", ">", ">=", "<="]
                terms.append(ast["terms"][i]["kind"])
            else:
                expn = evaluate_qac_expn(ast["terms"][i],parent,"scalar",args,kwargs)
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
                        return {"kind": "value_expn", "value": complex(0)}
                else:
                    all_values_known = False

        if all_values_known:
            return {"kind": "value_expn", "value": complex(1)}

        if mode == "block":
            return parent.block_boolean_cast({"kind": "boolean_expn", "terms":terms })
            pass

        if mode == "scalar":
            return {"kind": "boolean_expn", "terms":terms }

        assert False # unreachable

    # {"kind":"sum_expression", "terms": [<expn>, <expn>, <expn>], "loc":<loc> }
    if ast["kind"] == "sum_expression":

        # evaluate the terms
        values = []
        for term in ast["terms"]:
            values.append(evaluate_qac_expn(term,parent,mode,args,kwargs))

        # check if all the values are value_expns
        # also add up as many of these as you can now
        scalar = None
        all_values = True
        for value in values:
            if isinstance(value,dict) and value["kind"] == "value_expn":
                if scalar is None: scalar = 0j
                scalar += value["value"]
            else: all_values = False

        if all_values:
            return {"kind": "value_expn", "value": scalar}

        if mode == "block":
            casted_values = []

            if scalar is not None:
                casted_values.append(parent.block_scalar(scalar))

            firstshape = None
            for i in range(len(values)):
                value = values[i]
                if isinstance(value,dict):
                    if firstshape is None: firstshape = "Unit <- Unit"
                    elif firstshape != "Unit <- Unit":
                        raise SyntaxError("Shape of sum is '"+firstshape+"' but got a term of shape 'Unit <- Unit' "+error_at(ast["terms"][i]["loc"],args))

                    if value["kind"] == "value_expn": continue

                    assert value["kind"] == "named_register_expn"
                    casted_values.append(parent.block_cast(value["name"], value["key"]))
                    continue

                assert isinstance(value, QAC)
                thisshape = value.unnamed_shape()
                if firstshape is None: firstshape = thisshape
                elif firstshape != thisshape:
                    raise SyntaxError("Shape of sum is '"+thisshape+"' but got a term of shape 'Unit <- Unit' "+error_at(ast["terms"][i]["loc"],args))

                casted_values.append(value)

            return block_add(*casted_values)

        if mode == "scalar":
            terms = []

            if scalar is not None:
                terms.append({"kind": "value_expn", "value": scalar})

            for value in values:
                assert isinstance(value,dict)
                if value["kind"] == "value_expn": continue
                terms.append(value)

            return {"kind": "sum_expn", "terms":terms }

        assert False # unreachable



    # {"kind":"tensorproduct_expression", "terms": [<expn>, <expn>, <expn>], "loc":<loc> }
    if ast["kind"] == "tensorproduct_expression":
        assert mode == "block"

        scalar = 1
        terms = []
        for term in ast["terms"]:
            value = evaluate_qac_expn(term,parent,"block",args,kwargs)
            if isinstance(value,QAC):
                terms.append(value)
                continue

            assert isinstance(value,dict)
            if value["kind"] == "value_expn":
                scalar *= value["value"]
            else:
                assert value["kind"] == "named_register_expn"
                terms.append(parent.block_cast(value["name"], value["key"]))

        if len(terms) == 0:
            return {"kind": "value_expn", "value": scalar}

        if scalar != 1:
            terms.append(parent.block_scalar(scalar))

        out = terms[0]
        for i in range(1,len(terms)):
            out = block_tensor(out,terms[i])


        return out

    # {"kind":"product_expression", "terms": [<expn>, <expn>, <expn>], "loc":<loc> }
    if ast["kind"] == "product_expression":
        # evaluate the terms
        values = []
        for term in ast["terms"]:
            values.append(evaluate_qac_expn(term,parent,mode,args,kwargs))

        # check if all the values are value_expns
        # also add up as many of these as you can now
        scalar = None
        all_values = True
        for value in values:
            if isinstance(value,dict) and value["kind"] == "value_expn":
                if scalar is None: scalar = 1+0j
                scalar *= value["value"]
            else: all_values = False

        if all_values:
            return {"kind": "value_expn", "value": scalar}

        if mode == "block":
            out = None
            for i in range(len(values)):
                value = values[i]

                if isinstance(value,dict):
                    if value["kind"] == "value_expn":
                        term = parent.block_scalar(scalar)
                    else:
                        assert value["kind"] == "named_register_expn"
                        term = parent.block_cast(value["name"], value["key"])
                else:
                    term = value

                if out is None: out = term
                else:
                    if out.unnamed_rshape() != term.unnamed_lshape():
                        raise SyntaxError("Right shape '"+out.unnamed_rshape()+"' does not align with left shape '"+term.unnamed_lshape()+"' of term "+error_at(ast["terms"][i]["loc"],args))
                    out = block_mul(out, term)

            return out

        if mode == "scalar":

            terms = []

            if scalar is not None:
                terms.append({"kind": "value_expn", "value": scalar})

            for value in values:
                assert isinstance(value,dict)
                if value["kind"] == "value_expn": continue
                terms.append(value)

            return {"kind": "product_expn", "terms":terms }

        assert False # unreachable


    # {"kind":"division_expression", "dividend": <expn>, "divisor": <comp-expn>, "loc":<loc> }
    if ast["kind"] == "division_expression":

        divisor = evaluate_qac_expn(ast["divisor"],parent,"comp",args,kwargs)
        assert divisor["kind"] == "value_expn"
        if divisor["value"] == 0:
            raise ZeroDivisionError("Divisor "+str(divisor["value"])+" must be nonzero"+error_at(ast["loc"],args))

        if mode == "block":
            dividend = evaluate_qac_expn(ast["dividend"],parent,"block",args,kwargs)


            if isinstance(dividend,dict):
                if dividend["kind"] == "value_expn":
                    return {"kind": "value_expn", "value": dividend["value"]/divisor["value"]}
                else:
                    assert dividend["kind"] == "named_register_expn"
                    dividend = parent.block_cast(dividend["name"], dividend["key"])

            assert isinstance(dividend, QAC)

            return block_tensor(parent.block_scalar(1/divisor["value"]),dividend)

        if mode == "scalar":
            dividend = evaluate_qac_expn(ast["dividend"],parent,"scalar",args,kwargs)

            if dividend["kind"] == "value_expn":
                return {"kind": "value_expn", "value": dividend["value"]/divisor["value"]}

            return {"kind":"division_expn", "dividend": dividend, "divisor":divisor["value"]}

        assert False # unreachable



    # {"kind":"modulo_expression", "dividend": <expn>, "divisor": <comp-expn>, "loc":<loc> }
    if ast["kind"] == "modulo_expression":
        def complex_mod_int(complex_nr,integer):
            # Honestly, I don't know what generalization of modulo to complex numbers I should pick.
            # There are some sensible candidates but none of them are obvious or standard.
            # But I need a modulo operation, so I'm going with this simple thing for now - should be replaced later:
            # Insist the divisor is a positive integer, and shift the real part of the dividend into the range [0,divisor).
            # According to the wikipedia article on modulo, languages vary wildly in their implementation of this operation.
            # Sadly the python implementation won't work for me because I need it to be well defined for any complex dividend.
            out = complex_nr
            while out.real < 0: out += integer
            while out.real >= integer: out -= integer
            return out

        assert mode == "scalar"

        dividend = evaluate_qac_expn(ast["dividend"],parent,"scalar",args,kwargs)
        divisor = evaluate_qac_expn(ast["divisor"],parent,"comp",args,kwargs)

        assert divisor["kind"] == "value_expn"
        if int(divisor["value"].real) != divisor["value"] or divisor["value"] == 0:
            raise ValueError("Divisor "+str(divisor["value"])+" must be a positive nonzero integer"+error_at(ast["loc"],args))
        v = int(divisor["value"].real)

        if dividend["kind"] == "value_expn":
            return {"kind": "value_expn", "value": complex_mod_int(dividend["value"],v)}

        return {"kind":"modulo_expn", "dividend": dividend, "divisor":v}


    # {"kind":"exponent_expression", "base": <expn>, "exponent": <expn>/<index-expn>, "loc":<loc> }
    if ast["kind"] == "exponent_expression":

        if mode == "block":
            base = evaluate_qac_expn(ast["base"],parent,"block",args,kwargs)

            if isinstance(base, dict) and base["kind"] == "value_expn":
                exponent = evaluate_qac_expn(ast["exponent"],parent,"comp",args,kwargs)
                assert exponent["kind"] == "value_expn"
                return {"kind": "value_expn", "value": base["value"] ** exponent["value"]}

            exponent = evaluate_qac_expn(ast["exponent"],parent,"index",args,kwargs)

            if isinstance(base, dict):
                assert base["kind"] == "named_register_expn"
                base = parent.block_cast(base["name"], base["key"])

            if exponent["kind"] == "value_expn":
                return block_exponent(base,int(exponent["value"].real))

            if exponent["kind"] == "named_register_expn":
                assert exponent["key"] is None

                return block_exponent_index(base,exponent["name"])

            assert False # unreachable


        if mode == "scalar":
            base = evaluate_qac_expn(ast["base"],parent,"scalar",args,kwargs)
            exponent = evaluate_qac_expn(ast["exponent"],parent,"scalar",args,kwargs)

            if base["kind"] == "value_expn" and exponent["kind"] == "value_expn":
                return {"kind": "value_expn", "value": base["value"] ** exponent["value"]}

            return {"kind":"exponent_expn", "base": base, "divisor": exponent}

        assert False # unreachable


    # {"kind":"tensorexponent_expression", "base": <expn>, "exponent": <comp-expn>, "loc":<loc> }
    if ast["kind"] == "tensorexponent_expression":
        assert mode == "block"

        base = evaluate_qac_expn(ast["base"],parent,"block",args,kwargs)
        exponent = evaluate_qac_expn(ast["exponent"],parent,"comp",args,kwargs)

        assert exponent["kind"] == "value_expn"
        if exponent["value"] != int(exponent["value"].real) or int(exponent["value"].real) < 0:
            raise ValueError("Divisor "+str(exponent["value"])+" must be a non-negative integer"+error_at(ast["loc"],args))

        if isinstance(base,dict):
            if base["kind"] == "value_expn":
                return {"kind": "value_expn", "value": base["value"] ** int(exponent["value"].real) }
            else:
                assert base["kind"] == "named_register_expn"
                base = parent.block_cast(base["name"], base["key"])

        return block_tensorexponent(base, int(exponent["value"].real))

    assert False # should be unreachable


# extracts the value from an ast that is either a literal, identifier, or extern
def obtain_value(ast, check, args, kwargs):
    assert check in ["int","complex","blq"]

    if ast["kind"] == "literal":
        assert isinstance(ast["value"],str)

        if check == "blq":
            raise SyntaxError("Expected literal '"+str(ast["value"])+"' to be a block."+error_at(ast["loc"],args))

        if check == "int":
            try:
                return int(ast["value"])
            except:
                raise SyntaxError("Expected literal '"+str(ast["value"])+"' to be an integer."+error_at(ast["loc"],args))

        if check == "complex":
            try:
                return complex(ast["value"])
            except:
                raise SyntaxError("Expected literal '"+str(ast["value"])+"' to be a number."+error_at(ast["loc"],args))

        return dim


    if ast["kind"] == "identifier":
        key = ast["value"]
        if key not in kwargs:
            raise SyntaxError("No symbol '"+key+"' provided in keyword arguments."+error_at(ast["loc"],args))

        value = kwargs[key]
        if isinstance(value,complex):
            if check == "int":
                out = int(value.real)
                if out != value:
                    raise SyntaxError("Expected keyword extern '"+str(value)+"' to be an integer."+error_at(ast["loc"],args))
                return out
            if check == "complex": return value
            if check == "blq":
                raise SyntaxError("Expected keyword extern '"+str(value)+"' to be a block."+error_at(ast["loc"],args))
        else:
            if check == "int":
                raise SyntaxError("Expected keyword extern '"+str(value)+"' to be an integer."+error_at(ast["loc"],args))
            if check == "complex":
                raise SyntaxError("Expected keyword extern '"+str(value)+"' to be a number."+error_at(ast["loc"],args))
            if check == "blq":
                return value
        assert False # should be unreachable


    if ast["kind"] == "extern":
        if isinstance(ast["value"],complex):
            if check == "int":
                out = int(ast["value"].real)
                if out != ast["value"]:
                    raise SyntaxError("Expected extern '"+str(ast["value"])+"' to be an integer"+error_at(ast["loc"],args))
                return out
            if check == "complex": return ast["value"]
            if check == "blq":
                raise SyntaxError("Expected extern '"+str(ast["value"])+"' to be a block"+error_at(ast["loc"],args))
        else:
            if check == "int":
                raise SyntaxError("Expected extern '"+str(ast["value"])+"' to be an integer"+error_at(ast["loc"],args))
            if check == "complex":
                raise SyntaxError("Expected extern '"+str(ast["value"])+"' to be a number"+error_at(ast["loc"],args))
            if check == "blq":
                return ast["value"]
        assert False # should be unreachable

    assert False # should be unreachable



def build_block(in_decls,out_decls,instrs,args,kwargs):

    out = QAC()
    is_in_decls = True
    for decls in [in_decls, out_decls]:

        for ident in decls:
            name = ident["identifier"]["value"]

            if is_in_decls:
                assert name not in out.scope_inputs
            else:
               if name in out.scope_inputs: continue

            dim = obtain_value(ident["dim"],"int", args, kwargs)

            if ident["slots"] is None:
                out.scope_inputs[name] = Register(dim)
                out.scope_inputs[name].name_hints.append(name)
                out.scope_outputs[name] = out.scope_inputs[name]
            else:
                slots = obtain_value(ident["slots"],"int", args, kwargs)

                out.scope_inputs[name] = []
                out.scope_outputs[name] = []
                for i in range(slots):
                    r = Register(dim)
                    r.name_hints.append(name)
                    out.scope_inputs[name].append(r)
                    out.scope_outputs[name].append(r)

        is_in_decls = False

    def apply_instruction(instr,qac):
        # { "kind": "declare_instruction",  "identifier": <identifier>, "dim": <comp-expn>, "slots": <None or comp-expn>, "loc":<loc>}
        if instr["kind"] == "declare_instruction":
            name = instr["identifier"]["value"]
            dim = obtain_value(instr["dim"],"int", args, kwargs)
            if instr["slots"] is None:
                qac.declare(name, dim, slots=None)
            else:
                slots = obtain_value(instr["slots"],"int", args, kwargs)
                qac.declare(name, dim, slots=slots)
            return

        # { "kind": "discard_instruction",  "identifiers": [<identifier>], "loc":<loc> }
        if instr["kind"] == "discard_instruction":
            for ident in instr["identifiers"]:
                name = ident["value"]
                qac.discard(name)
            return

        # { "kind": "uncompute_instruction",  "identifiers": [<identifier>], "loc":<loc> }
        if instr["kind"] == "uncompute_instruction":
            raise SyntaxError("The uncompute instruction is not part of the current language. It may be added at some point. Invoked "+error_at(ast["loc"],args))
            return

        # { "kind": "scalar_instruction",  "expn": <block-expn>, "loc":<loc> }
        if instr["kind"] == "scalar_instruction":
            expn = evaluate_qac_expn(instr["expn"],qac,"block",args,kwargs)

            # convert to qac
            if isinstance(expn,dict):
                if expn["kind"] == "value_expn":
                    expn = qac.block_scalar(expn["value"])
                else:
                    assert expn["kind"] == "named_register_expn"
                    expn = qac.block_cast(expn["name"],expn["key"])

            assert isinstance(expn, QAC)

            expn_rshape = expn.unnamed_rshape()
            if expn_rshape != "Unit":
                raise SyntaxError("Block in scalar expression should be a scalar, but has right shape '"+expn_rshape+"'"+error_at(instr["loc"],args))

            expn_lshape = expn.unnamed_lshape()
            if expn_lshape != "Unit":
                raise SyntaxError("Block in scalar expression should be a scalar, but has left shape '"+expn_lshape+"'"+error_at(instr["loc"],args))

            qac.scalar_instr(expn)

            return



        # { "kind": "pass_instruction",  "loc":<loc> }
        if instr["kind"] == "pass_instruction":
            # Do nothing.
            return

        # { "kind": "repeat_instruction",  "expn": <index-expn>, "instructions",[], "loc":<loc> }
        if instr["kind"] == "repeat_instruction":

            child_qac = QAC()
            child_qac.parent = qac
            for child_instr in instr["instructions"]:
                apply_instruction(child_instr,child_qac)

            expn = evaluate_qac_expn(instr["expn"],qac,"index",args,kwargs)

            if expn["kind"] == "value_expn":
                count = int(expn["value"].real)
                assert count == expn["value"] and count >= 0
                qac.repeat_fixed(child_qac, count)
                return

            if expn["kind"] == "named_register_expn":
                assert expn["key"] is None
                qac.repeat_variable(child_qac, expn["name"])
                return

            assert False # unreachable


        # { "kind": "if_instruction",  "conditions": [{"condition":<scalar-expn>, "instructions",[]}], "else": [], "loc":<loc> }
        if instr["kind"] == "if_instruction":

            assert len(instr["conditions"]) > 0

            controls = []
            for cond in instr["conditions"]:
                expn = evaluate_qac_expn(cond["condition"],qac,"scalar",args,kwargs)

                child_qac = QAC()
                child_qac.parent = qac
                for child_instr in cond["instructions"]:
                    apply_instruction(child_instr,child_qac)

                controls.append({
                    "expn":expn,
                    "qac":child_qac
                })

            else_qac = None
            if len(instr["else"]) > 0:
                child_qac = QAC()
                child_qac.parent = qac
                for child_instr in instr["else"]:
                    apply_instruction(child_instr,child_qac)

            qac.compound_if_statement(controls,else_qac=else_qac)

            return


        # { "kind": "init_instruction", "targets":[{"target":<identifier>, "key":<index-expn> or None}], "expn":<block-expn>, "loc":<loc> }
        if instr["kind"] == "init_instruction":

            expn = evaluate_qac_expn(instr["expn"],qac,"block",args,kwargs)

            # convert to qac
            if isinstance(expn,dict):
                if expn["kind"] == "value_expn":
                    expn = qac.block_scalar(expn["value"])
                else:
                    assert expn["kind"] == "named_register_expn"
                    expn = qac.block_cast(expn["name"],expn["key"])

            targets = []
            target_shape = []
            for target in instr["targets"]:
                name = target["target"]["value"]
                if name not in qac.scope_outputs:
                    qac.promote(name)

                if target["key"] is None:
                    key = None
                else:
                    key = evaluate_qac_expn(target["key"],qac,"index",args,kwargs)

                    if key["kind"] == "value_expn":
                        idx = int(key["value"].real)
                        if idx < 0 or idx >= len(qac.scope_outputs[name]):
                            raise IndexError("Array index '"+str(key_int)+"' out of range"+error_at(key["loc"],args))
                        key = idx
                    else:
                        assert key["kind"] == "named_register_expn"
                        assert key["key"] is None
                        key = key["name"]
                        if name not in qac.scope_outputs:
                            qac.promote(key)

                targets.append({"name":name, "key":key})

                reg = qac.scope_outputs[name]
                if isinstance(reg,list):
                    if key is None:
                        for r in reg: target_shape.append(r.dim)
                    else:
                        target_shape.append(reg[0].dim)
                else:
                    target_shape.append(reg.dim)

            expn_rshape = expn.unnamed_rshape()
            if expn_rshape != "Unit":
                raise SyntaxError("Block in init expression should be a vector, but has right shape '"+expn_rshape+"'"+error_at(instr["loc"],args))

            target_shape = " @ ".join([str(x) for x in target_shape])
            expn_lshape = expn.unnamed_lshape()
            if target_shape != expn_lshape:
                raise SyntaxError("Shape of targets '"+target_shape+"' does not match left shape '"+expn_lshape+"' of block in init expression"+error_at(instr["loc"],args))

            qac.init_instr(targets,expn)
            return

        # { "kind": "assign_instruction", "target":<identifier>, "key":<index-expn> or None, "expn":<scalar-expn>, "undo":<scalar-expn> or None, "loc":<loc> }
        if instr["kind"] == "assign_instruction":

            expn = evaluate_qac_expn(instr["expn"],qac,"scalar",args,kwargs)

            name = instr["target"]["value"]
            if name not in qac.scope_outputs:
                qac.promote(name)

            assert name in qac.scope_outputs
            assert (instr["key"] is None) or (instr["undo"] is None)

            if instr["undo"] is not None:
                undo = evaluate_qac_expn(instr["undo"],qac,"scalar",args,kwargs)
                qac.assign(name,expn,undo=undo)
                return

            if instr["key"] is not None:
                key = evaluate_qac_expn(instr["key"],qac,"index",args,kwargs)

                if key["kind"] == "value_expn":
                    key_int = int(key["value"].real)

                    if key_int < 0 or key_int >= len(qac.scope_outputs[name]):
                        raise IndexError("Array index '"+str(key_int)+"' out of range"+error_at(key["loc"],args))

                    qac.assign_array_fixed(name,key_int,expn)
                    return

                if key["kind"] == "named_register_expn":
                    assert key["key"] is None
                    qac.assign_array_variable(name,key["name"],expn)
                    return

                assert False # unreachable

            qac.assign(name,expn,undo=None)

            return


        # { "kind": "increment/decrement_instruction", "target":<identifier>, "key":<index-expn> or None, "expn":<scalar-expn>, "loc":<loc> }
        if instr["kind"] in ["increment_instruction", "decrement_instruction"]:

            name = instr["target"]["value"]
            if name not in qac.scope_outputs:
                qac.promote(name)

            expn = evaluate_qac_expn(instr["expn"],qac,"scalar",args,kwargs)
            if instr["kind"] == "decrement_instruction":
                expn = {"kind": "negate_expn", "expn":expn }

            if instr["key"] is None:
                qac.increment(name,expn)
                return

            if instr["key"] is not None:
                key = evaluate_qac_expn(instr["key"],qac,"index",args,kwargs)

                if key["kind"] == "value_expn":
                    key_int = int(key["value"].real)

                    if key_int < 0 or key_int >= len(qac.scope_outputs[name]):
                        raise IndexError("Array index '"+str(key_int)+"' out of range"+error_at(key["loc"],args))

                    qac.increment_array_fixed(name,key_int,expn)
                    return

                if key["kind"] == "named_register_expn":
                    assert key["key"] is None
                    qac.increment_array_variable(name,key["name"],expn)
                    return

            assert False # unreachable

        assert False # unreachable


    for instr in instrs:
        apply_instruction(instr,out)

    assert len(out.unnamed_inputs) == 0
    assert len(out.unnamed_outputs) == 0

    # populate out.unnamed_inputs with the correct out.scope_inputs
    finished_input_decls = []
    for ident in in_decls:
        name = ident["identifier"]["value"]
        assert name in out.scope_inputs
        reg = out.scope_inputs[name]

        # check if dimension is right, then append to out.unnamed_inputs
        dim = obtain_value(ident["dim"],"int", args, kwargs)
        if ident["slots"] is None:
            assert isinstance(reg, Register)
            assert reg.dim == dim

            out.unnamed_inputs.append(reg)
        else:
            slots = obtain_value(ident["slots"],"int", args, kwargs)
            assert isinstance(reg, list)
            assert len(reg) == slots
            for i in range(slots):
                assert reg[i].dim == dim
                out.unnamed_inputs.append(reg[i])

        finished_input_decls.append(name)


    # declare all the output declarations that aren't inputs.
    declared_output_decls = []
    for ident in out_decls:
        name = ident["identifier"]["value"]
        assert name in out.scope_inputs
        reg = out.scope_inputs[name]

        # check if dimension is right
        dim = obtain_value(ident["dim"],"int", args, kwargs)
        if ident["slots"] is None:
            assert isinstance(reg, Register)
            assert reg.dim == dim
        else:
            slots = obtain_value(ident["slots"],"int", args, kwargs)
            assert isinstance(reg, list)
            assert len(reg) == slots
            for i in range(slots):
                assert reg[i].dim == dim

        # Was this already an input? If so, we are done.
        if name in finished_input_decls: continue

        # prepend some declaration statements
        if ident["slots"] is None:
            out.instrs.insert(0, {"kind":"qac_declare", "reg":reg, "dim":dim})
        else:
            for i in reversed(range(slots)):
                out.instrs.insert(0, {"kind":"qac_declare", "reg":reg[i], "dim":dim})

        declared_output_decls.append(name)



    # check that I got all of them
    for name in out.scope_inputs:
        assert name in finished_input_decls or name in declared_output_decls

    # put the output declarations into out.unnamed_outputs
    extracted_output_decls = []
    for ident in out_decls:
        name = ident["identifier"]["value"]

        assert name in out.scope_outputs
        reg = out.scope_outputs[name]

        extracted_output_decls.append(name)

        # check if dimension is right, append to unnamed_outputs
        dim = obtain_value(ident["dim"],"int", args, kwargs)
        if ident["slots"] is None:
            assert isinstance(reg, Register)
            assert reg.dim == dim
            out.unnamed_outputs.append(reg)
        else:
            slots = obtain_value(ident["slots"],"int", args, kwargs)
            assert isinstance(reg, list)
            assert len(reg) == slots
            for i in range(slots):
                assert reg[i].dim == dim
                out.unnamed_outputs.append(reg[i])

    # check that I got all of them
    for name in out.scope_outputs:
        assert name in extracted_output_decls

    # reset these, they are now all in opt.unnamed_in/outputs
    out.scope_inputs = {}
    out.scope_outputs = {}

    return out
