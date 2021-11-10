
from .tokens import error_at, merge_locations

from .debug import recur_print

def make_ast(tokens,args):
    # kept up to date by match
    # is an array so match can edit it
    lastloc = [tokens[0]["loc"]]

    # match one of several tokens by kind,
    #   remove that token from the tokens list,
    #   return the token that was matched
    # if None is one of the arguments and none of the token arguments match, returns None and leaves list intact
    # if None is not an argument and nothing matches this raises a SyntaxError
    def match(*targets):
        if len(tokens) == 0:
            if None not in targets:
                raise SyntaxError("Unexpected end of input, expected "+repr(targets))
            else:
                return None

        token = tokens[0]
        if token["kind"] not in targets:
            if None in targets: return None
            raise SyntaxError("Unexpected "+repr(token["kind"])+", expected "+repr(targets) + error_at(token["loc"],args))

        lastloc[0] = tokens[0]["loc"]
        return tokens.pop(0)

    ###### block declaration

    def parse_declaration():
        token = match("identifier")
        match(":")
        decl = {
            "identifier": token,
            "dim": match("extern","identifier","literal"),
            "slots": None,
        }
        decl["loc"] = merge_locations(token["loc"], decl["dim"]["loc"])
        if match("@@",None) is not None:
            decl["slots"] = match("extern","identifier","literal")
            decl["loc"] = merge_locations(token["loc"], decl["slots"]["loc"])
        return decl


    match("[")

    out_decls = []
    if match("<-",None) is None:
        while True:
            out_decls.append(parse_declaration())
            if match("<-","@")["kind"] == "<-":
                break

    in_decls = []
    if match("]",None) is None:
        while True:
            in_decls.append(parse_declaration())
            if match("]","@")["kind"] == "]":
                break


    # allow user to put a new line after the block decl, or not, up to them.
    match("\n",None)

    ######

    def parse_instruction_block():
        match(":")
        if match("\n",None) is None:
            instrs = [parse_instruction()]
        else:
            match("indent")
            instrs = []
            while True:
                if match("dedent",None) is not None:
                    break
                instrs.append(parse_instruction())
                match("\n")
        return instrs

    # doesn't eat up any new lines
    def parse_instruction():
        t = match("declare","discard", "uncompute", "scalar", "pass", "repeat", "if", "identifier")

        #    { "kind": "declare_instruction",  "identifier": <identifier>, "dim": <int>, "slots": <None or int>, "loc":<loc>}
        #    declare <identifier>:<regtype>
        if t["kind"] == "declare":
            decl = parse_declaration()
            return {
                "kind": "declare_instruction",
                "identifier": decl["identifier"],
                "dim": decl["dim"],
                "slots": decl["slots"],
                "loc": merge_locations(t["loc"], decl["loc"])
            }

        #    { "kind": "discard_instruction",  "identifiers": [<identifier>], "loc":<loc> }
        #    discard <identifier>,<identifier>
        if t["kind"] == "discard":
            idents = [match("identifier")]
            while match(",",None) is not None:
                idents.append(match("identifier"))

            return {
                "kind": "discard_instruction",
                "identifiers": idents,
                "loc": merge_locations(t["loc"], idents[-1]["loc"])
            }

        #    { "kind": "uncompute_instruction",  "identifier": <identifier>, "loc":<loc> }
        #    uncompute <identifier>
        if t["kind"] == "uncompute":
            idents = [match("identifier")]
            while match(",",None) is not None:
                idents.append(match("identifier"))

            return {
                "kind": "uncompute_instruction",
                "identifiers": idents,
                "loc": merge_locations(t["loc"], idents[-1]["loc"])
            }

        #    { "kind": "scalar_instruction",  "expn": <expn>, "loc":<loc> }
        #    scalar <expn>
        if t["kind"] == "scalar":
            expn = parse_expression()
            return {
                "kind": "scalar_instruction",
                "expn": expn,
                "loc": merge_locations(t["loc"], expn["loc"])
            }

        #    { "kind": "pass_instruction", "loc":<loc> }
        #    pass
        if t["kind"] == "pass":
            return {
                "kind": "pass_instruction",
                "loc": t["loc"]
            }


        # <expn> must be known at compile time
        #    { "kind": "repeat_instruction",  "expn": <expn>, "instructions",[], "loc":<loc> }
        #    repeat <expn>: <instrs>
        if t["kind"] == "repeat":
            expn = parse_expression()
            instrs = parse_instruction_block()
            return {
                "kind": "repeat_instruction",
                "expn": expn,
                "instructions": instrs,
                "loc": merge_locations(t["loc"], lastloc[0])
            }

        #    { "kind": "if_instruction",  "conditions": [{"condition":<expn>, "instructions",[]}], "else": [], "loc":<loc> }
        #    if <expn>: <instrs>
        #    elif <expn>:
        #    else: <instrs>
        if t["kind"] == "if":
            expn = parse_expression()
            instrs = parse_instruction_block()

            conditions = [{
                "condition": expn,
                "instructions": instrs,
            }]
            else_instrs = []

            while True:
                if not (len(tokens) >= 2 and tokens[0]["kind"] == "\n" and tokens[1]["kind"] in ["elif","else"]):
                    break
                match("\n")
                next_token = match("elif","else")

                if next_token["kind"] == "elif":
                    expn = parse_expression()
                    instrs = parse_instruction_block()
                    conditions.append({
                        "condition": expn,
                        "instructions": instrs
                    })
                else:
                    else_instrs = parse_instruction_block()

            return {
                "kind": "if_instruction",
                "conditions": conditions,
                "else": else_instrs,
                "loc": merge_locations(t["loc"], lastloc[0])
            }

        # "key" <expn> must be known at compile time or just be an identifier
        #    { "kind": "init_instruction", "targets":[{"target":<identifier>, "key":<expn> or None}], "expn":<expn>, "loc":<loc> }
        #    <ident> @ <ident>[<expn>] <- <expn>

        # "key" <expn> must be known at compile time or just be an identifier
        #    { "kind": "assign_instruction", "target":<identifier>, "key":<expn> or None, "expn":<expn>, "undo":<expn> or None, "loc":<loc> }
        #    <ident>[<expn>] = <expn> undo <expn>
        #    <ident>[<expn>] = <expn>
        #    <ident> = <expn> undo <expn>
        #    <ident> = <expn>
        if t["kind"] == "identifier":
            next_token = match("@","<-","[","+=","-=","=")

            key = None
            if next_token["kind"] == "[":
                key = parse_expression()
                match("]")
                next_token = match("@","<-","+=","-=","=")

            if next_token["kind"] == "<-":
                expn = parse_expression()
                return {
                    "kind": "init_instruction",
                    "targets": [{"target":t, "key":key}],
                    "expn": expn,
                    "loc": merge_locations(t["loc"], expn["loc"])
                }

            if next_token["kind"] == "@":
                targets = [{"target":t, "key":key}]

                while True:
                    next_target = {
                            "target": match("identifier"),
                            "key": None,
                        }

                    next_token = match("[", "@", "<-")

                    if next_token["kind"] == "[":
                        next_target["key"] = parse_expression()
                        match("]")
                        next_token = match("@","<-")

                    targets.append(next_target)

                    if next_token["kind"] == "<-": break

                expn = parse_expression()
                return {
                    "kind": "init_instruction",
                    "targets": targets,
                    "expn": expn,
                    "loc": merge_locations(t["loc"], expn["loc"])
                }

            expn = parse_expression()
            loc = merge_locations(t["loc"], expn["loc"])

            if next_token["kind"] == "+=":
                return {
                    "kind": "increment_instruction",
                    "target": t,
                    "expn": expn,
                    "key": key,
                    "loc": loc,
                }

            if next_token["kind"] == "-=":
                return {
                    "kind": "decrement_instruction",
                    "target": t,
                    "expn": expn,
                    "key": key,
                    "loc": loc,
                }

            undo = None

            if match("undo",None) != None:
                undo = parse_expression()
                loc = merge_locations(t["loc"], undo["loc"])

            return {
                "kind": "assign_instruction",
                "target": t,
                "expn": expn,
                "key": key,
                "undo": undo,
                "loc": loc,
            }


    def parse_expression():
        t = match("identifier", "literal", "extern", "|", 'ket(', "~", "(", "+", "-")

        # {"kind":"symbol_expression", "identifier": <identifier>, "key": <expn>/None, "loc":<loc> }
        #   <identifier>
        #   <identifier>[<expn>]
        if t["kind"] == "identifier":
            key = None
            arguments = []

            if match("[", None) is not None:
                key = parse_expression()
                match("]")

            out = {
                "kind": "symbol_expression",
                "identifier": t,
                "key": key,
                "loc": merge_locations(t["loc"], lastloc[0])
            }

        # {"kind":"scalar_expression", "value": <literal>, "loc":<loc> }
        if t["kind"] == "literal":
            out = {
                "kind": "scalar_expression",
                "value": complex(t["value"]),
                "loc": t["loc"]
            }

        # {"kind":"scalar_expression", "value": <literal>, "loc":<loc> }
        # {"kind":"block_expression", "value": <value>, "loc":<loc> }
        if t["kind"] == "extern":
            if isinstance(t["value"],complex):
                out = {
                    "kind": "scalar_expression",
                    "value": t["value"],
                    "loc": t["loc"]
                }
            else:
                out = {
                    "kind": "block_expression",
                    "value": t["value"],
                    "loc": t["loc"]
                }

        # {"kind":"consume_expression", "identifier": <identifier>, "loc":<loc> }
        #   |<ident>
        if t["kind"] == "|":
            out = {
                "kind": "consume_expression",
                "identifier": match("identifier"),
                "loc": merge_locations(t["loc"], lastloc[0])
            }

        # {"kind":"create_expression", "expn": <expn>, "dim":<expn> or None, "loc":<loc> }
        #   ket(<expn>, <expn>) or just ket(<expn>)
        if t["kind"] == "ket(":
            out = {
                "kind": "create_expression",
                "expn": parse_expression(),
                "dim": None,
            }
            if match(",", None) is not None:
                out["dim"] = match("extern","identifier","literal")
            match(')')
            out["loc"] = merge_locations(t["loc"], lastloc[0])

        # {"kind":"parenthetical_expression", "expn": <expn>, "loc":<loc> }
        #   (<expn>)
        if t["kind"] == "(":
            out = {
                "kind": "parenthetical_expression",
                "expn": parse_expression(),
            }
            match(")")
            out["loc"] = merge_locations(t["loc"], lastloc[0])



        # {"kind":"adjoint_expression", "expn": <expn>, "loc":<loc> }
        #   ~<expn>  adjoint
        # if t["kind"] == "~":
        #    out = {
        #        "kind": "adjoint_expression",
        #        "expn": parse_expression(),
        #        "loc": merge_locations(t["loc"], lastloc[0])
        #    }

        # removed
        #   +<expn>
        # if t["kind"] == "+":
        #    out = parse_expression()
        #    out["loc"] = merge_locations(t["loc"], out["loc"])

        # {"kind":"negate_expression", "expn": <expn>, "loc":<loc> }
        #   -<expn>
        # if t["kind"] == "-":
        #    out = {
        #        "kind": "negate_expression",
        #        "expn": parse_expression(),
        #    }
        #    out["loc"] = merge_locations(t["loc"], out["expn"]["loc"])

        # unary operators generate arithmetic expressions
        if t["kind"] in ["+","-","~"]:
            expn = parse_expression()
            if expn["kind"] == "arithmetic_expression":
                out = {
                    "kind": "arithmetic_expression",
                    "terms": [t] + expn["terms"],
                    "loc": merge_locations(t["loc"], expn["loc"])
                }
            else:
                out = {
                    "kind": "arithmetic_expression",
                    "terms": [t, expn["terms"]],
                    "loc": merge_locations(t["loc"], expn["loc"])
                }

        # {"kind":"arithmetic_expression", "terms": [mix of <expn> and <symbol_token>], "loc":<loc> }
        next_token = match("+", "-", "*", "**", "/", "%", "@", "@@",
                "==", "!=", ">", "<", ">=", "<=", None)

        if next_token is None: return out

        next_expn = parse_expression()

        if out["kind"] == "arithmetic_expression": begin = out["terms"]
        else: begin = [out]

        if next_expn["kind"] == "arithmetic_expression": end = next_expn["terms"]
        else: end = [next_expn]

        return {
            "kind": "arithmetic_expression",
            "terms": begin + [next_token] + end,
            "loc": merge_locations(out["loc"], next_expn["loc"])
        }

    ##############

    instrs = []
    while True:
        # check end of file
        instr = parse_instruction()
        instrs.append(instr)
        if match("\n", None) is None: break


    #######

    def process_arithmetic(ast):
        # handle anything that isnt an arithmetic expression
        if isinstance(ast,list):
            for i in range(len(ast)):
                ast[i] = process_arithmetic(ast[i])
            return ast

        if not isinstance(ast,dict):
            return ast

        if "kind" not in ast or ast["kind"] != "arithmetic_expression":
            for key in ast.keys():
                ast[key] = process_arithmetic(ast[key])
            return ast

        # handle arithmetic expressions
        assert ast["kind"] == "arithmetic_expression"

        # handle unary operators
        terms = []
        prv = None
        expn = True # alternating variable
        for term in ast["terms"]:
            if expn:
                if term["kind"] in ["+", "-", "~"]:
                    assert prv is None
                    prv = term
                else:
                    if prv is None or prv["kind"] == "+":
                        terms.append(term)
                    else:
                        term = {
                            "kind": "tmp",
                            "expn": term,
                            "loc": merge_locations(prv["loc"], term["loc"])
                        }

                        if prv["kind"] == "-":
                            term["kind"] = "negate_expression"
                        else:
                            assert prv["kind"] == "~"
                            term["kind"] = "adjoint_expression"

                        prv = None
                        terms.append(term)
                    expn = False
            else:
                terms.append(term)
                expn = True

        ast["terms"] = terms

        # recurse_first
        for i in range(len(ast["terms"])):
            if i % 2 == 1: continue
            ast["terms"][i] = process_arithmetic(ast["terms"][i])

        # regular associative: + - * / % @ **
        # boolean associative: == != > < >= <=
        # not associative: @@
        # EMDAS: @@ * / +-
        # x / y * z = x * (1/y) * z
        # boolean expression x != y < 3 == z can be left as is, highest precedence

        # replace - <expn> with + (negate <expn>)
        for i in range(len(ast["terms"])):
            if i % 2 == 0: continue
            if ast["terms"][i]["kind"] == "-":
                ast["terms"][i]["kind"] = "+"
                ast["terms"][i+1] = {
                    "kind": "negate_expression",
                    "expn": ast["terms"][i+1],
                    "loc": merge_locations(ast["terms"][i]["loc"], ast["terms"][i+1]["loc"])
                }

        # Need to be explicit about / and %
        #    A / B * C throws SyntaxError, need either (A / B) * C or A / (B * C), sim for %
        #    A / B + C and A / B @ C are allowed tho


        # Need to also be explicit about tensor exponent vs regular exponent

        # Determine lowest precedence of this expression
        # precedence order "exponent" "product (mod, div)" "tensor" "sum" "boolean"
        kinds = []
        kinds_lookup = {0:"boolean", 1:"sum", 2:"tensor", 3:"product", 4:"exponent"}
        for i in range(len(ast["terms"])):
            if i % 2 == 0: continue
            if ast["terms"][i]["kind"] == "+": kinds.append(1)
            elif ast["terms"][i]["kind"] == "@": kinds.append(2)
            elif ast["terms"][i]["kind"] in ["*", "%", "/"]: kinds.append(3)
            elif ast["terms"][i]["kind"] in ["**", "@@"]: kinds.append(4)
            else: kinds.append(0)

        kind = kinds_lookup[min(kinds)]

        def chop_terms(separators):

            terms = []
            tmp_terms = []

            for i in range(len(ast["terms"])):
                if i % 2 == 0:
                    tmp_terms.append(ast["terms"][i])
                else:
                    if ast["terms"][i]["kind"] not in separators:
                        tmp_terms.append(ast["terms"][i])
                    else:
                        if len(tmp_terms) == 1:
                            terms.append(tmp_terms[0])
                        else:
                            terms.append(process_arithmetic({
                                "kind": "arithmetic_expression",
                                "terms": tmp_terms,
                                "loc": merge_locations(tmp_terms[0]["loc"], tmp_terms[-1]["loc"])
                            }))
                        terms.append(ast["terms"][i])
                        tmp_terms = []

            if len(tmp_terms) == 1:
                terms.append(tmp_terms[0])
            else:
                terms.append(process_arithmetic({
                    "kind": "arithmetic_expression",
                    "terms": tmp_terms,
                    "loc": merge_locations(tmp_terms[0]["loc"], tmp_terms[-1]["loc"])
                }))

            return terms


        if kind == "boolean":
            return {
                "kind": "boolean_expression",
                "terms": chop_terms(["==","!=",">","<",">=","<="]),
                "loc": ast["loc"]
            }

        if kind == "sum":
            raw_terms = chop_terms(["+"])
            terms = []
            for i in range(len(raw_terms)):
                if i%2 == 0: terms.append(raw_terms[i])

            return {
                "kind": "sum_expression",
                "terms": terms,
                "loc": ast["loc"]
            }

        if kind == "tensor":
            raw_terms = chop_terms(["@"])
            terms = []
            for i in range(len(raw_terms)):
                if i%2 == 0: terms.append(raw_terms[i])

            return {
                "kind": "tensorproduct_expression",
                "terms": terms,
                "loc": ast["loc"]
            }

        if kind == "product":
            raw_terms = chop_terms(["*", "/", "%"])
            terms = []

            divide = None
            modulo = None

            for i in range(len(raw_terms)):
                if raw_terms[i]["kind"] in ["/", "%"]:
                    if i == len(raw_terms)-2:
                        if raw_terms[i]["kind"] == "/": divide = raw_terms[i+1]
                        if raw_terms[i]["kind"] == "%": modulo = raw_terms[i+1]
                        break
                    else:
                        raise SyntaxError("Ambiguous arithmetic expression at "+ast["loc"]+".\n"+\
                                          "Need to be explicit about parentheses with modulo and division: A / (B * C)  or  (A / B) * C.\n"+\
                                          "Can only be last operations of expression: A * B / C" + error_at(token["loc"],args))
                if i%2 == 0: terms.append(raw_terms[i])

            out = {
                "kind": "product_expression",
                "terms": terms,
                "loc": merge_locations(terms[0]["loc"], terms[-1]["loc"]),
            }

            if divide is not None:
                return {
                    "kind": "division_expression",
                    "dividend": out,
                    "divisor": divide,
                    "loc": ast["loc"]
                }
            if modulo is not None:
                return {
                    "kind": "modulo_expression",
                    "dividend": out,
                    "divisor": modulo,
                    "loc": ast["loc"]
                }
            return out

        if kind == "exponent":
            terms = chop_terms(["**", "@@"])
            if len(terms) > 3:
                raise SyntaxError("Ambiguous arithmetic expression at "+ast["loc"]+".\n"+\
                                  "Need to be explicit about parentheses with exponent ** and tensor exponent @@" + error_at(token["loc"],args))

            if terms[1]["kind"] == "**":
                return {
                    "kind": "exponent_expression",
                    "base": terms[0], "exponent": terms[2],
                    "loc": ast["loc"]
                }
            if terms[1]["kind"] == "@@":
                return {
                    "kind": "tensorexponent_expression",
                    "base": terms[0], "exponent": terms[2],
                    "loc": ast["loc"]
                }

        raise ValueError("Couldn't process arithmetic expression. (This should be unreachable)")



    ######

    instrs = process_arithmetic(instrs)

    return in_decls, out_decls, instrs
