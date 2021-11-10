

from .tokens import error_at

def scope_analysis(in_decls,out_decls,instrs,args,kwargs):

    scope = {}
    # { "x": {"loc" <loc of identifier that defined it>, "array": False} }

    controls = {}
    # { "x": <loc> }

    # enforces restrictions on expressions
    # returns rough_type, consume_symbs, linear_symbs
    # rough_type is either "scalar" or "block", arrays never make it through.
    # consume_symbs and linear symbs are lists of identifiers
    def analyze_expn(ast):

        if ast["kind"] == "symbol_expression":
            name = ast["identifier"]["value"]

            if name not in scope:
                if name not in kwargs:
                    raise SyntaxError("Undeclared identifier '"+name+"'"+error_at(ast["identifier"]["loc"],args))

                if isinstance(kwargs[name],complex):
                    return "scalar", [], []
                return "block", [], []

            if scope[name]["array"] and ast["key"] is None:
                raise SyntaxError("Symbol '"+name+"' is an array and is missing a key"+error_at(ast["loc"],args))

            if not scope[name]["array"] and ast["key"] is not None:
                raise SyntaxError("Symbol '"+name+"' is not an array but has a key"+error_at(ast["loc"],args))

            keySymb = []
            if ast["key"] is not None:
                keySymb = list_index_check(ast["key"])

            return "scalar", [], [ast["identifier"]]+keySymb


        if ast["kind"] == "scalar_expression":
            return "scalar", [], []

        if ast["kind"] == "block_expression":
            return "block", [], []

        if ast["kind"] == "consume_expression":
            name = ast["identifier"]["value"]

            if name not in scope:
                if name not in kwargs:
                    raise SyntaxError("Undeclared identifier '"+name+"'"+error_at(ast["identifier"]["loc"],args))

                raise SyntaxError("Can't consume keyword argument '"+name+"'"+error_at(ast["identifier"]["loc"],args))

            return "block", [ast["identifier"]], []


        if ast["kind"] == "create_expression":
            t, c, l = analyze_expn(ast["expn"])

            if t != "scalar":
                raise SyntaxError("Ket expected a scalar but got a block"+error_at(ast["expn"]["loc"],args))

            if len(c) > 0:
                raise SyntaxError("Can't consume '"+c[0]["value"]+"' in a create expression"+error_at(ast["loc"],args)\
                                    +"referenced"+error_at(c[0]["loc"]))

            if ast["dim"] is not None: decl_check(ast["dim"])

            return "block", [], l

        if ast["kind"] == "adjoint_expression":
            t, c, l = analyze_expn(ast["expn"])

            # Wrong: this can totally be a scalar, just take the complex conjugate.
            # if t != "block":
            #    raise SyntaxError("Adjoint expression expected a block but got a scalar"+error_at(ast["expn"]["loc"],args))

            if len(c) > 0:
                raise SyntaxError("Can't consume '"+c[0]["value"]+"' in an adjoint expression"+error_at(ast["loc"],args)\
                                    +"referenced"+error_at(c[0]["loc"],args))

            return "block", [], l

        if ast["kind"] == "parenthetical_expression":
            return analyze_expn(ast["expn"])

        if ast["kind"] == "negate_expression":
            return analyze_expn(ast["expn"])

        if ast["kind"] == "boolean_expression":
            l_out = []
            for i in range(len(ast["terms"])):
                if i % 2 == 0:
                    t, c, l = analyze_expn(ast["terms"][i])
                    if t != "scalar":
                        raise SyntaxError("Boolean expression expected a scalar but got a block"+error_at(ast["terms"][i]["loc"],args))
                    assert len(c) == 0
                    l_out += l

            return "scalar", [], l_out

        if ast["kind"] == "sum_expression":
            l_out = []
            t_out = "scalar"
            for term in ast["terms"]:
                t, c, l = analyze_expn(term)
                if t == "block": t_out = "block"
                if len(c) > 0:
                    raise SyntaxError("Can't consume '"+c[0]["value"]+"' in a sum expression"+error_at(ast["loc"],args)\
                                    +"referenced"+error_at(c[0]["loc"]))
                l_out += l

            return t_out, [], l_out

        if ast["kind"] == "tensorproduct_expression":
            c_out = []
            l_out = []
            for term in ast["terms"]:
                t, c, l = analyze_expn(term)
                c_out += c
                l_out += l

            return "block", c_out, l_out

        if ast["kind"] == "product_expression":
            l_out = []
            c_out = []
            t_out = "scalar"
            for term in ast["terms"]:
                t, c, l = analyze_expn(term)
                if t == "block": t_out = "block"
                c_out += c
                l_out += l

            return t_out, c_out, l_out

        if ast["kind"] == "division_expression":
            t1, c1, l1 = analyze_expn(ast["dividend"])
            t2, c2, l2 = analyze_expn(ast["divisor"])

            if t2 != "scalar":
                raise SyntaxError("Divisor should be a scalar but is a block"+error_at(ast["divisor"]["loc"],args))
            assert len(c2) == 0

            if len(l2) > 0:
                raise SyntaxError("Divisor"+error_at(ast["divisor"]["loc"],args)+\
                                    "must be known at compile time but references '"+l2[0]["value"]+"'"+error_at(l2[0]["loc"],args))

            return t1, c1, l1

        if ast["kind"] == "modulo_expression":
            t1, c1, l1 = analyze_expn(ast["dividend"])
            t2, c2, l2 = analyze_expn(ast["divisor"])

            if t1 != "scalar":
                raise SyntaxError("Dividend should be a scalar but is a block"+error_at(ast["dividend"]["loc"],args))

            if t2 != "scalar":
                raise SyntaxError("Divisor should be a scalar but is a block"+error_at(ast["divisor"]["loc"],args))

            assert len(c1) == 0
            assert len(c2) == 0

            if len(l2) > 0:
                raise SyntaxError("Divisor"+error_at(ast["divisor"]["loc"],args)+\
                                    "must be known at compile time but references '"+l2[0]["value"]+"'"+error_at(l2[0]["loc"],args))

            return "scalar", [], l1

        if ast["kind"] == "exponent_expression":
            t1, c1, l1 = analyze_expn(ast["base"])

            if t1 == "block":
                keySymb = list_index_check(ast["exponent"])

                return "block", c1, l1+keySymb
            else:
                assert t1 == "scalar"

                t2, c2, l2 = analyze_expn(ast["exponent"])
                if t2 != "scalar":
                    raise SyntaxError("Exponent should be a scalar but is a block"+error_at(ast["exponent"]["loc"],args))
                assert len(c1) == 0
                assert len(c2) == 0

                return "scalar", [], l1+l2


        if ast["kind"] == "tensorexponent_expression":
            t1, c1, l1 = analyze_expn(ast["base"])

            if t1 != "block":
                raise SyntaxError("Tensor base should be a block but is a scalar"+error_at(ast["base"]["loc"],args))

            keySymb = list_index_check(ast["exponent"])

            return "block", c1, l1+keySymb

        assert False # should be unreachable

    ###############################

    # a scalar.
    # known at compile time, or exactly an identifier
    # returns [] or [<identifier>]
    def list_index_check(ast):
        t, c, l = analyze_expn(ast)

        if t != "scalar":
            raise SyntaxError("Index or count should be a scalar but is a block"+error_at(ast["loc"],args))
        assert len(c) == 0

        if len(l) == 0: return []

        if ast["kind"] == "symbol_expression":
            if ast["key"] is not None:
                raise SyntaxError("Index or count cannot be an array element"+error_at(ast["loc"],args))
            return [ast["identifier"]]

        raise SyntaxError("Index or count must either be known at compile time or be a single symbol expression"+error_at(ast["loc"],args))

    # ast should evaluate to a dimension
    def decl_check(ast):
        if ast["kind"] == "extern":
            if not isinstance(ast["value"], complex):
                raise SyntaxError("Dimension should be a scalar but is a block"+error_at(ast["loc"],args))
            return

        if ast["kind"] == "identifier":
            name = ast["value"]
            if name in scope:
                raise SyntaxError("Can't use symbol '"+name+"'"+error_at(ast["identifier"]["loc"],args)+\
                                    "because value must be known at compile time")

            if name not in kwargs:
                raise SyntaxError("Unknown keyword argument '"+name+"'"+error_at(ast["identifier"]["loc"],args))
            return

        if ast["kind"] == "literal":
            return

        assert False # should be unreachable

    # can't consume a linear reference and also reference it
    # can't consume something being controlled
    # can only consume something once per expn
    def consume_check(c,l):
        c_dict = {}
        for ident in c:
            name = ident["value"]

            if name in controls:
                raise SyntaxError("Can't consume symbol '"+name+"'"+error_at(ident["loc"],args)+\
                                    "that is also being controlled"+error_at(controls[name],args))

            if name in c_dict:
                raise SyntaxError("Can't consume symbol '"+name+"' twice"+error_at(ident["loc"],args)+\
                                    "and"+error_at(c_dict[name],args))

            c_dict[name] = ident["loc"]

        for ident in l:
            if ident["value"] in c_dict:
                raise SyntaxError("Can't consume symbol '"+ident["value"]+"'"+error_at(c_dict[ident["value"]],args)+\
                                    "and also reference it"+error_at(ident["loc"],args))


    ###############################

    def analyze_instruction(ast):
        if ast["kind"] == "declare_instruction":
            name = ast["identifier"]["value"]

            if name in scope:
                raise SyntaxError("Can't override symbol '"+name+"'"+error_at(ast["loc"],args)+\
                                    "previously declared"+error_at(scope[name]["loc"],args))

            if name in kwargs:
                raise SyntaxError("Can't override a keyword argument with symbol '"+name+"'"+error_at(ast["loc"],args))

            decl_check(ast["dim"])
            if ast["slots"] is not None: decl_check(ast["slots"])

            scope[name] = {
                "loc": ast["identifier"]["loc"],
                "array": ast["slots"] is not None
            }

            return

        if ast["kind"] in ["discard_instruction", "uncompute_instruction"]:
            for ident in ast["identifiers"]:
                name = ident["value"]
                loc = ident["loc"]
                verb = ast["kind"].split("_")[0]

                if name in kwargs:
                    raise SyntaxError("Can't "+verb+" a keyword argument '"+name+"'"+error_at(ast["loc"],args))

                if name not in scope:
                    raise SyntaxError("Can't "+verb+" an undefined symbol '"+name+"'"+error_at(ast["loc"],args))

                if name in controls:
                    raise SyntaxError("Can't "+verb+" symbol '"+name+"'"+error_at(ast["loc"],args)+\
                                        "because it is being controlled"+error_at(controls[name],args))

                del scope[name]

            return

        if ast["kind"] == "scalar_instruction":
            _, c, l = analyze_expn(ast["expn"])
            consume_check(c,l)

            for ident in c:
                del scope[ident["value"]]

            return

        if ast["kind"] == "pass_instruction":
            return

        if ast["kind"] == "repeat_instruction":
            keySymb = list_index_check(ast["expn"])

            added_ctrl = False
            if len(keySymb) > 0:
                if keySymb[0]["value"] not in controls:
                    added_ctrl = True
                    controls[keySymb[0]["value"]] = keySymb[0]["loc"]

            scope_backup = {key:{"loc":scope[key]["loc"], "array":scope[key]["array"] } for key in scope.keys()}

            for instr in ast["instructions"]:
                analyze_instruction(instr)

            for key in scope_backup.keys():
                if key not in scope:
                    raise SyntaxError("Symbol '"+key+"' no longer present after a repeat block"+error_at(ast["loc"])+\
                            "declared"+error_at(scope_backup[key]["loc"]))

            for key in scope.keys():
                if key not in scope_backup:
                    raise SyntaxError("New symbol '"+key+"' present at the end of a repeat block"+error_at(ast["loc"])+\
                            "declared"+error_at(scope[key]["loc"]))

                if scope[key]["array"] != scope_backup[key]["array"]:
                    before = "" if scope_backup[key]["array"] else " not"
                    beforeLoc = error_at(scope_backup[key]["loc"])
                    after = "" if scope[key]["array"] else " not"
                    raise SyntaxError("Symbol '"+key+"' was"+before+" an array before repeat block, declared"+beforeLoc+\
                                        "but is"+after+" after"+error_at(scope[key]["loc"]))

            if added_ctrl:
                del controls[keySymb[0]["value"]]

            return

        if ast["kind"] == "if_instruction":

            new_controls = []
            for cond in ast["conditions"]:
                t, c, l = analyze_expn(cond["condition"])

                if t != "scalar":
                    raise SyntaxError("If condition should be a scalar but is a block"+error_at(cond["condition"],args))
                assert len(c) == 0

                for ident in l:
                    if ident["value"] not in scope and ident["value"] not in new_controls:
                        controls[ident["value"]] = ident["loc"]
                        new_controls.append(ident["value"])

            for instrs in [cond["instructions"] for cond in ast["conditions"]] + [ast["else"]]:

                scope_backup = {key:{"loc":scope[key]["loc"], "array":scope[key]["array"] } for key in scope.keys()}

                for instr in instrs:
                    analyze_instruction(instr)

                for key in scope_backup.keys():
                    if key not in scope:
                        raise SyntaxError("Symbol '"+key+"' no longer present after an if block,"+\
                                "declared"+error_at(scope_backup[key]["loc"]))

                for key in scope.keys():
                    if key not in scope_backup:
                        raise SyntaxError("New symbol '"+key+"' present at the end of an if block"+\
                                "declared"+error_at(scope[key]["loc"]))

                    if scope[key]["array"] != scope_backup[key]["array"]:
                        before = "" if scope_backup[key]["array"] else " not"
                        beforeLoc = error_at(scope_backup[key]["loc"])
                        after = "" if scope[key]["array"] else " not"
                        raise SyntaxError("Symbol '"+key+"' was"+before+" an array before if block, declared"+beforeLoc+\
                                            "but is"+after+" after"+error_at(scope[key]["loc"]))

            for ctrl in new_controls:
                del controls[ctrl]

            return


        if ast["kind"] == "init_instruction":

            prv_targets = {}
            key_linsymbs = []

            for target in ast["targets"]:
                name = target["target"]["value"]
                loc = target["target"]["loc"]

                if name in kwargs:
                    raise SyntaxError("Can't initialize a keyword argument '"+name+"'"+error_at(loc,args))

                if name not in scope:
                    raise SyntaxError("Can't initialize an undefined symbol '"+name+"'"+error_at(loc,args))

                if name in prv_targets:
                    raise SyntaxError("Can't initialize symbol '"+name+"' twice"+error_at(loc,args)+\
                                        "and"+error_at(prv_targets[name],args))
                prv_targets[name] = loc

                if name in controls:
                    raise SyntaxError("Can't initialize symbol '"+name+"'"+error_at(loc,args)+\
                                    "because it is being controlled"+error_at(controls[name],args))

                if scope[name]["array"] and target["key"] is None:
                    raise SyntaxError("Symbol '"+name+"' is an array and is missing a key"+error_at(loc,args))

                if not scope[name]["array"] and target["key"] is not None:
                    raise SyntaxError("Symbol '"+name+"' is not an array but has a key"+error_at(loc,args))

                if target["key"] is not None:
                    key_linsymbs += list_index_check(target["key"])

            for symb in key_linsymbs:
                name = symb["value"]
                if name in prv_targets:
                    raise SyntaxError("Can't use symbol '"+name+"' in a list index"+error_at(symb["loc"],args)+\
                                        "if it is also being initialized"+error_at(prv_targets[name],args))

            t, c, l = analyze_expn(ast["expn"])
            if t != "block":
                raise SyntaxError("Initialization expression should be a block but is a scalar"+error_at(ast["loc"],args))
            consume_check(c,l)


            for ident in c:
                if ident["value"] in prv_targets: continue
                del scope[ident["value"]]

            return

        if ast["kind"] in ["assign_instruction", "increment_instruction", "decrement_instruction"]:
            verb = ast["kind"].split("_")[0]
            if verb == "assign": verb += " to"

            name = ast["target"]["value"]
            loc = ast["target"]["loc"]

            if name in kwargs:
                raise SyntaxError("Can't "+verb+" a keyword argument '"+name+"'"+error_at(loc,args))

            if name not in scope:
                raise SyntaxError("Can't "+verb+" an undefined symbol '"+name+"'"+error_at(loc,args))

            if name in controls:
                raise SyntaxError("Can't "+verb+" symbol '"+name+"'"+error_at(loc,args)+\
                                "because it is being controlled"+error_at(controls[name],args))

            if scope[name]["array"] and ast["key"] is None:
                raise SyntaxError("Symbol '"+name+"' is an array and is missing a key"+error_at(loc,args))

            if not scope[name]["array"] and ast["key"] is not None:
                raise SyntaxError("Symbol '"+name+"' is not an array but has a key"+error_at(loc,args))

            if ast["key"] is not None:
                list_index_check(ast["key"])

            verb = ast["kind"].split("_")[0].capitalize()

            t,c,l = analyze_expn(ast["expn"])
            if t != "scalar":
                raise SyntaxError(verb.capitalize()+" expression should be a scalar but is a block"+error_at(ast["expn"],args))
            assert len(c) == 0

            if verb != "Assign":
                for ident in l:
                    if ident["value"] == name:
                        raise SyntaxError(verb.capitalize()+" expression can't refer to target symbol"+error_at(ident["loc"],args))
            else:
                if ast["key"] is not None and ast["undo"] is not None:
                    raise SyntaxError("Cant provide undo expression for array assignment"+error_at(ast["expn"],args))

                if ast["key"] is not None:
                    for ident in l:
                        if ident["value"] == name:
                            raise SyntaxError("Array assign expression can't refer to target symbol"+error_at(ident["loc"],args))

                if ast["undo"] is not None:
                    t,c,l = analyze_expn(ast["undo"])
                    if t != "scalar":
                        raise SyntaxError("Assign undo expression should be a scalar but is a block"+error_at(ast["expn"],args))
                    assert len(c) == 0


            return

        assert False # should be unreachable

    ############################

    for decl in in_decls:
        if decl["identifier"]["value"] in kwargs:
            thisloc = error_at(decl["identifier"]["loc"],args)
            raise SyntaxError("Can't override a keyword argument with symbol '"+decl["identifier"]["value"]+"'"+thisloc)

        if decl["identifier"]["value"] in scope:
            thisloc = error_at(decl["identifier"]["loc"],args)
            prvloc = error_at(scope[decl["identifier"]["value"]]["loc"],args)
            raise SyntaxError("Can't take symbol '"+decl["identifier"]["value"]+"' as input more than once"+thisloc\
                    +"previously used"+prvloc)

        decl_check(decl["dim"])
        if decl["slots"] is not None: decl_check(decl["slots"])

        scope[decl["identifier"]["value"]] = {
            "loc": decl["identifier"]["loc"],
            "array": decl["slots"] is not None
        }

    # temporary, to make sure we can't do this: [x:2 @ x:2 <- ]
    # only contains location
    out_decl_scope = {}

    for decl in out_decls:
        if decl["identifier"]["value"] in kwargs:
            thisloc = error_at(decl["identifier"]["loc"],args)
            raise SyntaxError("Can't override a keyword argument with symbol '"+decl["identifier"]["value"]+"' "+thisloc)

        if decl["identifier"]["value"] in out_decl_scope:
            thisloc = error_at(decl["identifier"]["loc"],args)
            prvloc = error_at(scope[decl["identifier"]["value"]]["loc"],args)
            raise SyntaxError("Can't take symbol '"+decl["identifier"]["value"]+"' as output more than once"+thisloc\
                    +"previously used"+prvloc)

        decl_check(decl["dim"])
        if decl["slots"] is not None: decl_check(decl["slots"])

        out_decl_scope[decl["identifier"]["value"]] = decl["identifier"]["loc"]
        if decl["identifier"]["value"] not in scope:
            scope[decl["identifier"]["value"]] = {
                "loc": decl["identifier"]["loc"],
                "array": decl["slots"] is not None
            }

    for instr in instrs:
        analyze_instruction(instr)

    for decl in out_decls:
        if decl["identifier"]["value"] not in scope:
            thisloc = error_at(decl["identifier"]["loc"],args)
            raise SyntaxError("Output symbol '"+decl["identifier"]["value"]+"' declared"+thisloc+"is not present at end of scope.")

    for var in scope:
        if var not in out_decl_scope:
            raise SyntaxError("Symbol '"+var+"' declared"+error_at(scope[var]["loc"],args)+"is still present at end of scope and is not an output.")





