

from .qac import Register, block_copy

def expose(qac, which=[]):
    qac = block_copy(qac)
    assert isinstance(which, list)
    for kind in which:
        assert kind in ["declare","discard","maxmixed","zero"]

    assert len(qac.scope_inputs) == 0
    assert len(qac.scope_outputs) == 0

    instrs, regs = expose_instrs(qac.instrs,is_if=False)

    qac.instrs = instrs

    if "declare" not in which:
        for reg in regs["declare"]:
            qac.instrs.insert(0,{
                "kind":"qac_declare", "reg":reg, "dim":reg.dim})
        del regs["declare"]
    else:
        qac.unnamed_inputs += regs["declare"]

    if "discard" not in which:
        for reg in regs["discard"]:
            qac.instrs.append({
                "kind":"qac_discard", "reg":reg})
        del regs["discard"]
    else:
        qac.unnamed_outputs += regs["discard"]

    if "maxmixed" not in which:
        for reg in regs["maxmixed"]:
            qac.instrs.insert(0,{
                "kind":"qac_maxmixed", "reg":reg, "dim":reg.dim})
        del regs["maxmixed"]
    else:
        qac.unnamed_inputs += regs["maxmixed"]

    if "zero" not in which:
        for reg in regs["zero"]:
            qac.instrs.append({
                "kind":"qac_zero", "reg":reg})
        del regs["zero"]
    else:
        qac.unnamed_outputs += regs["zero"]

    return qac, regs



def expose_instrs(instrs,is_if=True):

    out_instrs = []

    declare_regs = []
    discard_regs = []
    maxmixed_regs = []
    zero_regs = []

    case1_regs = []
    case2_partners = {}

    # If this is an if statement, then the scope must be preserved.
    # Thus, there are only two cases:

    # Case 1:
    #    a new register is created and destroyed within the instructions
    #    solution: just move the creation and destruction out.

    # Case 2:
    #    an input register is destroyed, and then later re-created
    #    solution: make a temporary register, created and destroyed outside.
    #              at point of destruction (or anywhere in between destruction
    #              and re-creation) swap the register with the temporary register.

    # If this is not an if statement then anything can happen.

    for instr in instrs:

        # {"kind":"qac_declare", "reg":<register>, "dim":<int>}
        if instr["kind"] == "qac_declare":

            if not is_if:
                declare_regs.append(instr["reg"])
                continue

            if instr["reg"] not in case2_partners:
                # case 1: this declared register must eventually be destroyed
                declare_regs.append(instr["reg"])
                case1_regs.append(instr["reg"])

            else:
                # case 2: this re-creates a previously destroyed register
                declare_regs.append(case2_partners[instr["reg"]])
                del case2_partners[instr["reg"]]
            continue

        # {"kind":"qac_discard", "reg":<register>}
        if instr["kind"] == "qac_discard":

            if not is_if:
                discard_regs.append(instr["reg"])
                continue

            if instr["reg"] in case1_regs:
                # case 1: this destroys a previously created register
                discard_regs.append(instr["reg"])
                case1_regs.remove(instr["reg"])

            else:
                # case 2:  this register must be created again later.
                tmpreg = Register(instr["reg"].dim)
                case2_partners[instr["reg"]] = tmpreg
                out_instrs.append({"kind":"qac_swap", "reg1":tmpreg, "reg2":instr["reg"]})
                discard_regs.append(tmpreg)
            continue

        # {"kind":"qac_maxmixed", "reg":<register>, "dim":<int>}
        if instr["kind"] == "qac_maxmixed":

            if not is_if:
                maxmixed_regs.append(instr["reg"])
                continue

            if instr["reg"] not in case2_partners:
                # case 1: this declared register must eventually be destroyed
                maxmixed_regs.append(instr["reg"])
                case1_regs.append(instr["reg"])

            else:
                # case 2: this re-creates a previously destroyed register
                maxmixed_regs.append(case2_partners[instr["reg"]])
                del case2_partners[instr["reg"]]
            continue

        # {"kind":"qac_zero", "reg":<register>}
        if instr["kind"] == "qac_zero":

            if not is_if:
                zero_regs.append(instr["reg"])
                continue

            if instr["reg"] in case1_regs:
                # case 1: this destroys a previously created register
                zero_regs.append(instr["reg"])
                case1_regs.remove(instr["reg"])

            else:
                # case 2:  this register must be created again later.
                tmpreg = Register(instr["reg"].dim)
                case2_partners[instr["reg"]] = tmpreg
                out_instrs.append({"kind":"qac_swap", "reg1":tmpreg, "reg2":instr["reg"]})
                zero_regs.append(tmpreg)

            continue

        # {"kind":"qac_if", "cond":<register>, "instructions":[<instrs>] }
        if instr["kind"] == "qac_if":

            if_instrs, regs = expose_instrs(instr["instructions"])

            declare_regs += regs["declare"]
            discard_regs += regs["discard"]
            maxmixed_regs += regs["maxmixed"]
            zero_regs += regs["zero"]

            out_instrs.append({
                "kind":"qac_if",
                "cond":instr["cond"],
                "instructions":if_instrs
                })

            continue

        # {"kind":"qac_rename", "source":<register>, "target": <register>}
        if instr["kind"] == "qac_rename":
            out_instrs.append(instr)

            if not is_if: continue

            if instr["source"] in case1_regs:
                case1_regs.remove(instr["source"])
                case1_regs.remove(instr["target"])

            if instr["source"] in case2_partners:
                case2_partners[instr["target"]] = case2_partners[instr["source"]]
                del case2_partners[instr["source"]]

            continue

        # {"kind":"qac_increment", "reg":<register>, "expn":<expn>}
        # {"kind":"qac_unitary", "reg":<register>, "mat":<matrix>}
        # {"kind":"qac_phase", "value":<complexnr>}
        # {"kind":"qac_swap", "reg1":<register>, "reg2": <register>}
        out_instrs.append(instr)


    assert len(case1_regs) == 0
    assert len(case2_partners) == 0

    return out_instrs, {
            "declare": declare_regs,
            "discard": discard_regs,
            "maxmixed": maxmixed_regs,
            "zero": zero_regs,
    }

