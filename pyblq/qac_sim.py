
import numpy as np
from .qac import Register, block_copy

# takes a qac with no inputs, returns the output state vector
def get_qac_state(qac):
    assert len(qac.scope_inputs) == 0
    assert len(qac.scope_outputs) == 0

    # purification simulation: state is stored as state vector
    # state = np.array of complex amplitudes for every nonzero entry
    # registers = {reg.id: np.array} where each np.array is the list of values
    # discards = list of np.array's corresponding to registers that have been discarded

    # possibilities may contain duplicates, so need to prune once in a while

    # cache performance relies on the idea that not too many registers are affected by a single instr
    # looping through all possibilities of a single register is fast,
    # but lookping through all possibilities of several registers is slow, since their arrays are far apart

    state = np.array([1]).astype(complex)
    registers = {}
    discards = []

    for reg in qac.unnamed_inputs:
        registers[reg.trace()] = np.array([0]).astype(int)

    def simulate_instruction(instr, state, registers, discards):
        # for debugging:
        # s = None
        # if "reg" in instr: s = str(instr["reg"].trace())
        # print(instr["kind"], s, registers.keys())
        print(instr["kind"])

        # {"kind":"qac_declare", "reg":<register>, "dim":<int>}
        if instr["kind"] == "qac_declare":
            tr = instr["reg"].trace()
            assert tr not in registers
            registers[tr] = np.zeros(len(state)).astype(int)
            return state, registers, discards

        # {"kind":"qac_discard", "reg":<register>}
        if instr["kind"] == "qac_discard":
            tr = instr["reg"].trace()
            arr = registers[tr]
            del registers[tr]

            # don't have to make garbage if all the values are the same
            # if not np.allclose(arr, arr[0]): discards.append(arr)
            # ^ actually, this doesn't work, because get_kraus will be feeding
            # this subroutine multiple circuits and comparing the results
            # leaving the garbage around is key to interpreting the results correctly

            discards.append(arr)

            return state, registers, discards

        # {"kind":"qac_maxmixed", "reg":<register>, "dim":<int>}
        if instr["kind"] == "qac_maxmixed":

            # make a bell state and discard half of it
            # the reason we do this this way is because we'd only like to
            # implement branching behavior once, namely in qac_unitary.

            dim = instr["reg"].dim
            tmp = Register(dim)
            omega = np.exp(2j*np.pi/dim)

            mat = []
            for i in range(dim):
                row = []
                for j in range(dim):
                    row.append(omega**(i*j) / np.sqrt(dim))
                mat.append(row)

            dummy_instrs = [
                    {"kind":"qac_declare", "reg":instr["reg"], "dim":dim},
                    {"kind":"qac_declare", "reg":tmp, "dim":dim},
                    {"kind":"qac_unitary", "reg":instr["reg"], "mat":mat},
                    {"kind":"qac_increment", "reg":tmp, "expn":{"kind": "register_expn", "register":instr["reg"]}},
                    {"kind":"qac_discard", "reg":tmp},
            ]

            for dummy_instr in dummy_instrs:
                state, registers, discards = simulate_instruction(dummy_instr, state, registers, discards)

            return state, registers, discards

        # {"kind":"qac_zero", "reg":<register>}
        if instr["kind"] == "qac_zero":
            arr = registers[instr["reg"].trace()]

            new_dim = 0
            for i in range(len(arr)):
                if arr[i] == 0: new_dim += 1

            new_state = np.zeros(new_dim).astype(complex)
            idx = 0
            for i in range(len(arr)):
                if arr[i] == 0:
                    new_state[idx] = state[i]
                    idx += 1

            new_registers = {}
            for reg in registers.keys():
                if reg == instr["reg"].trace(): continue
                new_arr = np.zeros(new_dim).astype(int)
                idx = 0
                for i in range(len(arr)):
                    if arr[i] == 0:
                        new_arr[idx] = registers[reg][i]
                        idx += 1
                new_registers[reg] = new_arr

            new_discards = []
            for old_arr in discards:
                new_arr = np.zeros(new_dim).astype(int)
                idx = 0
                for i in range(len(arr)):
                    if arr[i] == 0:
                        new_arr[idx] = old_arr[i]
                        idx += 1
                new_discards.append(new_arr)

            return new_state, new_registers, new_discards

        # {"kind":"qac_increment", "reg":<register>, "expn":<expn>}
        if instr["kind"] == "qac_increment":
            def eval_expn(idx, expn, registers):
                # {"kind": "register_expn", "register":<reg> }
                if expn["kind"] == "register_expn":
                    return complex(registers[expn["register"].trace()][idx])

                # {"kind": "value_expn", "value":5j}
                if expn["kind"] == "value_expn":
                    return complex(expn["value"])

                # {"kind": "sum_expn", "terms":[<linexp>] }
                if expn["kind"] == "sum_expn":
                    return sum([eval_expn(idx,sub_expn,registers) for sub_expn in expn["terms"]])

                # {"kind": "negate_expn", "expn":<linexp> }
                if expn["kind"] == "negate_expn":
                    return -eval_expn(idx,expn["expn"],registers)

                # {"kind": "adjoint_expn", "expn":<linexp> }
                if expn["kind"] == "adjoint_expn":
                    return eval_expn(idx,expn["expn"],registers).conj()

                # {"kind": "product_expn", "terms":[<linexp>] }
                if expn["kind"] == "product_expn":
                    return np.prod([eval_expn(idx,sub_expn,registers) for sub_expn in expn["terms"]])

                # {"kind": "division_expn", "dividend":<linexp>, "divisor":5j }
                if expn["kind"] == "division_expn":
                    return eval_expn(idx,expn["dividend"],registers) / expn["divisor"]

                # {"kind": "modulo_expn", "dividend":<linexp>, "divisor":5 }
                if expn["kind"] == "modulo_expn":
                    out = eval_expn(idx,expn["dividend"],registers)
                    assert isinstance(expn["divisor"],int)
                    assert expn["divisor"] > 0

                    # see comment in runtime.py
                    while out.real < 0: out += expn["divisor"]
                    while out.real >= expn["divisor"]: out -= expn["divisor"]

                    return out

                # {"kind": "boolean_expn", "terms":[<linexp>, <string>, <linexp>, <string>, ...] }
                # <string> is one of ==, !=, >, <, >=, <=
                if expn["kind"] == "boolean_expn":
                    terms = []
                    for i in range(len(expn["terms"])):
                        if i % 2 == 1:
                            assert expn["terms"][i] in ["==", "!=", "<", ">", ">=", "<="]
                            terms.append(expn["terms"][i])
                        else:
                            terms.append(eval_expn(idx,expn["terms"][i],registers))

                    for i in range(len(expn["terms"])):
                        if i % 2 == 0: continue

                        if terms[i] == "==": value = (terms[i-1] == terms[i+1])
                        elif terms[i] == "!=": value = (terms[i-1] != terms[i+1])
                        elif terms[i] == "<": value = (terms[i-1] < terms[i+1])
                        elif terms[i] == ">": value = (terms[i-1] > terms[i+1])
                        elif terms[i] == ">=": value = (terms[i-1] >= terms[i+1])
                        else:
                            assert terms[i] == "<="
                            value = (terms[i-1]["value"] <= terms[i+1]["value"])

                        if not value:
                            return complex(0)

                    return complex(1)

                assert False # unreachable

            arr = registers[instr["reg"].trace()]
            dim = instr["reg"].dim

            for i in range(len(state)):
                val = eval_expn(i, instr["expn"], registers)
                val = int(val.real)
                arr[i] = (arr[i] + val) % dim

            return state, registers, discards

        # {"kind":"qac_unitary", "reg":<register>, "mat":<matrix>}
        if instr["kind"] == "qac_unitary":
            arr = registers[instr["reg"].trace()]
            dim = instr["reg"].dim
            mat = instr["mat"]

            new_dim = len(state)*dim

            new_state = np.zeros(new_dim).astype(complex)
            idx = 0
            for i in range(len(arr)):
                val = arr[i]
                for j in range(dim):
                    new_state[idx] = state[i] * mat[j][val]
                    idx += 1

            new_registers = {}
            for reg in registers.keys():
                new_arr = np.zeros(new_dim).astype(int)
                idx = 0
                if reg == instr["reg"].trace():
                    for i in range(len(arr)):
                        for j in range(dim):
                            new_arr[idx] = j
                            idx += 1
                else:
                    for i in range(len(arr)):
                        for j in range(dim):
                            new_arr[idx] = registers[reg][i]
                            idx += 1
                new_registers[reg] = new_arr

            new_discards = []
            for old_arr in discards:
                new_arr = np.zeros(new_dim).astype(int)
                idx = 0
                for i in range(len(arr)):
                    for j in range(dim):
                        new_arr[idx] = old_arr[i]
                        idx += 1
                new_discards.append(new_arr)

            return new_state, new_registers, new_discards

        # {"kind":"qac_phase", "value":<complexnr>}
        if instr["kind"] == "qac_phase":
            for i in range(len(state)):
                state[i] *= complex(instr["value"])
            return state, registers, discards

        # {"kind":"qac_swap", "reg1":<register>, "reg2": <register>}
        if instr["kind"] == "qac_swap":
            tr1 = instr["reg1"].trace()
            tr2 = instr["reg2"].trace()

            registers[tr1], registers[tr2] = registers[tr2], registers[tr1]

            return state, registers, discards

        # {"kind":"qac_rename", "source":<register>, "target": <register>}
        if instr["kind"] == "qac_rename":
            src = instr["source"].trace()
            trg = instr["target"].trace()

            registers[trg] = registers[src]
            del registers[src]

            return state, registers, discards

        # {"kind":"qac_if", "cond":<register>, "instructions":[<instrs>] }
        if instr["kind"] == "qac_if":
            cond = instr["cond"].trace()
            arr = registers[cond]

            bad_dim = 0 # number of branches where we don't perform the instructions
            for i in range(len(arr)):
                if arr[i] == 0: bad_dim += 1

            bad_state = np.zeros(bad_dim).astype(complex)
            good_state = np.zeros(len(arr) - bad_dim).astype(complex)
            bad_idx, good_idx = 0,0
            for i in range(len(arr)):
                if arr[i] == 0:
                    bad_state[bad_idx] = state[i]
                    bad_idx += 1
                else:
                    good_state[good_idx] = state[i]
                    good_idx += 1

            bad_registers = {}
            good_registers = {}
            for reg in registers.keys():
                bad_arr = np.zeros(bad_dim).astype(int)
                good_arr = np.zeros(len(arr) - bad_dim).astype(int)
                bad_idx, good_idx = 0,0
                for i in range(len(arr)):
                    if arr[i] == 0:
                        bad_arr[bad_idx] = registers[reg][i]
                        bad_idx += 1
                    else:
                        good_arr[good_idx] = registers[reg][i]
                        good_idx += 1
                bad_registers[reg] = bad_arr
                good_registers[reg] = good_arr

            bad_discards = []
            good_discards = []
            for old_arr in discards:
                bad_arr = np.zeros(bad_dim).astype(int)
                good_arr = np.zeros(len(arr) - bad_dim).astype(int)
                bad_idx, good_idx = 0,0
                for i in range(len(arr)):
                    if arr[i] == 0:
                        bad_arr[bad_idx] = old_arr[i]
                        bad_idx += 1
                    else:
                        good_arr[good_idx] = old_arr[i]
                        good_idx += 1
                bad_discards.append(bad_arr)
                good_discards.append(good_arr)


            for i,sub_instr in enumerate(instr["instructions"]):
                good_state, good_registers, good_discards = simulate_instruction(sub_instr, good_state, good_registers, good_discards)

                if len(good_state) < 100: continue
                if sub_instr["kind"] not in ["qac_unitary", "qac_if"]: continue
                if i == len(instr["instructions"])-1: continue
                if instr["instructions"][i-1] in ["qac_unitary", "qac_if"]: continue

                good_state, good_registers, good_discards = prune(good_state, good_registers, good_discards)


            new_dim = len(good_state) + len(bad_state)

            new_state = np.zeros(new_dim).astype(complex)
            idx = 0
            for i in range(len(bad_state)):
                new_state[idx] = bad_state[i]
                idx += 1
            for i in range(len(good_state)):
                new_state[idx] = good_state[i]
                idx += 1

            for reg in good_registers.keys(): assert reg in bad_registers
            for reg in bad_registers.keys(): assert reg in good_registers

            new_registers = {}
            for reg in registers.keys():
                new_arr = np.zeros(new_dim).astype(int)
                idx = 0
                for i in range(len(bad_state)):
                    new_arr[idx] = bad_registers[reg][i]
                    idx += 1
                for i in range(len(good_state)):
                    new_arr[idx] = good_registers[reg][i]
                    idx += 1
                new_registers[reg] = new_arr

            assert len(good_discards) >= len(bad_discards)

            new_discards = []
            for i in range(len(good_discards)):
                new_arr = np.zeros(new_dim).astype(int)

                if i < len(bad_discards):
                    idx = 0
                    for j in range(len(bad_state)):
                        new_arr[idx] = bad_discards[i][j]
                        idx += 1
                    for j in range(len(good_state)):
                        new_arr[idx] = good_discards[i][j]
                        idx += 1
                else:
                    idx = len(bad_state)
                    for j in range(len(good_state)):
                        new_arr[idx] = good_discards[i][j]
                        idx += 1

                new_discards.append(new_arr)


            return new_state, new_registers, new_discards

        assert False # unreachable

    def prune(state, registers, discards):
        value_dict = {}

        for i in range(len(state)):
            if np.allclose(state[i], 0): continue

            key = []
            for reg in registers.keys():
                key.append(registers[reg][i])
            for disc in discards:
                key.append(disc[i])
            key = tuple(key)

            if key not in value_dict: value_dict[key] = 0j
            value_dict[key] += state[i]
            if np.allclose(value_dict[key],0):
                del value_dict[key]

        new_dim = len(value_dict)
        new_state = np.zeros(new_dim).astype(complex)
        new_registers = {}
        new_discards = []
        for reg in registers.keys():
            new_registers[reg] = np.zeros(new_dim).astype(int)
        for i in range(len(discards)):
            new_discards.append(np.zeros(new_dim).astype(int))

        for i,key in enumerate(value_dict.keys()):
            new_state[i] = value_dict[key]
            idx = 0
            for reg in registers.keys():
                new_registers[reg][i] = key[idx]
                idx += 1
            for j in range(len(discards)):
                new_discards[j][i] = key[idx]
                idx += 1

        # Can't do this, because Kraus operator calculation needs these
        # to_remove = []
        # for disc in new_discards:
        #    if np.allclose(disc,disc[0]): to_remove.append(disc)
        # for disc in to_remove:
        #    new_discards.remove(disc)

        return new_state, new_registers, new_discards


    for i,instr in enumerate(qac.instrs):
        state, registers, discards = simulate_instruction(instr, state, registers, discards)

        if len(state) < 100: continue
        if instr["kind"] not in ["qac_unitary", "qac_if"]: continue
        if i == len(qac.instrs)-1: continue
        if qac.instrs[i-1] in ["qac_unitary", "qac_if"]: continue

        state, registers, discards = prune(state, registers, discards)

    state, registers, discards = prune(state, registers, discards)

    return state, registers, discards


def sample_state(state,registers,discards):

    # state was pruned as a final step, so no duplicate branches exist.

    probabilities = np.zeros(len(state))
    for i in range(len(state)):
        probabilities[i] = np.abs(state[i])**2

    total = sum(probabilities)
    if total > 1:
        if np.allclose(total,1): probabilities /= total # small overshoot? fine. normalize.
        else: assert False # probability distribution is over-normalized.
    if np.random.random() > total: return None

    idx = np.random.choice(len(state), p=probabilities/total)

    out = []
    for reg in registers.keys():
        out.append(registers[reg][idx])

    return tuple(out)


def get_kraus(orig_qac):
    qac = block_copy(orig_qac)

    assert len(qac.scope_inputs) == 0
    assert len(qac.scope_outputs) == 0

    assert len(qac.unnamed_inputs) in [0,1]
    assert len(qac.unnamed_outputs) in [0,1]

    if len(qac.unnamed_outputs) == 1:
        out_reg = qac.unnamed_outputs[0]
        out_dim = out_reg.dim
    else:
        out_dim = 1

    # keys are tuples with the value of the discards
    out = {}

    if len(qac.unnamed_inputs) == 0:

        state, registers, discards = get_qac_state(qac)

        for j in range(len(state)):
            key = []
            for disc in discards: key.append(disc[j])
            key = tuple(key)

            if key not in out:
                out[key] = np.zeros((out_dim,1)).astype(complex)

            if len(qac.unnamed_outputs) == 0:
                out[key][0, 0] = state[j]
            else:
                val = registers[out_reg.trace()][j]
                out[key][val, 0] = state[j]

    else:
        in_reg = qac.unnamed_inputs[0]
        in_dim = in_reg.dim

        qac.unnamed_inputs = []
        qac.instrs = [
            {"kind":"qac_declare", "reg":in_reg, "dim":in_dim},
            {"kind":"qac_increment", "reg":in_reg, "expn": {"kind": "value_expn", "value": complex(0)} },
        ] + qac.instrs

        num_discards = None

        for i in range(in_dim):
            qac.instrs[1]["expn"]["value"] = complex(i)
            state, registers, discards = get_qac_state(qac)

            if num_discards is None: num_discards = len(discards)
            assert num_discards == len(discards)

            for j in range(len(state)):
                key = []
                for disc in discards: key.append(disc[j])
                key = tuple(key)

                if key not in out:
                    out[key] = np.zeros((out_dim,in_dim)).astype(complex)

                if len(qac.unnamed_outputs) == 0:
                    out[key][0, i] = state[j]
                else:
                    val = registers[out_reg.trace()][j]
                    out[key][val, i] = state[j]

    return [ qac.scale*mat for mat  in out.values()]




