
# December 2020 implementations of qac.symmetrize and block_add

# these were more ambitious in terms of optimization, but also look buggy,
# and the block_add one isn't even complete.

# maybe some code/ideas from here will come in handy to optimize the
# current implementation, which is simpler but also slower.


 def symmetrize(self):
        assert len(self.unnamed_inputs) == 0
        assert len(self.unnamed_outputs) == 0

        # assert qac.scope_input and qac.scope_output to have the same set of keys
        for key in self.scope_input.keys(): assert key in self.scope_output
        for key in self.scope_output.keys(): assert key in self.scope_input

        # Example 1:
        #   input = {"x":<reg1>, "y":<reg2>}
        #   output = {"x":<reg2>, "y":<reg1>}

        # solution:
        #   swap <reg1> <reg2>


        # Example 2:
        #   input = {"x":<reg1>}
        #   instrs:
        #       declare <reg2>
        #       <reg2> += <reg1>
        #       discard <reg1>
        #   output = {"x":<reg2>}

        # solution:
        #       declare <reg1>
        #       swap <reg1> <reg2>
        #       zero <reg2>

        # Step 1: make scope_output just be a relabeling of scope_input
        # that is, scope_output.values() == scope_input.values()

        # If scope_output[key] and scope_input[key] are never in scope at the same time,
        # then I can just do scope_output[key].substitute(scope_input[key]).
        # Otherwise I need re-declare scope_input[key], swap the data from scope_output[key],
        # and then zero scope_output[key].
        for key in self.scope_input.keys():

            reg = self.scope_input[key]
            if reg in self.scope_outputs.values():
                # Already an output. Don't have do do anything for this step.
                continue

            # Check if scope_output[key] and scope_input[key] are ever in scope at the same time, i.e. clash.
            # This clash-checker is not smart enough to find situations like
            # declare <reg1>
            # discard <reg2>
            # where technically <reg1> and <reg2> are in scope at the same time,
            # but you could just swap the order of these two instructions to resolve
            # the clash. But whatever. Best to try to not create situations like this
            # in init or assign. Certainly the user can't create this situation - what would they type?

            # If the output register is the inputs, then it has been in scope the whole time!
            # input = {"x": <reg1>, "y":<reg2>}
            # discard <reg1>
            # declare <reg3>
            # output = {"x": <reg2>, "y":<reg3>}
            # solution for x: rename <reg2> <reg1>
            # now output = {"x": <reg1>, "y":<reg3>}
            # then for y, swap <reg2> <reg3>; zero <reg2>
            # now output = {"x": <reg1>, "y":<reg2>}
            if self.scope_outputs[key] in self.scope_inputs.values():
                self.instrs.append({"kind":"qac_rename", "source":self.scope_outputs[key], "target":reg})
                scope.scope_outputs[key] = reg
                continue

            # Loop through the instructions in reverse, looking for declares.
            to_remove = None
            for instr in reversed(self.instrs):
                if instr["kind"] == "qac_declare":
                    if instr["reg"] == self.scope_outputs[key]:
                        # No clash: scope_outputs[key] was declared after scope_inputs[key] was removed.
                        self.scope_outputs[key].substitute(reg)
                        break

                if instr["kind"] in ["qac_discard", "qac_zero"]:
                    if instr["reg"] == reg:
                        # Clash: scope_outputs[key] was declared before scope_inputs[key] was removed.
                        self.instrs.append({"kind":"qac_rename", "source":self.scope_outputs[key], "target":reg})
                        scope.scope_outputs[key] = reg
                        break

                # can safely ignore qac_if because those instructions can't change the scope.

                if instr["kind"] == "qac_rename":
                    if instr["source"] == reg and instr["target"] == self.scope_outputs[key]:
                        # I can just remove this instruction now, and make the registers the same.
                        self.scope_outputs[key].substitute(reg)
                        to_remove = instr
                        break

                    if instr["source"] == reg:
                        # Clash: scope_outputs[key] was declared before (or at the same time as) scope_inputs[key] was removed.
                        self.instrs.append({"kind":"qac_rename", "source":self.scope_outputs[key], "target":reg})
                        scope.scope_outputs[key] = reg
                        break

                    if instr["target"] == self.scope_outputs[key]:
                        # No clash: scope_outputs[key] was declared after scope_inputs[key] was removed.
                        self.scope_outputs[key].substitute(reg)
                        break


            if to_remove is not None:
                self.instrs.remove(to_remove)

        # Step 2: Permute the data into the correct registers.
        while True:
            all_good = True
            for key in self.scope_input.keys():
                if self.scope_outputs[key] == self.scope_inputs[key]:
                    # Key is already correct.
                    continue

                # Greedily make this key correct by swapping in the register it wants.
                # The other_key is necessarily also wrong, so this operation can't make it more wrong.
                # Thus, the total number of correct keys will increase by at least one.
                # There are only so many incorrect keys, so this loop must eventually terminate.

                # The other_key holds the register that should be assigned to key.
                other_key = [k for k in self.scope_outputs if self.scope_outputs[k] == self.scope_inputs[key]][0]

                # Swap them.
                self.instrs.append({"kind":"qac_swap", "reg1":self.scope_outputs[other_key],
                                                         "reg2":self.scope_outputs[key] })
                self.scope_outputs[other_key], self.scope_outputs[key] = self.scope_outputs[key], self.scope_outputs[other_key]

                all_good = False
            if all_good: break



def block_add(*blocks):
    assert len(blocks) > 1

    b0 = blocks[0]

    out = QAC()
    out.parent = b0.parent

    # check some things
    for b in blocks:
        assert b.parent == b0.parent # all blocks must have the same parent.
        assert len(b.scope_inputs.keys()) == 0
        assert len(b.scope_outputs.keys()) == 0

        for i in range(len(b.unnamed_inputs)):
            assert b.unnamed_inputs[i].dim == b0.unnamed_inputs[i].dim

        for i in range(len(b.unnamed_outputs)):
            assert b.unnamed_outputs[i].dim == b0.unnamed_outputs[i].dim

    ##################################################
    # We will be creating lots of qac_if statements for the LCU circuit.
    # Each if_statement can't change the scope: it must have the same registers before and after.
    # The given blocks can have more inputs than outputs, or permute their registers, etc.
    # So we can't use the given blocks' registers as a starting point. We need
    # an official set of registers that are present at the start and end of each if statement,
    # and a template to map the given blocks onto them.

    official_registers = {}

    input_dim_counts = {}
    # count how many registers of each dimension we have in the inputs
    for reg in b0.unnamed_inputs:
        dim = reg.dim
        if dim not in final_dims: official_registers[dim] = []
        if dim not in input_dim_counts: input_dim_counts[dim] = 0
        input_dim_counts[dim] += 1

    # same for outputs
    output_dim_counts = {}
    for reg in b0.unnamed_outputs:
        dim = reg.dim
        if dim not in final_dims: official_registers[dim] = []
        if dim not in output_dim_counts: output_dim_counts[dim] = 0
        output_dim_counts[dim] += 1

    # create enough registers for both
    for dim in official_registers:
        if dim in input_dim_counts:
            while len(official_registers[dim]) < input_dim_counts[dim]:
                official_registers[dim].append(Register(dim))
        if dim in output_dim_counts:
            while len(official_registers[dim]) < output_dim_counts[dim]:
                official_registers[dim].append(Register(dim))

    # populate out.unnamed_inputs: these get mapped to the blocks[0].unnamed_inputs
    for reg in b0.unnamed_inputs:
        dim = reg.dim
        i = 0
        while official_registers[dim][i] in out.unnamed_inputs: i += 1
        out.unnamed_inputs.append(official_registers[dim][i])

    # populate out.unnamed_outputs: these get mapped to the blocks[0].unnamed_outputs
    for reg in b0.unnamed_outputs:
        dim = reg.dim
        i = 0
        while official_registers[dim][i] in out.unnamed_outputs: i += 1
        out.unnamed_outputs.append(official_registers[dim][i])

    # We need to make sure that all the official_registers are in scope
    # during the if statements, so we need to declare the ones not covered by
    # the out.unnamed_inputs.

    before_instrs = []
    for dim in official_registers:
        for reg in official_registers[dim]:
            if reg not in out.unnamed_inputs:
                before_instrs.append({"kind":"qac_declare", "reg":reg, "dim":reg.dim})

    # Finally, we can't have any official_registers lying around that
    # are not in out.unnamed_outputs. We need to zero these out.

    after_instrs = []
    for dim in official_registers:
        for reg in official_registers[dim]:
            if reg not in out.unnamed_outputs:
                after_instrs.append({"kind":"qac_zero", "reg":reg})

    ###############################################
    # Now we need to perform the substitutions on the blocks to line them up with
    # out.unnamed_in/outputs

    for b in blocks:
        # I can just blindly substitute the inputs.
        for i in range(len(b.unnamed_inputs)):
            b.unnamed_inputs[i].substitute(out.unnamed_inputs[i])

        # The outputs are not so easy.
        for i in range(len(b.unnamed_outputs)):
            block_out = b.unnamed_outputs[i]
            official_out = out.unnamed_outputs[i]

            # register was already lined up via input substitution
            if block_out == official_out:
                continue

            # indicates block_out has not been substituted yet
            block_available = block_out not in b.unnamed_inputs

            # indicates either that official_out has not been substituted yet,
            # or, if it has, it must have been zeroed or discarded by b.instrs.
            # Either way, indicates availability for substitution.
            official_available = official_out not in b.unnamed_outputs

            # Both block's output and official output were never substituted
            # with anything else so I am safe to merge them.
            if block_available and official_available:
                block_out.substitute(official_out)
                continue

            # For the following two cases this is a good example:
            # input=[<reg1>], declare <reg2>, output=[<reg2>,<reg1>]

            # Block output is free, but official output is not available
            # This is the first output in the example:
            #  b.unnamed_outputs[0] == <reg2> and out.unnamed_outputs[0] == <reg1>
            if block_available and not official_available:
                j = i
                while True:
                    k = [k for k in range(len(b.unnamed_outputs)) if b.unnamed_outputs[k] == out.unnamed_outputs[j]]

                    # TODO: how do we prevent swapping twice??

                    b.instrs.append({"kind":"qac_swap", "reg1":out.unnamed_outputs[k], "reg2":out.unnamed_outputs[j]})

                    # if out.unnamed_outputs[k] is not free for substitution, we need to go again
                    # otherwise we can substitute and be done.
                    if out.unnamed_outputs[k] in b.unnamed_outputs:
                        j = k
                    else:
                        block_out.substitute(out.unnamed_outputs[k])
                        break


            # Block output is an input, but official output is unbound
            # Second output of the example:
            #  b.unnamed_outputs[0] == <reg1> and out.unnamed_outputs[0] == <reg2>
            if not block_available and official_available:

                continue

            # Both block output and official output are bound.
            # Example: input=[<reg1>,<reg2>], output=[<reg2>,<reg1>]


    # TODO: something about this divide-and-conquer approach doesn't feel right.
    # I think there may be an angle from which this gets really simple.

    # Possibilility: generate the permutation that each block implements.
    # Then apply the permutation (not it's inverse, right?)

    # Or: read parser.py and look at the old solution.

    for b in blocks:

        while True:
            all_good = True
            for i in range(len(b.unnamed_outputs)):
                if i >= len(b.unnamed_inputs): continue

                if b.unnamed_inputs[i] == b.unnamed_outputs[i]:
                    continue

                all_good = False

                if b.unnamed_outputs[i] not in b.unnamed_inputs:
                    # these need to be dealt with via other swaps
                    continue

                j = [j for j in range(len(b.unnamed_inputs)) if b.unnamed_inputs[j] == b.unnamed_outputs[i]][0]
                b.instrs.append({"kind":"qac_swap", "reg1":b.unnamed_outputs[i], "reg2":b.unnamed_outputs[j]})
                b.unnamed_outputs[i], b.unnamed_outputs[j] = b.unnamed_outputs[j], b.unnamed_outputs[i]
                break

            # TODO: this implementation is a bit different from your symmetrize implementation above.
            # Compare the two. The symmetrize one looks simpler.
            if all_good: break





    ############################################################
    # Finally we can assemble the LCU circuit

    # TODO: synthesize control register
    out.scale = sum([b.scale for b in blocks])
    ctrl = Register(len(blocks))
    flag = Register(2)

    for i in range(blocks):
        instrs = []

        instrs += blocks[i].instrs



        # flip flag
        # if statement
        # unflip flag

    # example: an lcu of
    #   scalar ~ket(0) * |x
    #   scalar ~ket(1) * |x

    # solution: pull all the register manipulation to the end?
    # this is probably what would happen irl.

    # harder example:
    #    discard x
    #

    # build the LCU circuit


