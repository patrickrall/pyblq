
global register_count
register_count = 0
global register_refs
register_refs = {}
global register_dims
register_dims = {}

# a data structure such that:
#    identity unique upon initialization
#    can be merged with other registers
#    can can be an array or not: can specify dimension and slots

# TODO: should the Register perhaps know its user-given name, to make the final qaasm more readable?
# how does that gel with register merging? It'll make sense for a symmetrized merge, but not a block-algebra merge.

def Register():
    def __init__(self, dim, slots=None):
        assert int(dim) == dim
        assert dim > 1

        register_count += 1
        self.id = register_count
        register_refs[self.id] = None
        register_dims[self.id] = (dim,slots)

    def trace(self):
        out = self.id
        while register_refs[out] != None:
            out = register_refs[out]
        return out

    @property
    def dim(self):
        return register_dims[self.trace()][0]

    @property
    def slots(self):
        return register_dims[self.trace()][1]

    def __eq__(self,other):
        if not isinstance(other,Register): return False
        return self.trace() == other.trace()

    def substitute(self,other):
        assert isinstance(other,Register)
        assert self.dim == other.dim
        assert self.slots == other.slots
        target = other.trace()
        if target == self.trace(): return

        key = self.id
        while register_refs[key] != None:
            key = register_refs[key]
        register_refs[key] = target




###################################

# Where do the responsibilities of this class end and those of the runtime begin?
# Runtime should do:
#     parsing the ast.
#     pre-evaluation of expns
#     distinguish between block expns and value expns
# QAASM should do:
#     circuit synthesis
#     managing the scope

# scope ops:
#     key is removed because it was consumed
#     value is swapped out because of a relabeling
#     block is symmetrized: inputs must equal outputs (need to know scope before and after)
#     two blocks are matched: inputs=inputs, outputs=outputs (need to know scope before and after)
#        The whole {"x":1,...},["x"]  system worked pretty well for that.


# QAASM data structure keeps track of both reg objects and their names in the scope.
# Blq objects just keep track of the reg objects.

class QAASM():

    def __init__(self,parent=None):
        self.scale = 1

        self.instrs = []

        # Key idea: I can lazily swap/relabel registers by manipulating the self.outputs dictionary.
        # Only when I need to symmetrize or align do I need to actually implement a permutation using qaasm.

        self.inputs  = {}
        self.outputs  = {}

        # Both None if unshaped. Both are lists if shaped. Check via self.shaped.
        # Needs to be unshaped in order to add instructions.
        # Needs to be shaped in order to do algebra.
        self.lshape = None   # an ordered subset of self.inputs.keys(), those not in ordering are the consumed registers
        self.rshape = None   # an ordering on self.outputs.keys()

        # There seems to be a difference in the needs of the methods:
        #     if, repeat, increment, scalar, init, assign
        #        all only really care about IF a variable is in scope, not about scope order of target block
        #     add, multiply, adjoint, tensorproduct
        #        do care about scope order
        # When is scope order determined?

        # types of blocks
        #   blocks in if and repeat statements: dont care about scope order at all
        #   ket() expn, consume expn, block cast, Blq's: can just make scope order correct upon init

        assert isinstance(parent,QAASM)
        self.parent = parent

        # Expressions can refer to keys in parent scope. Only if a register is declared/discarded/zero'd
        # or permuted in scope must ot be an output.

        # Can make something an explicit output by promoting it.
        # Should promotion make it a consume or an input?

    @property
    def shaped(self):
        if self.lshape is None:
            assert self.rshape is None
            return False
        assert isinstance(self.lshape,list)
        assert isinstance(self.rshape,list)
        return True

    # get key from parent scope
    def parentGet(self,key):
        if self.parent is None:
            raise KeyError()
        if key in self.parent.outputs:
            return self.parent.outputs[key]
        return self.parent[key]

    # check if parent has key
    def parentHas(self,key):
        if self.parent is None: return False
        if key in self.parent.outputs: return True
        return key in self.parent


    def promote(self, name):
        assert self.lshape is None and self.rshape is None

        assert self.parentHas(name)

        # check that 'name' was never in scope
        assert name not in self.inputs
        for instr in self.instrs:
            if instr["kind"] == "nqaasm_declare":
                assert instr["name"] != name
        assert name not in self.outputs

        prvreg = self.parentGet(name)
        reg = Register(prvreg.dim, slots=prvreg.slots)

        self.inputs[name] = reg
        self.outputs[name] = reg


    # named-qaasm aka nqaasm
    # its unclear to me that this is really that different
    # uses string register names rather than reg objects
    # except for declare which includes both. Regobj can be an array.
    # {"kind":"nqaasm_declare", "reg":<regobj>, "name":<name>}
    # {"kind":"nqaasm_discard", "name":<name>}
    # {"kind":"nqaasm_zero", "name":<name>}
    # {"kind":"nqaasm_increment", "name":<name>, "expn":<expn>}
    # {"kind":"nqaasm_unitary", "name":<name>, "mat":<matrix>}
    # {"kind":"nqaasm_phase", "value":<complexnr>}
    # {"kind":"nqaasm_swap", "name1":<name>, "name2":<name>}
    # {"kind":"nqaasm_if", "name":<register>, "instructions":[<instrs>] }

    def declare(self, name, dim, slots=None):
        assert self.lshape is None and self.rshape is None
        assert name not in self.outputs
        reg = Register(dim,slots=slots)
        self.instrs.append({"kind":"nqaasm_declare", "name":name, "reg":reg})

    def discard(self, name):
        assert self.lshape is None and self.rshape is None
        if name not in self.outputs: self.promote(name)
        assert name in self.outputs
        self.instrs.append({"kind":"qaasm_discard", "name":name})
        del self.outputs[name]

    # zero

    # Boring:
    # unitary
    # phase

    def increment(self, reg, expn):
        # if reg is not in scope, it has to be in parent scope, and needs to be promoted.

        # assert expn's regs are either in parent scope or in current scope and have the right shape

        # perhaps all the array decompiling does is make all indexes integers rather than variables

        def process_expn(expn):
            if expn["kind"] == "register_expn":
                if expn["key"] is None:
                    pass
                if isinstance(expn["key"],int):
                    pass
                if isinstance(expn["key"],str):
                    pass
            # recurse
            pass

        process_expn(expn)

        pass


    def symmetrize(self):
        # assert qaasm.input and qaasm.scope need to have the same set of keys
        for key in qaasm.input.keys():
            if qaasm.input[key] == qaasm.scope[key]:
                continue

            # check if there is any point in time when both qaasm.input[key] and qaasm.output[key]
            # are in scope. If yes, need to do a swap.

        pass


    def if(self, reg, qaasm):
        # reg is either in scope or in parent scope.
        # assert qaasm.lshape == qaasm.rshape == []

        assert qaasm.parent = self

        # qaasm.inputs need to be in self.scope. Promote if needed.
        qaasm.symmetrize()

        for key in qaasm.input.keys():
            if key not in self.scope: self.promote(key)
            self.scope[key].substitute(qaasm.input[key])



    def repeat(self, qaasm, count):
        # same as if, basically.
        pass



    ###################

    def scalar_instr(self, qaasm):
        # how to tell the runtime how the scope changed?
        # qaasm.rshape == qaasm.lshape == []

        assert qaasm.parent = self

        # promote any qaasm.inputs if needed, and wire them up
        for key in qaasm.input.keys():
            if key not in self.scope: self.promote(key)
            self.scope[key].substitute(qaasm.input[key])

        # delete any consumed variables
        for key in self.scope.keys():
            if key in qaasm.input and key not in qaasm.scope:
                del qaasm.scope[key]

        assert len(qaasm.scope.keys()) == 0

        self.scale *= qaasm.scale

        for instr in qaasm.instrs:
            self.instrs.append(instr)


    def init_instr(self, targets, qaasm):
        assert qaasm.parent = self
        assert len(qaasm.rshape) == 0


        for key in qaasm.scope: assert key in qaasm.lshape # is this always true anyway?
        # for key in qaasm.lshape: assert key in qaasm.scope # this should be true anyway
        assert len(targets) = len(qaasm.lshape)

        # promote any qaasm.inputs if needed, and wire them up
        for key in qaasm.input.keys():
            if key not in self.scope: self.promote(key)
            self.scope[key].substitute(qaasm.input[key])

        # delete any consumed variables
        for key in self.scope.keys():
            if key in qaasm.input and key not in qaasm.scope:
                del qaasm.scope[key]

        for i in range(len(targets)):
            target = targets[i]
            key = qaasm.lshape[i]
            reg = qaasm.scope[key]
            assert

        pass


    def assign_instr(self, reg, expn):
        pass

    def assign_undo_instr(self, reg, expn, undo):
        pass

    def assign_array_instr(self, key, regs, expn):
        pass

    def assign_array_undo_instr(self, key, regs, expn, undo):
        pass

    ############################


# The difference between nqaasm and regular qaasm:
#   - nqaasm knows what names the user has given to the variables.
#   - nqaasm can implement user-level permutations and relabelings without actually generating instructions
#   - nqaasm can't really be obtained from matrix literals or create expressions. (this is a problem!)
#   - If nqaasm is serialized all the labeling information is lost. It can't be deserialized.
#        - Need support for temporary names in nqaasm, which is the very problem registers are supposed to solve.

# "nqaasm_unnamed_declare"?
# have register objects hold on to their user-level names? That merges nqaasm with qaasm, but gets rid of permutation facility.
# if swap is a qaasm instruction, then can't the swap overhead be reduced in post?


# Idea: instructions are inherently named. algebra is inherently unnamed
# consume, create, cast are sort-of the boundary between named and unnamed.

# three types: referenced / scoped / unnamed

# blocks as inputs to if statements can't have any unnamed registers.
# algebraic blocks can't have any scoped registers as output.
# user-level blocks can't have any scoped registers and referenced registers

# what registers are what is determined by their presence in the bookkeeping dictionaries
#  not by qaasm. Qaasm only knows about registers.

# should qaasm support arrays, just with fixed indices?


###################

# Proposal
#   QAASM blocks are unshaped, and instructions can be appended to them
#   Blq blocks are shaped and instruction immutable - can only be manipulated via block algebra

# problems with this proposal:
#   Blq objects still need to be able to refer to things in scope, and are thus still nqaasm. Different from userspace blqs.
#

# Three types blocks:
#   Unshaped QAASM. Basically a bag of instructions. Can add instructions, can't do algebra.
#   Shaped QAASM. Knows how to refer to parent scope. Only mutable through algebra.
#   Userspace Blocks. Doesn't know anything about scope.

# Question: why do userspace blocks and shaped qaasm need to be different?
#   It still seems userspace blocks are just a restricted version of shaped qaasm.
#   Especially if I need to convert back and forth between the two in order to do anything.
# Similarities and differences:
#   They both use reg objects.
#   Named qaasm vs regular qaasm. Named qaasm is optimized for a named scope.
#       The whole register/scope system is still somewhat unprincipled.
#   Userspace blocks don't know about parent scope, or scope at all.
#   Open: can userspace blocks permute through relabeling?

# should userspace blocks use reg objects? Yes.
#     if no: need to interconvert a lot
#     if yes: lots of spare reg objects floating around.
#   Motivation for yes: blocks exist to be manipulated. are usually not static.
# no such things as arrays in userspace blocks
# userspace blocks can't refer to things in scope, shaped QAASM can
# userspace blocks can't consume, shaped QAASM can

# Choice:
# userspace rshape,lshape are [<reg>,<reg>]
# userspace block shape should match declaration order.
#   declare x: 2
#   declare y: 3
#   -> should have lshape [2,3]

# Userspace block:
#     rshape = [<reg>,<reg>], lshape is by declaration order.
#     I believe this prevents swapping by relabeling. Is that what I want?
#     If userspace blocks have swapping by relabeling, then permutations automatically cancel.

# example:
#   rshape = [<reg1>,<reg2>]   # lshape = [<reg1>,<reg2>]
#   declare <reg3>             # lshape = [<reg1>,<reg2>,<reg3>]
#   <reg3> += <reg1>
#   <reg1> -= <reg3>
#   zero <reg1>                # lshape = [<reg2>,<reg3>]

# Question: make swapping a primitive?
#   yes, can do this via algebra, but is more inefficient.
#   Helps give hints to any future qaasm compilers.




    # these should all return shaped QAASM blocks
    def block_create(parent, expn, dim):
        pass

    def block_consume(parent, name):
        pass

    def block_cast(parent, name):
        pass

    def block_userspace(parent, blq):
        pass



###############################

def block_add(*blocks):
    scales = [b.scale for b in blocks]
    out = QAASM()
    out.scale = sum(scales)
    col1 = [(s/out.scale)**(1/2) for s in scales]
    mat = [] # make matrix with col1 as the first column

    # substitution business


def block_mul(b1,b2):
    pass

def block_tensor(b1,b2):
    pass

def block_adjoint(b1,b2):
    pass

