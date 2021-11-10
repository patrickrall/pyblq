
global register_count
register_count = 0
def reg():
    register_count += 1
    return "r_"+str(register_count)



class QAASM():
    def __init__(self):
        self.scale = 1
        self.instrs = []

        self.inputs = {}

        # how to do if/repeat statements?
        # how to do block algebra?
        # how to do consumes?
        # how to do arrays? <- solve later

        # how to do shaping?
        # maintain these externally when sequencing
        # only get used for block algebra
        # self.rshape = []
        # self.lshape = []
        # maybe this handles consumes too?
        #  |x               is like  inputs = {"x":2}, rshape = [], scope={"x":2}, lshape=["x"]
        #  ~ket(0) * |x     is like  inputs = {"x":2}, rshape = [], scope={}, lshape=[]

    @property
    def scope(self):
        scope = {}
        for key in self.inputs.keys():
            scope[key] = self.inputs[key]

        for instr in self.instrs:
            if instr["kind"] == "qaasm_declare":
                assert instr["reg"] not in scope
                scope[instr["reg"]] = instr["dim"]

            if instr["kind"] in ["qaasm_discard", "qaasm_zero"]:
                assert instr["reg"] in scope
                del scope[instr["reg"]]

        return scope

    def input(self, reg, dim):
        # TODO check that reg is never in scope
        self.inputs[reg] = dim

    def declare(self, reg, dim):
        assert reg not in self.scope:
        self.instrs.append({"kind":"qaasm_declare", "reg":reg, "dim":dim})

    def discard(self, reg, dim):
        if reg not in self.scope:
            self.input(reg,dim)
        self.instrs.append({"kind":"qaasm_discard", "reg":reg})

    def zero(self, reg, dim):
        if reg not in self.scope:
            self.input(reg,dim)
        self.instrs.append({"kind":"qaasm_zero", "reg":reg})

    # Boring:
    # unitary
    # phase

    def increment(self, reg, expn):
        # Maybe this spot should handle array references?
        # That way process_qaasm_expn's only responsibility is to parse the ast and try to pre-evaluate
        # and resolve {"kind":"array_expn", "regs":[...], "key":...} or something

        scope = self.scope
        assert reg in scope

        # put any dependencies of expn not in scope into inputs
        def process_expn(expn):
            # make sure references are in scope, collect array references
            # can't just add to inputs willy-nilly - what if that clashes?
            return expn

        self.instrs.append({"kind":"qaasm_increment", "reg":reg, "expn":expn})

    def if(self, reg, qaasm):
        # reg in scope
        # anything in qaasm.inputs not in self.inputs gets self.input'ed
        # qaasm.inputs should match qaasm.scope
        pass


    ######################

    def scalar_instr(self, qaasm):
        # anything in qaasm.inputs not in self.inputs gets self.input'ed
        # qaasm.lshape == qaasm.rshape == []
        # just append the instructions
        # and multiply by scale
        pass

    def init_instr(self, regs, qaasm):
        # anything in qaasm.inputs not in self.inputs gets self.input'ed
        # len(qaasm.lshape) == len(regs)
        # qaasm.rshape == []

        # multiply by scale

        # each linear symbol is either referenced or consumed
        # x,y <- |z @ |x
        #     inputs={"x":2,"z":2}, rshape=[], scope={"x":2,"z":2}, lshape=["z","x"]

        # appended:
        #   <nothing>
        # need to add:
        #

        # scope to regs
        # z -> x
        # x -> y

        # x <- x * ket(0)
        #     inputs={"x":2}, rshape=[], scope={"auto_1":2, "x":2}, lshape=["auto_1"]

        # lshape to regs
        # "auto_1" -> "x"
        # "x" -> discarded

        # appended:
        #   declare auto_1:2
        #   if x:
        #       declare tmp:2
        #       tmp += 1
        #       zero tmp
        # need to add:
        #   discard x
        #   declare x
        #   x += auto_1
        #   auto_1 -= x
        #   zero auto_1

        # just append the instructions
        # for each item in qaasm.scope, discard if needed
        # then determine permutation on remaining ones. implement permutation.

        # anything in regs can't be
        # append the instructions and multiply by scale
        # for reg in regs:
        #
        pass

    def assign_instr(self, reg, expn):
        # put any dependencies of expn not in scope into inputs

        # if qaasm.inputs doesn't contain reg:
        #    discard reg
        #    declare reg
        #    reg += expn

        # if it does:
        #    declare tmp
        #    tmp += expn
        #    discard reg
        #    declare reg
        #    reg += tmp
        #    tmp -= reg
        #    zero tmp
        pass

    def assign_undo_instr(self, reg, expn, undo):
        # e.g., say I know x < dim/2
        # x = x*2 undo x/2

        # put any dependencies of expn not in scope into inputs
        # put any dependencies of undo not in scope into inputs

        # expn and undo better depend on reg
        # maybe double check this and throw a nice error message

        # declare tmp
        # tmp += expn
        # x -= undo(but with tmp as x)
        # zero x
        # declare x
        # x += tmp
        # tmp -= x
        # zero tmp

        pass

    ########## The big problem with this approach is that it adds a bunch of needless swaps
    # some instructions create swaps that are really needed
    #  x @ y <- |y @ |x
    # others don't
    #  x <- ket(0)
    #  complex:
    #      declare auto_1
    #      discard x
    #      declare x
    #      x += auto_1
    #      auto_1 -= x
    #      zero auto_1
    #  simple:
    #      discard x
    #      declare x
    ##########
    #  Is this something for the optimizer?
    #  Like, is this just as hard to figure out now as it is to figure out later?
    # I have a decent sense for how hard it is now, but how hard is it later?

    ####### TODO: perhaps fix with a relabel instruction!!!
    # need to prevent 'active relabel':
    # if cond:
    #    relabel x tmp
    #    relabel y x
    #    relabel tmp y

    # other ideas: substitute



def block_add(*blocks):
    scales = [b.scale for b in blocks]
    out = QAASM()
    out.scale = sum(scales)
    col1 = [(s/out.scale)**(1/2) for s in scales]
    mat = [] # make matrix with col1 as the first column

    # declare ctrl : len(blocks)
    # unitary ctrl mat
    # if ctrl == 0: blocks[0].instructions
    # if ctrl == 1: blocks[1].instructions
    # if ctrl == 2: blocks[2].instructions
    # unitary ctrl mat.adjoint
    # zero ctrl

    # blocks need to have the same shape
    # if blocks have regs in different order/names, need to permute to standardize
    # need to also come up with other names referenced registers
    # guaranteed no consumes: if it's an input but not an output, it had to be in rshape

    # e.g. {"x":2, "z":2}, ["x"], {"y":2, "z":2}, ["y"]
    # +    {"x":2, "y":2}, ["y"], {"z":2, "x":2}, ["z"]

    # this sure is a lot of headache, that wouldn't be present if all the registers were already unique.
    # if all the regs start out unique then all you need to do is merge them.

    # alternative:
    #   references = {"x":<reg>, "y":<reg> }
    #   consumes = {"z":<reg>}
    #   inputs = [<reg>,<reg>]
    #   outputs = [<reg>,<reg>]

    pass

def block_mul(b1,b2):
    pass

def block_tensor(b1,b2):
    pass

def block_adjoint(b1,b2):
    pass


# types of QAASM:
#  if/repeat blocks     don't have a shape, interact with an external scope.  Can't consume
#  intermediate blocks  have shape,         interact with an external scope.  Can consume
#  Blq blocks           have shape,         have their own scope only         Can't consume


