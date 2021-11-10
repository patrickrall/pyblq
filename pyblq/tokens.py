
import string

#################
# location parsing

# how to make this work with this argument system?
# new lines must be explicit
# Blq("""[x:2 -> x:2]
# x += """, my_block, """
# x += 1
# """)

# format:
# line,col
# line,col-col
# line,col-line,col

def make_location(l1,c1,l2,c2):
    s = ""
    if l1 == l2:
        s += str(l1)+","+str(c1)
        if c1 != c2:
            s += "-"+str(c2)
    else:
        s += str(l1)+","+str(c1)
        s += "-"+str(l2)+","+str(c2)
    return s

def parse_location(s):
    if s.count(",") == 2:
        start, end = s.split("-")
        l1,c1 = start.split(",")
        l2,c2 = end.split(",")
        return int(l1), int(c1), int(l2), int(c2)
    else:
        assert s.count(",") == 1
        line, s = s.split(",")

        if s.count("-") == 0:
            return int(line), int(s), int(line), int(s)
        else:
            assert s.count("-") == 1
            c1, c2 = s.split("-")
            return int(line), int(c1), int(line), int(c2)

def merge_locations(loc1,loc2):
    l11,c11,l12,c12 = parse_location(loc1)
    l21,c21,l22,c22 = parse_location(loc2)

    assert l11 < l21 or c11 <= c21 # bounds
    assert l12 < l22 or c12 <= c22 # checks

    return make_location(l11, c11, l22, c22)

def print_location(loc, args, pad=3):
    FAIL = '\033[91m'
    ENDC = '\033[0m'

    l1,c1,l2,c2 = parse_location(loc)

    s = ""
    for arg in args:
        if isinstance(arg,str): s += arg
        else: s += "?"
    data = s.split("\n")

    # newline characters show up as the line number being highlighted
    if c1 == len(data[l1-1])+1:
        l1 += 1
        c1 = 0
    if c2 == len(data[l2-1])+1:
        l2 += 1
        c2 = 0

    start = max(1,l1-pad)
    end = min(len(data),l2+pad)

    data = data[start-1:end]

    s = ""
    i = start
    for line in data:
        line = "%" + line # line number character

        if i == l1:
            line = line[:c1]+FAIL+line[c1:]

        if i == l2:
            if i == l1:
                line = line[:c2+len(FAIL)+1] + ENDC + line[c2+len(FAIL)+1:]
            else:
                line = line[:c2+1] + ENDC + line[c2+1:]

        s += line.replace("%", str(i).zfill(len(str(end)))+" ", 1) + "\n"
        i += 1

    return s

def error_at(loc,args):
    return " at "+loc+"\n"+print_location(loc,args)

#################
# generating a list of tokens

"""
identifier
literal
extern
indent, dedent
static:
  declare uncompute discard repeat undo scalar pass
  if elif else
  ( ) [ ] <- \n : | ~ = ket( , += -=
  + - * ** / % @ @@ == != > < >= <=
comment (# comments only)
"""


def tokenize(args):

    # line, col, index in args, index in the string
    lcij = {0:  (1, 1, 0, 0) }  # a dictionary so I change it in a subroutine
    tokens = []


    def charsleft():
        return lcij[0][2] < len(args)

    def advance(n):
        line, col, idx, jdx = lcij[0]
        for i in range(n):
            if not isinstance(args[idx],str):
                idx += 1
                jdx = 0
                col += 1
                continue
            if args[idx][jdx] == "\n":
                line += 1
                col = 0
            col += 1
            jdx += 1
            if jdx >= len(args[idx]):
                idx += 1
                jdx = 0
        lcij[0] = (line,col,idx,jdx)

    def match(s):
        _,_,idx,jdx = lcij[0]

        for c in s:
            if not idx < len(args): return False
            if not isinstance(args[idx],str): return False
            if c != args[idx][jdx]: return False
            jdx += 1
            if jdx >= len(args[idx]):
                idx += 1
                jdx = 0
        return True

    def current():
        _,_,idx,jdx = lcij[0]
        if not isinstance(args[idx],str): return args[idx]
        else: return args[idx][jdx]

    def matchExtern():
        _,_,idx,jdx = lcij[0]
        return not isinstance(args[idx],str)

    # call this before advancing if you want it to be accurate.
    def loc(n):
        line, col, i, j = lcij[0]
        return make_location(line, col, line, col+n-1)

    # specifically for new lines - also parses subsequent indent
    def advanceNl():
        assert match("\n")
        s = "\n"
        while True:
            if match(s+" "): s += " "
            elif match(s+"\t"): s += "\t"
            else: break
        tokens.append({ "kind": "\n", "loc": loc(len(s)), "indent": s })
        advance(len(s))

    ################ Parsing loop

    linecomment = False
    while charsleft():

        ############## Comments
        if linecomment:
            if match("\n"):
                linecomment = False
                advanceNl()
            else:
                advance(1)
            continue

        if match("#"):
            linecomment = True
            continue

        ############################ Statics

        if match("\n"):
            advanceNl()
            continue

        if match(" ") or match("\t"):
            advance(1)
            continue

        statics = ["declare", "uncompute", "discard", "repeat", "undo", "scalar", "pass",
                   "if", "elif", "else", "(", ")", "[", "]", "<-", ":", "|", "~", "=" ,"ket(" ,",",
                   "+=","-=","+", ":", "-", "*", "**", "/", "%", "@", "@@", "==", "!=", ">", "<", ">=", "<="]

        # take longest one that matches
        found = False
        for s in reversed(sorted(statics)):
            if not match(s): continue

            tokens.append({
                "kind": s,
                "loc": loc(len(s))
            })

            advance(len(s))
            found = True
            break

        if found:
            continue


        ############################## Literals / Identifiers / Extern

        if matchExtern():
            _,_,idx,jdx = lcij[0]
            tokens.append({ "kind": "extern", "loc": loc(1), "value": args[idx] })
            advance(1)
            continue

        numbers = string.digits
        letters = string.ascii_letters
        lowercase = string.ascii_lowercase

        if current() in numbers+".":
            # parse number
            has_dot = current() == "."
            has_j = False

            s = current()
            start_loc = loc(1)
            end_loc = loc(1)
            advance(1)

            while charsleft() and current() in numbers+".j":
                if has_dot:
                    raise SyntaxError("Unexpected character 'j'" + error_at(loc(1),args))

                if current == ".":
                    if has_dot:
                        raise SyntaxError("Unexpected character '.'" + error_at(loc(1),args))
                    has_dot = True
                if current == "j":
                    has_j = True

                s += current()
                end_loc = loc(1)
                advance(1)

            tokens.append({
                "kind": "literal",
                "value": s,
                "loc": merge_locations(start_loc,end_loc)
            })

        elif current() in letters:

            s = current()
            start_loc = loc(1)
            end_loc = loc(1)
            advance(1)

            while charsleft() and current() in numbers+letters+"_":
                s += current()
                end_loc = loc(1)
                advance(1)

            tokens.append({
                "kind": "identifier",
                "value": s,
                "loc": merge_locations(start_loc,end_loc)
            })

        else:
            raise SyntaxError("Unexpected character '"+current()+"'"+error_at(loc(1),args))

    ###### Post processing indents

    # Honestly, the indententation and new line handling is pretty kludgey right now.
    # Don't tweak it or else it will break, since it interacts with make_ast.
    # Best develop a good idea for how it SHOULD work and then reimplement from scratch.

    outTokens = []

    indentLevel = None
    prvNl = False
    for i in range(len(tokens)):
        token = tokens[i]

        if token["kind"] != "\n":
            outTokens.append(token)
            prvNl = False
            continue

        # skip empty lines
        if prvNl:
            continue
        # if i+1 == len(tokens) or tokens[i+1]["kind"] == "\n":
        #     continue
        prvNl = True

        l1,c1,_,_ = parse_location(token["loc"])
        nlLoc = make_location(l1,c1,l1,c1)
        indentLoc = make_location(l1+1,1,l1+1,len(token["indent"])) # indents are on next line
        dedentLoc =  make_location(l1,c1+1,l1,c1+len(token["indent"])-1) # dedents are on this line

        if indentLevel is None:
            indentLevel = [token["indent"]]
            outTokens.append({"kind":"\n", "loc": nlLoc })
            continue

        # indent
        if indentLevel[-1] in token["indent"] and\
                token["indent"].index(indentLevel[-1]) == 0\
                and len(token["indent"]) > len(indentLevel[-1]):
            indentLevel.append(token["indent"])
            outTokens.append({"kind":"\n", "loc": nlLoc })
            outTokens.append({"kind":"indent", "loc": indentLoc })
            continue

        good = False
        for i in range(len(indentLevel)):
            if token["indent"] == indentLevel[-(i+1)]:
                outTokens.append({"kind":"\n", "loc": nlLoc })
                for j in range(i):
                    outTokens.append({"kind":"dedent", "loc": dedentLoc })
                    # duplicate the new line token so that every instruction ends on a new line
                    outTokens.append({"kind":"\n", "loc": nlLoc })
                    indentLevel.pop()
                good = True
                break

        if not good:
            raise SyntaxError("Indentation does not match any previous indentation level"+error_at(loc(0),args))


    # adjust final indentation level
    finalLoc = outTokens[-1]["loc"]
    if outTokens[-1]["kind"] != "\n":
        outTokens.append({"kind":"\n", "loc": finalLoc })

    if indentLevel is not None:
        while len(indentLevel) > 1:
            outTokens.append({"kind":"dedent", "loc": finalLoc })
            outTokens.append({"kind":"\n", "loc": finalLoc })
            indentLevel.pop()

    if outTokens[-1]["kind"] == "\n": outTokens.pop()

    return outTokens
