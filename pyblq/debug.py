
from .tokens import print_location

# this needs args in order to pretty print the locations of any symbols it finds
def recur_print(t, args, d=None):

    if d is None:
        print(recur_print(t,args,d=0))
        return
    tab = "  "

    if d > 100: return "recursion limit \n"

    if not isinstance(t,dict):
        if isinstance(t,list):
            if len(t) == 0: return "[]\n"
            s = "[\n"
            for elem in t:
                s += tab*(d) + recur_print(elem,args, d+1)
            s += tab*d + "]\n"
            return s
        else:
            return repr(t) + "\n"

    if "kind" in t:
        if t["kind"] == "identifier": return t["value"]+"\n"
        if t["kind"] == "literal": return t["value"]+"\n"

        s = t["kind"]
    else:
        s = "<dict>"

    if "loc" in t:
        s += " at " + t["loc"]
        s += "\n"+tab*d
        s += print_location(t["loc"],args,pad=0).replace("\n","\n"+tab*d)[:-len(tab)*d]
    else:
        s += "\n"

    for key in t.keys():
        if key == "kind": continue
        if key == "loc": continue
        s += tab*d + '"'+key+'": ' + recur_print(t[key],args,d+1)

    return s
