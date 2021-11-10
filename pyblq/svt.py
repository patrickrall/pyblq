
import numpy as np
from numpy.polynomial import Chebyshev as Cheb
from .qsp import quantum_signal_processing, sanitize_polynomial
from .qac import *
from .expose import expose

def svt(qac, P, allow_scale=True):

    assert len(qac.scope_inputs) == 0
    assert len(qac.scope_outputs) == 0

    # convert to chebyshev polynomial
    if not isinstance(P, Cheb):
        P = P.convert(kind=Cheb)

    # adjust scale
    if qac.scale != 1:
        P = P(Cheb([0,qac.scale]))
        qac.scale = 1

    P = sanitize_polynomial(P)
    d = len(P.coef)-1


    out_qac = QAC()

    qac, regs = expose(qac, which=["declare","discard","maxmixed","zero"])

    out_qac.unnamed_inputs = []
    out_qac.unnamed_outputs = []
    for reg in qac.unnamed_inputs:
        if reg not in regs["declare"]+regs["maxmixed"]:
            out_qac.unnamed_inputs.append(reg)
            if d % 2 == 0:
                out_qac.unnamed_outputs.append(reg)
    if d % 2 == 1:
        for reg in qac.unnamed_outputs:
            if reg not in regs["discard"]+regs["zero"]:
                out_qac.unnamed_outputs.append(reg)

    # convert maxmixed into bell pairs
    maxmixed_partners = []
    for reg in regs["maxmixed"]:
        qft = np.zeros((reg.dim,reg.dim)).astype(complex)
        for i in range(reg.dim):
            for j in range(reg.dim):
                qft[i,j] = np.exp(i*j*2j*np.pi/reg.dim)/np.sqrt(reg.dim)

        partner = Register(reg.dim)
        partner.name_hints += reg.name_hints + ["maxmixed_partner"]
        qac.instrs = [
            # {"kind":"qac_declare", "reg":partner, "dim":partner.dim},
            {"kind":"qac_unitary", "reg":partner, "mat":qft.tolist()},
            {"kind":"qac_increment", "reg":reg, "expn":
                    {"kind":"register_expn", "register":partner}}
                ] + qac.instrs
        qac.unnamed_inputs.append(partner)
        qac.unnamed_outputs.append(partner)
        maxmixed_partners.append(partner)

    # run qsp
    scale, phase, phis, completion_mode = quantum_signal_processing(P,allow_scale=allow_scale)
    out_qac.scale = scale

    # if in FG mode need to take LCU of positive phases and negative phases
    # to cancel out imaginary part of the implemented polynomial
    # see 1806.01838 Corollary 18
    if completion_mode == "FG":
        lcu_reg = Register(2)
        lcu_reg.name_hints.append("svt_lcu_reg")
        out_qac.instrs.append({"kind":"qac_declare", "reg":lcu_reg, "dim":2})
        out_qac.instrs.append({"kind":"qac_unitary", "reg":lcu_reg, "mat":
            [[1/np.sqrt(2),1/np.sqrt(2)],[1/np.sqrt(2),-1/np.sqrt(2)]]})

        out_qac.instrs.append({"kind":"qac_unitary", "reg":lcu_reg, "mat":
            [[phase,0],[0,1/phase]]})
    else:
        out_qac.instrs.append({"kind":"qac_phase", "value": phase })


    flag = Register(2)
    flag.name_hints.append("svt_flag")
    out_qac.instrs.append({"kind":"qac_declare", "reg":flag, "dim":2})

    # initialize ancillae
    for reg in regs["declare"] + regs["maxmixed"] + maxmixed_partners:
        out_qac.instrs.append({"kind":"qac_declare", "reg":reg, "dim":reg.dim})

    # convenience function for checking if a bunch of stuff is zero
    def phase_terms(regs):
        out = []
        for reg in regs:
            out.append({"kind":"register_expn", "register":reg})
            out.append("==")
        out.append({"kind":"value_expn", "value":0j})
        return out

    Udagger = block_adjoint(block_copy(qac))
    for i in range(len(qac.unnamed_outputs)):
        Udagger.unnamed_inputs[i].substitute(qac.unnamed_outputs[i])
    for i in range(len(qac.unnamed_inputs)):
        Udagger.unnamed_outputs[i].substitute(qac.unnamed_inputs[i])

    # assemble the SVT circuit
    for i in range(d):
        if i % 2 == 0:
            out_qac.instrs += qac.instrs

            if i != d-1:
                out_qac.instrs.append({"kind":"qac_increment", "reg":flag,
                    "expn":{"kind":"boolean_expn", "terms":phase_terms(regs["zero"])}})

                if completion_mode == "FG":
                    out_qac.instrs.append({"kind":"qac_unitary", "reg":flag, "mat":
                        [[np.exp(-1j*phis[i]),0],[0,np.exp(1j*phis[i])]]})
                    out_qac.instrs.append({"kind":"qac_if", "cond":lcu_reg,
                        "instructions":[
                            {"kind":"qac_unitary", "reg":flag, "mat":
                                [[np.exp(2j*phis[i]),0],[0,np.exp(-2j*phis[i])]]}
                        ]})

                else:
                    out_qac.instrs.append({"kind":"qac_unitary", "reg":flag, "mat":
                        [[np.exp(-1j*phis[i]),0],[0,np.exp(1j*phis[i])]]})

                out_qac.instrs.append({"kind":"qac_increment", "reg":flag,
                    "expn":{"kind":"boolean_expn", "terms":phase_terms(regs["zero"])}})

        else:
            # U dagger
            out_qac.instrs += Udagger.instrs

            if i != d-1:
                out_qac.instrs.append({"kind":"qac_increment", "reg":flag,
                    "expn":{"kind":"boolean_expn", "terms":phase_terms(regs["declare"] + regs["maxmixed"] + maxmixed_partners)}})
                if completion_mode == "FG":
                    out_qac.instrs.append({"kind":"qac_unitary", "reg":flag, "mat":
                        [[np.exp(-1j*phis[i]),0],[0,np.exp(1j*phis[i])]]})
                    out_qac.instrs.append({"kind":"qac_if", "cond":lcu_reg,
                        "instructions":[
                            {"kind":"qac_unitary", "reg":flag, "mat":
                                [[np.exp(2j*phis[i]),0],[0,np.exp(-2j*phis[i])]]}
                        ]})
                else:
                    out_qac.instrs.append({"kind":"qac_unitary", "reg":flag, "mat":
                        [[np.exp(-1j*phis[i]),0],[0,np.exp(1j*phis[i])]]})
                out_qac.instrs.append({"kind":"qac_increment", "reg":flag,
                    "expn":{"kind":"boolean_expn", "terms":phase_terms(regs["declare"] + regs["maxmixed"] + maxmixed_partners)}})

    # final postselection
    if d % 2 == 0:
        for reg in regs["declare"] + regs["maxmixed"] + maxmixed_partners:
            out_qac.instrs.append({"kind":"qac_zero", "reg":reg})

    else:
        for reg in regs["discard"] + maxmixed_partners:
            out_qac.instrs.append({"kind":"qac_discard", "reg":reg})

        for reg in regs["zero"]:
            out_qac.instrs.append({"kind":"qac_zero", "reg":reg})

    if completion_mode == "FG":
        out_qac.instrs.append({"kind":"qac_unitary", "reg":lcu_reg, "mat":
            [[1/np.sqrt(2),1/np.sqrt(2)],[1/np.sqrt(2),-1/np.sqrt(2)]]})
        out_qac.instrs.append({"kind":"qac_zero", "reg":lcu_reg})

    out_qac.instrs.append({"kind":"qac_zero", "reg":flag})


    return out_qac


