
import numpy as np
import os
import ctypes

def C_call(co_act_exp, co_act, co_rep_exp, co_rep, g, gene_disc_exp) :
    # Chargement la bibliothèque partagée contenant la fonction C, /!\ AJUSTER LE CHEMIN
    lib = ctypes.CDLL(os.path.abspath("./utils/comblicorn.c"))

    # definition of argument types for the C function
    lib.combnLicorn.argtypes=[
        ctypes.POINTER(ctypes.c_int),  # co_act_exp
        ctypes.POINTER(ctypes.c_int),  # ncoacts
        ctypes.POINTER(ctypes.c_int),  # co_rep_exp
        ctypes.POINTER(ctypes.c_int),  # ncoreps
        ctypes.POINTER(ctypes.c_int),  # gexp
        ctypes.POINTER(ctypes.c_int),  # nsamples
        ctypes.POINTER(ctypes.c_double),  # mae
        ctypes.POINTER(ctypes.c_int),  # iact
        ctypes.POINTER(ctypes.c_int),  # irep
    ]

    # definition of the values
    co_act_exp=co_act_exp.T
    nco_acts=len(co_act)
    co_rep_exp=co_rep_exp.T
    ncoreps=len(co_rep)
    g_exp=gene_disc_exp.iloc[g] 
    nsamples=len(g_exp)
    mae=np.full(nco_acts * ncoreps, -1)
    iact=np.full(nco_acts * ncoreps, -1)
    irep=np.full(nco_acts * ncoreps, -1)

    # C_call with right arguments
    lib.combnLicorn(
        co_act_exp.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        nco_acts.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        co_rep_exp.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        ncoreps.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        g_exp.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        nsamples.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        mae.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        iact.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        irep.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
    )

    # mae, iact, irep give the results of C function
    return (mae, iact, irep)