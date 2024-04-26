
import numpy as np
import ctypes
    
'''on appelle le bout de programme en C, A VERIFIER
x= .C("combnLicorn",
as.integer(t(coactexp)),as.integer(length(coact)),#expression and number of coactivators
as.integer(t(corepexp)),as.integer(length(corep)),#expression and number of corepressors
as.integer(geneDiscExp[g,]),      as.integer(ncol(geneDiscExp)),    #expression of gene, number of samples
as.double(rep(-1,(length(coact))*(length(corep)))),# vector to store MAE results
as.integer(rep(-1,(length(coact))*(length(corep)))),# vector to store index of coactivator
as.integer(rep(-1,(length(coact))*(length(corep))))# vector to store index of corepressor
)
'''
def C_call(coactexp, coact, corepexp, corep, g, geneDiscExp):
    # Chargement la bibliothèque partagée contenant la fonction C, /!\ AJUSTER LE CHEMIN
    lib = ctypes.CDLL("./library.so")

    # Définition des types des arguments de la fonction C
    lib.combnLicorn.argtypes = [
        ctypes.POINTER(ctypes.c_int),  # coactexp
        ctypes.POINTER(ctypes.c_int),  # ncoacts
        ctypes.POINTER(ctypes.c_int),  # corepexp
        ctypes.POINTER(ctypes.c_int),  # ncoreps
        ctypes.POINTER(ctypes.c_int),  # gexp
        ctypes.POINTER(ctypes.c_int),  # nsamples
        ctypes.POINTER(ctypes.c_double),  # mae
        ctypes.POINTER(ctypes.c_int),  # iact
        ctypes.POINTER(ctypes.c_int),  # irep
    ]

    # Définition des valeurs des arguments
    coactexp = coactexp.T
    ncoacts = len(coact)
    corepexp = corepexp.T
    ncoreps = len(corep)
    gexp = geneDiscExp.iloc[g] 
    nsamples = len(gexp)
    mae = np.full(ncoacts * ncoreps, -1)
    iact = np.full(ncoacts * ncoreps, -1)
    irep = np.full(ncoacts * ncoreps, -1)

    # Appel de la fonction C avec les arguments appropriés
    lib.combnLicorn(
        coactexp.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        ncoacts.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        corepexp.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        ncoreps.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        gexp.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        nsamples.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        mae.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        iact.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        irep.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
    )

    # result, iact, irep contiennent les résultats retournés par la fonction C
    return (mae, iact, irep)