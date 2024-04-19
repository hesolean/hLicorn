
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
def C_call(coactexp,corepexp,geneDiscExp):
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
        ctypes.POINTER(ctypes.c_double),  # result
        ctypes.POINTER(ctypes.c_int),  # iact
        ctypes.POINTER(ctypes.c_int),  # irep
    ]

    # Définition des valeurs des arguments
    coactexp = np.array(coactexp)
    ncoacts = np.array([len(coactexp)], dtype=np.int32)
    corepexp = np.array(corepexp)
    ncoreps = np.array([len(corepexp)], dtype=np.int32)
    gexp = np.array(geneDiscExp) 
    nsamples = np.array([len(gexp)], dtype=np.int32)
    result = np.zeros((ncoacts[0], ncoreps[0]), dtype=np.float64)
    iact = np.zeros((ncoacts[0], ncoreps[0]), dtype=np.int32)
    irep = np.zeros((ncoacts[0], ncoreps[0]), dtype=np.int32)

    # Appel de la fonction C avec les arguments appropriés
    lib.combnLicorn(
        coactexp.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        ncoacts.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        corepexp.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        ncoreps.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        gexp.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        nsamples.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        result.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        iact.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        irep.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
    )

    # result, iact, irep contiennent les résultats retournés par la fonction C
    return (result, iact, irep)