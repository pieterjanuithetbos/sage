from sage.all import *
attach("../sage-files/generate_macaulay.sage")
def wedge_attack(p):
    """
    input:
        p a polynomial sequence

    """
    # build Macaulay matrix of degree o
    p.Macaulay_matrix()
    # compute the kernel vector
    # build the oil space from the kernel vector
    v = m.right_kernel()