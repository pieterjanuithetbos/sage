def KeyGen(n, m, o, q):
    '''
    returns a tuple (pk, sk) containing the public key and the secret key.
    The secret is again a tuple (FF, M), where
            - FF is the central map and
            - M an invertible linear transformation.
    '''
    F = GF(q)
    FF = []
    # pick the quadratic equations at random
    for i in range(m):
        FF.append(generate_f(n, o, q))
    # print("FF[2] \n", FF[2])
    # choose an invertible linear transform F_n -> F_n
    M = GL(n, F).random_element().matrix()
    sk = (FF, M)
    # construct the public key
    PP = []
    for i in range(len(FF)):
        PP.append(M.transpose()*FF[i]*M)
    # print("PP[2] \n", PP[2])
    pk = PP
    # print(type(pk), type(sk))
    return pk, sk