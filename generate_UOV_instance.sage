import time 

def generate_f(n, o, q):
    '''
    generates a random quadratic equation over F_q in n variables. 
    The coefficient matrix use deglex order
    '''
    F = GF(q)
    
    R = PolynomialRing(F, n, 't', order='deglex')
    
    # construct a quadratic polynomial over F_q, such that every term has a vinegar variable
    f_k = sum(R.gens()[i]*sum(F.random_element()*R.gens()[j] for j in range(n)) for i in range(o, n)) 
    size = len(R.gens())
    gram = Matrix(F, size, size)
    for i in range(size):
        quadratic_terms = {R.gens()[i]: 2}
        gram[i, i] = f_k.coefficient(quadratic_terms)
        for j in range(i+1, n):
            other_terms = {R.gens()[i]: 1, R.gens()[j]: 1}
            gram[i, j] = f_k.coefficient(other_terms)
    
    # Other possibility: use the QuadraticForm and construct the matrix of the polar form
    # drawback: doesn't make sense in characteristic two
    
    x = vector(R.gens())
    assert (x * gram * x == f_k)
    
    return gram

def KeyGen(n, m, o, q):
    '''
    returns a tuple (pk, sk) containing the public key and the secret key.
    The public key is a list of m equations.
    The secret is again a tuple (FF, M), where
            - FF is the central map and
            - M an invertible linear transformation.
    '''
    F = GF(q)
    FF, pk = [], []
    # pick the quadratic equations at random
    for i in range(m):
        FF.append(generate_f(n, o, q))
        
    # choose an invertible linear transform F_n -> F_n
    M = GL(n, F).random_element().matrix()
    sk = (FF, M)
    
    # construct the public key
    for i in range(len(FF)):
        pk.append(M.transpose()*FF[i]*M)

    return pk, sk

def sign_uov(n, m, o, q, sk, t):

    '''
    returns the UOV signature as a vector of length n
    '''
    assert (len(t) == m) 
    F = GF(q)
    R = PolynomialRing(F, n, 'x')
    x = vector(R.gens())

    # we skip the hashing part
    
    i = 1
    while true:
        try:
            sol_temp = solve_FF(n, m, o, q, sk, t, F, R, x)
            break
        except ValueError:
            i+=1
            # system has no solutions
            print(f'Hélas ... Too much vinegar. Poging {i}')
            
    # create s'
    sol = [sol_temp[i] for i in range(o)]
    sol.extend([x[i] for i in range(o, n)])
    
    # create s
    signature = sk[1].inverse()*vector(sol)

    return signature

def solve_FF(n, m, o, q, sk, t, F, R, x):
    
    # pick random values for the vinegar variables
    for i in range(n-o):
        x[i+o] = F.random_element()

    # solve the system F(s') = t (m linear equations in o variables, m <= o)
    FF, M = sk
    system = [x*FF[i]*x - t[i] for i in range(m)]
    A = matrix([ [eq.coefficient(v) for v in x[0:o]] for eq in system])
    b = vector([-eq.constant_coefficient() for eq in system])
    sol = A.solve_right(b)
    
    return sol

def verify(n, m, o, q, PP, signature):
    F = GF(q)
    R = PolynomialRing(F, n, 'x')
    x = vector(R.gens())
    t_bar = []
    t_bar.append(signature*PP[i]*signature for i in range(len(PP)))
    return t_bar

def generate_UOV_variables(n, m, o, q):
    v = n - o
    F = GF(q)

    # generate public and private keys
    pk, sk = KeyGen(n, m, o, q)
    FF, T = sk
    P = PolynomialRing(F, n, 't', order='deglex')
    gens_vector = vector(P.gens())
    oil = gens_vector[0:o]
    vinegar = gens_vector[o:]
    scrabbled = T*gens_vector

    O_bar = block_matrix([
        [identity_matrix(F, o), zero_matrix(F, o, n-o)]
        ])

    # compute the oil space 
    O = (T.inverse()*O_bar.transpose()).transpose()
    V_bar = block_matrix([
        [zero_matrix(F, v, o), identity_matrix(F, v)]
        ])
    V = (T.inverse()*V_bar.transpose()).transpose()

    return pk, sk, P, oil, vinegar, O, V, gens_vector