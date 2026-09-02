def generate_f(n, o, q):
    '''
    generates a random quadratic equation over F_q in n variables.
    q has to be odd.
    '''
    F = GF(q)
    # define the polynomial ring where the variables live
    R = PolynomialRing(F, n, 't')
    # construct a quadratic polynomial over F_q, such that every term has a vinegar variable
    f_k = sum(R.gens()[i]*sum(F.random_element()*R.gens()[j] for j in range(n)) for i in range(o, n)) 
    # print(type(f_k))
    # print(f_k.coefficients())
    size = len(R.gens())
    gram = Matrix(F, size, size)
    for i in range(size):
        quadratic_terms = {R.gens()[i]: 2}
        gram[i, i] = f_k.coefficient(quadratic_terms)
        for j in range(i+1, n):
            other_terms = {R.gens()[i]: 1, R.gens()[j]: 1}
            gram[i, j] = f_k.coefficient(other_terms)
    

    # Q = QuadraticForm(f_k);
    # # gets matrix of polar form
    # gram = Q.matrix() * F(2)**(-1);
    
    # further checks
    x = vector(R.gens())
    assert (x * gram * x == f_k)
    return gram