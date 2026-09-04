from sage.all import *
load("key_gen.sage")
# SET VARIABLES
n, m, o, q = 5, 4, 2, 2

def generate_UOV_variables(n, m, o, q):
    v = n - o
    F = GF(q)
    pk, sk = KeyGen(n, m, o, q)
    FF, T = sk
    P = PolynomialRing(F, n, 't', order='degrevlex')
    gens_vector = vector(P.gens())
    oil = gens_vector[0:o]
    vinegar = gens_vector[o:]
    # print('oil and vinegar \t', oil, vinegar)
    scrabbled = T*gens_vector
    # print('scrabbled variables \t',scrabbled)

    O_bar = block_matrix([
        [identity_matrix(F, o), zero_matrix(F, o, n-o)]
        ])

    O = (T.inverse()*O_bar.transpose()).transpose()
    V_bar = block_matrix([
        [zero_matrix(F, v, o), identity_matrix(F, v)]
        ])
    V = (T.inverse()*V_bar.transpose()).transpose()
    print('O \n', O)
    return pk, sk, P, oil, vinegar, O, V, gens_vector

def experiment_1_sampling_vectors(n, m, o, q):
    pk, sk, P, oil, vinegar, O, V, gens_vector = generate_UOV_variables(n, m, o, q)
    F = GF(q)
    v = n - o
    # EXPERIMENT 1
    test_it = q^n
    counter_O, counter_V, counter_random = 0, 0, 0
    counter_V2 = 0
    counter_random_2 = 0

    # convert the public key to polynomials
    pk_poly = []
    for coef_matrix in pk:
        polynomial = P(gens_vector*coef_matrix*gens_vector)
        pk_poly.append(polynomial)
    for i in range(test_it):
        random_O_vector = vector([F.random_element() for i in range(o)])*O
        random_V_vector = vector([F.random_element() for i in range(v)])*V
        random_vector = vector([F.random_element() for i in range(n)])
        if all((poly(list(random_O_vector)) == 0) for poly in pk_poly):
            counter_O += 1
        if all((poly(list(random_V_vector)) == 0) for poly in pk_poly):
            counter_V += 1
        if all((poly(list(random_vector)) == 0) for poly in pk_poly):
            counter_random += 1
    print('oil vectors  \t \t', round(counter_O/test_it, 2))
    print('vinegar vectors \t', round(counter_V/test_it, 2))
    print('random vectors \t \t', round(counter_random/test_it, 2))

    # print('random O vector \t', random_O_vector, 'evaluation \t',  p1(list(random_O_vector)) == p2(list(random_O_vector)) == 0)
    # print('random V vector \t', random_V_vector, 'evaluation \t',  p1(list(random_V_vector)) == p2(list(random_V_vector)) == 0)
    # print('random vector \t \t', random_vector, 'evaluation \t',  p1(list(random_vector)) == p2(list(random_vector)) == 0

experiment_1_sampling_vectors(n, m, o, q)
