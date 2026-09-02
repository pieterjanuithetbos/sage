from sage.all import *

def generate_UOV_variables():
    # test
    n, m, o, q = 5, 4, 2, 2
    v = n - o
    F = GF(q);

    rank_T = 5
    T = matrix(F, rank_T, rank_T, [
        [1, 0, 1, 0, 3],
        [3, 0, 2, 0, 1],
        [2, 3, 0, 1, 2],
        [1, 0, 2, 3, 1],
        [0, 2, 3, 5, 1],
        ])
    assert rank(T) == rank_T

    P = PolynomialRing(F, n, 'x', order='deglex')
    gens_vector = vector(P.gens())
    oil = gens_vector[0:2]
    vinegar = gens_vector[2:]
    print('oil and vinegar \t', oil, vinegar)
    scrabbled = T*gens_vector
    print('scrabbled variables \t',scrabbled)
    im = [scrabbled[i] for i in range(len(scrabbled))]
    phi = P.hom(im)
    print(phi)

def generate_quadratic_forms():
    F_1 = [
        0, 0, 0, 1, 0,
           0, 0, 0, 0,
              1, 0, 0,
                 0, 0,
                    0
    ]
    F_2 = [
        0, 0, 0, 0, 0,
           0, 1, 1, 0,
              0, 1, 0,
                 0, 0,
                    1
    ]
    F_3 = [
        0, 0, 0, 1, 0,
           0, 0, 0, 0,
              1, 0, 1,
                 0, 0,
                    1
    ]
    F_4 = [
        0, 0, 0, 1, 1,
           0, 0, 0, 1,
              1, 1, 1,
                 0, 1,
                    1
    ]
    f1 = QuadraticForm(F, 5, F_1)
    f2 = QuadraticForm(F, 5, F_2)
    f3 = QuadraticForm(F, 5, F_3)
    f4 = QuadraticForm(F, 5, F_4)

    print('f1 \t', f1.polynomial())
    p1 = phi(f1.polynomial())
    print('p1 \t', p1)
    print('f2 \t', f2.polynomial())
    p2 = phi(f2.polynomial())
    p3 = phi(f3.polynomial())
    p4 = phi(f4.polynomial())
    p = [p1, p2, p3, p4]
    p_sequence = Sequence(p)

    print('p2 \t', p2) # checked manually, ok!
    print("type of p's", type(p2))


def generate_spaces():
    O_bar = matrix(F, o, n, [
       [1, 0, 0, 0, 0],
       [0, 1, 0, 0, 0],
                   ])
    O = (T.inverse()*O_bar.transpose()).transpose()

    V_bar = matrix(F, v, n, [
        [0, 0, 1, 0, 0],
        [0, 0, 0, 1, 0],
        [0, 0, 0, 0, 1],
        ])
    V = (T.inverse()*V_bar.transpose()).transpose()

    print('O \n', O)


def test():
    test_it = 1000
    counter_O, counter_V, counter_random = 0, 0, 0
    counter_V2 = 0
    counter_random_2 = 0

    for i in range(test_it):
        random_O_vector = vector([F.random_element() for i in range(o)])*O
        random_V_vector = vector([F.random_element() for i in range(v)])*V
        random_vector = vector([F.random_element() for i in range(n)])
        if all((poly(list(random_O_vector)) == 0) for poly in p):
            counter_O += 1
        if all((poly(list(random_V_vector)) == 0) for poly in p):
            counter_V += 1
        if all((poly(list(random_vector)) == 0) for poly in p):
            counter_random += 1
    print('oil vectors  \t \t', round(counter_O/test_it, 2))
    print('vinegar vectors \t', round(counter_V/test_it, 2))
    print('random vectors \t \t', round(counter_random/test_it, 2))

    # print('random O vector \t', random_O_vector, 'evaluation \t',  p1(list(random_O_vector)) == p2(list(random_O_vector)) == 0)
    # print('random V vector \t', random_V_vector, 'evaluation \t',  p1(list(random_V_vector)) == p2(list(random_V_vector)) == 0)
    # print('random vector \t \t', random_vector, 'evaluation \t',  p1(list(random_vector)) == p2(list(random_vector)) == 0
