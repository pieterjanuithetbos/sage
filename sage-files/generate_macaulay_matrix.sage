from sage.all import *
load("generate_f.sage")
load("generate_col_index.sage")
def generate_M(n, q, p, d, order_string):
    """
    input:
        q (int) characteristic of underlying field
        n (int) number of variables present
        p (list) of (QuadraticForm) a system of homogeneous quadratic equations,
        d (int) the desired degree of the output matrix
        order_string (String) the order in which the momomials should be indexed

    returns at tuple with:
        1. a macauley matrix of degree d, with (n d) columns and m(n d-2) rows
        2. a list representing the row indices (in the order specified by order_string)
        3. a list representing the columns indices (in degrevlex order)
    whose columns are index in lexicographical order
    """
    F = GF(q)
    R = PolynomialRing(F, n, 'x', order=order_string)
    x = R.gens()
    x_vec = vector(x)
    m = len(p)
    I = R.ideal([el**2 for el in x])
    Q = R.quotient(I)
    nb_rows_per_eq = Combinations(n, d-2).cardinality()
    # convert the system of equations to F/<x^2, ..., x_n^2>
    for eq in p:
        eq = eq - diagonal_matrix(eq.diagonal())
    # calculate dimensions of resulting matrix
    nb_cols = Combinations(n, d).cardinality()
    nb_rows = m*nb_rows_per_eq
    # print("De matrix heeft dimensies: ", nb_rows, nb_cols)
    M = matrix(F, nb_rows, nb_cols)

    # create row index
    row_index = []
    # for every polynomial equation
    for i in range(m):
        mono = R.monomials_of_degree(d-2)
        for i in range(len(mono)):
            current_element = mono[len(mono) - 1 - i]
            if current_element not in I:
                row_index.append((current_element, i))
    # print("de rij-index-lijst: \n \n", row_index)
    # print('nb_rows: \t', nb_rows, 'het huidige aantal: ', len(row_index))

    assert len(row_index) == nb_rows
    # create column index in lex order
    col_index = create_col_index(x, n, d)
    # print('col_index: \t', len(col_index), 'nb_cols: \t', nb_cols)
    assert len(col_index) == nb_cols
    # print("de column-index-lijst: \t", col_index)
    # fill the matrix
    for i in range(nb_rows):
        for j in range(nb_cols):
            poly_index = i // nb_rows_per_eq
            temp_poly = row_index[i][0]*(x_vec*p[poly_index]*x_vec)
            M[i, j] = temp_poly.monomial_coefficient(col_index[j])
    return M, row_index, col_index

m = 2
n = 5
o = 1
q = 2
d = 4
p = []
order = 'degrevlex'
for i in range(m):
    p.append(generate_f(n, o, q))
print(generate_M(n, q, p, d, order))