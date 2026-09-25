import random 
import math
import itertools
import numpy as np

# copied from Ran's implementation
def Upper(M):
    assert M.nrows() == M.ncols()
    n = M.nrows()
    F = M.base_ring()
    return matrix(F, n, n, lambda i, j: M[i, j] if i == j else (M[i,j] + M[j,i] if i < j else F(0)))

# own work
def from_gram_matrix_to_poly_sequence(pk, P):
    res_list = []
    x_vec = vector(P.gens())
    for mat in pk:
        res_list.append(x_vec*mat*x_vec)
    return Sequence(res_list)

def create_col_index(x, n, d, col_temp=None, degree_index=0):
    """
    generate all the 'boolean' monomials of total degree d, out of n variables, in deglex ordering
    """
    if col_temp is None:
        col_temp = 1
    if n == d:
        res = 1
        for el in x:
            res *= el
        return res
    col_index = []
    
    # base case: the total degree of the momomial is d
    if len(col_temp.factor()) == d:
        return [col_temp]
        
    # extend the solutions by multiplying with any monomial with higher index
    ext_sol = []
    for extension_index in range(degree_index, min(n - d + degree_index + 1, n)): # TODO understand the index issue
        temp_sol = col_temp*x[extension_index]
        ext_sol.append((temp_sol, extension_index))
        
    # call creat_col_index recursively
    for sol in ext_sol:
        res = create_col_index(x, n, d, sol[0], sol[1]+1)
        if [res] != []:
            col_index += res

    return col_index


def generate_M(n, q, p, d, order_string, educational_implementation=True, prune = False, can_print = False):
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
        3. a list representing the columns indices (in deglex order)
    whose columns are index in lexicographical order
    """
    F = GF(q)
    R = PolynomialRing(F, n, 'x', order=order_string)
    x = R.gens()
    x_vec = vector(x)
    o = d
    v = n - o
    m = len(p)
    I = R.ideal([el**2 for el in x])
    Q = QuotientRing(R, I)

    # convert the system of equations to F/<x^2, ..., x_n^2>
    polynomials = from_gram_matrix_to_poly_sequence(p, Q)
    for i in range(len(p)):
       p[i] = p[i] - diagonal_matrix(p[i].diagonal())
    qs = [Upper(p[i]) for i in range(len(p))]
    # alternatively, consider as elements of quotient ring
    p_in_Q = Sequence([Q(polynomials[i]) for i in range(len(polynomials))])

    """
    using the built-in seems not possible, because it requires a function is_homogeneous() to be defined on the 
    structure where the polynomials live in. Unfortunately, that function is not defined in the above defined 
    quotient ring Q
    
    tested the code: M_check, row_index_check, col_index_check = p_in_Q.macaulay_matrix(o-2, homogeneous=True, variables=None, return_indices=True)
    """
    
    full = set(range(n))
    time_lookup_start = time.time()
    lookup_I = {I: tuple(sorted(full-set(I))) for I in itertools.combinations(range(n), v)}
    print("Time to calculate the lookup dict: \t \t", time.time() - time_lookup_start)
    # instead, we use Ran's implementation
    time_idx_rows = time.time()
    idx = {I: i for i, I in enumerate(itertools.combinations(range(n), o))}
    rows = [(J, k) for J in itertools.combinations(range(n), v + 2) for k in range(m)]
    print("Time to calculate the idx and row dicts \t", time.time() - time_idx_rows)

    # prune
    if prune:
        kept_rows = set(random.sample(range(m * math.comb(n, v + 2)), math.comb(n, v)))
        rows = [r for i, r in enumerate(rows) if i in kept_rows]

    entries = dict()
    for r, (J, k) in enumerate(rows):
        entries.update(
            {(r, idx[lookup_I[J[:x] + J[x+1:y] + J[y+1:]]]) : qs[k][J[x]][J[y]] for x in range(v+2) for y in range(x+1, v+2)})

    M_check = matrix(F, len(rows), math.comb(n, v), entries)

    # build col index
    time_col_index_start = time.time()
    col_index_check = []
    for comb in itertools.combinations(range(n), o):
        temp = 1
        for k in comb:
            temp *= x[k]
        col_index_check.append(temp) 
    print("Time to compute the col_index: \t \t", time.time() - time_col_index_start)

    if educational_implementation:
        # calculate dimensions of resulting matrix
        nb_cols = math.comb(n, d)
        nb_rows_per_eq = math.comb(n, d-2) 
        nb_rows = m*nb_rows_per_eq 
    
        # create row index
        row_index = []
        for i in range(m):
            mono = R.monomials_of_degree(d-2)
            for j in range(len(mono)):
                current_element = mono[len(mono) - j - 1] # append the elements in reverse order
                if current_element not in I:
                    row_index.append((current_element, i))
        assert len(row_index) == nb_rows
        
        # create column index in lex order
        col_index = create_col_index(x, n, d)
        assert len(col_index) == nb_cols
    
        # prune the matrix
        if prune:
            is_rectangular = nb_cols > nb_rows
            if is_rectangular:
                selected_rows = random.sample(range(nb_rows), nb_cols)
                new_row_index = []
                for i in range(nb_rows):
                    new_row_index[i] = row_index[selected_rows[i]] 
                row_index = new_row_index
                assert len(new_row_index) == nb_cols
            
        # fill the matrix
        M = matrix(F, len(row_index), nb_cols)
        for i in range(nb_rows):
            for j in range(nb_cols):
                poly_index = i // nb_rows_per_eq
                temp_poly = row_index[i][0]*(x_vec*p[poly_index]*x_vec)
                M[i, j] = temp_poly.monomial_coefficient(col_index[j])

        # for debugging purposes only: check if result of educational implementation matches built-in function
        M_check_equals_M = (M_check.rref() == M.rref())
        print(f"Is the result of the educational implementation the same as that of the built-in? {'Yes' if M_check_equals_M else 'No'}")
        # if not M_check_equals_M:
            # print(f"M of educational implementation: \n{M.rref()}")
            # print(f"M of Ran'simplementation: \n{M_check.rref()}")
            # print("-------------"*5)
            # print(f"col_index of the educational implementation \t {col_index}")
            # print(f"col_index of Ran's implementation \t {col_index_check}")
            # print("-------------"*5)
            # print(f"rows in Ran's implementation: \t {rows}")
        
    # if not educational_implementation
    else:
        M, row_index, col_index = M_check, rows, col_index_check
        

    if can_print: # nb_rows can differ from len(row_index) if prune is True (idem for nb_cols)       
        print("row_index: \t", row_index)
        print('nb_rows: \t', nb_rows) 
        print('col_index: \t', len(col_index), 'nb_cols: \t', nb_cols)
        print("col_index: \t", col_index)

    return M, row_index, col_index



