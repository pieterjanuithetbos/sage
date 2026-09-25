import math
import itertools
import time

def retrieve_kernel(n, q, pk, o, P, x_vec, prune, educational_implementation=True):
    """
    input:
        pk    (Sequence) the public key a sequence of polynomial equations
        x_vec (vector) a vector containing the generators of P

    output: 
        v (vector) kernel vector of M
        col_index (list) contains (boolean) monomials of degree d, in deglex order
        M (matrix), the macaulay matrix of degree o of pk
        dim (integer) the dimension of the kernel 
    """
    timings = []

    print("prune? \t",  prune)
    
    # build Macaulay matrix of degree o
    time_generate_M_start = time.time()
    M, row_index, col_index = generate_M(n, q, pk, o, 'deglex', educational_implementation, prune)
    time_generate_M_delta = time.time() - time_generate_M_start 
    timings.append(("generate Macaulay matrix", time_generate_M_delta))
    print("time to generate Macaulay matrix: \t", time_generate_M_delta)
    print("M has dimensions: \t", M.nrows(), M.ncols())

    # compute the kernel vector
    time_compute_kernel_start = time.time()
    v = M.right_kernel()
    time_compute_kernel_delta = time.time() - time_compute_kernel_start
    timings.append(("compute the kernel of M", time_compute_kernel_delta))
    print("time compute the kernel of M \t \t", time_compute_kernel_delta)
    v_vector = v.basis()

    dim = v.dimension()


    return v_vector[0], col_index, M, dim, timings


def build_oil_space(col_index, kernel_vector, n, o, q, educational_implementation, can_print = False):
    v = n-o

    # convert the monomial labeling of the columns to a corresponding labeling in the wedge basis 
    col_index_wedge_basis = []
    for i in range(len(col_index)):
        monomial = col_index[i]
        complement = [(1 - exp) % 2 for exp in monomial.exponents()[0]]
        # TODO check ordering of the monomials: is it not reverse?
        col_index_wedge_basis.append(complement)

    # iterating over the above created indices, select the first subset K of indices with |K| = v, 
    # such that the coefficient of K in col_index_wedge_basis is 1
    found = False
    v_ind, o_ind, i = [], [], 0
    while not found and i < len(kernel_vector):
        current_wedge_index = col_index_wedge_basis[i]
        if sum(current_wedge_index) == v and kernel_vector[i] != 0:
            vinegar_coefficient = kernel_vector[i]

            # create two lists, containing the indices of the oil and vinegar variables, respectively
            for k in range(n):
                if current_wedge_index[k] == 1:
                    v_ind.append(k)
                else:
                    o_ind.append(k)
            found = True
        i += 1
    
    # create a list comparing the col_index with its wedge basis correspondents (for debugging purposes)
    temp = []
    for i in range(len(col_index)):
        temp.append((col_index[i],  col_index_wedge_basis[i]))
        
    # filter basis elements e_J satisfying J = {i} \union {o+1, ..., n} \ {j} for 1 <= i <= o < j <= n
    res = Matrix(GF(q), o, n)
    for k in range(len(col_index_wedge_basis)):
        current_exp = col_index_wedge_basis[k]
        
        # check if only one oil vector is included in the wedge)
        one_i = sum([current_exp[oo] for oo in o_ind]) == 1
        # check if only one vinegar vector is not included in the wedge
        all_but_one_j = sum([current_exp[vv] for vv in v_ind]) == v - 1
        # check if all the vinegar vectors are included
        all_vinegar = sum([current_exp[vv] for vv in v_ind]) == v
        
        # fill the matrix
        if (one_i & all_but_one_j):
            i = [current_exp[oo] for oo in o_ind].index(1)
            j = [current_exp[vv] for vv in v_ind].index(0)
            res[i, v_ind[j]] = kernel_vector[k]
            if can_print:
                print("(i, j) \t", f"({i}, {j})")
                print("o \t", o, "\t", current_exp, "one_i \t", one_i, "all_but_one_j \t", all_but_one_j)
                print(res)
        elif all_vinegar:
            pass # TODO can this case happen oly once?
    
    # add a diagonal matrix
    for i in range(o):
        res[i, o_ind[i]] = vinegar_coefficient
    
    if can_print:
        print("vinegar coefficient \t", vinegar_coefficient)
        print("de index van vinegar variables \t \t", v_ind)
        print("de index van oil variables \t \t", o_ind)
        # print("de waarde van de kernel vector voor v_ind \t", kernel_vector[col_index_wedge_basis.index(v_ind)])
        # print("both_indices: \t",temp)

    return res


def wedge_attack(pk, n, q, o, prune, educational_implementation=True, can_print = True):
    timings = []
    P = PolynomialRing(GF(q), n, 'x', order='deglex')
    x_vec = vector(P.gens())

    # retrieve the kernel
    v, column_index, M, dim, timings_retrieve_kernel = retrieve_kernel(n, q, pk, o, P, x_vec, prune, educational_implementation)
    timings.extend(timings_retrieve_kernel)

    # build the oil space from the kernel vector
    time_build_oil_space_start = time.time()
    res = build_oil_space(column_index, v, n, o, q, False)
    time_build_oil_space_delta  = time.time() - time_build_oil_space_start
    timings.append(("build the oil space from the kernel", time_build_oil_space_delta))
    print("time to build the oil space \t \t", time_build_oil_space_delta)

    # print statements
    if can_print:
        print(f"retrieved oil space in rref: \n{res.rref()}")
        # print("Macaulay matrix \t", M)

    return res, v, M, dim, timings

def predict_rank(n, m, o):
    """
    return the o-th term of the Hilbert series of a system of m quadratic equations in n variables over F/<x^2, ..., x_n^2>
    """
    v = n - o
    s_python = max(0, 
        sum((-1)**j * math.comb(m + j - 1, j) * math.comb(v + o, v + 2*j) 
            for j in range(o//2 + 1)))
    
    # TODO why is this implementation not equivalent?
    # s = max(0, sum(
    #     (-1)**i*binomial(m+i-1, i) * binomial(v+o, v+2*i), 
    #     "i", 0, floor(o/2)))
    
    return s_python
    
def is_admissible(n, m, o):
    """
    checks if the two conditions mentioned in Ran's paper are satisfied 
    """
    v = n-o
    cond_1 = (min((o-1)/2, 2)*m > v) and (v > o)
    s = predict_rank(n, m, o)
    cond_2 = s <= 1
    return cond_1 and cond_2
