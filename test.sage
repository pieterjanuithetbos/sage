import time
load("wedge_attack.sage")
load("generate_macaulay.sage")
load("generate_UOV_instance.sage")


"""
generates a UOV instance and performs a wedge attack. 
Returns information about the created UOV instance and the result of the wedge attack.
"""
def test_wedge_attack(q, v, o, m, prune=False, educational_implementation=True, can_print=False):
    pk, sk, P, oil, vinegar, O, V, gens_vector = generate_UOV_variables(n, m, o, q)
    res, v, M, dim, timings = wedge_attack(pk, n, q, o, prune, educational_implementation, can_print)
    return ((pk, sk, P, oil, vinegar, O, V, gens_vector), (res, v, M, dim, timings))

# test_wedge_attack(q, v, o, m)


"""
generate UOV instance and perform the wedge attack on it.
"""
def wedge_once(n, m, o, q, prune, educational_implementation, can_print=False):
    timings = []
    time_generate_UOV_variables_start = time.time()
    pk, sk, P, oil, vinegar, O, V, gens_vector = generate_UOV_variables(n, m, o, q)
    time_generate_UOV_variables_delta = time.time() - time_generate_UOV_variables_start
    timings.append(("building a UOV instance", time_generate_UOV_variables_delta))

    x_vec = vector(P.gens())
    retrieved_space, kernel_vector, M, dim, timings_wedge_attack = wedge_attack(pk, n, q, o, prune, educational_implementation, can_print=False)
  
    timings.extend(timings_wedge_attack)
    succeeded = O.rref() == retrieved_space.rref()
    
    if can_print:
        if dim != 1:
            print(f"De kernel heeft rang {dim}")
        # print("vinegar_coefficient: \t", vinegar_coefficient)
        print(f"O in reduced echelon form: \n{O.rref()}")
        assert O.rank() == o
        print("the Macaulay matrix has rank: \t", M.rank(), f" and dimensions:\t{M.nrows()}, {M.ncols()}")
        # print("the kernel vector: \t", kernel_vector[1]) # prints only the second element. We expect the first element to be the zero matrix
        print("-----------"*7)
        print(f"the predicted rank of the kernel is: \t {predict_rank(n, m, o)}, its observed rank is: \t {dim}")
        print("-----------"*7)
        print(f"the retrieved oil space in rref: \t \n{retrieved_space.rref()}")
        print("-----------"*7)
        if succeeded:
            print("O.rref == retrieved_space.rref(), attack succeeded.")
        print("-----------"*7)
        # for timing in timings:
        #     print(timing[0], "\t", timing[1])
            
    return dim, O.rref(), retrieved_space.rref(), timings, (M.nrows(), M.ncols()), succeeded

"""
Perform the wedge attack (wedge_once) multiple times
"""
def wedge_party(q, v, o, m, test_it, prune, educational_implementation, can_print=False):
    n = v + o
    if is_admissible(n, m, o):
        counter_rank = 0
        counter_too_high = 0
        counter_attack = 0
        for i in range(test_it):
            kernel_dim, O, retrieved, timings, dimensions, succeeded = wedge_once(n, m, o, q, prune, educational_implementation, can_print)
            if exp == 1:
                counter_rank += 1
            if exp > 1:
                counter_too_high += 1
            if O == retrieved:
                counter_attack += 1
            # print("rank prediction", predict_rank(n, m, o))
            if can_print:
                print(f"iteratie {i}")
                print("aandeel van de testen waar de rang gelijk is aan 1: \t", round(counter_rank/test_it, 2))
                print("aandeel van de testen waar de rang groter is dan 1: \t", round(counter_too_high/test_it, 2))
                print("aandeel van de testen waar O.rref() == res.rref(): \t", round(counter_attack/test_it, 2))
     
    else: 
        print("parameters not admissible")
        if not min((o-1)/2, 2)*m > v:
            print("not min((o-1)/2, 2)*m > v")
        if not v > o:
            print("not v > o")
        if predict_rank(n, m, o) > 1:
            print(f"{predict_rank(n, m, o)} > 1")
    return kernel_dim, dimensions, succeeded
    
"""
call wedge_once with the given parameters
"""
def test(educational_implementation, prune, can_print):
    q = 16
    v, o, m = 8, 6, 8
    test_it = 1
    n = v + o

    start = time.time()
    dim, M_dim, succeeded = wedge_party(q, v, o, m, test_it, prune, educational_implementation, can_print)
    stop = time.time()
    # test = wedge_once(n, m, o, q, educational_implementation, can_print, prune)
    print(f"q: \t {q} \nn: \t {n} \no: \t {o} \nm: \t {m}")
    print("---------------"*6)
    print(f"Attack took \t\t\t\t\t {stop-start} \t seconds")
    print(f"M has dimensions: \t\t\t\t {M_dim[0], M_dim[1]}")
    print(f"Kernel dimension: \t\t\t\t {dim}")
    print(f"O.rref() == retrieved_O.rref(): \t\t {succeeded}")

test(False, True, False)
