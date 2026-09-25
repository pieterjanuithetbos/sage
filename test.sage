import time
load("wedge_attack.sage")
load("generate_macaulay.sage")
load("generate_UOV_instance.sage")


def test_wedge_attack(q, v, o, m, prune=False, educational_implementation=True, can_print=False):
    pk, sk, P, oil, vinegar, O, V, gens_vector = generate_UOV_variables(n, m, o, q)
    res, v, M, dim, timings = wedge_attack(pk, n, q, o, prune, educational_implementation, can_print)

q = 2 
v, o, m = 4, 3, 6
n = o + v

test_wedge_attack(q, v, o, m)


"""
generate UOV instance and perform the wedge attack on it
"""
def wedge_once(n, m, o, q, prune, educational_implementation, can_print=False):
    timings = []

    # generate UOV instance
    time_generate_UOV_variables_start = time.time()
    pk, sk, P, oil, vinegar, O, V, gens_vector = generate_UOV_variables(n, m, o, q)
    time_generate_UOV_variables_delta = time.time() - time_generate_UOV_variables_start
    timings.append(("building a UOV instance", time_generate_UOV_variables_delta))

    x_vec = vector(P.gens())
    retrieved_space, kernel_vector, M, dim, timings_wedge_attack = wedge_attack(pk, n, q, o, prune, educational_implementation, can_print=False)
  
    timings.extend(timings_wedge_attack)
    
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
        if O.rref() == retrieved_space.rref():
            print("O.rref == retrieved_space.rref(), attack succeeded.")
        print("-----------"*7)
        # for timing in timings:
        #     print(timing[0], "\t", timing[1])
            
    return dim, O.rref(), retrieved_space.rref(), timings

"""
Perform the wedge attack multiple times
"""
def wedge_party(q, v, o, m, test_it, prune, educational_implementation, can_print=False):
    n = v + o
    if is_admissible(n, m, o):
        counter_rank = 0
        counter_too_high = 0
        counter_attack = 0
        for i in range(test_it):
            exp, O, retrieved, timings = wedge_once(n, m, o, q, prune, educational_implementation, can_print)
            if exp == 1:
                counter_rank += 1
            if exp > 1:
                counter_too_high += 1
            if O == retrieved:
                counter_attack += 1
            # print("rank prediction", predict_rank(n, m, o))
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

"""
call wedge_once with the given parameters
"""
def test(educational_implementation, prune, can_print):
    q = 16
    v, o, m = 8, 6, 8
    test_it = 1
    n = v + o

    wedge_party(q, v, o, m, test_it, prune, educational_implementation, can_print)
    # test = wedge_once(n, m, o, q, educational_implementation, can_print, prune)

test(False, True, True)