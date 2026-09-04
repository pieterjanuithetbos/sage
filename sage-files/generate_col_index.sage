def create_col_index(x, n, d, col_temp=None, degree_index=0):
    if col_temp is None:
        col_temp = 1
    if n == d:
        res = 1
        for el in x:
            res *= el
        return res
    col_index = []
    # base case
    if len(list(col_temp.factor())) == d:
        # print(col_temp)
        return [col_temp]
    for possible_start_index in range(degree_index, n - d + degree_index):
        # extend the solutions
        ext_sol = []
        for extension_index in range(possible_start_index, min(n - d + degree_index + 1, n)): # TODO understand the index issue
            temp_sol = col_temp*x[extension_index]
            ext_sol.append((temp_sol, extension_index))
        # continue backtracking
        for sol in ext_sol:
            res = create_col_index(x, len(x), d, sol[0], sol[1]+1)
            # print(res)
            if [res] != []:
                col_index += res
            # print(possible_start_index, col_index)
    return col_index

q = 9
n = 5
F = GF(q)
d = 4
R = PolynomialRing(F, n, 'x')
# print(R.monomials_of_degree(2))
print(create_col_index(R.gens(), n, d))