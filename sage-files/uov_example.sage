from sage.all import *
load("sign.sage")
load("key_gen.sage")
load("verify.sage")
# testsetup
n, m, o, q = 8, 4, 4, 9
kg = KeyGen(n, m, o, q)
pk, sk = kg[0], kg[1]
t = [GF(q).random_element() for i in range(m)]
s = sign_uov(n, m, o, q, sk, t)
v = verify(n, m, o, q, pk, s)
# TEST
print('te ondertekenen:', t)
print('response:', v)
print(v == t)