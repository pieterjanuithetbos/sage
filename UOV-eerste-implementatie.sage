{
 "cells": [
  {
   "cell_type": "code",
   "execution_count": 1,
   "id": "7ebb4ebc-0ce1-4cbf-9089-974ec1a96ffc",
   "metadata": {
    "scrolled": true
   },
   "outputs": [
    {
     "data": {
      "text/plain": [
       "[       0        0        0        0 2*z2 + 1       z2   z2 + 2     2*z2]\n",
       "[       0        0        0        0 2*z2 + 2       z2   z2 + 2     2*z2]\n",
       "[       0        0        0        0 2*z2 + 1        0   z2 + 2        0]\n",
       "[       0        0        0        0   z2 + 2        1       z2 2*z2 + 2]\n",
       "[2*z2 + 1 2*z2 + 2 2*z2 + 1   z2 + 2 2*z2 + 1       z2   z2 + 1        0]\n",
       "[      z2       z2        0        1       z2   z2 + 2       z2       z2]\n",
       "[  z2 + 2   z2 + 2   z2 + 2       z2   z2 + 1       z2        1        0]\n",
       "[    2*z2     2*z2        0 2*z2 + 2        0       z2        0        2]"
      ]
     },
     "execution_count": 1,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "def generate_f(n, o, q):\n",
    "    '''\n",
    "    generates a random quadratic equation over F_q in n variables.\n",
    "    q has to be odd.\n",
    "    '''\n",
    "    F = GF(q)\n",
    "    # define the polynomial ring where the variables live\n",
    "    R = PolynomialRing(F, n, 't')\n",
    "    # construct a quadratic polynomial over F_q, such that every term has a vinegar variable\n",
    "    f_k = sum(R.gens()[i]*sum(F.random_element()*R.gens()[j] for j in range(n)) for i in range(o, n)) \n",
    "    Q = QuadraticForm(f_k);\n",
    "    # gets matrix of polar form\n",
    "    gram = Q.matrix() * F(2)**(-1);\n",
    "    \n",
    "    # further checks\n",
    "    det = Q.Gram_det()\n",
    "    assert det == gram.det()\n",
    "    x = vector(R.gens())\n",
    "    assert (x * gram * x == f_k)\n",
    "    return gram\n",
    "\n",
    "generate_f(8, 4, 9)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 2,
   "id": "8921d8f3-80ab-42aa-b96e-6b493df14e16",
   "metadata": {},
   "outputs": [],
   "source": [
    "def KeyGen(n, m, o, q):\n",
    "    '''\n",
    "    returns a tuple (pk, sk) containing the public key and the secret key.\n",
    "    The secret is again a tuple (FF, M), where \n",
    "            - FF is the central map and \n",
    "            - M an invertible linear transformation.\n",
    "    '''\n",
    "    F = GF(q)\n",
    "    FF = []\n",
    "    # pick the quadratic equations at random\n",
    "    for i in range(m):\n",
    "        FF.append(generate_f(n, o, q))\n",
    "    # print(\"FF[2] \\n\", FF[2])\n",
    "    # choose an invertible linear transform F_n -> F_n \n",
    "    M = GL(n, F).random_element().matrix()\n",
    "    sk = (FF, M)\n",
    "    # construct the public key\n",
    "    PP = []\n",
    "    for i in range(len(FF)):\n",
    "        PP.append(M.transpose()*FF[i]*M)\n",
    "    # print(\"PP[2] \\n\", PP[2])\n",
    "    pk = PP\n",
    "    return pk, sk"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 3,
   "id": "d7436d99-023d-4583-8a8f-d5c68fd93ccd",
   "metadata": {},
   "outputs": [],
   "source": [
    "def solve_FF(n, m, o, q, sk, t, F, R, x):\n",
    "    # pick random values for the 'vinegar' variables\n",
    "    for i in range(n-o):\n",
    "        # print(i+o-1)\n",
    "        x[i+o] = F.random_element()\n",
    "    # print(x)\n",
    "    # solve the system F(s') = t (m linear equations in o variables, m <= o)\n",
    "    FF, M = sk\n",
    "    system = [x*FF[i]*x - t[i] for i in range(m)]\n",
    "    A = matrix([ [eq.coefficient(v) for v in x[0:o]] for eq in system])\n",
    "    # TODO try with and without max\n",
    "    b = vector([-eq.constant_coefficient() for eq in system])\n",
    "    # print(\"het rechterlid\", b)\n",
    "    sol = A.solve_right(b)\n",
    "    return sol"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 4,
   "id": "d09061fa-611a-4c6a-8543-134a91933563",
   "metadata": {},
   "outputs": [],
   "source": [
    "def sign(n, m, o, q, sk, t):\n",
    "    \n",
    "    '''\n",
    "    returns the UOV signature as a vector of length n\n",
    "    '''\n",
    "    assert (len(t) == m)\n",
    "    F = GF(q)\n",
    "    R = PolynomialRing(F, n, 'x')\n",
    "    x = vector(R.gens())\n",
    "    \n",
    "    # we skip the hashing part\n",
    "    i = 1\n",
    "    while true:\n",
    "        try:\n",
    "            sol_temp = solve_FF(n, m, o, q, sk, t, F, R, x)\n",
    "            # print(\"sol\", sol_temp)\n",
    "            break\n",
    "        except ValueError:\n",
    "            i+=1\n",
    "            print(f'Hélas ... Too much vinegar. Poging {i}')\n",
    "    # print(\"de inverse transformatie\", sk[1].inverse())\n",
    "    # print('element', sol_temp[2])\n",
    "    # create s'\n",
    "    sol = [sol_temp[i] for i in range(o)]\n",
    "    sol.extend([x[i] for i in range(o, n)])\n",
    "    assert (len(sol) == n)\n",
    "    # print(sol)\n",
    "    signature = sk[1].inverse()*vector(sol)\n",
    "    return signature"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 5,
   "id": "5ab51774-7d9e-40ab-bea7-319a3a89156e",
   "metadata": {},
   "outputs": [],
   "source": [
    "def verify(n, m, o, q, PP, signature):\n",
    "    F = GF(q)\n",
    "    R = PolynomialRing(F, n, 'x')\n",
    "    x = vector(R.gens())\n",
    "    t_bar = []\n",
    "    for i in range(len(PP)):\n",
    "        t_bar.append(signature*PP[i]*signature)\n",
    "    return t_bar\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 6,
   "id": "c134962a-ebfe-4baa-8446-a2d4deb17d5a",
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "te ondertekenen: [z2, 2*z2, 2, 2*z2 + 2]\n",
      "response: [z2, (-z2), -1, (-z2 - 1)]\n"
     ]
    }
   ],
   "source": [
    "# testsetup\n",
    "n, m, o, q = 8, 4, 4, 9\n",
    "kg = KeyGen(n, m, o, q)\n",
    "pk, sk = kg[0], kg[1]\n",
    "t = [GF(q).random_element() for i in range(m)]\n",
    "s = sign(n, m, o, q, sk, t)\n",
    "v = verify(n, m, o, q, pk, s)\n",
    "# TEST\n",
    "print('te ondertekenen:', t)\n",
    "print('response:', v)"
   ]
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "SageMath 10.7",
   "language": "sage",
   "name": "sagemath"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.11.15"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 5
}
