""" sample generator """
import os
from random import randint
from Crypto.Util.number import getPrime

COUNT = 5

M = 3
LOG_N = 4

if not os.path.exists("multivar_sample"):
    os.makedirs("multivar_sample")

def recur_form(sh:list, _prime:int)->str:
    """ recursively generate a nice-looking values """
    if len(sh) == 0:
        return ""
    if len(sh) == 1:
        return " ".join([str(randint(1, _prime - 1)) for _ in range(sh[0] + 1)])
    return "".join([recur_form(sh[1:], _prime) + "\n" for _ in range(sh[0] + 1)])

for i in range(COUNT):
    prime = getPrime(LOG_N)
    shape = [randint(prime >> 1, prime - 2) for _ in range(M)]
    coeff = recur_form(shape, prime)
    with open(f"multivar_sample/input{i+1}.txt", "w", encoding="utf-8") as fil:
        fil.write(f"0 0\n{M}\n{prime}\n{' '.join([str(sh) for sh in shape])}\n\n{coeff}")
