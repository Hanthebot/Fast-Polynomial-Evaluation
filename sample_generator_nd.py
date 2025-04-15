""" sample generator """
import os
from random import randint
from Crypto.Util.number import getPrime

COUNT = 5

PRINT_COEFF = 0
PRINT_RESULT = 1
VERIFY_RESULT = 1
M = 3
LOG_N = 6
SAMPLE_RANGE = [10, 30]

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
    modulo = randint(2, getPrime(LOG_N) * 2)
    shape = [randint(modulo >> 1, modulo) for _ in range(M)]
    num_samples = randint(SAMPLE_RANGE[0], SAMPLE_RANGE[1])
    eval_points = "\n".join([" ".join([str(randint(1, modulo - 1)) for _ in range(M)]) for _ in range(num_samples)])
    coeff = recur_form(shape, modulo)
    with open(f"multivar_sample/input{i+1}.txt", "w", encoding="utf-8") as fil:
        fil.write(f"{PRINT_COEFF} {PRINT_RESULT} {VERIFY_RESULT}\n{M}\n{modulo}\n{' '.join([str(sh) for sh in shape])}\n{num_samples}\n{eval_points}\n\n{coeff}")
