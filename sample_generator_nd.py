""" sample generator """
import os
from random import randint

COUNT = 10

M = 5
LOG_N = 4
MIN_DEGREE = 2
MAX_DEGREE = 15
MAX_VAL = 1 << LOG_N

if not os.path.exists("multivar_sample"):
    os.makedirs("multivar_sample")

def recur_form(sh:list)->str:
    """ recursively generate a nice-looking values """
    if len(sh) == 0:
        return ""
    if len(sh) == 1:
        return " ".join([str(randint(1, MAX_VAL)) for _ in range(sh[0] + 1)])
    return "".join([recur_form(sh[1:]) + "\n" for _ in range(sh[0])])

for i in range(COUNT):
    shape = [randint(MIN_DEGREE, MAX_DEGREE) for _ in range(M)]
    str_form = recur_form(shape)
    with open(f"multivar_sample/input{i+1}.txt", "w", encoding="utf-8") as fil:
        fil.write(f"0 0\n{M}\n{LOG_N}\n\n{' '.join([str(sh) for sh in shape])}\n{str_form}")
