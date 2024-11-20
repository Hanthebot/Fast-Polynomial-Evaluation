""" sample generator """
import os
from random import randint

COUNT = 10
LOG_N = 20
MIN_DEGREE = 1000
MAX_DEGREE = 3000
MAX_VAL = (1 << LOG_N) + 1

if not os.path.exists("sample"):
    os.makedirs("sample")

for i in range(COUNT):
    degree = randint(MIN_DEGREE, MAX_DEGREE)
    l = [str(randint(1, MAX_VAL)) for _ in range(degree + 1)]
    with open(f"sample/input{i+1}.txt", "w", encoding="utf-8") as fil:
        fil.write(f"0 0\n{LOG_N}\n{degree}\n{' '.join(l)}")
