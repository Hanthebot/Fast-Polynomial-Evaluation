""" sample generator """
from random import randint

COUNT = 10
MIN_DEGREE = 1000
MAX_DEGREE = 3000
MAX_VAL = (1 << 16) + 1

for i in range(COUNT):
    degree = randint(MIN_DEGREE, MAX_DEGREE)
    l = [str(randint(1, MAX_VAL)) for _ in range(degree + 1)]
    with open(f"sample/input{i+1}.txt", "w", encoding="utf-8") as fil:
        fil.write(f"{degree}\n{' '.join(l)}")
