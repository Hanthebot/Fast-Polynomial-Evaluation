""" sample generator """
import os
import sys
from random import randint
from Crypto.Util.number import getPrime

COUNT = 5

PRINT_COEFF = 0
PRINT_RESULT = 1
VERIFY_RESULT = 1
M = 3
LOG_N = 5
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

if __name__ == "__main__":
    PRINT_COEFF = int(input("Enter PRINT_COEFF: "))
    PRINT_RESULT = int(input("Enter PRINT_RESULT: "))
    VERIFY_RESULT = int(input("Enter VERIFY_RESULT: "))
    COUNT = int(input("Enter no. of input to generate: "))
    M = int(input("Enter M: "))
    SAMPLE_RANGE[0] = int(input("Enter lower bound on no. of points: "))
    SAMPLE_RANGE[1] = int(input("Enter upper bound on no. of points: "))
    custom_modulo = int(input("Enter custom modulo: "))
    custom_shape = []
    print(f"Enter custom shape ({M} times): ")
    for i in range(M):
        custom_shape.append(int(input(f"Dimension [{i+1}]: ")))
    prefix = input("Enter prefix for input files: ")

    for i in range(COUNT):
        modulo = custom_modulo if (custom_modulo != 0) else randint(2, getPrime(LOG_N) * 2)
        shape = custom_shape if (custom_shape[0] != 0) else [randint(max(1, modulo >> 2), modulo) for _ in range(M)]
        num_samples = randint(SAMPLE_RANGE[0], SAMPLE_RANGE[1])
        eval_points = "\n".join([" ".join([str(randint(1, modulo - 1)) for _ in range(M)]) for _ in range(num_samples)])
        coeff = recur_form(shape, modulo)
        with open(f"multivar_sample/input_{prefix}_{i+1}.txt", "w", encoding="utf-8") as fil:
            fil.write(f"{PRINT_COEFF} {PRINT_RESULT} {VERIFY_RESULT}\n{M}\n{modulo}\n{' '.join([str(sh) for sh in shape])}\n{num_samples}\n{eval_points}\n\n{coeff}")
