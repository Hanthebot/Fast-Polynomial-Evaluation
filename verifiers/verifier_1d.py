""" Python verifier for polynomial computations in 1D. """

def evaluate_1d(_x: int, _coeff: list, _p: int)->int:
    """ emulator """
    num = 0
    for ind, c in enumerate(_coeff):
        num += c * pow(_x, ind, _p)
        if num >= _p:
            num %= _p
    return num

def parse_coeff_1d(_i: str)->tuple:
    """ parse coefficients """
    with open(f"sample/input{_i}.txt", "r", encoding="utf-8") as fil:
        _lines = fil.readlines()
        _logn = int(_lines[1])
        _p = (1 << _logn) + 1
        return (_p, [int(i) for i in _lines[3].split(" ")])
    return None

if __name__ == "__main__":
    i = input("Enter file name: sample/input{???}.txt: ")
    p, coeff = parse_coeff_1d(i)
    while True:
        try:
            x = int(input("Enter x: "))
            print(evaluate_1d(x, coeff, p))
        except EOFError:
            break
