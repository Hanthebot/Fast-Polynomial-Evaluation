""" Python verifier for polynomial computations in ND. """

def evaluate_nd(xyz: list, _coeff: list, _p: int, val: int = 1)->int:
    """ emulator """
    num = 0
    if not isinstance(_coeff[0], list):
        for ind, c in enumerate(_coeff):
            num += c * pow(xyz[0], ind, _p) * val
            if num >= _p:
                num %= _p
        return num
    for ind, c in enumerate(_coeff):
        num += evaluate_nd(xyz[1:], c, _p, val * pow(xyz[0], ind, _p))
        if num >= _p:
            num %= _p
    return num

def parse_coeff_recur(_coeff: list, _lines: list, _shape: list, _ind: list)->None:
    """ parse coefficients recursively """
    if len(_shape) == 0:
        return
    if len(_shape) == 1:
        while _lines[_ind[0]] in ['', ' ', '\n']:
            _ind[0] += 1
        _nums =_lines[_ind[0]].split(" ")
        for _ix in range(_shape[0] + 1):
            _i = _nums[_ix]
            if _i not in [' ', '\n', '']:
                _coeff.append(int(_i))
        _ind[0] += 1
        return
    for _i in range(_shape[0] + 1):
        _coeff.append([])
        parse_coeff_recur(_coeff[_i], _lines, _shape[1:], _ind)

def parse_coeff_nd(_i: str)->tuple:
    """ parse coefficients """
    with open(f"multivar_sample/input{_i}.txt", "r", encoding="utf-8") as fil:
        _lines = fil.readlines()
        _m = int(_lines[1])
        _prime = int(_lines[2])
        _shape = [int(i) for i in _lines[3].split(" ") if i not in ['\n', '', ' ']]
        _samples  = int(_lines[4])
        _coeff = []
        parse_coeff_recur(_coeff, _lines[6 + _samples:], _shape, [0])
        return (_prime, _m, _coeff)
    return None

if __name__ == "__main__":
    i = input("Enter file name: multivar_sample/input{???}.txt: ")
    p, m, coeff = parse_coeff_nd(i)
    while True:
        try:
            vals = [int(i) for i in input(f"Enter {m} coordinates: ").split(" ")]
            if vals[0] == vals[1] == -1:
                print(coeff)
            print(evaluate_nd(vals, coeff, p))
        except EOFError:
            break
