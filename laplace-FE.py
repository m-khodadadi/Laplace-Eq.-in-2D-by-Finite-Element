import csv
import numpy


def zeros(n, m=None):
    if m is not None:
        return [[0] * m for i in range(n)]
    return [0] * n


def multiply(A, B):
    rows = len(A)
    c = len(A[0])
    r = len(B)
    cols = len(B[0])
    if c != r:
        raise Exception("Cannot multiply the two matrices. Incorrect dimensions.")
    result = []
    for ii in range(rows):
        row = []
        for jj in range(cols):
            cel = 0
            for kk in range(c):
                cel += A[ii][kk] * B[kk][jj]
            row.append(cel)
        result.append(row)
    return result


def c_multiply(c, B):
    result = []
    for i in range(len(B)):
        row = []
        for j in range(len(B[i])):
            row.append(c * B[i][j])
        result.append(row)
    return result


def sum_matrix(A, B):
    ra = len(A)
    ca = len(A[0])
    rb = len(B)
    cb = len(B[0])
    if ca != cb and ra != rb:
        raise Exception("invalid arguments!")
    result = zeros(ra, ca)
    for i in range(ra):
        for j in range(ca):
            result[i][j] = A[i][j] + B[i][j]
    return result


def transpose(A):
    if len(A) == 0:
        raise Exception("Invalid argument!")
    return [list(i) for i in zip(*A)]


kx = 1
ky = 1
with open('coor.txt') as f:
    coor = [[int(float(c.strip())) for c in list(line)] for line in csv.reader(f)]
with open('noc.txt') as ff:
    noc = [[int(float(c.strip())) for c in list(line)] for line in csv.reader(ff)]
[r.pop(0) for r in coor]
[r.pop(0) for r in noc]
f = lambda x=None: - x ** 2 + 10 * x
bc = list(range(1, 12)) + [12, 22, 23, 33, 34, 44, 45, 55] + list(range(56, 67))
tbc = [zeros(19) + [f(r[0]) for r in coor[55:66]]]
nn = len(coor)
ne = len(noc)
ka = zeros(nn, nn)
for k in range(ne):
    xel = [coor[c - 1][0] for c in noc[k]]
    yel = [coor[c - 1][1] for c in noc[k]]
    detj = (xel[0] - xel[2]) * (yel[1] - yel[2]) - (xel[1] - xel[2]) * (yel[0] - yel[2])
    dN_dx = [[(yel[1] - yel[2]) / detj, (yel[2] - yel[0]) / detj, (yel[0] - yel[1]) / detj]]
    dN_dy = [[(xel[2] - xel[1]) / detj, (xel[0] - xel[2]) / detj, (xel[1] - xel[0]) / detj]]
    ke = c_multiply(detj * 0.5, sum_matrix(c_multiply(kx, multiply(transpose(dN_dx), dN_dx)),
                                           c_multiply(ky, multiply(transpose(dN_dy), dN_dy))))
    q = [c - 1 for c in noc[k]]
    for i in range(3):
        for j in range(3):
            ka[q[i]][q[j]] += ke[i][j]
for i in bc:
    for j in range(len(ka[i - 1])):
        ka[i - 1][j] = 0
fa = [[r[j - 1] for j in bc] for r in ka]
fa = multiply(c_multiply(-1, fa), transpose(tbc))
for r in ka:
    for j in bc:
        r[j - 1] = 0
for i, j in zip(bc, range(len(bc))):
    fa[i - 1][0] = tbc[0][j]
for i in range(len(bc)):
    ka[bc[i] - 1][bc[i] - 1] = 1
T = numpy.linalg.solve(ka, fa)
x = [r[0] for r in coor]
y = [r[1] for r in coor]
ro = "  {0:03d}    {1:03d}    {2:07.4f}"
print("  x      y      T")
for i in range(len(coor)):
    print(ro.format(coor[i][0], coor[i][1], T[i][0]))
