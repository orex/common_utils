#/usr/bin/env python3

import random
import math

def cmb(l):
    sm = 0
    bt = 1
    for k in l:
        bt = bt * math.factorial(k)
        sm = sm + k
    return math.factorial(sm) // bt

wb = open("data_combinatorics/ctest.dat", "w")

for i in range(1000000):
    ln = random.randint(2, 5)
    ls = []
    for i in range(ln):
        ls.append(random.randint(2, 30))
    ref = cmb(ls)
    if(ref > 2**90):
        continue
    if(ref >= 2**63-1):
        ref = -1
    lss = ", ".join([str(x) for x in ls])
    wb.write("{:d}, {}\n".format(ref, lss))