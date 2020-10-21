#!/usr/bin/env python
import numpy as np
import sys

def csvdiff(tol, file1, file2):
    data1 = np.genfromtxt(file1, delimiter=',', skip_header=1)
    data2 = np.genfromtxt(file2, delimiter=',', skip_header=1)
    diff = abs(data2 - data1)
    if (diff > tol).any():
        print(diff)

def main():
    assert(len(sys.argv) == 4)
    tol = float(sys.argv[1])
    file1 = sys.argv[2]
    file2 = sys.argv[3]
    csvdiff(tol, file1, file2)

if __name__ == "__main__":
    main()
