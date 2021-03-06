import numpy as np
import time
import tests
import scipy as sc
from IPython.display import display
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import sys
import neumann

def get_arguments():
    with open("arguments.txt", "r") as f:
        file_str = f.read()
        args = file_str.split(" ")
        u, v, alpha = float(args[0]), float(args[1]), float(args[2])
        return u, v, alpha

def main():
    sys.stdout = open("output.txt", "w+")

    u, v, alpha = get_arguments()
    
    amount = 1000
    seq = [neumann.generate(u, v) for x in range(amount)]
    #with open("normal_{}_{}_{}.txt".format(len(seq), u, v), "w") as f:
    #    f.write(str(seq))
    tests.smirnov(seq, u, v, alpha)
    tests.chisqr_test(seq, alpha, u, v)


if __name__ == "__main__":
    main()