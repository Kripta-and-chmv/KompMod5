import max_value
import numpy as np
import time
import tests
import scipy as sc
from IPython.display import display
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import sys

def get_arguments():
    with open("arguments.txt", "r") as f:
        file_str = f.read()
        args = file_str.split(" ")
        v, u, alpha = float(args[0]), float(args[1]), float(args[2])
        return v, u, alpha



def main():
    sys.stdout = open("output.txt", "w+")

    v, u, alpha = get_arguments()
    testing_seq(v, u, alpha, 50)
    testing_seq(v, u, alpha, 200)
    testing_seq(v, u, alpha, 1000)

    show_prob_density_and_function(v, u)

if __name__ == "__main__":
    main()