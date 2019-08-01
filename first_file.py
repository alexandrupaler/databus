import operator as op
from functools import reduce

import math

def ncr(n, r):
    """
    Code from:  https://stackoverflow.com/questions/4941753/is-there-a-math-ncr-function-in-python
    """
    r = min(r, n-r)
    numer = reduce(op.mul, range(n, n-r, -1), 1)
    denom = reduce(op.mul, range(1, r+1), 1)

    return numer / denom

def compute_prob_error_parity(err, nr_bits):
    """
    Computes the sum of all probabilities of odd weight errors which can flip
    a sequence of nr bits encoding a single parity bit
    1, 3, 5, 7,..., nr - 1 (if nr is even)

    :return: Total probability of all odd weight errors
    """

    sum = 0.0
    max_err = nr_bits
    for nr_err in range(max_err):
        if nr_err%2 == 1:
            # is this an odd number?
            sum += ncr(max_err, nr_err) * math.pow(1 - err, nr_bits - nr_err) * math.pow(err, nr_err)

    return sum

def compute_prob_fail_decode_repetition(err, nr_bits):
    """
        Assume a repetition code of length nr_bits
        Errors which are not correctable have a weight larger than (nr-1)/2
        -> This means that the correct majority cannot be computed back

        This function computes the probability that a majority from a repetition code will be wrongly read
    """

    sum = 0.0
    max_err = nr_bits // 2
    for i in range(max_err, nr_bits):
            sum += ncr(max_err, i) * math.pow(1 - err, nr_bits - i) * math.pow(err, i)

    return sum

def main():

    # number of potential error locations
    # for the moment, I am not discussing the if this is a square or a line.
    nr_error_locations = 1000

    # per qubit (location) physical error rate
    phys_err = 0.001

    # targetted reliability of logical measurement
    targeterr = math.pow(10, -12)
    print("targeted error probability of logical bus measurement", targeterr)

    # a single bus measurement encodes a parity
    # what is the probability that the bus holds the wrong parity after it has been measured?
    trial_err = compute_prob_error_parity(phys_err, nr_error_locations)
    print("a bus measurement will fail with", trial_err)

    currenterr = 1
    measurement_rep = 1
    while currenterr > targeterr:
        measurement_rep += 2
        currenterr = compute_prob_fail_decode_repetition(trial_err, measurement_rep)

    print("repeat N=", measurement_rep, "for an error of", currenterr)


if __name__ == "__main__":
    main()