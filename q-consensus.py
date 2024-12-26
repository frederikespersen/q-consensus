#!/usr/bin/env python3
"""
    Q-CONSENSUS
    --------------------------------------------------------------------------------

    Tools for calculating Q-scores of a consensus sequence
    made from an alignment of reads

    // Frederik Espersen Knudsen, December 2024

    --------------------------------------------------------------------------------
"""

# ················································································· #

import numpy as np

# ················································································· #

def q_consensus(R: np.array,
                Q: np.array,
                A: np.array = np.array(['A', 'C', 'G', 'T'])) -> tuple[np.array, np.array]:
    """
    Takes an alignment of reads and their position-wise Q-scores,
    returns the maximum likelihood sequence and its consensus Q-score.

    Has an optional parameter to specify the library of possible bases in reads.

    Assumes that there are no prior base biases for both reads and true sequence,
    and that reads are conditionally independent given the true sequence.

    TODO: Handle alignment gaps

    :param R: Aligned reads as an LxN ``np.array`` of bases, where L is the length of alignment and N is the number of reads.
    :param Q: The Q-scores for the called bases in R  as an LxN ``np.array``; thus must match the dimensions of R.
    :param A: The possible nucleotides that can be called for a read as an LxN ``np.array``. Defaults to {A,C,G,T}.
    :return: A tuple 's, q', where s is the consensus sequence and q is the associated Q-scores.
    """
    # Check inputs
    assert (np.unique(A) == A).all(), "The provided alphabet 'A' must not contain duplicate values!"
    assert set(np.unique(R)) <= set(A), "Read alignment 'R' must only contain letters found in provided alphabet 'A'!"
    assert R.shape == Q.shape, "'R' (aligned read bases) and 'Q' (aligned read Q-scores) must have the same dimensions!"
    assert (Q > 0).all(), "All Q-scores in 'Q' must be greater than zero!"

    # Determining sequence length and read number
    L, N = Q.shape

    # Calculating read error rates
    E = 10**(-Q/10)

    # Setting Read-Sequence mask
    S = np.repeat(A[None,:], L, axis=0)
    M = S[:,None,:] == R[:,:,None]

    # Calculating base probabilities per read
    P = (M) * (1 - E[:,:,None]) + (1 - M) * (E[:,:,None] / (A.shape[0] - 1))

    # Bayesian normalization with uniform priors
    P = P.prod(axis=1) / P.prod(axis=1).sum(axis=1)[:,None]

    # Determining maximum likelihood sequence
    s = S[np.arange(L), P.argmax(axis=1)]

    # Calculating Q-score for maximum likelihood sequence
    p = P.max(axis=1)
    q = np.round(-10 * np.log10(1 - p), 0)

    return s, q

# ················································································· #

if __name__ == '__main__':
    R = np.array([['A', 'A', 'C', 'A', 'T', 'G', 'A', 'G'],
                  ['C', 'C', 'C', 'T', 'T', 'A', 'A', 'G']])
    Q = np.array([[20, 20, 20, 15, 13, 20, 22, 19],
                  [10, 20, 20, 15, 13, 20, 22, 19]])
    print(q_consensus(R, Q))
