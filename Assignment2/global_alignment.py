"""Global Alignment using Needleman-Wunsch Algorithm."""

import numpy as np


def initiate_matrix(first_seq, second_seq, gap_penalty_value):
    taller_seq = max(len(first_seq), len(second_seq))
    shorter_seq = min(len(first_seq), len(second_seq))
    alignment_matrix = np.zeros((shorter_seq + 1, taller_seq + 1))
    # Fill the first Row
    for i in range(taller_seq):
        alignment_matrix[0, i + 1] = (i + 1) * gap_penalty_value
    # Fill the first column
    for i in range(shorter_seq):
        alignment_matrix[i + 1, 0] = (i + 1) * gap_penalty_value
    return alignment_matrix


# Returns the max score and its position
def get_max_score(scores):
    scores = np.array(scores)
    max_scores = max(scores)
    indices = np.where((scores == max_scores))
    indices = np.array(indices)

    return max(scores), indices


def check_alignment_score(first_seq, sec_seq, max_score, match, mismatch, gap):
    score = 0
    for i in range(len(first_seq)):
        y = first_seq[i]
        x = sec_seq[i]
        if first_seq[i] == '-' or sec_seq[i] == '-':
            score = score + gap
        elif first_seq[i] == sec_seq[i]:
            score = score + match
        elif first_seq[i] != sec_seq[i]:
            score = score + mismatch
    if score == max_score:
        check = "Correct_possible_Alignment"
    else:
        check = "Incorrect_Alignment"
    return score, check


def trace_back(path_matrix, first_seq, sec_seq, match, mismatch, gap):
    still_in = True
    new_first_seq = []
    new_sec_seq = []
    I = len(first_seq) - 1
    J = len(sec_seq) - 1
    i = path_matrix.shape[0] - 1
    j = path_matrix.shape[1] - 1
    while still_in:
        if path_matrix[i, j, 1] == 'D':
            new_first_seq.append(first_seq[I])
            new_sec_seq.append(sec_seq[J])
            I = I - 1
            J = J - 1
            i = i - 1
            j = j - 1
        elif path_matrix[i, j, 1] == 'F':
            new_first_seq.append('-')
            new_sec_seq.append(sec_seq[J])
            J = J - 1
            i = i - 1
        elif path_matrix[i, j, 1] == 'S':
            new_first_seq.append(first_seq[I])
            new_sec_seq.append('-')
            I = I - 1
            j = j - 1
        if j == 0:
            still_in = False
    return new_first_seq[::-1], new_sec_seq[::-1]


def global_alignment(first_seq, second_seq, match_penalty_value, mismatch_penalty_value, gap_penalty_value):
    """ First step is to fill the matrix with preset initial values. """
    alignment_matrix = initiate_matrix(first_seq, second_seq, gap_penalty_value)
    path_matrix = np.zeros((alignment_matrix.shape[0], alignment_matrix.shape[1], 3), dtype=str)
    scores = []
    """ Second step is to apply get the max score method, then assign it to the current cell. """
    for i in range(1, alignment_matrix.shape[0]):
        for j in range(1, alignment_matrix.shape[1]):
            row_score = alignment_matrix[i, j - 1] + gap_penalty_value
            column_score = alignment_matrix[i - 1, j] + gap_penalty_value
            if second_seq[i - 1] == first_seq[j - 1]:
                diagonal_score = alignment_matrix[i - 1, j - 1] + match_penalty_value
            else:
                diagonal_score = alignment_matrix[i - 1, j - 1] + mismatch_penalty_value
            scores.append(row_score)
            scores.append(column_score)
            scores.append(diagonal_score)
            alignment_matrix[i, j], ex_cell = get_max_score(scores)
            scores.clear()
            for I in range(ex_cell.size):
                if ex_cell[0][I] == 0:
                    path_matrix[i, j, 1] = "S"
                elif ex_cell[0][I] == 1:
                    path_matrix[i, j, 1] = "F"
                elif ex_cell[0][I] == 2:
                    path_matrix[i, j, 1] = "D"

    max_score = alignment_matrix[i, j]
    """ Third step is to trace back."""
    f, s = trace_back(path_matrix, first_seq, second_seq, match_penalty_value, mismatch_penalty_value,
                      gap_penalty_value)
    """Last step is to check the max score with the aligned sequences score. """
    new_s, check = check_alignment_score(f, s, max_score, match_penalty_value, mismatch_penalty_value,
                                         gap_penalty_value)
    return f, s, new_s, check


# Testing
a = ['G', 'A', 'A', 'T', 'T', 'C', 'A', 'G', 'T', 'T', 'A']
b = ['G', 'G', 'A', 'T', 'C', 'G', 'A']
c = ['C', 'T', 'A', 'T', 'T', 'G', 'A', 'C', 'G', 'T', 'A', 'A', 'C', 'A', 'T']
d = ['C', 'T', 'A', 'T', 'T', 'G', 'A', 'A', 'C', 'A', 'T']
gap_penalty = -3
Match = 4
Mismatch = -1
first_seq, sec_seq, Max_score, trace_back_check = global_alignment(c, d, Match, Mismatch, gap_penalty)
first_seq_1, sec_seq_1, Max_score_1, trace_back_check_1 = global_alignment(a, b, 5, -3, -4)
""" Checking."""
print(first_seq)
print(sec_seq)
print(Max_score, trace_back_check)
print(first_seq_1)
print(sec_seq_1)
print(Max_score_1, trace_back_check_1)
