import csv
from numba import njit, typed
import numpy as np
from typing import Any, Tuple, List

#### PARAMETERS FOR CONSTRAINTS
MINSUPPERCENT = 0
MINSUP = 1
EPSILON = 2
GAP = 3
TIMEOUTSECONDS = 4
LEVEL = 5
PRINTF = 6
PRINTL = 7
FORGETTABLE = 8

#### PARAMETERS FOR RELATIONS
NONE = -1
EQUALS = 0
LEFTMATCHES = 1
RIGHTMATCHES = 2
CONTAINS = 3
OVERLAPS = 4
MEETS = 5
FOLLOWS = 6

R = {EQUALS, LEFTMATCHES, RIGHTMATCHES, CONTAINS, OVERLAPS, MEETS, FOLLOWS}

# PART 1: BASIC LOAD AND PREPROCESSING
#
# load_file, preprocess
#

def load_file(filename: str) -> np.ndarray:
    with open(filename, "r") as f:
        reader = csv.reader(f, delimiter =' ')
        database = np.array(list(reader), dtype=np.int64)
    return database

@njit
def preprocess(database: np.ndarray) -> Tuple[int, int, int, float, float, List[Any], List[Any], List[Any]]:
    """
    Function for preprocessing data into the shape that the algorithm takes

    :param database: a loaded event interval sequence database

    :return tseq: size of the database
    :retuen tdis: number of unique event labels
    :return tintv: total number of intervals in the database
    :return aintv: average number of intervals per e-sequence
    :return avgtime: average timespan (temporal length of e-sequence)
    :return eseqdb: a processed event sequence database
    :return unique_labels: unique event labels per e-sequence
    :return initial_support: initial support of each event label
    """

    eseq_size = database[:, 0].max() + 1
    distinct_events = np.unique(database[:, 1])

    # Each e-sequence has different size so it cannot be numpy array
    eseqdb = typed.List()

    for i in range(eseq_size):
        eseq = typed.List()
        for row in database[database[:, 0] == i]:
            eseq.append((row[1], row[2], row[3]))
        eseqdb.append(eseq)
    #eseqdb = [[(row[1], row[2], row[3]) for row in database[database[:, 0] == i]] for i in range(eseq_size)]

    unique_labels = typed.List()

    for eseq in eseqdb:
        eseq_label = typed.List()
        for interval in eseq:
            eseq_label.append(interval[0])
        unique_labels.append(set(eseq_label))
    #unique_labels = [set([interval[0] for interval in eseq]) for eseq in eseqdb]
    timelength = [np.max(database[database[:, 0] == i][:, -1]) for i in range(eseq_size)]
    timelength = np.array(timelength)
    
    tseq = len(eseqdb)
    tdis = len(distinct_events)
    tintv = len(database)
    aintv = len(database) / len(np.unique(database[:, 0]))
    avgtime = timelength.sum() / len(timelength)

    # From 0 to np.max() + 1
    initial_support = np.zeros(np.max(database[:, 1]) + 1)

    for eseq in unique_labels:
        for label in eseq:
            #if label not in initialSupport:
            #    initialSupport[label] = 0
            initial_support[label] += 1

    return tseq, tdis, tintv, aintv, avgtime, eseqdb, unique_labels, initial_support


# PART 2: INTERVAL HANDLING
#
@njit
def get_interval_duration(interval: Tuple[int, int, int]) -> int:
    # (event label, start time, end time)
    return interval[2] - interval[1]

@njit
def hash_interval(interval: Tuple[int, int, int]) -> int:
    # for python-supported hash function (for set operations)
    return hash((interval[0], interval[1], interval[2]))

@njit
def compare_intervals(one, two) -> bool:
    # our ordering is based on start -> end -> label
    if one[1] == two[1]:
        if one[2] == two[2]:
            return one[0] < two[0]
        else:
            return one[2] < two[2]
    else:
        return one[1] < two[1]

@njit
def getRelation(A, B, constraints):
    relation = NONE

    epsilon = constraints[EPSILON]
    gap = constraints[GAP]

    if abs(B[1] - A[1]) <= epsilon:
        if abs(B[2] - A[2]) <= epsilon:
            relation = EQUALS
        elif B[2] - A[2] > epsilon:
            relation = LEFTMATCHES
    elif abs(B[2] - A[2]) <= epsilon and B[1] - A[1] > epsilon:
        relation = RIGHTMATCHES
    elif B[1] - A[1] > epsilon and A[2] - B[2] > epsilon:
        relation = CONTAINS
    elif A[2] - B[1] > epsilon and B[1] - A[1] > epsilon:
        relation = OVERLAPS
    elif abs(B[1] - A[2]) <= epsilon:
        relation = MEETS
    elif B[1] - A[2] > epsilon and (gap == 0 or B[1] - A[2] < gap):
        relation = FOLLOWS

    return relation

@njit
def remove_intervals_from_database(eseqdb, initial_support, min_seq):
    """
    This function removes all the intervals in the database that do not satisfy mininum support requirement.

    :param eseqdb: event sequence database
    :return eseqdb: event sequence database without intervals not meeting the requirement.
    """
    frequent_events = set()
    eseqdb_new = typed.List()
    for eseq in eseqdb:
        eseq_new = typed.List()
        for interval in eseq:
            if initial_support[interval[0]] >= min_seq:
                eseq_new.append(interval)
                frequent_events.add(interval[0])
        eseqdb_new.append(eseq_new)
    
    #eseqdb_new = [[interval for interval in List(eseq) if initial_support[interval[0]] >= min_seq] for eseq in List(eseqdb)]
    return eseqdb_new, frequent_events


#PART 3: Export, Import, Parameter settings

@njit
def make_constraints(argv, database):
    """
    
    0 - minSupPercent
    1 - minSup
    2 - epsilon
    3 - gap
    4 - timeoutseconds
    """
    constraints = np.zeros(9)
    constraints[MINSUPPERCENT] = float(argv[MINSUPPERCENT]) if (len(argv) > 0) else 0
    constraints[MINSUP] = constraints[MINSUPPERCENT] * len(database)
    constraints[EPSILON] = float(argv[1]) if (len(argv) > 1) else 0
    constraints[GAP] = float(argv[2]) if (len(argv) > 2) else np.inf
    if constraints[GAP] == -1:
        constraints[GAP] = np.inf
    constraints[TIMEOUTSECONDS] = int(argv[3]) if (len(argv) > 3) else 10000
    if constraints[TIMEOUTSECONDS] == -1:
        constraints[TIMEOUTSECONDS] = np.inf
    constraints[LEVEL] = float(argv[4]) if (len(argv) > 4) else np.inf
    if constraints[LEVEL] == -1:
        constraints[LEVEL] = np.inf
    constraints[PRINTF] = True if ((len(argv) > 5) and argv[5] == True) else False
    constraints[PRINTL] = True if ((len(argv) > 6) and argv[6] == True) else False
    constraints[FORGETTABLE] = True if ((len(argv) > 7) and argv[7] == True) else False
    return constraints
