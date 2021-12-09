import csv
import sys
from time import process_time
from numba import typed, njit, types
from . import utils

def return_numba_tuple_type(size):
    return types.UniTuple(types.int64, size)

tuple_triple = return_numba_tuple_type(3)
listtype = types.ListType(tuple_triple)
dicttype_4 = types.DictType(types.int64, listtype)
dicttype_3 = types.DictType(types.int64, dicttype_4)
dicttype_2 = types.DictType(types.int64, dicttype_3)
dicttype = types.DictType(types.int64, dicttype_2)
listtype = types.ListType(tuple_triple)

@njit
def moveFollowsToEarlierIntervals(self, z, sk, k):
    fft = min(sk.end, z[1])
    znew = (z[0] + (sk,), fft)
    return znew

@njit
def extendZ(self, z, ek, Rnew, Si, Rprev, Ak, k):
    candidates = self.getZTableSecondTable((z[0][-1].label, ek), Rnew, Si, z[0][-1])
    if candidates == False:
        return
    foundRelation = None

    # Trivial, just "follow"
    # firstNewRelat = Rprev + "b" * len(z[0])

    for sk in candidates:
        # z[-1] is always the value of z.EarliestEndTime conceptually
        if sk.start - z[-1] >= self.constraints["gap"]:
            continue

        # initialize current relation set with previously defined follows

        # This is an implementation of earlier intervals
        # This string contain indices of earlierIntervals
        Rprevrst = Rprev[-(k - 2) :]
        Rk = []

        # indice set
        for idx, si in enumerate(z[0][:-1]):
            # print(len(z[0][:-1]), k, idx)
            foundRelation = False
            # We have some overlapped events to be checked. We check it with frequentRelations
            if (si.label, sk.label) not in self.FL[2]:
                break

            # This is an implementation of earlier intervals
            # indices of earlier intervals are marked as 'b', so that we can just skip it
            if Rprevrst[idx] == "b":
                Rk.append("b")  # automatic follow for current relation set
                foundRelation = True
                # if we find one, we do not need more iteration for this pair
                continue

            for Rcand in self.FL[2][(si.label, sk.label)]:
                # laterevent and new event
                self.comparisoncount += 1
                relatArr = self.getZTableSecondTable((si.label, sk.label), Rcand, Si, si)

                if relatArr != False and sk in relatArr:  # hash operation in -> O(1)
                    Rk.append(Rcand)
                    foundRelation = True
                    # if we find one, we do not need more iteration for this pair
                    break
            if foundRelation == False:
                break
        if foundRelation == False:
            continue

        znew = self.moveFollowsToEarlierIntervals(z, sk, k)

        # new final relationship
        Rk = Rprev + "".join(Rk) + Rnew

        if Rk not in Ak:
            Ak[Rk] = {Si: []}
        elif Si not in Ak[Rk]:
            Ak[Rk][Si] = []

        Ak[Rk][Si].append(znew)

@njit
def initZ(FL, constraints, s1, s2, ek, Rnew, Si, Rprev, Ak, comparison_count):
    candidates = getZTableSecondTable((s2.label, ek), Rnew, Si, s2)

    if candidates == False:
        return

    for s3 in candidates:
        # gap skipping
        if s3.start - s1.end >= constraints[utils.GAP]:
            continue
        intervals = (s1, s2, s3)
        Rk = Rprev
        pair = (s1.label, s3.label)
        # We have some overlapped events to be checked. We check it with frequentRelations
        foundRelation = False

        for Rcand in FL[2][pair]:
            # laterevent and new event
            comparison_count += 1
            relatArr = getZTableSecondTable(pair, Rcand, Si, s1)
            if relatArr != False and s3 in relatArr:  # hash operation in -> O(1)
                Rk += Rcand
                foundRelation = True
                # if we find one, we do not need more iteration for this pair
                break
        if foundRelation == False:
            continue

        Rk += Rnew

        if Rk not in Ak:
            Ak[Rk] = {Si: []}
        elif Si not in Ak[Rk]:
            Ak[Rk][Si] = []

        z = (intervals, min(s1.end, s2.end, s3.end))

        Ak[Rk][Si].append(z)
        
    
# O(1) checking of frequency -> second level is hashed
@njit
def isFrequent(fp, relation, constraints):
    return len(fp[relation].keys()) >= constraints[utils.MINSUP]

# get following sequence (relation, label, and sequence id) -> it will have multiple sequences
@njit
def getZArrangement(FL, k, E, R, Si):
    return FL[k][E][R][Si]

@njit
def getZTableSecondTable(FL, E, R, Si, s1):
    try:
        return FL[2][E][R][Si][s1]
    except:
        return False


# growZ function is affected by mode D
@njit
def growZ(Eprev, Rprev, k, FL, database, unique_labels, frequent_events, total_frequency, comparison_count, constraints, timeout):
    if k not in FL:
        FL[k] = typed.Dict.empty(types.int64, dicttype_2)

    prevSeq = FL[k - 1][hash(Eprev)][Rprev].keys()

    for ek in frequent_events:
        # ex: if newLabel = ABC, then lastNewCand = BC
        Ek = Eprev + (ek,)

        if Ek[-2:] not in FL[2] or ((Eprev[0], ek) not in FL[2]):
            continue

        if Ek not in FL[k]:
            FL[k][Ek] = {}
        Ak = {}

        for Rnew in FL[2][Ek[-2:]]:
            # only care about previous bitmap + new bitmap events - to get sequences
            newSeq = FL[2][Ek[-2:]][Rnew].keys()

            for Si in set(prevSeq).intersection(set(newSeq)):

                if k == 3:
                    for s1 in FL[2][Eprev][Rprev][Si]:
                        for s2 in FL[2][Eprev][Rprev][Si][s1]:
                            #FL, constraints, s1, s2, ek, Rnew, Si, Rprev, Ak, comparison_count
                            initZ(FL, constraints, s1, s2, ek, Rnew, Si, Rprev, Ak, comparison_count)
                else:
                    for z in getZArrangement(k - 1, Eprev, Rprev, Si):
                        extendZ(z, ek, Rnew, Si, Rprev, Ak, k)

        for Rk in Ak:
            if isFrequent(Ak, Rk) == True:
                FL[k][Ek][Rk] = Ak[Rk]
                total_frequency += 1
                growZ(Ek, Rk, k + 1, FL, database, unique_labels, frequent_events, total_frequency, constraints, timeout)
                if constraints[utils.FORGETTABLE] == True:
                    FL[k][Ek][Rk] = len(FL[k][Ek][Rk])
        if len(FL[k][Ek]) == 0:
            del FL[k][Ek]

@njit
def createZTable(database, FL, constraints, total_frequency):
    # each e-sequence id to generate next
    event_list = typed.List()
    for S in database:
        # iterate every event pairs in the same sequence
        for s1 in range(len(S)):
            for s2 in range(len(S)):
                # we keep the order not to make duplication
                if utils.compare_intervals(S[s1], S[s2]):
                    R2 = utils.getRelation(S[s1], S[s2], constraints)
                    if R2 != utils.NONE:
                        E2_unhashed = (S[s1][0], S[s2][0])
                        E2 = hash(E2_unhashed)
                        Ss1 = hash(S[s1])

                        # initialization
                        # F parts: frequent arrangements
                        # event pair hash table
                        if E2 not in FL[2]:
                            event_list.append(E2_unhashed)
                            FL[2][E2] = typed.Dict.empty(types.int64, dicttype_3)
                            FL[2][E2][R2] = typed.Dict.empty(types.int64, dicttype_4)
                            FL[2][E2][R2][s1] = typed.Dict.empty(types.int64, listtype)
                            FL[2][E2][R2][s1][Ss1] = typed.List.empty_list(tuple_triple)
                        # relation part of frequent arrangements
                        # relation hash table
                        elif R2 not in FL[2][E2]:
                            FL[2][E2][R2] = typed.Dict.empty(types.int64, dicttype_4)
                            FL[2][E2][R2][s1] = typed.Dict.empty(types.int64, listtype)
                            FL[2][E2][R2][s1][hash(S[s1])] = typed.List.empty_list(tuple_triple)
                        # L parts: sequence location
                        # e-sequence hash table
                        elif s1 not in FL[2][E2][R2]:
                            
                            FL[2][E2][R2][s1] = typed.Dict.empty(types.int64, listtype)
                            FL[2][E2][R2][s1][hash(S[s1])] = typed.List.empty_list(tuple_triple)
                        # first interval hash table
                        elif hash(S[s1]) not in FL[2][E2][R2][s1]:
                            FL[2][E2][R2][s1][hash(S[s1])] = typed.List.empty_list(tuple_triple)

                        # second interval hash table
                        # adding all addresses of s2 having a relation with s1
                        FL[2][E2][R2][s1][hash(S[s1])].append(S[s2])
    
    event_list_pruned = typed.List()
    for E2_unhashed in event_list:
        E2 = hash(E2_unhashed)
        for R2 in list(FL[2][E2]):
            if len(FL[2][E2][R2]) < constraints[utils.MINSUP]:
                del FL[2][E2][R2]
            else:
                total_frequency += 1

        if len(FL[2][E2]) == 0:
            del FL[2][E2]
        else:
            event_list_pruned.append(E2_unhashed)

    
    return event_list_pruned

def zminer(eseqdb, unique_labels, frequent_events, constraints):
    # F (arrangements information - level, event, relation) and 
    # L (locations - sequence, interval, candidate) are practially saved together.
    # F is 0-2 dimension of each level, L is 3-5 dimension for Z-Table and 2-3 dimension for Z-Arrangements
    # F is used as a key of locations as stated in the paper.
    FL = typed.Dict.empty(types.int64, dicttype)
    FL[2] = typed.Dict.empty(types.int64, dicttype_2)

    #self.mode = mode
    #self.oldFL = oldFL
    comparison_count = 0
    total_frequency = 0

    print("########## Z-MINER ##########")
    print("1-1. MINIMUM SUPPORT:", constraints[utils.MINSUP])
    print("1-2. EPSILON CONSTRAINT:", constraints[utils.EPSILON])
    print("1-3. GAP CONSTRAINT:", constraints[utils.GAP])
    print("1-4. TIMEOUT:", constraints[utils.TIMEOUTSECONDS])
    print("1-5. LEVEL:", constraints[utils.LEVEL])
    print("1-6. FORGETTABLE:", constraints[utils.FORGETTABLE])
    print("2. NUMBER OF E-SEQUENCES:", len(eseqdb))

    t1 = process_time()
    # Create Z-Table with pruned database and minSup
    # pruned database (D') and minSup are global variables of the class
    event_list = createZTable(eseqdb, FL, constraints, total_frequency)
    timeout = t1 + constraints[utils.TIMEOUTSECONDS]

    FL[3] = typed.Dict.empty(types.int64, dicttype_2)
    
    #def growZ(Eprev, Rprev, k, FL, database, unique_labels, frequent_events, total_frequency, comparison_count, constraints, timeout):
    # run the iteration growZ to grow and check frequent arrangements iteratively
    for E2 in event_list:
        for R2 in FL[2][hash(E2)]:
            print(E2)
            #Eprev, Rprev, k, FL, database, total_frequency, constraints, timeout
            growZ(E2, R2, 3, FL, eseqdb, unique_labels, frequent_events, total_frequency, comparison_count, constraints, timeout)

    if process_time() > timeout:
        print("TIMED OUT")
        return 0, 0, process_time() - t1, True, None

    t2 = process_time()

    print("3. TOTAL COMPARISON COUNTS:", comparison_count)
    print("4. TOTAL FREQUENT ARRANGEMENTS:", total_frequency)
    print("5. TOTAL TIME CONSUMED:", t2 - t1)

    return comparison_count, total_frequency, t2 - t1, False, FL


def run(dataname, argv):
    # Arguments
    # ==========================================================================
    # datename: File name, file should be located into data directory
    # [0]: Minimum support threshold: it is a percentage form. default is 0 (retrieve all).
    # [1]: epsilon constraint: default is 0, meaning we do not allow flexibility
    # [2]: gap constraint: default is infinity, meaning we do not have any gap to be considered
    # [3]: timeout seconds of the algorithm. default is 10000
    # [4]: printF: export frequent arrangements to a csv file
    # [5]: printL: export frequent arrangements and their locations to a csv file
    # [6]: forgettable: do not save any frequent arrangements to save memory - experiment only
    # ==========================================================================

    filename = str(dataname).split("/")[-1].split(".")[0]
    
    print("TEST WITH", filename, "DATASET")
    
    # Preprocess
    database = utils.load_file(dataname)
    tseq, tdis, tintv, aintv, avgtime, eseqdb, unique_labels, initial_support = utils.preprocess(database)
    constraints = utils.make_constraints(argv, eseqdb) 

    # Final e-sequence database is created 
    eseqdb, frequent_events = utils.remove_intervals_from_database(eseqdb, initial_support, constraints[utils.MINSUP])

    print("##### PREPROCESSING #####")
    print("TOTAL SEQUENCE:", tseq)
    print("TOTAL DISTINCT EVENTS:", tdis)
    print("TOTAL INTERVALS:", tintv)
    print("AVERAGE INTERVAL PER SEQ:", aintv)
    print("AVERAGE TIMESPAN:", avgtime)
    print("TEST WITH", dataname, "DATASET")

    # run z-miner with a processed algorithm
    count, freq, timedelta, timeout, FL = zminer(eseqdb, unique_labels, frequent_events, constraints)
