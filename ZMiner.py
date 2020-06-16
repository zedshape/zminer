import sys
import csv
from utils import *
from time import process_time

# ZMiner
# We did not use any external library for the fare comparison
class ZMiner:
    def __init__(self, database, constraints, forgettable=False, mode=None, oldFL=None):
        # F (arrangements information) and L (locations) are practially saved together.
        # F is 1-2 dimension of each level, L is 3-5 dimension for Z-Table and 2-3 dimension for Z-Arrangements
        # F is used as a key of locations as stated in the paper.
        self.FL = {}
        self.mode = mode
        self.oldFL = oldFL
        self.comparisoncount = 0
        self.totalfrequency = 0
        self.constraints = constraints
        self.database = database
        self.forgettable = forgettable

    def pruneWithMinsup(self):
        copiedEvents = self.database.initialSupport.copy()
        # remove event below threshold
        for label, support in self.database.initialSupport.items():
            if support < self.constraints["minSup"]:
                del copiedEvents[label]

        self.database.initialSupport = copiedEvents
        for seq in self.database.sequences:
            prunedSequences = []
            for event in seq.sequences:
                if (
                    event.label in self.database.initialSupport
                    and self.database.initialSupport[event.label]
                    >= self.constraints["minSup"]
                ):
                    prunedSequences.append(event)
            seq.sequences = prunedSequences
        return

    def createZTable(self):
        self.FL[2] = {}
        # each e-sequence id to generate next
        for S in self.database.sequences:
            # iterate every event pairs in the same sequence
            for s1 in S.sequences:
                for s2 in S.sequences:
                    # we keep the order not to make duplication
                    if s1 < s2:
                        R2 = getRelation(s1, s2, self.constraints)
                        if R2 != None:

                            E2 = (s1.label, s2.label)
                            # initialization
                            # F parts: frequent arrangements
                            # event pair hash table
                            if E2 not in self.FL[2]:
                                self.FL[2][E2] = {R2: {S.id: {s1: set()}}}
                            # relation part of frequent arrangements
                            # relation hash table
                            elif R2 not in self.FL[2][E2]:
                                self.FL[2][E2][R2] = {S.id: {s1: set()}}
                            # L parts: sequence location
                            # e-sequence hash table
                            elif S.id not in self.FL[2][E2][R2]:
                                self.FL[2][E2][R2][S.id] = {s1: set()}
                            # first interval hash table
                            elif s1 not in self.FL[2][E2][R2][S.id]:
                                self.FL[2][E2][R2][S.id][s1] = set()
                            # second interval hash table
                            # adding all addresses of s2 having a relation with s1
                            self.FL[2][E2][R2][S.id][s1].add(s2)

        for E2 in list(self.FL[2]):
            for R2 in list(self.FL[2][E2]):
                if len(self.FL[2][E2][R2]) < self.constraints["minSup"]:
                    del self.FL[2][E2][R2]
                else:
                    self.totalfrequency += 1

            if len(self.FL[2][E2]) == 0:
                del self.FL[2][E2]

    def initZ(self, s1, s2, ek, Rnew, Si, Rprev, Ak):
        candidates = self.getZTableSecondTable((s2.label, ek), Rnew, Si, s2)

        if candidates == False:
            return

        for s3 in candidates:
            # gap skipping
            if s3.start - s1.end >= self.constraints["gap"]:
                continue
            intervals = (s1, s2, s3)
            Rk = Rprev
            pair = (s1.label, s3.label)
            # We have some overlapped events to be checked. We check it with frequentRelations
            foundRelation = False

            for Rcand in self.FL[2][pair]:
                # laterevent and new event
                self.comparisoncount += 1
                relatArr = self.getZTableSecondTable(pair, Rcand, Si, s1)
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

            """
            lE = []
            eE = []

            for idx, event in enumerate(intervals[:-1]):
                # sequences[-1].start = lastStartTime
                if intervals[-1].start - event.end <= self.constraints["epsilon"]:
                    lE.append((event, idx))
                else:
                    eE.append(event)
                    #Here we only save position to fix follow relation
                    ch[idx] = 'b'
                    #eE.append(idx)
                fft = min(event.end, fft)
            """
            # We no longer explicitly save separate data structure, but separate them by indices
            # for time and memory efficiency

            #z = (tuple(eE), tuple(lE), intervals[-1], fft, ch)

            z = (intervals, min(s1.end, s2.end, s3.end))

            Ak[Rk][Si].append(z)

    def moveFollowsToEarlierIntervals(self, z, sk, k):
        """

        Here the commented codes are when we actually develop the data structure introduced in the paper
        In practice, we only save the indices for earlier intervals for memory and time efficiency

        lE = []
        eE = list(z[0])
        fft = min(sk.end, z[3])
        ch = z[-1] + [False]

        # lastStartTime = 0
        if sk.start - z[2].end <= self.constraints["epsilon"]:
            # Still later
            #Here we only save position to fix follow relation
            lE.append((z[2], k-2))
        else:
            # converted to eE
            eE.append(z[2])
            ch[k-2] = 'b'

        for event in z[1]:
            if sk.start - event[0].end <= self.constraints["epsilon"]:
                lE.append(event)
            else:
                eE.append(event[0])
                ch[event[1]] = 'b'
        """

        fft = min(sk.end, z[1])
        znew = (z[0]+(sk,), fft)
        return znew

    def extendZ(self, z, ek, Rnew, Si, Rprev, Ak, k):
        candidates = self.getZTableSecondTable((z[0][-1].label, ek), Rnew, Si, z[0][-1])
        if candidates == False:
            return
        foundRelation = None

        # Trivial, just "follow"
        #firstNewRelat = Rprev + "b" * len(z[0])

        for sk in candidates:
            #z[-1] is always the value of z.EarliestEndTime conceptually
            if sk.start - z[-1] >= self.constraints["gap"]:
                continue

            # initialize current relation set with previously defined follows

            #This is an implementation of earlier intervals
            #This string contain indices of earlierIntervals
            Rprevrst = Rprev[-(k-2):]
            Rk = []

            #indice set
            for idx, si in enumerate(z[0][:-1]):
                #print(len(z[0][:-1]), k, idx)
                foundRelation = False
                # We have some overlapped events to be checked. We check it with frequentRelations
                if (si.label, sk.label) not in self.FL[2]:
                    break

                #This is an implementation of earlier intervals
                #indices of earlier intervals are marked as 'b', so that we can just skip it
                if Rprevrst[idx] == 'b':
                    Rk.append('b') # automatic follow for current relation set
                    foundRelation = True
                    # if we find one, we do not need more iteration for this pair
                    continue

                for Rcand in self.FL[2][(si.label, sk.label)]:
                    # laterevent and new event
                    self.comparisoncount += 1
                    relatArr = self.getZTableSecondTable(
                        (si.label, sk.label), Rcand, Si, si
                    )

                    if (
                        relatArr != False and sk in relatArr
                    ):  # hash operation in -> O(1)
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

    #growZ function is affected by mode D
    def growZ(self, Eprev, Rprev, k):
        if self.timeout < process_time():
            return
        if k not in self.FL:
            self.FL[k] = {}

        prevSeq = self.FL[k - 1][Eprev][Rprev].keys()
        for ek in self.database.initialSupport.keys():
            # ex: if newLabel = ABC, then lastNewCand = BC
            Ek = Eprev + (ek,)

            if Ek[-2:] not in self.FL[2] or ((Eprev[0], ek) not in self.FL[2]):
                continue

            if Ek not in self.FL[k]:
                self.FL[k][Ek] = {}
            Ak = {}

            for Rnew in self.FL[2][Ek[-2:]]:
                # only care about previous bitmap + new bitmap events
                newSeq = self.FL[2][Ek[-2:]][Rnew].keys()

                for Si in set(prevSeq).intersection(set(newSeq)):

                    if k == 3:
                        for s1 in self.FL[2][Eprev][Rprev][Si]:
                            for s2 in self.FL[2][Eprev][Rprev][Si][s1]:
                                self.initZ(s1, s2, ek, Rnew, Si, Rprev, Ak)
                    else:
                        for z in self.getZArrangement(k - 1, Eprev, Rprev, Si):
                            self.extendZ(z, ek, Rnew, Si, Rprev, Ak, k)

            for Rk in Ak:
                if self.isFrequent(Ak, Rk) == True:
                    self.FL[k][Ek][Rk] = Ak[Rk]
                    self.totalfrequency += 1
                    self.growZ(Ek, Rk, k + 1)
                    if self.forgettable == True:
                        self.FL[k][Ek][Rk] = len(self.FL[k][Ek][Rk])
            if len(self.FL[k][Ek]) == 0:
                del self.FL[k][Ek]

    # O(1) checking of frequency -> second level is hashed
    def isFrequent(self, fp, relation):
        return len(fp[relation].keys()) >= self.constraints["minSup"]

    # get following sequence (relation, label, and sequence id) -> it will have multiple sequences
    def getZArrangement(self, k, E, R, Si):
        return self.FL[k][E][R][Si]

    def getZTableSecondTable(self, E, R, Si, s1):
        try:
            return self.FL[2][E][R][Si][s1]
        except:
            return False

    # Algorithm 1: Z-Miner
    def ZMiner(self):
        print("########## Z-MINER ##########")
        print("1-1. MINIMUM SUPPORT:", self.constraints["minSup"])
        print("1-2. EPSILON CONSTRAINT:", self.constraints["epsilon"])
        print("1-3. GAP CONSTRAINT:", self.constraints["gap"])
        print("1-4. TIMEOUT:", self.constraints["timeoutseconds"])
        print("1-5. LEVEL:", self.constraints["level"])

        print("2. NUMBER OF E-SEQUENCES:", len(self.database.sequences))

        t1 = process_time()
        # step 1: prune with duration constraint Cd
        # pseudo code line 1

        # we have frequent event labels information when we initialize ZMiner
        # Database and Sigma_freq are global variables of the class
        if self.constraints["minSup"] != 0:
            self.pruneWithMinsup()

        # Create Z-Table with pruned database and minSup
        # pruned database (D') and minSup are global variables of the class
        self.createZTable()
        self.timeout = t1 + self.constraints["timeoutseconds"]

        self.FL[3] = {}
        # run the iteration growZ to grow and check frequent arrangements iteratively
        for pair in self.FL[2]:
            for Rpair in self.FL[2][pair]:
                self.growZ(pair, Rpair, 3)

        if process_time() > self.timeout:
            print("TIMED OUT")
            return 0, 0, process_time() - t1, True

        t2 = process_time()

        print("3. TOTAL COMPARISON COUNTS:", self.comparisoncount)
        print("4. TOTAL FREQUENT ARRANGEMENTS:", self.totalfrequency)
        print("5. TOTAL TIME CONSUMED:", t2 - t1)

        return self.comparisoncount, self.totalfrequency, t2 - t1, False, self.FL


def main(argv):
    # Arguments
    # ==========================================================================
    # [0]: File name, file should be located into data directory
    # [1]: Minimum support threshold: it is a percentage form. default is 0 (retrieve all).
    # [2]: epsilon constraint: default is 0, meaning we do not allow flexibility
    # [3]: gap constraint: default is infinity, meaning we do not have any gap to be considered
    # [4]: timeout seconds of the algorithm. default is 10000
    # [5]: printF: export frequent arrangements to a csv file
    # [6]: printL: export frequent arrangements and their locations to a csv file
    # ==========================================================================

    # FILE LOAD
    dataname = argv[0]
    filename = str(dataname).split("/")[-1].split(".")[0]
    print("TEST WITH", dataname, "DATASET")
    tseq, tdis, tintv, aintv, avgtime, dataset = preprocess(dataname)
    print("TOTAL SEQUENCE:", tseq)
    print("TOTAL DISTINCT EVENTS:", tdis)
    print("TOTAL INTERVALS:", tintv)
    print("AVERAGE INTERVAL PER SEQ:", aintv)
    print("AVERAGE TIMESPAN:", avgtime)
    print("TEST WITH", dataname, "DATASET")
    database = Database(dataset)
    constraints = makeConstraints(argv[1:], dataset)
    algorithm = ZMiner(database, constraints)
    count, freq, timedelta, timeout, FL = algorithm.ZMiner()
    if constraints["printF"] == True:
        print("EXPORTING F ...")
        resultname = exportF(filename, FL, constraints)
        print("CHECK THE FILE: ", resultname)
    if constraints["printL"] == True:
        print("EXPORTING L ...")
        resultname = exportL(filename, FL, constraints)
        print("CHECK THE FILE: ", resultname)


if __name__ == "__main__":
    main(sys.argv[1:])
