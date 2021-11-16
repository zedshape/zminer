# DEPENDENCY: NUMPY
import sys

import numpy as np


def ZGenerator(meanE, sigE, meanI, sigI, meanH, sigH, meanL, sigL, vLen):
    database = ""

    # for each sequence,
    for i in range(vLen):
        # define the horizontal length of this sequence
        hLen = np.int(np.random.normal(meanH, sigH))
        # define the number of event and its label by picking it with replacement
        numEvt = np.int(np.random.normal(meanI, sigI))
        # select event... We do not accept the replacement to test everything correctly

        while True:
            if numEvt > 0:
                break
            else:
                numEvt = np.int(np.random.normal(meanI, sigI))

        localEvents = np.round(np.random.normal(meanE, sigE, numEvt)).astype(int)
        # for each events to put
        for event in localEvents:
            # define its size first
            # resampling if it does not meet the criteria ( 0 <= size <= length of seq)
            size = np.int(np.random.normal(meanL, sigL))
            while True:
                size = np.int(np.random.normal(meanL, sigL))
                if size > 0 and size <= hLen:
                    break
            # print(hLen, size)

            starttime = 0  # if size == hLen -> only possibility is zero
            if size != hLen:
                # find the place to put the event uniformly (it has to have enough space to put)
                # resampling if it does not meet the creteria ( 0 <= starttime)
                starttime = np.random.choice(hLen - size, 1)[0]
                while starttime < 0:
                    starttime = np.random.choice(hLen - size, 1)[0]

            database += str(i) + " " + str(event) + " " + str(starttime) + " " + str(starttime + size) + "\n"
    return database


def makeConstraints(argv):
    meanE = float(argv[0]) if (len(argv) > 0) else 0
    stdE = float(argv[1]) if (len(argv) > 1) else 0
    meanI = float(argv[2]) if (len(argv) > 2) else 0
    stdI = float(argv[3]) if (len(argv) > 3) else 0
    meanH = float(argv[4]) if (len(argv) > 4) else 0
    stdH = float(argv[5]) if (len(argv) > 5) else 0
    meanL = float(argv[6]) if (len(argv) > 6) else 0
    stdL = float(argv[7]) if (len(argv) > 7) else 0
    vLen = int(argv[8]) if (len(argv) > 8) else 0
    return meanE, stdE, meanI, stdI, meanH, stdH, meanL, stdL, vLen


def main(argv):
    # Arguments
    # ==========================================================================
    # [0]: meanE: mean of the distribution of event labels
    # [1]: stdE: standard deviation of the distribution of event labels
    # [2]: meanI: mean of the distribution of the number of event intervals
    # [3]: stdI: standard deviation of the distribution of the number of event intervals
    # [4]: meanH: mean of the distribution of e-sequence length
    # [5]: stdH: standard deviation of the distribution of e-sequence length
    # [6]: meanL: mean of the distribution of interval length
    # [7]: stdL: standard deviation of the distribution of interval length
    # [8]: vlen(l): size of e-sequence database
    # ==========================================================================

    print("########## Z-GENERATOR ##########")
    meanE, stdE, meanI, stdI, meanH, stdH, meanL, stdL, vLen = makeConstraints(argv)
    data = ZGenerator(meanE, stdE, meanI, stdI, meanH, stdH, meanL, stdL, vLen)
    filename = (
        "DATA_"
        + str(meanE)
        + "_"
        + str(stdE)
        + "_"
        + str(meanI)
        + "_"
        + str(stdI)
        + "_"
        + str(meanH)
        + "_"
        + str(stdH)
        + "_"
        + str(meanL)
        + "_"
        + str(stdL)
        + "_"
        + str(vLen)
    )
    text_file = open(filename + ".txt", "w")
    text_file.write(data)
    text_file.close()
    print("PARAMETERS: ", meanE, stdE, meanI, stdI, meanH, stdH, meanL, stdL, vLen)
    print("CHECK THE FILE: " + filename + ".txt")
    return


if __name__ == "__main__":
    main(sys.argv[1:])
