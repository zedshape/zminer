# Z-Miner

Hello! **Z-Miner** is a novel algorithm for solving the temporal arrangement mining problem with fast speed and comsuming less memory.  Here we offer our new algorithm **Z-Miner** with some more tools we used in the paper submission and for the future usage and research. In this README, we describe how to use each tool here.

## Z-Miner

**Z-Miner** is run on a single file called `ZMiner.py` . For the research purpose, we did **not** use any external library that can speed up **Z-Miner** more, for the fare comparison. This file accepts few parameters as follows:
>  `python ZMiner.py filename minSup epsilon gap timeout level printF printL` 

Here is brief information about the parameters we need.
- [0]: File name
- [1]: Minimum support threshold: it is a percentage form. default is 0 (retrieve all).
- [2]: epsilon constraint: default is 0, meaning we do not allow flexibility
- [3]: gap constraint: default is infinity, meaning we do not have any gap to be considered
- [4]: timeout seconds of the algorithm. default is 10000
- [5]: maximum level (number of events) of arrangements. default is infinity
- [6]: printF: This option enables exporting frequent arrangements to a csv file. The file will contain {event set, relation set} for each frequent arrangement. The file will have the following form:
  `F_filename_support_epsilon_gap_timeout.csv` 
- [7]: printL: Unlike other algorithms, we can track all the information about the location in database where the frequent arrangement happens while managing small usage of the memory space. This option enables exporting frequent arrangements and their locations to a csv file. The file will contain {event set, relation set, sequence id, intervals} and will have all locations having a specific {event set, relation set}. The file will have the following form:
`L_filename_support_epsilon_gap_timeout.csv` 

Here is one example,
>  `python ZMiner.py data/ASL_BU_1.txt 0.01 10 10 300 5 True True` 

This code above will run **ZMiner** on ASL_BU_1 dataset, with 0.01 minimum support, epsilon 10, gap 10, timeout 300 sec, and level 5. It will also save all frequent arrangements information and corresponding locations. The saved files will have the following form in this case:
>  `F_ASL_BU_1_0.01_10_10_300.csv`  
>  `L_ASL_BU_1_0.01_10_10_300.csv` 

If you just want to run this algorithm with the basic options, you can just run it with filename and a minimum support! Note that you should also mention the directories from the place you run the code.

>  `python ZMiner.py data/ASL_BU_1.txt 0.01 ` 
>  
Which means we run **Z-Miner** with a minimum support 0.01, without epsilon, without a gap, and without printing any frequent arrangement information within the time limit 10,000.


## Z-Generator

**Z-Generator** is used to generate synthetic e-sequence datasets used in our experiment. It receives nine parameters and generates the datasets having the distribution properties that the user inputs. Here is a brief explanation about each parameter. Detailed information can be found in the paper.

- [0]: meanE: mean of the distribution of event labels
- [1]: stdE: standard deviation of the distribution of event labels
- [2]: meanI: mean of the distribution of the number of event intervals
- [3]: stdI: standard deviation of the distribution of the number of event intervals
- [4]: meanH: mean of the distribution of e-sequence length
- [5]: stdH: standard deviation of the distribution of e-sequence length
- [6]: meanL: mean of the distribution of interval length
- [7]: stdL: standard deviation of the distribution of interval length
- [8]: vlen(l): size of e-sequence database

**Note** that **Z-Generator** has a library dependency on `numpy` to generate distributions!
The way to generate a dataset is to simply input all the parameters in order described above. For example, if we want to generate a dataset similar to `HSPAN200` explained in the paper, we can run the generator as follows:

>  `python ZGenerator.py 100 2 30 3 200 20 20 6 100` 

Then it will create the dataset and save it as the following name:
>  `DATA_100.0_2.0_30.0_3.0_200.0_20.0_20.0_6.0_100.txt`

## Dependencies
### Z-Miner
For the fair and reproducible experiment **Z-Miner** did not use any external libraries and even built-in data structures that can help with speed-up in the procedure. To save file and check the time, we used the following built-in libraries.
> `sys`: to recieve arguments from console  
> `csv`: to save the result file (*F* and *L*) into csv  
> `time`: to check the runtime

As they are built-in libraries in pure Python distribution, **Z-Miner** can be executed with pure Python 3.7 without any additional installation.

### Z-Generator
**Z-Generator** used one external library `numpy` to create normal distributions. It needs to be installed before running **Z-Generator**. We also used `sys` with the same purpose as **Z-Miner**.
> `sys`: to recieve arguments from console  
> `numpy`: to make normal distributions with received parameters

