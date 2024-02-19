# MGP prediction #

### How to install
---
It works on Linux operating system and has been tested on Ubuntu 16.04.
The runtimes below are generated using a computer with the specs (32 GB RAM, Intel(R) Core(TM) i7-4790 CPU @ 3.60 GHz) and internet of speed around 80 Mbit/s.

1. Clone the repository
```
git clone https://bitbucket.org/kaistsystemsbiology/mgp_prediction.git
```
2. Create and activate a conda environment
```
cd mgp_prediction
conda env create -f environment.yaml
source activate MGPprediction
```
which should install in about 2 minutes.

---
### Implementation
1. Create .csv file for flux (e.g.,exampleFlux.csv) and mutation (e.g.,exampleMutation.csv). Please referred to the format of each file in exampleInput directory.

2. Run MGP prediction. The following arguments are required: -o/--output_dir, -f/--flux, -mut/--mutation
	
	#### Example
	```
	mkdir output
	python runMGPprediction.py -o ./output -f ./exampleInput/exampleFlux.csv -mut ./exampleInput/exampleMutation.csv
	```
    The output file "predictedMGPs.csv" will be generated in output directory after computation.
    
    The run time for this source code on example input is about 90 seconds.
    
3. The output file is in .csv format and consists of four columns (i.e., Target gene, Target metabolite and its MetaNetX ID, and Target pathway) which represent information of predicted MGPs.

---
### Reference
Publication will be added here.

