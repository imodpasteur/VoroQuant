# High resolution chromosome analysis

High resolution chromosome analysis is a python program allowing to analyse simulated and real data. Requirements and examples to run the program are given below.

## Requirements


The program is written and has been tested in python (version 2.7). The program runs with a standard computer with 64Gb of memory (you may be able to process some of the data if your computer has only 32Gb).


We recommend the use of conda in order to install libraries and to run our scripts. 


## Installation

Installation requires typically 10 minutes.
First, download this github repository. 
Then, open a terminal and go to the downloaded directory.
A '.yml' file will allow you to load all necessary libraries. 
To create the conda environment, run this command line:
`conda env create -f voroquant.yml`
Then, to load the conda environment, run this command line:
`source activate voroquant`





## Data

Simulated and real data can be downloaded here: [link]







## Simulations analysis

### How to run


To run analysis on simulation data, first, load the conda environment: `source activate voroquant`. 

A list of files with simulated chromosomes is available in `./data/simulations/` folder. The simulated chromosome length is 248Mb. The program will split the chromosome in multiple parts (using positions saved in `./randomCutting`) folder to simulate the variation of chromosome lengths in the genome and to simulate sister chromatin exchange. 

To process one file, please run: `python simulationAnalysis.py file_to_process.txt`, where `file_to_process.txt` argument correspond to the path of the file to process.

For each split denoted by `[N]`, the program outputs four files:
* `Coord_X_timepoint_XX.txt_6split=[N]_GT.csv`: The localization table of the simulated chromosome splitted.
* `Coord_X_timepoint_XX.txt_6split=[N]_GT_withNoise.csv`: The localization table with additive localization errors and background.
* `Coord_X_timepoint_XX.txt_6_split=[N]_segmentation.csv`: The localization table after segmentation (voronoi-based segmentation)
* `Coord_X_timepoint_XX.txt_6split=[N].result`: A text file containing 4 values: {surface (nm^2), volume (nm^3), surface/volume normalized (in [0,1]), gyration radius (nm^2)}



### Demo

The total processing time is around 5 minutes. It needs 64Gb of memory for this example. If the computer has only 32Gb of memory, the program will stop before the end and results will be uncomplete.

To run analysis on the chromosome at path `./data/simulations/0_cohesin_Mb/Coord_1_timepoint_12.txt` (0 loop/Mb), run the following commands: 
* `source activate voroquant`. 
* `python simulationAnalysis.py ./data/simulations/0_cohesin_Mb/Coord_1_timepoint_12.txt`. 
Then, the following 6 files are created:
* `Coord_1_timepoint_12.txt_6times=1split=0_rank=2_cluster1.result`
* `Coord_1_timepoint_12.txt_6times=1split=1_rank=2_cluster1.result`
* ...




## Real data analysis


### How to run

To run analysis on real data, first, load the conda environment: `source activate voroquant`. 

In `./data/realData` folder, we provide localization tables with and without Auxin. Each localization table is associated to a `crop` file containing the coordinates of a regin of interest in the background (manually selected). Local background density is used to segment chromosome.

To do the analysis on the one chromosome, please, run : `python realDataAnalysis.py file_to_process.csv`, where `file_to_process.csv` argument correspond to the path of the file to process.


After processing, the program outputs four files:
* `2a.csv_rank=2_cluster1.csv`: The localization table after segmentation (voronoi-based segmentation)
* `2a.csv.result`: A text file containing 4 values: {surface (nm^2), volume (nm^3), surface/volume normalized (in [0,1]), gyration radius (nm^2)}


### Demo

The total processing time is around 1 minute. It needs 32Gb of memory for this example but 64Gb of memory is necessary to process the full dataset.

To run analysis on the file at path `./data/realData/-Auxin/2a.csv`, please run: 
`source activate voroquant`. 
`python realDataAnalysis.py ./data/realData/-Auxin/2a.csv`. 
Then, the following file is created `2a.csv_rank=2_cluster1.result` containing quantifications (density, gyration, ...). Also, the results shown in the terminal is:
* `('gyration radius:', 1.51, 'um')`
* `('volume concav alpha=100:', 5503423073, 'nm^3')`
* `('volume concav alpha=500:', 9593397527, 'nm^3')`
* `('SMOOTHNESS (%) :', 57.3)`
 



## Localization tables

Each localization table saved as a `.csv` file can be opened and rendered using [ZOLA-3D](https://github.com/imodpasteur/ZOLA-3D) Fiji plugin. 



