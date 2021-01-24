# BactPipes: Bioinformatics Pipeline

One Paragraph of project description goes here

## Getting Started

1. Download this repository
2. Ensure the requirements for the pipeline code itself(listed further down) are installed
2. Download the and install the following requisite programs. Currently we do not have an installer for packages, so you must do so manually.
    - prokka
    - orthofinder
    - TandemRepeatsFinder
    - kSNP3
    - scoary
    - BayesTraitsV3

Following are the current list of requirements for the pipeline. Currently you must install these yourself, we do hope to change this in the near future.
   #### Pipeline code requirements
    ##### Python Requirements
    1. Python 3
    2. Modules:
        * OS, shutil, time, sys, argparse, subprocess, multiprocessing, pandas, re, itertools
    ##### R Requirements
    1. Libraries:
        * ape, argparse, dplyr, phangorn, reshape2, stringr
   #### Exterior program requirements
    ##### Tandem Repeat Finder
    ##### RVDMiner
    ##### DisTAL
    ##### Prokka
    ##### Orthofinder
    ##### KSNP3 Parse
    ##### Scoary
    ##### Bayes Traits V3

## Calling

### Basics 
To call the pipeline please be inside of the overarching directory, currently named BioPipeline. Please use Python 3. From the command line you can call it as follows:

    python BactPipes_v1.0.py <Name of Project> <Number of Processors> [-B] <path/to/fasta/directory> [-S] <path/to/your.csv>


 ##### Notes for Calling
 1. '<>' Designates a file or directory call. Use either the relative or absolute path, do not include the '<>'
 2. '[]' Designates an "optional" flag. Call as shown without '[]', and the requisite following path. If neither is called, the pipeline will not run.
 3. <Name of project> must be a currently non-existant directory
 4. <Number of Processors> needs to be an integer value. Based on current busy server being used, if the value is less than or equal to 10, it will use that value for all processes. If it is greater than 10 it will use half of the value for each individual process
 6. [-P] is an optional argument to call first section of pipeline. Needs to be supplied with <path/to/fasta/directory>, which      needs to be the directory path, either a relative path to the current working directory or an absolute path. Currently        set to use .fa, .fna, .fasta, and .fas files. Case doesn't matter, only checks the file extension
 7. [-S] is an optional argument to call latter section of pipeline. Needs to be supplied with <path/to/your.csv>, which needs      to be a csv file. 
    * Initial checks will not verify the contents of this csv, but for proper processing this csv needs to have column names that match the rows of the created Bound Matrix. 
    * The rows of the Bound Matrix are based off names of the FASTA files given to the program. 
    * Periods in names of genomes will be removed in the Bound Matrix.
    * Provided csv file will be a presence/absence matrix in the following comma delimited format:
    
                            Name,Presence
                            X_albilineans__FIJ080,0
                            X_albilineans__HVO082,0
                            X_albilineans__LKA070,0
                            X_alfalfae__GEVRose07,0
                            X_arboricola__CFBP6827,0
                            X_arboricola__CFBP7634,0
                            X_arboricola__MEDVA37,0
                            X_arboricola_arracaciae_CFBP7407,0
                            X_arboricola_celebensis_NCPPB1630,0
                            X_arboricola_celebensis_NCPPB1832,0
                            X_arboricola_corylina_CFBP1159,0
                            X_arboricola_corylina_CFBP2565,1


## Deployment

Add additional notes about how to deploy this on a live system

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Authors

* **Alvaro L Perez-Quintero** - *Design and function* - Post-Doc Position, CSU Dept. of Bioagricultural Sciences and Pest Management
    [Research Gate Profile](https://www.researchgate.net/profile/Alvaro_L_Perez-Quintero)

* **Rex Steele**  - *Coding* - CSU CS Junior

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details
