# BactPipes: Bioinformatics Pipeline
Pipeline project for CSU Plant Science Department
Coded by Rex Steele, CSU CS Junior
Project overseen by Alvaro Pérez

To call the pipeline please be inside of the overarching directory, currently named BioPipeline. Then from the command line
you can call it as follows (not including <>):

    python IdunsPipeline_V0.2.py <Name of Project> <path/to/FASTA/Files/> <Number of Processors> <-s> <path/to/your.csv>


Name of project must be a currently non-existant directory in the current working directory of the program folder

Path to FASTA files needs to be the directory path relative to the current working directory. Currently set to use .fa, .fna, .fasta, and .fas files. Case doesn't matter.

Processors used, if set to 10 or less, will be set that number for all processes that allow setting of processors.
If greater than 10, it will be divided by two, then assigned to relevant processes

Scoary is an optional argument, not required. Needs a csv file, must end with .csv, case doesn't matter. CSV column names must match what is given by boundmatrix,
which will be based of the names of the FASTA files given to the program. Due to R script usage, periods in names will be removed when it comes to the boundmatrix.
