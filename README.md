# BactPipes: Bioinformatics Pipeline
Pipeline project for CSU Plant Science Department
Coded by Rex Steele, CSU CS Junior
Project overseen by Alvaro Pérez

### Basics 
To call the pipeline please be inside of the overarching directory, currently named BioPipeline. Then from the command line
you can call it as follows:

    python BactPuipes_v0.7.py <Name of Project> <Number of Processors> [-P] <path/to/fasta/directory> [-S] <path/to/your.csv>


 ##### Notes for Calling
 1. '<>' Designates a file or directory call. Use either the relative or absolute path, do not include the '<>'
 2. '[]' Designates an "optional" flag. Call as shown without '[]', and the requisite following path
 3. <Name of project> must be a currently non-existant directory
 4. <Number of Processors> needs to be an integer value. Based on current busy server being used, if the value is less than or      equal to 10, it will use that value for all processes. If it is greater than 10 it will use half of the value for each      individual process
 6. [-P] is an optional argument to call first section of pipeline. Needs to be supplied with <path/to/fasta/directory>, which      needs to be the directory path, either a relative path to the current working directory or an absolute path. Currently        set to use .fa, .fna, .fasta, and .fas files. Case doesn't matter, only checks the file extension
 6. [-S] is an optional argument to call latter section of pipeline. Needs to be supplied with <path/to/your.csv>, which needs      to be a csv file. 
        ⋅⋅*Initial checks will not verify the contents of this csv, but for proper processing this csv needs to have column names that match the rows of the created Bound Matrix. 
        ⋅⋅*The rows of the Bound Matrix are based off names of the FASTA files given to the program. 
        ⋅⋅*Periods in names of genomes will be removed in the Bound Matrix.

Scoary is an optional argument, not required. Needs a csv file, must end with .csv, case doesn't matter. CSV column names must match what is given by boundmatrix,
which will be based of the names of the FASTA files given to the program. Due to R script usage, periods in names will be removed when it comes to the boundmatrix.
