#Bioinformatics Pipeline
import os #Necessary for directory checks
import time #To breakup processes and make output readable
import sys #Allows exiting of code in case of irreconcilable error

#Variables needed globally for Pipeline
projName = None
pipePath = None

#Verify directory, get name of project
def pipeStart():
    dirPath = os.getcwd()
    projName = raw_input("Enter name of Project here: ")
    pipePath = dirPath + "/" + projName
    print pipePath
    while os.path.isdir(pipePath) == True:
        projName = raw_input("This directory already exists and is invalid. Please enter a valid \
                    project name for a new directory, or use CTRL + C to exit out: ")
        pipePath = dirPath + "/" + projName
    else:
        print "Creating directory" + pipePath + "..."
        time.sleep(2)
        try:
            os.mkdir(pipePath)
        except OSError:
            print("Creation of the directory %s failed. Exiting..." % pipePath)
            sys.exit()
        else:
            print("Succesfully created the directory %s." % pipePath)


pipeStart()
