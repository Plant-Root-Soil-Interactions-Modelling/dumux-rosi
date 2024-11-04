#!/usr/bin/env python3

"""
TO USE ONLY ON THE JÜLICH SUPERCOMPUTER
One click install script for dumux
created by the dumux team, adapted by Mona Giraud, Dirk Helmrich
"""
import os
import sys
import shutil
import subprocess


def show_message(message):
    print("*" * 120) 
    print(message)
    print("*" * 120)


def run_command(command):
    with open("../installDumuxRosi.log", "a") as log:
        popen = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        for line in popen.stdout:
            log.write(line)
            print(line, end='')
        for line in popen.stderr:
            log.write(line)
            print(line, end='')
        popen.stdout.close()
        popen.stderr.close()
        return_code = popen.wait()
        if return_code:
            print("\n")
            message = "\n    (Error) The command {} returned with non-zero exit code\n".format(command)
            message += "\n    If you can't fix the problem yourself consider reporting your issue\n"
            message += "    on the mailing list (dumux@listserv.uni-stuttgart.de) and attach the file 'installdumux.log'\n"
            show_message(message)
            sys.exit(1)


def git_clone(url, branch=None):
    clone = ["git", "clone"]
    if branch:
        clone += ["-b", branch]
    result = run_command(command=[*clone, url])


# clear the log file
open('installdumux_JSC.log', 'w').close()

show_message("TO USE ONLY ON THE JÜLICH SUPERCOMPUTER \nDo not forget to check for loaded modules and environments")
show_message("""
run the following commands in the terminal:
module --force purge
module load Stages/2024
module load GCC
module load ParaStationMPI
module load Python 
module load OpenCV
module load CMake
module load SciPy-Stack
module load VTK

""")

#################################################################
#################################################################
## (1/3) Check some prerequistes
#################################################################
#################################################################
# programs = ['wget', 'git', 'gcc', 'g++', 'cmake', 'pkg-config','clang', 'gfortran','python3.8'] 
#

# # check some prerequistes
error = []
# for program in programs:
#     try:
#         output2 = subprocess.run([program, "--version"], capture_output=True)
#     except FileNotFoundError:
#         error.append(program)
        
        
programs = ['Java', 'Boost', 'VTK','SciPy-Stack','ParaStationMPI/5.4.9-1'] 
show_message("(1/3) (a) Checking cluster prerequistes: " + " ".join(programs) + "...")
for program in programs:
    output = subprocess.run(["module","avail", program], capture_output=True,shell=True)
    if ('No module(s) or extension(s) found!' in str(output)):
        error.append(program)
        
if len(error) > 0:
    print("Program(s) {0} has/have not been found. try running module load {0}".format(" ".join(error)))
    raise Exception('import modules')


import pip

# check some prerequistes
modules = ['numpy', 'scipy', 'matplotlib', 'mpi4py', 'ipython'] 
show_message("(1/3) (b) Checking python prerequistes: " + " ".join(modules) + "...")

for module in modules :
    output = subprocess.run(["python3" ,"-m","pip", "show", module], capture_output=True)
    #subprocess.run(["conda", "install", module]) 
    if "WARNING: Package(s) not found" in str(output) :
        subprocess.run(["pip3", "install", module]) 
        #raise Exception("Some python library "+str(output)[29:]+" are not found, consider loading SciPy-Stack")

#for mymodule in modules:
#	subprocess.run(["pip3", "install", mymodule]) 
      
show_message("(1/3) Step completed. All prerequistes found.")

#################################################################
#################################################################
## (2/3) Clone modules
#################################################################
#################################################################

url= 'https://raw.githubusercontent.com/Plant-Root-Soil-Interactions-Modelling/dumux-rosi/master/installdumux.py'
subprocess.run(["wget", url, "-P", "."], check=True)

subprocess.run(["python3","installdumux.py"], check=True)

os.chdir("./dumux")

subprocess.run(["python3","dumux/bin/installexternal.py","ug","spgrid","foamgrid"], check=True)

# dumux-rosi
if not os.path.exists("dumux-rosi"):
    subprocess.run(['git', 'clone', '--depth','1', '-b', 'master', 'https://github.com/Plant-Root-Soil-Interactions-Modelling/dumux-rosi.git'])
else:
    print("-- Skip cloning dumux-rosi because the folder already exists.")

# CPlantBox
if not os.path.exists("CPlantBox"):
    subprocess.run(['git', 'clone', '--depth','1', '-b', 'master', 'https://github.com/Plant-Root-Soil-Interactions-Modelling/CPlantBox.git'])
else:
    print("-- Skip cloning CPlantBox because the folder already exists.")

#################################################################
#################################################################
## (3/3) Configure and build
#################################################################
#################################################################
show_message("(3/3) Configure and build dune modules and dumux using dunecontrol. This may take several minutes...")

os.chdir("CPlantBox")


subprocess.run(['git', 'submodule', 'update',  '--recursive', '--init'])
subprocess.run(['cmake', '.']) 
subprocess.run(['make'])  
os.chdir("..")

# run dunecontrol
subprocess.run(["./dune-common/bin/dunecontrol", "--opts=dumux-rosi/cmake.opts", "all"])

print("(3/3) Step completed. Succesfully configured and built CPlantBox, dune and dumux.")


print("to test installation, run \n cd dumux/dumux-rosi/python/coupled \n python3 example7b_coupling.py")
 
