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
    print("*" * 80) 
    print(message)
    print("*" * 80)


def run_command(command):
    with open("../installdumux.log", "a") as log:
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
open('installdumux.log', 'w').close()

show_message("TO USE ONLY ON THE JÜLICH SUPERCOMPUTER \nDo not forget to check for loaded modules and environments")

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
        subprocess.run(["conda", "install", module]) 
        #raise Exception("Some python library "+str(output)[29:]+" are not found, consider loading SciPy-Stack")

#for mymodule in modules:
#	subprocess.run(["pip3", "install", mymodule]) 
      
show_message("(1/3) Step completed. All prerequistes found.")

#################################################################
#################################################################
## (2/3) Clone modules
#################################################################
#################################################################
# make a new folder containing everything
os.makedirs("./DUMUX", exist_ok=True)
os.chdir("DUMUX")

show_message("(2/3) Cloning repositories. This may take a while. Make sure to be connected to the internet...")

dune_version=2.7
dumux_version=3.0
# the core modules
for module in ['common', 'geometry', 'grid', 'localfunctions', 'istl' ]:
    if not os.path.exists("dune-{}".format(module)):
        git_clone('https://gitlab.dune-project.org/core/dune-{}.git'.format(module), "releases/{}".format(dune_version))
    else:
        print("-- Skip cloning dune-{} because the folder already exists.".format(module))
        

# the extensions
for module in ['foamgrid','grid-glue', 'alugrid', 'spgrid' ]:
    if not os.path.exists("dune-{}".format(module)):
        git_clone('https://gitlab.dune-project.org/extensions/dune-{}.git'.format(module), "releases/{}".format(dune_version))
    else:
        print("-- Skip cloning dune-{} because the folder already exists.".format(module))
    
# staging    
for module in ['uggrid' ]:
    if not os.path.exists("dune-{}".format(module)):
        git_clone('https://gitlab.dune-project.org/staging/dune-{}.git'.format(module), "releases/{}".format(dune_version))
    else:
        print("-- Skip cloning dune-{} because the folder already exists.".format(module))

# pybind
for module in ['pybindxi' ]:
    if not os.path.exists("dune-{}".format(module)):
        git_clone('https://github.com/dune-community/dune-{}.git'.format(module), branch = 'master')
    else:
        print("-- Skip cloning dune-{} because the folder already exists.".format(module))

# dumux and course
if not os.path.exists("dumux"):
    git_clone('https://git.iws.uni-stuttgart.de/dumux-repositories/dumux.git', "releases/{}".format(dumux_version))
else:
    print("-- Skip cloning dumux because the folder already exists.")

#if not os.path.exists("dumux-course"):
#    git_clone('https://git.iws.uni-stuttgart.de/dumux-repositories/dumux-course.git', "releases/{}".format(dumux_version))
#else:
 #   print("-- Skip cloning dumux-course because the folder already exists.")


# dumux-rosi
if not os.path.exists("dumux-rosi"):
    git_clone('https://github.com/Plant-Root-Soil-Interactions-Modelling/dumux-rosi.git', branch = 'master')
else:
    print("-- Skip cloning dumux-rosi because the folder already exists.")

# CPlantBox
if not os.path.exists("CPlantBox"):
    git_clone('https://github.com/Plant-Root-Soil-Interactions-Modelling/CPlantBox.git', branch = 'master')
    os.chdir("CPlantBox")
    subprocess.run(['cmake', '.']) 
    subprocess.run(['make']) 
    os.chdir("..")
else:
    print("-- Skip cloning CPlantBox because the folder already exists.")


show_message("(2/3) Step completed. All repositories have been cloned into a containing folder.")

#################################################################
#################################################################
## (3/3) Configure and build
#################################################################
#################################################################
show_message("(3/3) Configure and build dune modules and dumux using dunecontrol. This may take several minutes...")


#copy cplantbox so file to dumux-rosi folder
for f_name in os.listdir("CPlantBox"):
    if f_name.startswith('plantbox.cpython') and f_name.endswith('.so'):
        shutil.copyfile("CPlantBox/{0}".format(f_name), "dumux-rosi/{0}".format(f_name))
    

# run dunecontrol
if not os.path.isfile("cmake.opts"):
    shutil.move("dumux/cmake.opts", "cmake.opts")
else:
    print("-- The file cmake.opts already exists. The existing file will be used to configure dumux.")


subprocess.check_output(["./dune-common/bin/dunecontrol", "--opts=cmake.opts", "all"])
#./dune-common/bin/dunecontrol --opts=cmake.opts all
show_message("(3/3) Step completed. Succesfully configured and built CPlantBox, dune and dumux.")


show_message("to test installation, run n\ cd DUMUX/dumux-rosi/python/coupled \n python3 example7b_coupling.py")
 
