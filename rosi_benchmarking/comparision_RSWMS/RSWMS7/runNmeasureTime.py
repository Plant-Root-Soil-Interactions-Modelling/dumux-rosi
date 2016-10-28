import datetime
import subprocess
from subprocess import call

a = datetime.datetime.now()
call(["./rswms"])
b = datetime.datetime.now()
print("")
print("SIMULATION TIME: ")
print(b-a)
