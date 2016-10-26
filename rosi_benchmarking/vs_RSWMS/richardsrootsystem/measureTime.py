import datetime
import subprocess
from subprocess import call

a = datetime.datetime.now()
call(["./test_rosi_vs_rswms"])
b = datetime.datetime.now()
print(b-a)
