import datetime
import subprocess
from subprocess import call

a = datetime.datetime.now()
call(["./rswms"])
b = datetime.datetime.now()
print(b-a)
