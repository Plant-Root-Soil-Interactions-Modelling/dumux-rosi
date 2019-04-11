import py_rootbox as rb

rootsystem = rb.RootSystem()

# Open plant and root parameter from a file
name = "singleroot"  # "Anagallis_femina_Leitner_2010"
rootsystem.openFile(name, "../modelparameter/")

# Initialize
rootsystem.initialize()

# Simulate
rootsystem.simulate(50, True)

# Export final result (as vtp)
rootsystem.write("../results/singleroot.vtp")

print("done!")
