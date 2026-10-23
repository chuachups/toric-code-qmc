# This is a list of jobs, to be manually edited by the user. Each job is stored as a dictionary with the following keys:
# 1. Parameters: L, beta, J_1, J_2, basis
# 2. Run number
# 3. Task: "thermalize" or "simulate"
# 4. Number of sweeps
# 5. Optional for thermalization warm start: L, beta, J_1, J_2, basis, run number of the previous run to use as a warm start

jobs = [
    {
        "L" : 2,
        "beta" : 1,
        "J_1" : 1,
        "J_2" : 1,
        "basis" : "Z",

        "run" : 1,

        "task" : "simulate",
        "num_sweeps" : 1000000
    }
]

# Example: cold start thermalize, simulate, and warm start thermalize respectively
# jobs = [
#   {
#       "L" : 2,
#       "beta" : 1,
#       "J_1" : 1,
#       "J_2" : 1,
#       "basis" : "Z",
#
#       "run" : 0,
#       
#       "task" : "thermalize",
#       "num_sweeps" : 100000
#   },
#   {
#       "L" : 2,
#       "beta" : 1,
#       "J_1" : 1,
#       "J_2" : 1,
#       "basis" : "Z",
#
#       "run" : 0,
#       
#       "task" : "simulate",
#       "num_sweeps" : 500000
#   },
#   {
#       "L" : 2,
#       "beta" : 1,
#       "J_1" : 1,
#       "J_2" : 1,
#       "basis" : "X",
#
#       "run" : 0,
#       
#       
#       "task" : "thermalize",
#       "num_sweeps" : 100000,
#
#       "previous": {
#           "L" : 2,
#           "beta" : 1,
#           "J_1" : 1,
#           "J_2" : 1,
#           "basis" : "Z",
#           "run" : 0
#       }
#   },
# ]