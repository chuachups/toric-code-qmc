using Random
using HDF5
using Serialization
using Dates
using ArgParse

include("SSE.jl")
include("thermalize.jl")
include("simulate.jl")
include("write_file.jl")
include("run_job.jl")

# This function reads the command line for inputs and passes it to run_job()

function parse_command_line()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--filename"
            required = true
        "--L"
            arg_type = Int
            required = true
        "--beta"
            arg_type = Float64
            required = true
        "--J1"
            arg_type = Float64
            required = true
        "--J2"
            arg_type = Float64
            required = true
        "--C1"
            arg_type = Float64
            required = true
        "--C2"
            arg_type = Float64
            required = true
        "--C3"
            arg_type = Float64
            required = true
        "--C4"
            arg_type = Float64
            required = true
        "--basis"
            required = true
        "--task"
            required = true
        "--num-sweeps"
            arg_type = Int
            required = true

        "--prev-file"
            default = ""

        "--state"
            arg_type = UInt64   # Eventually change this
            default = 0

        "--M"
            arg_type = Int
            default = 0

        "--use-prev-rng"
            action = :store_true

        "--seed"
            arg_type = Int
            default = 0
    end

    return parse_args(s)
end

# Get input
args = parse_command_line()

# Pass to run_job()
run_job(
    args["filename"],
    args["L"],
    args["beta"],
    args["J1"],
    args["J2"],
    args["C1"],
    args["C2"],
    args["C3"],
    args["C4"],
    args["basis"],
    args["task"],
    args["num-sweeps"];
    prev_file=args["prev-file"],
    state=args["state"],
    M=args["M"],
    use_prev_rng=args["use-prev-rng"],
    seed=args["seed"],
)