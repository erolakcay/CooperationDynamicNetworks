#!/usr/bin/env julia

### Code for Akcay, 2017 bioRxiv.
### Executable file. Type "./main.jl --help" in shell for usage

# the functions required to run the simulations are defined in MainFunctions.jl
include("MainFunctions.jl")

using ArgParse
using JLD

## Parameter object
type NetworkParam

    pn::Float64
    pr::Float64
    netsize::Int64
    generations::Int64
    b::Float64
    c::Float64
    d::Float64
    mu::Float64
    evollink::Bool
    mulink::Float64
    sigmapn::Float64
    sigmapr::Float64
    clink::Float64
    retint::Int64
    funnoevollink::String
    funevollink::String
    delta::Float64
    payfun::String
    replicates::Int64

    # Constructor. Takes keyword arguments to make it easier to call with command line arguments
    NetworkParam(;pn::Float64=0.5, pr::Float64=0.1, netsize::Int64=100, generations::Int64=100, b::Float64=1.0, c::Float64=0.5, d::Float64=0.0, mu::Float64=0.01, evollink::Bool=false, mulink::Float64=0.0, sigmapn::Float64=0.05, sigmapr::Float64=0.01,clink::Float64=0.0,retint::Int64=0, funnoevollink::String="Coauthor", funevollink::String="Coauthor", delta::Float64=1.0, payfun::String="Lin", replicates::Int64=1) = new(pn, pr, netsize, generations, b, c, d, mu, evollink, mulink, sigmapn, sigmapr, clink, retint, funnoevollink, funevollink, delta, payfun, replicates)

end

### Defining the parameters, default values, and help text.
function parse_commandline()
  s = ArgParseSettings()

  @add_arg_table s begin
    "--pn"
        help = "probability of inheriting connections, if evollink=true, gives the initial mean value"
        arg_type = Float64
        default = 0.5
    "--pr"
        help = "probability of random connections, if evollink=true, gives the initial mean value"
        arg_type = Float64
        default = 0.01
    "--netsize"
        help = "network size"
        arg_type = Int64
        default = 100
    "--generations"
        help = "number of generations to iterate network. Total time-steps given by n*gen"
        arg_type = Int64
        default = 100
    "--b"
        help = "benefit going out from cooperator"
        arg_type = Float64
        default = 1.0
    "--c"
        help = "cost of cooperation"
        arg_type = Float64
        default = 0.5
    "--d"
        help = "synergism coefficient"
        arg_type = Float64
        default = 0.0
    "--mu"
        help = "mutation rate for cooperator locus"
        arg_type = Float64
        default = 0.01
    "--evollink"
        help = "Whether linking probabilities are evolving or not"
        arg_type = Bool
        default = false
    "--mulink"
        help = "mutation rate for linking probabilities"
        arg_type = Float64
        default = 0.0
    "--sigmapn"
        help = "standard deviation for mutations in pn"
        arg_type = Float64
        default = 0.05
    "--sigmapr"
        help = "standard deviation for mutations in pr"
        arg_type = Float64
        default = 0.01
    "--clink"
        help = "cost per connection"
        arg_type = Float64
        default = 0.0
    "--file"
        help = "name of the file to save things in"
        arg_type = String
        default = "output"
    "--retint"
        help = "Interval at which the simulation saves output, in terms of number of death/birth events. Default = 0 reverts to saving things at intervals equal to network size"
        arg_type = Int64
        default = 0
    "--funnoevollink"
        help = "function to run for no evolution of linking probabilities: either networkitCoauthor or networkitPD, default is networkitCoauthor"
        arg_type = String
        default = "Coauthor"
    "--funevollink"
        help = "function to run for evolving linking probabilities: either networkitCoauthorEvolLink or networkitPDEvolLink, default is networkitCoauthorEvolLink"
        arg_type= String
        default = "Coauthor"
    "--payfun"
        help="which payoff function to use. Options Lin (1+delta*pi) or Exp ((1+delta)^pi), default is Lin"
        arg_type=String
        default= "Lin"
    "--delta"
        help = "The strength of selection (how much the total payoffs are weighted relative to the baseline payoff)"
        arg_type = Float64
        default = 1.0
    "--replicates"
        help = "number of replicates to run for the same parameter value"
        arg_type = Int64
        default = 1
    "--saveall"
        help = "whether to save all trajectories, or calculate the mean across all replicates"
        arg_type = Bool
        default = true
  end

  return parse_args(s)
end

# Function to parse input parameters and call up the main simulation function
function main()
  parsed_args = parse_commandline()
  sim_args = copy(parsed_args)
  delete!(sim_args, "file")
  delete!(sim_args, "saveall")

  sim_params = Dict()
  for (arg, val) in sim_args
      sim_params[parse(arg)] = val
  end
  params = NetworkParam(;sim_params...)

  if params.funnoevollink == "PD"
      funnoevollinkin = networkitPD
  elseif params.funnoevollink == "Coauthor"
      funnoevollinkin = networkitCoauthor
  end

  if params.funevollink == "PD"
      funevollinkin = networkitPDEvolLink
  elseif params.funevollink == "Coauthor"
      funevollinkin = networkitCoauthorEvolLink
  end

  if params.payfun == "Lin"
      payfunin = linpay
  elseif params.payfun == "Exp"
      payfunin = exppay
  end

  if parsed_args["saveall"] == true
      for i in 1:params.replicates
          typehist, pnhist, prhist, degreehist, payoffhist = runSim(;pn=params.pn, pr=params.pr, netsize=params.netsize, generations=params.generations, b=params.b, c=params.c, d=params.d, mu=params.mu, evollink=params.evollink, mulink=params.mulink, sigmapn=params.sigmapn, sigmapr=params.sigmapr, clink=params.clink, retint=params.retint, funnoevollink=funnoevollinkin, funevollink=funevollinkin, delta=params.delta,payfun=payfunin)

          file = ismatch(r"\.jl", parsed_args["file"]) ? parsed_args["file"] : parsed_args["file"]*"-"*lpad(string(i), length(digits(params.replicates)), "0")*".jld"

          save(file, "params", params, "typehist", typehist, "pn", pnhist, "pr", prhist, "degree", degreehist, "payoff", payoffhist)
      end
  elseif parsed_args["saveall"] == false
      typehistAve, pnhistAve, prhistAve, degreehistAve, payoffhistAve = (zeros(params.generations), zeros(params.generations), zeros(params.generations), zeros(params.generations), zeros(params.generations))
      for i in 1:params.replicates
          typehist, pnhist, prhist, degreehist, payoffhist = runSim(;pn=params.pn, pr=params.pr, netsize=params.netsize, generations=params.generations, b=params.b, c=params.c, d=params.d, mu=params.mu, evollink=params.evollink, mulink=params.mulink, sigmapn=params.sigmapn, sigmapr=params.sigmapr, clink=params.clink, retint=params.retint, funnoevollink=funnoevollinkin, funevollink=funevollinkin,delta=params.delta,payfun=payfunin)
          typehistAve += typehist
          pnhistAve += pnhist
          prhistAve += prhist
          degreehistAve += degreehist
          payoffhistAve += payoffhist
      end
      typehist = typehistAve/params.replicates
      pnhist = pnhistAve/params.replicates
      prhist = prhistAve/params.replicates
      degreehist = degreehistAve/params.replicates
      payoffhist = payoffhistAve/params.replicates

      file = ismatch(r"\.jl", parsed_args["file"]) ? parsed_args["file"] : parsed_args["file"]*".jld"

      save(file, "params", params, "typehist", typehist, "pn", pnhist, "pr", prhist, "degree", degreehist, "payoff", payoffhist)
  end

end

# Call up the main function.
main()
