# packages -
using CSV
using DataFrames
using LinearAlgebra
using Plots
using Distributions
using DifferentialEquations
using GlobalSensitivity
using Tables
using NumericalIntegration
using Optim
using Colors
using JLD2
using BSTModelKit
using Interpolations

# paths -
_PATH_TO_MODEL = joinpath(pwd(),"model")
_PATH_TO_DATA = joinpath(pwd(),"data")
_PATH_TO_TMP = joinpath(pwd(),"tmp")
_PATH_TO_ACTUAL_ENSEMBLE = joinpath(pwd(),"actual_ensemble_s_system")

# code -
include("src/Compute.jl")
include("src/Learn.jl")