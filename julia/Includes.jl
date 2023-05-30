using DayAheadEmbedding
using Dates
using SQLite
using DataFrames, DataFramesMeta
using Plots
using Measures
using StatsPlots
using ProgressMeter
using Distributions
using Embeddings
using StatsBase
using LaTeXStrings, Latexify
using Chain
using ForwardDiff
using Interpolations
using BivariateCopulas
using DiscretizedCopulas
using Memoize, Random

import Base:
    Threads
    length

import SQLite:
    DB,
    load!

import SQLite.DBInterface:
    execute

function get_data(date::Date, hour::Int, side::String)
    @assert 1 ≤ hour ≤ 24
    @assert side ∈ ["Buy","Sell"]

    db = DB("./julia/TestData/DayAhead.db")

    execute(db, "SELECT * FROM '$(string(date)*"_"*string(hour)*"_"*side)'") |> DataFrame
end