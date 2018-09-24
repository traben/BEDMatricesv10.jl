__precompile__()

"""
    module BEDMatrices1.0

A collection of tools for efficiently reading or indexing into .bed
files for calculation or other manipulations. Functionality primarily
revolves around the `immutable` `BEDMatrix`; see `BEDMatrix`
documentation for more details.

"""
module BEDMatrices1.0
using LinearAlgebra.apxy!

module Consts
include("constants.jl")
end

const NA_byte = Consts.NA_byte

include("bedmatrix.jl")
include("coltools.jl")
include("gwas.jl")

export BEDintomatrix,
    BEDintomatrix!,
    BEDMatrix,
    path,
    rownames,
    colnames,
    getflips,
    setflips!,
    NArep,
    getquartermap,
    hasNAs,
    countNAs,
    BEDdot,
    column_dist,
    column_sum,
    column_dot,
    column_NAsup_dot,
    column_sumabs2,
    column_norm,
    column_mean_and_std,
    column_mean_and_var,
    GWAS,
    gwas_writecsv

end
