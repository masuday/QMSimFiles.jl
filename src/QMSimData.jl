module QMSimData

using Printf
using Distributions
using HDF5

# chrosomome map
struct QMSimChromosomeMap
   nLoci::Int
   nSNP::Int
   nQTL::Int
   seqQTL::Vector{Int}
   pos::Vector{Float64}
   maxAllele::Int
   naQTL::Vector{Int}
   effQTL::Matrix{Float64}
end

# a set of chromosome maps
struct QMSimMap
   nchr::Int
   totalSNP::Int
   totalQTL::Int
   chr::Vector{QMSimChromosomeMap}
end

# genotypes in a chromosome
struct QMSimChromosomeGenome
   gp::Vector{Int8}
   gm::Vector{Int8}
end

# genotypes of an individual
struct QMSimIndividualGenome
   tbv::Float64
   chr::Vector{QMSimChromosomeGenome}
end

# genotypes
struct QMSimPopulationGenome
   map::QMSimMap
   individual::Vector{QMSimIndividualGenome}
end

import Base: ==
export QMSimChromosomeMap, QMSimMap, QMSimChromosomeGenome, QMSimIndividualGenome, QMSimPopulationGenome
export read_maps, read_genotypes, get_true_breeding_value, read_qmsim_data, 
       save_qmsim_data_hdf5, read_qmsim_map_hdf5, read_qmsim_data_hdf5,
       read_qmsim_individual_hdf5, add_qmsim_individual_hdf5, create_qmsim_map_hdf5,
       export_qmsim_individual_blupf90
export mating

include("io.jl")
include("hdf5io.jl")
include("tool.jl")
include("mating.jl")

end
