module QMSimFiles

using Printf
using Random
using Statistics
using Distributions
using HDF5

# chrosomome map
struct QMSimChromosomeMap
   nLoci::Int
   nSNP::Int
   nQTL::Int
   seqQTL::Vector{Int}   # sequential number for QTL (length = nLoci)
   pos::Vector{Float64}  # physical position (length = nLoci)
   maxAllele::Int        # max. number of QTL alleles
   naQTL::Vector{Int}    # number of QTL alleles at each locus (length = nQTL)
   # effQTL[allele,seq] or effQTL[allele,seq,tr]
   effQTL::Union{Matrix{Float64},Array{Float64,3}} 
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
   tbv::Union{Float64,Vector{Float64}}
   chr::Vector{QMSimChromosomeGenome}
end

# genotypes
struct QMSimPopulationGenome
   map::QMSimMap
   individual::Vector{QMSimIndividualGenome}
end

# allele frequency on a chromosome
# SNP and QTL merged into a single matrix
struct QMSimChromosomeAlleleFrequency
   freq::Matrix{Float64}   # max(maxAllele,2) by nLoci
end

# a set of chromosome allele frequency
struct QMSimAlleleFrequency
   map::QMSimMap
   chr::Vector{QMSimChromosomeAlleleFrequency}
end

import Base: ==
export QMSimChromosomeMap, QMSimMap, QMSimChromosomeGenome, QMSimIndividualGenome, QMSimPopulationGenome
export read_maps, read_genotypes, get_true_breeding_value, read_qmsim_data, read_freq,
       save_qmsim_data_hdf5, read_qmsim_map_hdf5, read_qmsim_data_hdf5, get_qmsim_genotype_size,
       read_qmsim_individual_hdf5, add_qmsim_individual_hdf5, create_qmsim_map_hdf5,
       export_qmsim_individual_blupf90,
       is_mtmap, generate_effQTL!
export mating

include("io.jl")
include("hdf5io.jl")
include("tool.jl")
include("mating.jl")
include("update.jl")

end
