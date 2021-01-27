# functions for HDF5
"""
    save_qmsim_data_hdf5(pop,hdf5file)

Save a genome structure to a file (HDF5), which will replace the old file.
"""
function save_qmsim_data_hdf5(pop,hdf5file)
   map = pop.map
   h5open(hdf5file,"w") do fid
      # map
      create_group(fid, "map")
      write(fid, "map/nchr", map.nchr)
      write(fid, "map/totalSNP", map.totalSNP)
      write(fid, "map/totalQTL", map.totalQTL)
      create_group(fid, "map/chr")
      for i=1:map.nchr
         n = @sprintf("%d",i)
         create_group(fid, "map/chr/"*n)
         write(fid, "map/chr/"*n*"/nLoci", map.chr[i].nLoci)
         write(fid, "map/chr/"*n*"/nSNP", map.chr[i].nSNP)
         write(fid, "map/chr/"*n*"/nQTL", map.chr[i].nQTL)
         write(fid, "map/chr/"*n*"/seqQTL", map.chr[i].seqQTL)
         write(fid, "map/chr/"*n*"/pos", map.chr[i].pos)
         write(fid, "map/chr/"*n*"/maxAllele", map.chr[i].maxAllele)
         write(fid, "map/chr/"*n*"/naQTL", map.chr[i].naQTL)
         write(fid, "map/chr/"*n*"/effQTL", map.chr[i].effQTL)
      end

      # genotypes
      create_group(fid, "genome")
      gr = fid["genome"]
      packedsz = get_packed_genome_size(pop.map)
      write(gr, "packedsz", packedsz)
      packed = zeros(Int8,packedsz)
      nind = length(pop.individual)
      write(gr, "nind", nind)
      dset = create_dataset(gr, "genotypes", datatype(Int8), ((packedsz,nind),(packedsz,-1)), chunk=(packedsz,1))
      for k=1:nind
         pack_genome!(pop.map, pop.individual[k],packed)
         dset[1:packedsz,k] = packed
      end
   end
   return nothing
end

"""
    add_qmsim_individual_hdf5(map,hdf5file)

Adds a genotyped individual to pre-created HDF5 file.
"""
function add_qmsim_individual_hdf5(map,hdf5file,individual)
   h5open(hdf5file,"r+") do fid
      gr = fid["genome"]
      packedsz = get_packed_genome_size(map)
      nind = read(gr, "nind")
      dset = gr["genotypes"]
      (dims,max_dims) = HDF5.get_extent_dims(dset)
      if nind>=dims[2]
         # doubled size
         HDF5.set_extent_dims(dset, (packedsz,dims[2]*2))
      end
      packed = zeros(Int8,packedsz)
      pack_genome!(map, individual, packed)
      dset[1:packedsz,nind+1] = packed
      q = gr["nind"]
      write(gr["nind"], nind+1)
   end
end

"""
    pop = read_qmsim_map_hdf5(hdf5file)

Read a map from `hdf5file`, which has been created by `save_qmsim_data`.
"""
function read_qmsim_map_hdf5(hdf5file)
   nchr = 0
   totalSNP = 0
   totalQTL = 0
   h5open(hdf5file,"r") do fid
      gr = fid["map"]
      nchr = read(gr,"nchr")
      totalSNP = read(gr,"totalSNP")
      totalQTL = read(gr,"totalQTL")
   end
   chr = Vector{QMSimChromosomeMap}(undef,nchr)
   h5open(hdf5file,"r") do fid
      gr = fid["map"]
      for i=1:nchr
         n = @sprintf("%d",i)
         cgr = fid["map/chr/"*n]
         nLoci = read(cgr,"nLoci")
         nSNP = read(cgr,"nSNP")
         nQTL = read(cgr,"nQTL")
         seqQTL = read(cgr,"seqQTL")
         pos = read(cgr,"pos")
         maxAllele = read(cgr,"maxAllele")
         naQTL = read(cgr,"naQTL")
         effQTL = read(cgr,"effQTL")
         chr[i] = QMSimChromosomeMap(nLoci,nSNP,nQTL,seqQTL,pos,maxAllele,naQTL,effQTL)
      end
   end
   return QMSimMap(nchr,totalSNP,totalQTL,chr)
end

"""
    genotype = read_qmsim_individual_hdf5(map,hdf5file,id)

Read a genotype of particular individual (`id`) from `hdf5file`, which has been created by `save_qmsim_data`.
"""
function read_qmsim_individual_hdf5(gmap,hdf5file,idx)
   chromosome_set = read_qmsim_chromosome_set_hdf5(gmap,hdf5file,idx)
   tbv = get_true_breeding_value(gmap,chromosome_set)
   return QMSimIndividualGenome(tbv,chromosome_set)
end

"""
    genotype = read_qmsim_chromosome_set_hdf5(map,hdf5file,id)

Read a set of chromosomes of particular individual (`id`) from `hdf5file`, which has been created by `save_qmsim_data`.
"""
function read_qmsim_chromosome_set_hdf5(gmap,hdf5file,idx)
   chromosome_set = Vector{QMSimChromosomeGenome}()
   h5open(hdf5file,"r") do fid
      gr = fid["genome"]
      #packedsz = read(gr, "packedsz")
      #nind = read(gr, "nind")
      dset = gr["genotypes"]
      packed = dset[1:end,idx]
      chromosome_set = unpack_genome(gmap,packed)
   end
   return chromosome_set
end

"""
    genome = read_qmsim_all_genotypes_hdf5(map,hdf5file)

Read all genotypes from `hdf5file`, which has been created by `save_qmsim_data`.
"""
function read_qmsim_all_genotypes_hdf5(gmap,hdf5file)
   # number of genotypes in the file
   nind = 0
   h5open(hdf5file,"r") do fid
      gr = fid["genome"]
      nind = read(gr, "nind")
   end

   # read the genotypes
   individual = Vector{QMSimIndividualGenome}(undef, nind)
   h5open(hdf5file,"r") do fid
      gr = fid["genome"]
      packedsz = read(gr, "packedsz")
      nind = read(gr, "nind")
      dset = gr["genotypes"]
      packed = zeros(Int8,packedsz)
      for i=1:nind
         packed .= dset[1:end,i]
         chromosome_set = unpack_genome(gmap,packed)
         tbv = get_true_breeding_value(gmap,chromosome_set)
         individual[i] = QMSimIndividualGenome(tbv, chromosome_set)
      end
   end
   return individual
end

"""
    data = read_qmsim_data_hdf5(hdf5file)

Read the enture data from `hdf5file`, which has been created by `save_qmsim_data`.
"""
function read_qmsim_data_hdf5(hdf5file)
   gmap = read_qmsim_map_hdf5(hdf5file)
   individual = QMSimData.read_qmsim_all_genotypes_hdf5(gmap,hdf5file)
   return QMSimPopulationGenome(gmap,individual)
end




function get_packed_genome_size(gmap)
   return 2*(gmap.totalSNP + gmap.totalQTL)
end

function pack_genome(gmap::QMSimMap,g::QMSimIndividualGenome)
   sz = get_packed_genome_size(gmap)
   packed = zeros(Int8,sz)
   pack_genome!(gmap,g,packed)
   return packed
end

function pack_genome(gmap::QMSimMap,chromosome_set::Vector{QMSimChromosomeGenome})
   sz = get_packed_genome_size(gmap)
   packed = zeros(Int8,sz)
   pack_genome!(gmap,chromosome_set,packed)
   return packed
end

function pack_genome!(gmap::QMSimMap,g::QMSimIndividualGenome,packed::Vector{Int8})
   pack_genome!(gmap,g.chr,packed)
   return nothing
end

function pack_genome!(gmap::QMSimMap,chr::Vector{QMSimChromosomeGenome},packed::Vector{Int8})
   fst = 1
   lst = 0
   for i=1:gmap.nchr
      fst = lst + 1
      lst = fst + gmap.chr[i].nLoci - 1
      packed[fst:lst] = chr[i].gp[1:end]
      fst = lst + 1
      lst = fst + gmap.chr[i].nLoci - 1
      packed[fst:lst] = chr[i].gm[1:end]
   end
   return nothing
end

function unpack_genome(gmap::QMSimMap,packed::Vector{Int8})
   chromosome_set = generate_chromosome_set(gmap)
   fst = 1
   lst = 0
   for i=1:gmap.nchr
      fst = lst + 1
      lst = fst + gmap.chr[i].nLoci - 1
      chromosome_set[i].gp[1:end] = packed[fst:lst]
      fst = lst + 1
      lst = fst + gmap.chr[i].nLoci - 1
      chromosome_set[i].gm[1:end] = packed[fst:lst]
   end
   return chromosome_set
end
