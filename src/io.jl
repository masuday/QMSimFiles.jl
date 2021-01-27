# functions for io

"""
    map = read_maps(snpmapfile,qtlmapfile,qtleffectfile)

Read three map and effect files, and create a single map object, `QMSimMap`.
"""
function read_maps(snpmapfile,qtlmapfile,qtleffectfile)
   # number of chromosomes
   maxChrM = get_number_of_chromosomes(snpmapfile)
   maxChrQ = get_number_of_chromosomes(qtlmapfile)
   maxChr = max(maxChrM,maxChrQ)
   chrMap = Vector{QMSimChromosomeMap}()

   # number of markers/QTLs
   nSNP = get_number_of_loci(snpmapfile,maxChr)
   nQTL = get_number_of_loci(qtlmapfile,maxChr)
   nLoci = nSNP .+ nQTL

   # number of QTL allele
   maxna,maxAllele = get_number_of_QTL_allele(qtleffectfile,maxChr)

   # allocate mutable chromosome maps for temporary use
   maps = Vector{QMSimChromosomeMap}(undef,maxChr)
   @inbounds for i in 1:maxChr
      na = maxAllele[i]
      maps[i] = QMSimChromosomeMap(nLoci[i],nSNP[i],nQTL[i],zeros(Int,nLoci[i]),zeros(Float64,nLoci[i]),na,zeros(Int,nQTL[i]),zeros(Float64,na,nQTL[i]))
   end

   # read positions
   load_SNP_QTL_maps!(snpmapfile,qtlmapfile,maps)

   # sort positions to merge SNP and QTL loci
   for m in maps
      sort_by_position!(m.pos, m.seqQTL)
   end

   # QTL effect
   load_QTL_effect!(qtleffectfile,maps)

   # a set of maps
   return QMSimMap(maxChr,sum(nSNP),sum(nQTL),maps)
end

function get_number_of_chromosomes(mapfile)
   maxChr = 0
   open(mapfile,"r") do io
      line = readline(io)
      while !eof(io)
         line = readline(io)
         items = split(line)
         c = tryparse(Int,items[2])
         maxChr = max(maxChr,c)
      end
   end
   return maxChr
end

function get_number_of_loci(mapfile,maxChr::Int)
   nLoci = zeros(Int,maxChr)
   open(mapfile,"r") do io
      line = readline(io)
      while !eof(io)
         line = readline(io)
         items = split(line)
         c = tryparse(Int,items[2])
         nLoci[c] = nLoci[c] + 1
      end
   end
   return nLoci
end

function load_SNP_QTL_maps!(snpmapfile,qtlmapfile,maps)
   nChr = length(maps)
   nLoci = zeros(Int,nChr)
   open(snpmapfile,"r") do io
      line = readline(io)
      while !eof(io)
         line = readline(io)
         items = split(line)
         c = tryparse(Int,items[2])
         pos = tryparse(Float64,items[3])
         nLoci[c] = nLoci[c] + 1
         maps[c].pos[nLoci[c]] = pos
      end
   end
   SeqNo = 0
   open(qtlmapfile,"r") do io
      line = readline(io)
      while !eof(io)
         line = readline(io)
         items = split(line)
         c = tryparse(Int,items[2])
         pos = tryparse(Float64,items[3])
         nLoci[c] = nLoci[c] + 1
         maps[c].pos[nLoci[c]] = pos
         SeqNo = SeqNo + 1
         maps[c].seqQTL[nLoci[c]] = SeqNo
      end
   end
   return nothing
end

function sort_by_position!(pos,seqQTL)
   perm = sortperm(pos)
   pos .= pos[perm]
   seqQTL .= seqQTL[perm]
   return nothing
end

function get_number_of_QTL_allele(qtleffectfile,maxChr)
   maxna = 0
   maxAllele = zeros(Int,maxChr)
   open(qtleffectfile,"r") do io
      line = readline(io)
      while !eof(io)
         line = readline(io)
         na = Int(round((length(line)-13)/14))
         maxna = max(maxna,na)

         items = split(line)
         c = tryparse(Int,items[2])
         maxAllele[c] = max(maxAllele[c],na)
      end
   end
   return maxna,maxAllele
end

function load_QTL_effect!(qtleffectfile,maps)
   maxChr = length(maps)
   nQTL = zeros(Int,maxChr)
   open(qtleffectfile,"r") do io
      line = readline(io)
      while !eof(io)
         line = readline(io)
         items = split(line)
         c = tryparse(Int,items[2])
         na = Int(round((length(line)-13)/14))
         nQTL[c] = nQTL[c] + 1
         # <--13--><--5--><--9--><--5--><--9--><--5--><--9-->
         for i=1:na
            fst = 13 + (i-1)*14 + 5 + 1
            lst = fst + 8
            maps[c].effQTL[i,nQTL[c]] = parse(Float64,strip(line[fst:lst]))
         end
         maps[c].naQTL[nQTL[c]] = na
      end
   end
end

"""
    g = read_genotypes(snpfile,qtlfile,map)

Read a set of genotypes from files given a map structure.
"""
function read_genotypes(snpfile,qtlfile,gmap)
   if !isfile(snpfile); throw(ArgumentError("file $(snpfile) not found")); end
   if !isfile(qtlfile); throw(ArgumentError("file $(qtlfile) not found")); end

   # temporary storage
   snp = zeros(Int8, gmap.totalSNP*2) 
   qtl = zeros(Int8, gmap.totalQTL*2)

   # number of animals
   nanimsnp = countlines(snpfile)-1
   nanimqtl = countlines(qtlfile)-1
   if nanimsnp != nanimqtl
      throw(ArgumentError("different number of individuals in snp and qtl files"))
   end
   if nanimsnp==0
      throw(ArgumentError("no individuals in the file"))
   end

   # initialization
   g = Vector{Vector{QMSimChromosomeGenome}}(undef,nanimsnp)

   # qtl data for individual k
   k = 0
   use_snpcode = is_snpcode_file(snpfile,gmap.totalSNP)
   if use_snpcode
      open(snpfile,"r") do io
         line = readline(io)
         while !eof(io)
            k = k + 1
            line = readline(io)
            text_to_snpcode!(line,gmap.totalSNP,snp)
            # snp code to haplotypes
            chromosome_set = generate_chromosome_set(gmap)
            convert_snpcode_to_haplotype!(snp,gmap,chromosome_set)
            g[k] = chromosome_set
         end
      end
   else
      open(snpfile,"r") do io
         line = readline(io)
         while !eof(io)
            k = k + 1
            line = readline(io)
            text_to_code!(line,gmap.totalSNP,snp)
            # snp code to haplotypes
            chromosome_set = generate_chromosome_set(gmap)
            convert_snpdata_to_haplotype!(snp,gmap,chromosome_set)
            g[k] = chromosome_set
         end
      end
   end

   # qtl data for individual k
   k = 0
   open(qtlfile,"r") do io
      line = readline(io)
      while !eof(io)
         line = readline(io)
         text_to_code!(line,gmap.totalQTL,qtl)
         # qtl code to haplotypes
         k = k + 1
         convert_qtldata_to_haplotype!(qtl,gmap,g[k])
      end
   end

   # true breeding value added to genome
   genotypes = Vector{QMSimIndividualGenome}(undef,nanimsnp)
   for i in 1:nanimsnp
      tbv = get_true_breeding_value(gmap,g[i])
      genotypes[i] = QMSimIndividualGenome(tbv,g[i])
   end

   return genotypes
end

function generate_chromosome_set(gmap)
   chromosome_set = Vector{QMSimChromosomeGenome}(undef, gmap.nchr)
   for i in 1:gmap.nchr
      chromosome_set[i] = QMSimChromosomeGenome(zeros(Int8,gmap.chr[i].nLoci),zeros(Int8,gmap.chr[i].nLoci))
   end
   return chromosome_set
end

function text_to_code!(line,totalSNPQTL,snpqtl)
   # <---7---><-1-><-[1]-><-1-><-[1]->
   for i in 1:totalSNPQTL*2
      k = 7 + (i-1)*2 + 1 + 1
      snpqtl[i] = parse(Int8,line[k:k])
   end
   return nothing
end

function text_to_snpcode!(line,totalSNPQTL,snp)
   # <---7---><-[1]-><-[1]->
   for i in 1:totalSNPQTL
      k = 7 + i
      snp[i] = parse(Int8,line[k:k])
   end
   return nothing
end

function is_snpcode_file(snpfile,totalSNP)
   n = 0
   open(snpfile,"r") do io
      line = readline(io)
      while !eof(io)
         line = readline(io)
         n = length(line)
         break
      end
   end
   if n<1
      error("file '$(snpfile)' not valid")
   elseif n > 7 + totalSNP
      return false
   else
      return true
   end
end

# paternal and maternal haplotype
function convert_snpcode_to_haplotype!(snp,gmap,chromosome_set)
   # decode
   loc = 0
   @inbounds for i in 1:gmap.nchr
      nLoci = gmap.chr[i].nLoci
      for j in 1:nLoci
         # non QTL
         if gmap.chr[i].seqQTL[j]==0
            loc = loc + 1
            if snp[loc] == 0
               chromosome_set[i].gp[j] = 0
               chromosome_set[i].gm[j] = 0
            elseif snp[loc]==2
               chromosome_set[i].gp[j] = 1
               chromosome_set[i].gm[j] = 1
            elseif snp[loc]==3
               chromosome_set[i].gp[j] = 0
               chromosome_set[i].gm[j] = 1
            elseif snp[loc]==4
               chromosome_set[i].gp[j] = 1
               chromosome_set[i].gm[j] = 0                  
            else
               error("illegal SNP code: $(snp[loc])")
            end                  
         end
      end
   end
   return nothing
end

function convert_snpdata_to_haplotype!(snp,gmap,chromosome_set)
   # decode
   loc = 0
   @inbounds for i in 1:gmap.nchr
      nLoci = gmap.chr[i].nLoci
      for j in 1:nLoci
         # non QTL
         if gmap.chr[i].seqQTL[j]==0
            loc = loc + 1
            chromosome_set[i].gp[j] = snp[loc]-1
            loc = loc + 1
            chromosome_set[i].gm[j] = snp[loc]-1
         end
      end
   end
   return nothing
end

function convert_qtldata_to_haplotype!(qtl,gmap,chromosome_set)
   # decode
   loc = 0
   @inbounds for i in 1:gmap.nchr
      nLoci = gmap.chr[i].nLoci
      for j in 1:nLoci
         # QTL
         if gmap.chr[i].seqQTL[j]>0
            loc = loc + 1
            chromosome_set[i].gp[j] = qtl[loc]
            loc = loc + 1
            chromosome_set[i].gm[j] = qtl[loc]
         end
      end
   end
   return nothing
end

function get_true_breeding_value(gmap::QMSimMap,chromosome_set::Vector{QMSimChromosomeGenome})
   tbv = 0.0
   for i in 1:gmap.nchr
      qlocus = 0
      nLoci = gmap.chr[i].nLoci
      for j in 1:nLoci
         # QTL
         if gmap.chr[i].seqQTL[j]>0
            qlocus = qlocus + 1
            tbv = tbv + gmap.chr[i].effQTL[chromosome_set[i].gp[j],qlocus]
            tbv = tbv + gmap.chr[i].effQTL[chromosome_set[i].gm[j],qlocus]
         end
      end
   end
   return tbv
end

function get_true_breeding_value(gmap::QMSimMap,individual_genome::QMSimIndividualGenome)
   return get_true_breeding_value(gmap,individual_genome.chr)
end

"""
    g = read_qmsim_data(snpmapfile,qtlmapfile,qtleffectfile,snpfile,qtlfile)

Read map, QTL effect, and genotype files to give a unified structure, `g::QMSimPopulationGenome`.
"""
function read_qmsim_data(snpmapfile,qtlmapfile,qtleffectfile,snpfile,qtlfile)
   gmap = read_maps(snpmapfile,qtlmapfile,qtleffectfile)
   genotypes = read_genotypes(snpfile,qtlfile,gmap)
   return QMSimPopulationGenome(gmap, genotypes)
end
