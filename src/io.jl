# functions for io

"""
    map = read_maps(snpmapfile,qtlmapfile,qtleffectfile; ntr=0)

Read three map and effect files, and create a single map object, `QMSimMap`.
With `ntr>0`, it returns the map object with QTL effects of `ntr` traits.
In this case, the QTL effects will be loaded for trait 1; zero QTL effects will be assigned to the other traits.
"""
function read_maps(snpmapfile,qtlmapfile,qtleffectfile; ntr=1)
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
      if ntr>1
         maps[i] = QMSimChromosomeMap(ntr,nLoci[i],nSNP[i],nQTL[i],zeros(Int,nLoci[i]),zeros(Float64,nLoci[i]),na,zeros(Int,nQTL[i]),zeros(Float64,na,nQTL[i],ntr))
      else
         maps[i] = QMSimChromosomeMap(1,nLoci[i],nSNP[i],nQTL[i],zeros(Int,nLoci[i]),zeros(Float64,nLoci[i]),na,zeros(Int,nQTL[i]),zeros(Float64,na,nQTL[i]))
      end
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
   # current restriction
   if maxna>9
      throw(ArgumentError("no support: 10 or more allele at QTL"))
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
            #maps[c].effQTL[i,nQTL[c]] = parse(Float64,strip(line[fst:lst]))
            parse_QTL_effect!(i,nQTL[c],maps[c].effQTL,strip(line[fst:lst]))
         end
         maps[c].naQTL[nQTL[c]] = na
      end
   end
end

function parse_QTL_effect!(allele,pos,effQTL::Matrix{Float64},str)
   effQTL[allele,pos] = parse(Float64,str)
end
function parse_QTL_effect!(allele,pos,effQTL::Array{Float64,3},str; tr=1)
   effQTL[allele,pos,tr] = parse(Float64,str)
end

"""
    af = read_freq(snpfreqfile,qtlfreqfile,map)

Read a set of allele frequencies from files given a map structure.
"""
function read_freq(snpfreqfile,qtlfreqfile,gmap)
   if !isfile(snpfreqfile); throw(ArgumentError("file $(snpfreqfile) not found")); end
   if !isfile(qtlfreqfile); throw(ArgumentError("file $(qtlfreqfile) not found")); end

   snpfreq = read_raw_snp_freq(snpfreqfile,gmap)
   qtlfreq = read_raw_qtl_freq(qtlfreqfile,gmap)

   af = Vector{QMSimChromosomeAlleleFrequency}(undef,gmap.nchr)
   for i in 1:gmap.nchr
      na = max(gmap.chr[i].maxAllele,2)
      af[i] = QMSimChromosomeAlleleFrequency(zeros(Float64,na,gmap.chr[i].nLoci))
   end

   snpIdx = Vector{Vector{Int64}}()
   qtlIdx = Vector{Vector{Int64}}()
   for i in 1:gmap.nchr
      push!(snpIdx, findall(gmap.chr[i].seqQTL .== 0))
      push!(qtlIdx, findall(gmap.chr[i].seqQTL .>  0))
   end

   # loading
   ns = 0
   nq = 0
   for i in 1:gmap.nchr
      # snp
      for k in 1:gmap.chr[i].nSNP
         ns = ns + 1
         af[i].freq[1,snpIdx[i][k]] = snpfreq[ns]
         af[i].freq[2,snpIdx[i][k]] = 1.0 - snpfreq[ns]
      end
      # qtl
      for k in 1:gmap.chr[i].nQTL
         nq = nq + 1
         for l in 1:length(qtlfreq[nq])
            af[i].freq[l,qtlIdx[i][k]] = qtlfreq[nq][l]
         end
      end
   end
   return QMSimAlleleFrequency(gmap,af)
end

function read_raw_snp_freq(snpfreqfile,gmap)
   if !isfile(snpfreqfile); throw(ArgumentError("file $(snpfreqfile) not found")); end
   n = 0
   snpfreq = zeros(Float64,gmap.totalSNP)
   open(snpfreqfile,"r") do io
      line = readline(io)
      while !eof(io)
         # <--10--><--2--><--4--><--5--><--1--><--8.6-->
         line = readline(io)
         n = n + 1
         allele = parse(Int,strip(line[10+2+4+1:10+2+4+5]))
         fst = 10+2+4+5+1 + 1
         lst = fst + 7
         snpfreq[n] = parse(Float64,strip(line[fst:lst]))
         if allele==2
            snpfreq[n] = 0.0
         end
      end
   end
   return snpfreq
end

function read_raw_qtl_freq(qtlfreqfile,gmap)
   if !isfile(qtlfreqfile); throw(ArgumentError("file $(qtlfreqfile) not found")); end
   n = 0
   qtlfreq = Vector{Vector{Float64}}()
   open(qtlfreqfile,"r") do io
      line = readline(io)
      while !eof(io)
         # <--10--><--2--><--4--><--11--><--1-->(<--4--><--1--><--8.6-->)
         # |--------------28-------------------|(|---------13----------|)
         line = readline(io)
         n = n + 1
         ni = Int(round((length(line)-28)/13))
         na = tryparse(Int,line[28+(ni-1)*13+1:28+(ni-1)*13+4])
         freq = zeros(na)
         for i in 1:ni
            allele = tryparse(Int,line[(28+(i-1)*13+1):(28+(i-1)*13+4)])
            fst = 28+(i-1)*13 + 6
            lst = 28+(i-1)*13 + 13
            freq[allele] = parse(Float64,strip(line[fst:lst]))
         end
         push!(qtlfreq,freq)
      end
   end
   return qtlfreq
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
      offset = get_leading_width_snpcode(snpfile)
      open(snpfile,"r") do io
         line = readline(io)
         while !eof(io)
            k = k + 1
            line = readline(io)
            text_to_snpcode!(line,gmap.totalSNP,snp,offset)
            # snp code to haplotypes
            chromosome_set = generate_chromosome_set(gmap)
            convert_snpcode_to_haplotype!(snp,gmap,chromosome_set)
            g[k] = chromosome_set
         end
      end
   else
      offset = get_leading_width_standard(snpfile)
      open(snpfile,"r") do io
         line = readline(io)
         while !eof(io)
            k = k + 1
            line = readline(io)
            text_to_code!(line,gmap.totalSNP,snp,offset)
            # snp code to haplotypes
            chromosome_set = generate_chromosome_set(gmap)
            convert_snpdata_to_haplotype!(snp,gmap,chromosome_set)
            g[k] = chromosome_set
         end
      end
   end

   # qtl data for individual k
   k = 0
   offset = get_leading_width_standard(qtlfile)
   open(qtlfile,"r") do io
      line = readline(io)
      while !eof(io)
         line = readline(io)
         text_to_code!(line,gmap.totalQTL,qtl,offset)
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

function text_to_code!(line,totalSNPQTL,snpqtl,offset)
   # <---7/8---><-1-><-[1]-><-1-><-[1]->
   for i in 1:totalSNPQTL*2
      k = offset + (i-1)*2 + 1 + 1
      snpqtl[i] = parse(Int8,line[k:k])
   end
   return nothing
end

# offset = 7 for v1, or 8 for v2
function text_to_snpcode!(line,totalSNPQTL,snp,offset)
   # <---7/8---><-[1]-><-[1]->
   for i in 1:totalSNPQTL
      k = offset + i
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
   elseif n > 8 + totalSNP
      return false
   else
      return true
   end
end

function get_leading_width_snpcode(snpfile)
   open(snpfile,"r") do io
      line = readline(io)
      while !eof(io)
         line = readline(io)
         if line[8]==' '
            return 8   # v2
         else
            return 7   # v1
         end
      end
   end
end

function get_leading_width_standard(file)
   open(file,"r") do io
      line = readline(io)
      while !eof(io)
         line = readline(io)
         if line[9]==' '
            return 8   # v2
         else
            return 7   # v1
         end
      end
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
   if is_mtmap(gmap)
      return _get_true_breeding_value_mt(gmap,chromosome_set)
   else
      return _get_true_breeding_value(gmap,chromosome_set)
   end
end

function get_true_breeding_value(gmap::QMSimMap,individual_genome::QMSimIndividualGenome)
   return get_true_breeding_value(gmap,individual_genome.chr)
end

function _get_true_breeding_value(gmap::QMSimMap,chromosome_set::Vector{QMSimChromosomeGenome})
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

function _get_true_breeding_value_mt(gmap::QMSimMap,chromosome_set::Vector{QMSimChromosomeGenome})
   ntr = size(gmap.chr[end].effQTL,3)
   tbv = zeros(ntr)
   for t in 1:ntr
      for i in 1:gmap.nchr
         qlocus = 0
         nLoci = gmap.chr[i].nLoci
         for j in 1:nLoci
            # QTL
            if gmap.chr[i].seqQTL[j]>0
               qlocus = qlocus + 1
               tbv[t] = tbv[t] + gmap.chr[i].effQTL[chromosome_set[i].gp[j],qlocus,t]
               tbv[t] = tbv[t] + gmap.chr[i].effQTL[chromosome_set[i].gm[j],qlocus,t]
            end
         end
      end
   end
   return tbv
end


"""
    g = read_qmsim_data(snpmapfile,qtlmapfile,qtleffectfile,snpfile,qtlfile; ntr=1)

Read map, QTL effect, and genotype files to give a unified structure, `g::QMSimPopulationGenome`.
"""
function read_qmsim_data(snpmapfile,qtlmapfile,qtleffectfile,snpfile,qtlfile; ntr=1)
   gmap = read_maps(snpmapfile,qtlmapfile,qtleffectfile, ntr=ntr)
   genotypes = read_genotypes(snpfile,qtlfile,gmap)
   return QMSimPopulationGenome(gmap, genotypes)
end

"""
    export_qmsim_individual_blupf90(map,id,individual,snpfile; append=true)

Export individual genotypes to a text file, `snpfile` with an integer code `id`, using a format supported by BLUPF90.
If the file does not exist, the function creates a new file.
If there is the same file name as `snpfile`, the funtion appends the data to existing file by default.
If you want to recreate the file, please remove the file before calling this function, or use `append=false`.
"""
function export_qmsim_individual_blupf90(gmap,id,individual,snpfile; append=true)
   # existing file
   if !append && isfile(snpfile)
      rm(snpfile,force=true)
   end
   # append
   nSNP = gmap.totalSNP
   open(snpfile,append=true) do io
      packed = pack_markrs(gmap,individual) .+ UInt8(48)
      write(io,@sprintf("%-8d",id)*String(packed)*"\n")
   end
end

# for BLUPF90 output
function pack_markrs(gmap::QMSimMap,individual::QMSimIndividualGenome)
   # storage
   packed = Vector{UInt8}(undef,gmap.totalSNP)
   fst = 1
   lst = 0
   for i in 1:gmap.nchr
      fst = lst + 1
      lst = fst + gmap.chr[i].nSNP - 1
      packed[fst:lst] = individual.chr[i].gp[ gmap.chr[i].seqQTL .==0 ]
      packed[fst:lst] .+= individual.chr[i].gm[ gmap.chr[i].seqQTL .==0 ]
   end
   return packed
end
