# functions for mating

function mate_genotypes(gmap,sire,dam)

end

# chelen: cM
function number_of_crosses(chrlen)
   # to Morgan
   lambda = chrlen/100.0
   return rand(Poisson(lambda))
end

# chelen: cM
function location_of_crossover(chrlen,pos,ncross)
   nloci = length(pos)
   loccross = zeros(Int,ncross)
   ntrials = 0
   while(true)
      ntrials = ntrials +1
      poscross = sort(rand(Uniform(0,chrlen),ncross))
      loccross .= 0
      k = 1
      for i=1:nloci
         if poscross[k]<=pos[i]
            loccross[k] = i
            if k==ncross
              break
            end
            k = k + 1
         end
      end
      if sum(loccross .> 0)==ncross
         break
      end
      if ntrials>10
         # too many trials: no crossover
         @warn "too many trials to find crossovers"
         return (0,Vector{Tuple{Int,Int}}(),Vector{Float64}())
      end
   end
   # which chromosomes
   gameteid = Vector{Tuple{Int,Int}}(undef,ncross)
   for i=1:ncross
      whichp = ifelse(rand()<0.5, 1, 2)
      whichm = ifelse(rand()<0.5, 3, 4)
      gameteid[i] = (whichp,whichm)
   end
   return adjusted_location_of_crossover(gameteid,loccross)
end

# removing double cross etc
function adjusted_location_of_crossover(gameteid,loccross)
   # effective location
   ncross = length(loccross)
   effective = zeros(Int,ncross)
   n = zeros(Int,4)
   prev = zeros(Int,4)
   which = zeros(Int,ncross)
   k = 0
   # remove double cross
   adj_ncross = 0
   for i=1:ncross
      # (1,3)->1, (1,4)->2, (2,3)->3, (2,4)->4
      which[i] = ifelse(gameteid[i][1]==1,gameteid[i][1]+gameteid[i][2]-3,gameteid[i][1]+gameteid[i][2]-2)
      n[which[i]] = n[which[i]] + 1
      if n[which[i]] == 2
         n[which[i]] = 0
         which[prev[which[i]]] = 0
         prev[which[i]] = 0
         which[i] = 0
         adj_ncross = adj_ncross - 1
      else
         prev[which[i]] = i
         adj_ncross = adj_ncross + 1
      end
   end

   # adjusted crossover
   if adj_ncross < ncross
      new_gameteid = Array{Tuple{Int,Int}}(undef,adj_ncross)
      new_loccross = zeros(Int, adj_ncross)
      k = 0
      for i=1:ncross
         if which[i]>0
            k = k +1
            new_gameteid[k] = gameteid[i]
            new_loccross[k] = loccross[i]
         end
      end
      return adj_ncross, new_gameteid, new_loccross
   else
      return ncross, gameteid, loccross
   end
end


# sample a gamete from an individual's paternal and maternal haplotypes
#   gamete1: copy of gp ----------
#   gamete2: copy of gp ----------
#   gamete3: copy of gm ==========
#   gamete4: copy of gm ==========
# crossover occurring between (1,3), (1,4), (2,3), or (2,4)
function generate_gamete(chrlen,pos,gp,gm)
   nloci = length(pos)
   ncross = number_of_crosses(pos[end])
   chosen = rand(1:4)
   if ncross==0 || nloci<ncross+1
      # no crossover; either gp or gm
      if chosen<=2
         return copy(gp)
      else
         return copy(gm)
      end
   else
      adj_ncross,gameteid,loccross = location_of_crossover(chrlen,pos,ncross)
      return generate_gamete_with_crossover(gameteid,loccross,chosen,gp,gm)
   end
end

function generate_gamete_with_crossover(gameteid,loccross,chosen,gp,gm)
   # make a copy
   gp_copy1 = copy(gp)
   gp_copy2 = copy(gp)
   gm_copy1 = copy(gm)
   gm_copy2 = copy(gm)
   tmp = copy(gp)
   tmp .= 0
   # pointers
   gam = [gp_copy1, gp_copy2, gm_copy1, gm_copy2]
   # simulate crossover
   ncross = length(loccross)
   for i=1:ncross
      # which chromosome
      (whichp,whichm) = gameteid[i]
      # swap genotypes between gam[whichp][loc:nloci] and gam[whichm][loc:nloci]
      loc = loccross[i]
      tmp[loc:end] = gam[whichp][loc:end]
      gam[whichp][loc:end] = gam[whichm][loc:end]
      gam[whichm][loc:end] = tmp[loc:end]
   end
   # choose the gamete
   return gam[chosen]
end

"""
    progeny = mating(map,sire,dam)

Given parents `sire` and `dam` as `QMSimIndividualGenome` and a map structure `QMSimMap`, this function produces one progeny as `QMSimIndividualGenome`.
"""
function mating(gmap::QMSimMap,sire::QMSimIndividualGenome,dam::QMSimIndividualGenome)
   chromosome_set = generate_chromosome_set(gmap)
   for i in 1:gmap.nchr
      chrlen = gmap.chr[i].pos[end]
      chromosome_set[i].gp .= generate_gamete(chrlen, gmap.chr[i].pos, sire.chr[i].gp, sire.chr[i].gm)
      chromosome_set[i].gm .= generate_gamete(chrlen, gmap.chr[i].pos, dam.chr[i].gp, dam.chr[i].gm)
   end
   tbv = get_true_breeding_value(gmap,chromosome_set)
   return QMSimIndividualGenome(tbv,chromosome_set)
end
