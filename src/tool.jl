# tools

function ==(a::QMSimChromosomeGenome,b::QMSimChromosomeGenome)
   if length(a.gp) != length(b.gp); return false; end
   if length(a.gm) != length(b.gm); return false; end
   rp = all( a.gp .== b.gp)
   rm = all( a.gm .== b.gm)
   if !rp || !rm
      return false
   end
   return true
end

function ==(a::Vector{QMSimChromosomeGenome},b::Vector{QMSimChromosomeGenome})
   nchr_a = length(a)
   nchr_b = length(b)
   if nchr_a != nchr_b
      return false
   end
   for i in 1:nchr_a
      ret = (a[i] == b[i])
      if !ret
         return false
      end
   end
   return true
end

function ==(a::QMSimIndividualGenome,b::QMSimIndividualGenome)
   if !isapprox(a.tbv,b.tbv)
      return false
   end
   return a.chr == b.chr
end

function ==(m1::QMSimMap,m2::QMSimMap)
   if m1.nchr != m2.nchr; return false; end
   if m1.totalSNP != m2.totalSNP; return false; end
   if m1.totalQTL != m2.totalQTL; return false; end
   try
      if length(m1.chr) != length(m2.chr); return false; end
      for i in 1:length(m1.chr)
         if m1.chr[i].nLoci != m2.chr[i].nLoci; return false; end
         if m1.chr[i].nSNP != m2.chr[i].nSNP; return false; end
         if m1.chr[i].nQTL != m2.chr[i].nQTL; return false; end
         if !all( m1.chr[i].seqQTL .== m2.chr[i].seqQTL); return false; end
         if !all( m1.chr[i].pos .== m2.chr[i].pos); return false; end
         if m1.chr[i].maxAllele != m2.chr[i].maxAllele; return false; end
         if !all( m1.chr[i].naQTL .== m2.chr[i].naQTL); return false; end
         if !all( m1.chr[i].effQTL .== m2.chr[i].effQTL); return false; end
      end
   catch e
      return false
   end
   return true
end

function ==(a::QMSimPopulationGenome,b::QMSimPopulationGenome)
   ret = (a.map == b.map)
   if !ret
      return false
   end
   try
      if length(a.individual) != length(b.individual)
         return false
      end
      for i in 1:length(a.individual)
         ret = (a.individual[i]==b.individual[i])
         if !ret
            return false
         end
      end
   catch e
      return false
   end
   return true   
end
