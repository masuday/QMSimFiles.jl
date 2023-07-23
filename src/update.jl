# update maps

"""
    generate_effQTL!(gmap::QMSimMap, af::QMSimAlleleFrequency, dist_or_func, args...; expvar=1.0)

Generate and adjust QTL effects given allele frequency `af` and update the map `gmap`.
A distribution `dist` (defined in the `Distribution` package) or a function `func` shoud be provided to generate QTL effects optionally accompanied with arguments `args...`.
The generated effects will be scaled so that the expected genetic variance equals to `expvar`.
The function supports a map file with multiple traits.
"""
function generate_effQTL!(gmap::QMSimMap, af::QMSimAlleleFrequency, dist_or_func, args...; expvar=1.0)
   rand_effQTL!(gmap, dist_or_func, args...)
   obsvar = var_effQTL(gmap, af)
   scale = sqrt.(expvar ./ obsvar)
   adjust_effQTL!(gmap, af, scale)
   return nothing
end

function generate_effQTL!(chrmap::QMSimChromosomeMap, chraf::QMSimChromosomeAlleleFrequency, dist_or_func, args...; expvar=1.0)
   rand_effQTL!(chrmap, dist_or_func, args...)
   obsvar = var_effQTL(chrmap, chraf)
   scale = sqrt.(expvar ./ obsvar)
   adjust_effQTL!(chrmap, chraf, scale)
   return nothing
end

function generate_effQTL!(effQTL, naQTL, freq, dist_or_func, args...; expvar=1.0)
   rand_effQTL!(effQTL, dist_or_func, args...)
   obsvar = var_effQTL(effQTL, naQTL, freq)
   scale = sqrt.(expvar ./ obsvar)
   adjust_effQTL!(effQTL, naQTL, freq, scale)
   return nothing
end

"""
    rand_effQTL!(gmap::QMSimMap, dist)

Generate QTL effects following `dist` (defined by functions in `Distributions`) and update the map `gmap`.
"""
function rand_effQTL!(gmap::QMSimMap, dist::T) where {T<:Distribution}
   for c in eachindex(gmap.chr)
       rand_effQTL!(gmap.chr[c],dist)
   end
   return nothing
end

function rand_effQTL!(chrmap::QMSimChromosomeMap, dist::T) where {T<:Distribution}
   rand_effQTL!(chrmap.effQTL, dist)
   return nothing
end

function rand_effQTL!(effQTL, dist::T) where {T<:Distribution}
   for i in eachindex(effQTL)
      effQTL[i] = rand(dist)
   end
   return nothing
end

function rand_effQTL!(gmap::QMSimMap, func::Function, args...)
   for c in eachindex(gmap.chr)
       rand_effQTL!(gmap.chr[c],func,args...)
   end
   return nothing
end

function rand_effQTL!(chrmap::QMSimChromosomeMap, func::Function, args...)
   rand_effQTL!(chrmap.effQTL, func, args...)
   return nothing
end

function rand_effQTL!(effQTL::Matrix{Float64}, func::Function, args...)
   for i in eachindex(effQTL)
      effQTL[i] = func(args...)
   end
   return nothing
end

function rand_effQTL!(effQTL::Array{Float64,3}, func::Function, args...)
   (na, nQTL, ntr) = size(effQTL)
   for j in 1:nQTL
      for i in 1:na
        val = func(args...)
        for k in 1:ntr
            effQTL[i,j,k] = val[k]
         end
      end
   end
   return nothing
end
 
"""
    adjust_effQTL!(gmap::QMSimMap, af::QMSimAlleleFrequency, scale)
"""
function adjust_effQTL!(gmap::QMSimMap, af::QMSimAlleleFrequency, scale)
   for c in eachindex(gmap.chr)
       adjust_effQTL!(gmap.chr[c], af.chr[c], scale)
   end
   return nothing
end

function adjust_effQTL!(chrmap::QMSimChromosomeMap, chraf::QMSimChromosomeAlleleFrequency, scale)
   freq = collect_qtlfreq(chrmap, chraf)
   adjust_effQTL!(chrmap.effQTL, chrmap.naQTL, freq, scale)
   return nothing
end

function adjust_effQTL!(effQTL, naQTL, freq, scale)
   nQTL = size(effQTL,2)
   for i in 1:nQTL
      na = naQTL[i]
      avgeff = freq[1:na,i]' * effQTL[1:na,i]
      effQTL[1:na,i] .= (effQTL[1:na,i] .- avgeff)*scale
      effQTL[1:end,i] .= effQTL[1:end,i] .* (freq[1:end,i] .> 0.0)
   end
   return nothing
end

function adjust_effQTL!(effQTL::Array{Float64,3}, naQTL, freq, scale::Vector{Float64})
   ntr = size(effQTL,3)
   for i in 1:ntr
      adjust_effQTL!(view(effQTL,:,:,i), naQTL, freq, scale[i])
   end
   return nothing
end


"""
    var_effQTL(effQTL, naQTL, freq)
"""
function var_effQTL(gmap::QMSimMap, af::QMSimAlleleFrequency)
   if is_mtmap(gmap)
      ntr = size(gmap.chr[end].effQTL,3)
      v = zeros(ntr)
      for c in eachindex(gmap.chr)
         v .= v .+ var_effQTL(gmap.chr[c], af.chr[c])
      end
   else
      v = 0.0
      for c in eachindex(gmap.chr)
         v = v + var_effQTL(gmap.chr[c], af.chr[c])
      end
   end
   return v
end
 
function var_effQTL(chrmap::QMSimChromosomeMap, chraf::QMSimChromosomeAlleleFrequency)
   freq = collect_qtlfreq(chrmap, chraf)
   return var_effQTL(chrmap.effQTL, chrmap.naQTL, freq)
end

function var_effQTL(effQTL::Matrix{Float64}, naQTL, freq)
   v = 0.0
   nQTL = size(effQTL,2)
   for i in 1:nQTL
      na = naQTL[i]
      avgeff = freq[1:na,i]' * effQTL[1:na,i]
      vareff = freq[1:na,i]' * (effQTL[1:na,i] .- avgeff).^2
      if vareff > 0.0
         v = v + vareff
      end
   end
   return v
end

function var_effQTL(effQTL::Array{Float64,3}, naQTL, freq)
   ntr = size(effQTL,3)
   v = zeros(ntr)
   for i in 1:ntr
      v[i] = v[i] + var_effQTL(effQTL[:,:,i], naQTL, freq)
   end
   return v
end

"""
    qtlfreq = collect_qtlfreq(chrmap::QMSimChromosomeMap, chraf::QMSimChromosomeAlleleFrequency)
"""
function collect_qtlfreq(chrmap::QMSimChromosomeMap, chraf::QMSimChromosomeAlleleFrequency)
   freq = zeros(chrmap.maxAllele,chrmap.nQTL)
   n = 0
   for j in 1:chrmap.nLoci
      if chrmap.seqQTL[j]>0
         n = n + 1
         freq[:,n] .= chraf.freq[:,j]
      end
   end
   return freq
end

