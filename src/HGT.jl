module HGT

using MashTun: MinHashSketch, minhash, mashdistance
using Bio.Seq
using DataStructures: SortedSet

immutable SubSketch
    sketch::MinHashSketch
    contig::String
    pos::Int
end

type GenomeSketch
    sketches::Dict{String,Vector{SubSketch}}
    name::String

    function GenomeSketch(n::String)
        new(Dict{String,Vector{MinHashSketch}}(), n)
    end
end

b1s = GenomeSketch("Brachybacterium alimentarium 738.10")
b2s = GenomeSketch("Brachybacterium alimentarium 862.8")
cs = GenomeSketch("Corynebacterium sp JB4")
as = GenomeSketch("Alkalibacterium kapii FAM208 38")

Base.push!(g::GenomeSketch, s::SubSketch) = push!(g.sketches, s)

function subsketches(seq::BioSequence, n::String, k::Int, s::Int, subseqlen::Int)
    contiglen = length(seq)
    function _it()
        for i in 1:Int(round(subseqlen / 2)):contiglen
            if i + subseqlen <= contiglen
                sketch = SubSketch(minhash(seq[i:i+subseqlen], k, s), n, i)
                produce(sketch)
            end
        end
    end
    Task(_it)
end

function getsubsketches{T<:BioSequence}(seq::FASTAReader{T}, n::String, k::Int, s::Int, subseqlen::Int)
    gs = GenomeSketch(n)
    for contig in seq
        gs.sketches[contig.name] = []
        for subsketch in subsketches(contig.seq, String(contig.name), k, s, subseqlen)
            push!(gs.sketches[contig.name], subsketch)
        end
    end
    return gs
end

function comparesubsketches(a::GenomeSketch, b::GenomeSketch)
    akeys = Dict(y => x for (x, y) in enumerate(sort([k for k in keys(a.sketches)])))
    bkeys = Dict(y => x for (x, y) in enumerate(sort([k for k in keys(b.sketches)])))
    function _it()
        for acontig in keys(a.sketches), y in keys(b.sketches)
            (isempty(a.sketches[x]) || isempty(b.sketches[y])) && continue
            produce(mashdistance(a.sketches[x], b.sketches[y]), akeys[x], akeys[y])
        end
    end
    Task(_it)
end

end # module HGT
