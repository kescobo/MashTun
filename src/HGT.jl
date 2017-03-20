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

function getsubsketches(seq::FASTA.Reader, n::String, k::Int, s::Int, subseqlen::Int)
    gs = GenomeSketch(n)
    for contig in seq
        gs.sketches[identifier(contig)] = []
        for subsketch in subsketches(sequence(contig), identifier(contig), k, s, subseqlen)
            push!(gs.sketches[identifier(contig)], subsketch)
        end
    end
    return gs
end

function comparesubsketches(a::GenomeSketch, b::GenomeSketch)
    akeys = Dict(y => x for (x, y) in enumerate(sort([k for k in keys(a.sketches)])))
    bkeys = Dict(y => x for (x, y) in enumerate(sort([k for k in keys(b.sketches)])))
    function _it()
        for acontig in keys(a.sketches), bcontig in keys(b.sketches)
            (isempty(a.sketches[acontig]) || isempty(b.sketches[bcontig])) && continue
            produce(mashdistance(a.sketches[acontig], b.sketches[bcontig]), akeys[x], akeys[y])
        end
    end
    Task(_it)
end

function subsketchdistances(a::MinHashSketch, b::GenomeSketch)
    mingenomedist = 1
    for contig in keys(b.sketches)
        isempty(b.sketches[contig]) && continue
        mincontigdist = 1
        for bsub in b.sketches[contig]
            d = mashdistance(a, bsub.sketch)
            if d < mincontigdist
                mincontigdist = d
            end
        end

        if mincontigdist < mingenomedist
            mingenomedist = mincontigdist
        end
    end
    return mingenomedist
end
