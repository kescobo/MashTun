module MashTun

import Base.length,
       Base.size,
       Base.start,
       Base.next,
       Base.done

export minhash,
       jaccarddist,
       mashdist,
       MASHSketch

using Bio.Seq
using DataStructures


type MASHSketch
    sketch::Vector{UInt64}
    kmersize::Int

    function MASHSketch(sketch::Vector, kmersize::Int)
        length(sketch) > 0 || error("Sketch cannot be empty")
        kmersize > 0 || error("Kmersize must be greater than 0")
        new(sketch, kmersize)
    end
end

function length(s::MASHSketch)
    return length(s.sketch)
end

function size(s::MASHSketch)
    return (length(s), s.kmersize)
end

function start(s::MASHSketch)
    return start(s.sketch)
end

function next(s::MASHSketch, state)
    return next(s.sketch, state)
end

function done(s::MASHSketch, state)
    return done(s.sketch, state)
end


function rchash(kmer::Kmer)
    k = minimum((kmer, reverse_complement(kmer)))
    return hash(k)
end


function kmerminhash(seq::BioSequence, kmerset, kmerhashes::Vector{UInt64}, k::Int, s::Int)
    typeof(kmerset) <: Set || typeof(kmerset) <: SortedSet || error("Kmerset must be a `Set` or `SortedSet`")
    length(kmerhashes) <= s || error("Kmerhashes cannot be larger than the set size")

    for kmer in each(DNAKmer{k}, seq)
        if length(kmerhashes) == 0
            if length(kmerset) < s
                push!(kmerset, hash(kmer[2]))
            elseif length(kmerset) == s
                kmerset = SortedSet(kmerset)
                for i in kmerset
                    push!(kmerhashes, pop!(kmerset))
                end
            end
        end

        if length(kmerhashes) == s
            h = hash(kmer[2])
            if h < kmerhashes[end]
                i = searchsortedlast(kmerhashes, h)
                if i == 0 && i != kmerhashes[1]
                    pop!(kmerhashes)
                    unshift!(kmerhashes, h)
                elseif h != kmerhashes[i]
                    pop!(kmerhashes)
                    insert!(kmerhashes, i+1, h)
                end
            end
        end

    end
    return (kmerset, kmerhashes)
end


function minhash(seq::BioSequence, k::Int, s::Int)
    kmerset = Set{UInt64}()
    kmerhashes = Vector{UInt64}()
    kmerset, kmerhashes = kmerminhash(seq, kmerset, kmerhashes, k, s)
    return MASHSketch(kmerhashes, k)
end

function minhash{T<:BioSequence}(seqs::Vector{T}, k::Int, s::Int)
    kmerset = Set{UInt64}()
    kmerhashes = Vector{UInt64}()

    for seq in seqs
        kmerset, kmerhashes = kmerminhash(seq, kmerset, kmerhashes, k, s)
    end
    return MASHSketch(kmerhashes, k)
end

function minhash{T<:BioSequence}(seqs::FASTAReader{T}, k::Int, s::Int)
    kmerset = Set{UInt64}()
    kmerhashes = Vector{UInt64}()
    for seq in seqs
        kmerset, kmerhashes = kmerminhash(seq.seq, kmerset, kmerhashes, k, s)
    end
    return MASHSketch(kmerhashes, k)
end


function jaccarddist(sketch1::MASHSketch, sketch2::MASHSketch)
    sketch1.kmersize == sketch2.kmersize || error("sketches must have same kmer length")

    d = length(setdiff(sketch1.sketch, sketch2.sketch))
    l = length(sketch1)
    return (l-d) / (l+ d)
end


function mashdist(k::Int, j::Float64)
    return 1/k * log(2j / (1+j))
end

function mashdist(sketch1::MASHSketch, sketch2::MASHSketch)
    length(sketch1) == length(sketch2) || error("sketches must be the same size")
    j = jaccarddist(sketch1, sketch2)
    k = sketch1.kmersize
    return mashdist(k, j)
end


end # module MashTun
