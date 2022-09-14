module Compat

export readinput

const HEADER = r"\s*(\d+)[\s,]*(\d+)[\s,]*(\d+)[\s,]*(\d+)"
const REGEX2 = r"P\s*=\s*-?\d*\.?\d*\s*V\s*=(\s*\d*\.?\d*)\s*E\s*=\s*(-?\d*\.?\d*)"i

function readinput(file)
    str = open(file, "r") do io
        read(io, String)
    end
    lines = Iterators.Stateful(split(str, r"\R"; keepempty=false))
    nv, nq, nm, n = ntuple(Returns(nothing), 4)
    for line in lines
        line = lstrip(line)
        if !isempty(line) || !startswith(line, '#')
            if match(HEADER, line) !== nothing
                nv, nq, nm, n = parse.(Int, match(HEADER, line).captures)
                break
            end
        end
    end
    if any(x === nothing for x in (nv, nq, nm, n))
        error("At least one of 'nv', 'nq', 'np', 'n' is not found in file!")
    end
    volumes = Vector{Float64}(undef, nv)
    energies = Vector{Float64}(undef, nv)  # Static
    frequencies = Array{Float64,3}(undef, nv, nq, nm)
    weights = Vector{Float64}(undef, nq)
    i = 1
    for line in lines
        line = strip(line)
        if isempty(line)
            continue
        end
        if '=' in line
            m = match(REGEX2, line)
            if m !== nothing
                volumes[i], energies[i] = map(Base.Fix1(parse, Float64), m.captures)
                i += 1  # Meet 1 volume
                j = 1
            end
            continue  # Next line
        end
        sp = split(line, r"[ \t]"; keepempty=false)
        if length(sp) == 3
            for k in 1:nm  # Note `k` is the index of mode
                line = first(iterate(lines))
                frequencies[i - 1, j, k] = parse(Float64, line)
            end
            j += 1
            continue
        end
        if occursin("weight", lowercase(line))
            break
        end
    end
    if i != nv + 1
        error("the number of volumes detected is not equal to what specified in head!")
    end
    j = 1
    for line in lines
        str = split(line, r"[ \t]"; keepempty=false)[end]
        weights[j] = parse(Float64, str)
        j += 1
    end
    return (n=n, v=volumes, e=energies, f=frequencies, w=weights)
end

end
