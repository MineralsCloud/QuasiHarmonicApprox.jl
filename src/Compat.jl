module Compat

export readinput

const HEADER = r"\s*(\d+)[\s,]*(\d+)[\s,]*(\d+)[\s,]*(\d+)"
const REGEX2 = r"P\s*=\s*-?\d*\.?\d*\s*V\s*=(\s*\d*\.?\d*)\s*E\s*=\s*(-?\d*\.?\d*)"i

function readinput(file)
    open(file, "r") do io
        nv, nq, nm, n = ntuple(Returns(nothing), 4)
        for line in eachline(io)
            line = lstrip(line)
            if !isempty(line) || !startswith(line, '#')
                if match(HEADER, line) !== nothing
                    nv, nq, nm, n = parse.(Int, match(HEADER, line).captures)
                    mark(io)
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
        i, isqpoint = 1, false
        reset(io)
        for line in eachline(io)
            line = strip(line)
            if isempty(line)
                continue
            end
            if occursin("weight", lowercase(line))
                mark(io)
                break
            end
            if '=' in line
                m = match(REGEX2, line)
                if m !== nothing
                    volumes[i], energies[i] = parse.(Float64, m.captures)
                    i += 1  # Meet 1 volume
                    j = 0
                end
                continue  # Next line
            end
            if isqpoint
                for k in 1:nm  # Note `k` is the index of mode
                    frequencies[i - 1, j, k] = parse(Float64, line)
                    line = readline(io)
                end
                isqpoint = false
                continue
            end
            sp = split(line, r"[ \t]"; keepempty=false)
            if length(sp) == 3
                j += 1
                isqpoint = true
                continue
            end
        end
        if i != nv + 1
            error("the number of volumes detected is not equal to what specified in head!")
        end
        j = 1
        reset(io)
        for line in eachline(io)
            str = split(line, r"[ \t]"; keepempty=false)[end]
            weights[j] = parse(Float64, str)
            j += 1
        end
        return (n=n, v=volumes, e=energies, f=frequencies, w=weights)
    end
end

end
