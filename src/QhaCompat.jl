module QhaCompat

export readinputfile

function readinputfile(file)
    str = open(file, "r") do io
        read(io, String)
    end
    regex0 = r"\s*(\d+)[\s,]*(\d+)[\s,]*(\d+)[\s,]*(\d+)"
    regex1 = r"P\s*=\s*-?\d*\.?\d*\s*V\s*=(\s*\d*\.?\d*)\s*E\s*=\s*(-?\d*\.?\d*)"i
    for line in split(str, r"\R")
        if !startswith(line, r"\s*#")
            if match(regex0, line) !== nothing
                nv, nq, nm, _ = parse.(Int, match(regex0, line).captures)
            end
        end
    end
    i, j = 0, 0
    for line in Iterators.filter(!isempty, gen)
        if '=' in line
            m = match(regex1, line)
            if m !== nothing
                volumes[i], static_energies[i] = match.groups()
                i += 1
                j = 0
            continue
            end
        end
        sp = line.split()
        if length(sp) == 3
            for k in range(modes_per_q_point_amount)  # Note `k` is the index of mode, like `j`, not count like `i`.
                line = next(gen)
                frequencies[i - 1, j, k] = line

            j += 1
            continue
            end
        end

        if "weight" in line.lower()
            break
        end
    end
end

end
