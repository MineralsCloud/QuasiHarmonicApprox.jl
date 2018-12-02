"""
# module Sampling

- Julia version: 1.0.2
- Author: qz
- Date: 2018-12-02

# Examples

```jldoctest
julia>
```
"""
module Sampling

export sample_brillouin_zone

function sample_brillouin_zone(q_weights::Vector{Float64}, quantity::Matrix{Float64})
    length(q_weights) == size(quantity, 2) && throw(DimensionMismatch)
    all(q_weights .>= 0) || throw(DomainError("All the values of the weights should be greater than 0!"))

    q_weights /= sum(q_weights)
    return quantity * q_weights
end

end