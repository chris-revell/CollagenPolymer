#
#  AdaptTimestep.jl
#  CollagenPolymer
#
#  Created by Christopher Revell on 06/07/2020.
#
#

module AdaptTimestep

# Import Julia packages
using LinearAlgebra
using StaticArrays

# Import local program modules


function adaptTimestep(F,W,magsF,magsW,σ)

    magsF .= norm.(F)
    maxF = maximum(magsF)

    magsW .= norm.(W)
    maxW = maximum(magsW)

    Δt = min(σ^2/(4.0*maxW^2),σ/(2.0*maxF))

    return Δt

end

export adaptTimestep

end
