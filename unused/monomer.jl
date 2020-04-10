#
#  monomer.jl
#  collagen-model
#
#  Created by Christopher Revell on 30/03/2020.
#
#

module Monomer

mutable struct monomer
  label::Int64
  Ndomains::Int64
  domainLength::Float64
  pos::Array{Float64,1}
  v::Array{Float64,1}
end

export monomer

end
