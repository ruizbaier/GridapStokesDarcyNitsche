module GridapStokesDarcyNitsche
  using Gridap
  using Gridap.FESpaces
  using Gridap.ReferenceFEs
  using Printf
  #using Krylov

  include("DiscreteModels.jl")

  export generate_model_unit_square_biot_stokes
  export generate_model_unit_square
end
