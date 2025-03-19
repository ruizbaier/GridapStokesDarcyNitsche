module StokesDarcyNitscheAccuracyTestsTHRT
  using Gridap
  using GridapStokesDarcyNitsche
  import Gridap: ∇
  using Printf

  (Gridap.CellData.cross)(a::Gridap.CellData.CellField,
  b::Gridap.CellData.SkeletonPair{<:Gridap.CellData.CellField}) = Operation(cross)(a,b)

  Id  = TensorValue(1,0,0,1)
  const μ   = 1.0
  const α   = 1.0
  const κ   = 1.0
  const δ   = 1.0

  # manufactured solutions
  uS_ex(x)  = VectorValue(-cos(π*x[1])*sin(π*x[2]), sin(π*x[1])*cos(π*x[2]))
  pS_ex(x)  = cos(π*x[1])*cos(π*x[2])
  uD_ex(x)  = VectorValue(x[1]*x[2]+x[1], -x[2]^2*x[1])
  pD_ex(x)  = sin(π*x[1])*sin(π*x[2])

  function solve_stokes_darcy_TH_RT(model)

    labels = get_face_labeling(model)
    tags   = Gridap.Geometry.get_face_tag(labels,num_dims(model))

    # Method of manufactured solutions (RHS)

    σS_ex(x) = 2*μ*ε(uS_ex)(x)- pS_ex(x)*Id
    fS(x) = -(∇⋅σS_ex)(x)
    gS(x) = (∇⋅uS_ex)(x)
    fD(x) = μ/κ*uD_ex(x) + ∇(pD_ex)(x)
    gD(x) = (∇⋅uD_ex)(x)

    # Reference FEs
    order = 2 
    reffe_uS = ReferenceFE(lagrangian,VectorValue{2,Float64},order)
    reffe_pS = ReferenceFE(lagrangian,Float64,order-1)
    reffe_uD = ReferenceFE(raviart_thomas,Float64,order-2)
    reffe_pD = ReferenceFE(lagrangian,Float64,order-2)

    Ω   = Interior(model)
    Ω_S = Interior(model,tags="fluid")
    Ω_D = Interior(model,tags="porous")

    # Boundary triangulations and outer unit normals
    Γ_S = BoundaryTriangulation(model,tags="wallED")
    n_ΓS = get_normal_vector(Γ_S)
    Γ_D = BoundaryTriangulation(model,tags="wallPD")
    n_ΓD = get_normal_vector(Γ_D)

    # Triangulation and normal on interface (pointing outwards with respect to Ω_S)
    Σ = InterfaceTriangulation(Ω_S,Ω_D)
    n_Σ  = get_normal_vector(Σ)

    # Numerical integration
    degree = 2*order
    dΩ = Measure(Ω,degree)
    dΩ_S = Measure(Ω_S,degree)
    dΩ_D = Measure(Ω_D,degree)

    bdegree = 2*order
    dΓ_S = Measure(Γ_S,bdegree)
    dΓ_D = Measure(Γ_D,bdegree)

    idegree = 2*order
    dΣ = Measure(Σ, idegree)

    h_e_Σ = CellField(get_array(∫(1)dΣ)  , Σ  )

    # FE spaces
    Vs  = TestFESpace(Ω_S, reffe_uS, conformity=:H1, dirichlet_tags=["wallED"])
    Vd  = TestFESpace(Ω_D, reffe_uD, conformity=:HDiv, dirichlet_tags=["wallPD"])

    Qs = TestFESpace(Ω_S, reffe_pS , conformity=:C0)
    Qd = TestFESpace(Ω_D, reffe_pD, conformity=:L2)

    L_= ConstantFESpace(Ω)
    L = TrialFESpace(L_) 

    Us = TrialFESpace(Vs,[uS_ex])
    Ud = TrialFESpace(Vd,[uD_ex])
    Ps = TrialFESpace(Qs)
    Pd = TrialFESpace(Qd)

    Y = MultiFieldFESpace([Vs,Vd,Qs,Qd,L_])
    X = MultiFieldFESpace([Us,Ud,Ps,Pd,L])

    # Galerkin and Nitsche terms 
    a11(u,v) = ∫((2.0*μ)*(ε(u)⊙ε(v)))dΩ_S + 
                 ∫((μ*α/sqrt(κ))*(u.⁺×n_Σ.⁺)*(v.⁺×n_Σ.⁺))dΣ +
                 ∫((δ/h_e_Σ)*(u.⁺⋅n_Σ.⁺)*(v.⁺⋅n_Σ.⁺))dΣ 
                     
    a12(v,u) = ∫((δ/h_e_Σ)*(u.⁻⋅n_Σ.⁻)*(v.⁺⋅n_Σ.⁺))dΣ

    b13(v,q) = ∫(-(∇⋅v)*q)dΩ_S

    b14(v,q) = ∫(q.⁻*(v.⁺⋅n_Σ.⁺))dΣ
    
    b24(v,q) = ∫(-(∇⋅v)*q)dΩ_D + ∫(q.⁻*(v.⁻⋅n_Σ.⁻))dΣ
    
    a22(u,v) = ∫((μ/κ)*(u⋅v))dΩ_D + ∫((δ/h_e_Σ)*(u.⁻⋅n_Σ.⁻)*(v.⁻⋅n_Σ.⁻))dΣ

    b34(ψ,q) = ∫(ψ*q)dΩ_D           

    lhs((uS,uD,pS,pD,φ),(vS,vD,qS,qD,ψ)) = 
    a11(uS,vS) + a12(vS,uD) + b13(vS,pS) + b14(vS,pD)             +
    a12(uS,vD) + a22(uD,vD)              + b24(vD,pD) + b34(φ,qD) +
    b13(uS,qS)                                                    +
    b14(uS,qD) + b24(uD,qD)                                       +
                 b34(ψ,pD)

    lvS(v) = ∫(v⋅fS)dΩ_S + ∫(((σS_ex⋅n_Σ.⁺)⋅n_Σ.⁺)*(v.⁺⋅n_Σ.⁺))dΣ + 
             ∫((-(μ*α/sqrt(κ))*(uS_ex×n_Σ.⁺)-(σS_ex⋅n_Σ.⁺)×n_Σ.⁺)*(v.⁺×n_Σ.⁺))dΣ +
             ∫((δ/h_e_Σ)*(uS_ex⋅n_Σ.⁺ + uD_ex⋅n_Σ.⁻)*(v.⁺⋅n_Σ.⁺))dΣ

    lvD(v) = ∫(v⋅fD)dΩ_D + ∫((δ/h_e_Σ)*(uS_ex⋅n_Σ.⁺ + uD_ex⋅n_Σ.⁻)*(v.⁻⋅n_Σ.⁻))dΣ 

    lqS(q) = ∫(-q*gS)dΩ_S 

    lqD(q) = ∫(-q*gD)dΩ_D 
    
    lψ(ψ)  = ∫(ψ*pD_ex)dΩ_D
    
    rhs((vS,vD,qS,qD,ψ)) = lvS(vS) + lvD(vD) + lqS(qS) + lqD(qD) + lψ(ψ)

    # Build affine FE operator
    op = AffineFEOperator(lhs,rhs,X,Y) 
    uSh, uDh, pSh, pDh, ψh = solve(op)


    euS=uS_ex-uSh
    epS=pS_ex-pSh
    euD=uD_ex-uDh
    epD=pD_ex-pDh

    error_uS = sqrt(sum(∫(euS⋅euS + ε(euS)⊙ε(euS))*dΩ_S))
    error_pS = sqrt(sum(∫(epS*epS)*dΩ_S))
    error_uD = sqrt(sum(∫(euD⋅euD +(∇⋅(euD))*(∇⋅(euD)))*dΩ_D))
    error_pD = sqrt(sum(∫(epD*epD)*dΩ_D))

    error_uS, error_uD, error_pS, error_pD, Gridap.FESpaces.num_free_dofs(X) 
  end

  function convergence_test(nkmax)

    euS   = Float64[]
    euD   = Float64[]
    epS   = Float64[]
    epD   = Float64[]
    
    ruS   = Float64[]
    ruD   = Float64[]
    rpS   = Float64[]
    rpD = Float64[]
    hh   = Float64[]
    nn   = Int[]

    push!(ruS,0.)
    push!(ruD,0.)
    push!(rpS,0.)
    push!(rpD,0.)

    for nk in 1:nkmax
        println("******** Refinement step: $nk")
        model=generate_model_unit_square_biot_stokes(nk;
                                                         simplexify_model=true,
                                                         bcs_type=:full_dirichlet)
        push!(hh,sqrt(2)/2^nk)
        error_uS,error_uD,error_pS,error_pD,ndofs=solve_stokes_darcy_TH_RT(model)
        push!(nn,ndofs)
        println("******** Total DoFs: ", nn[nk])
        push!(euS,error_uS)
        push!(euD,error_uD)
        push!(epS,error_pS)
        push!(epD,error_pD)
        if nk>1
            push!(ruS, log(euS[nk]/euS[nk-1])/log(hh[nk]/hh[nk-1]))
            push!(ruD, log(euD[nk]/euD[nk-1])/log(hh[nk]/hh[nk-1]))
            push!(rpS, log(epS[nk]/epS[nk-1])/log(hh[nk]/hh[nk-1]))
            push!(rpD, log(epD[nk]/epD[nk-1])/log(hh[nk]/hh[nk-1]))
        end
    end
    println("======================================================================================================")
    println("   DoF  &    h   &   e(uS)   &  r(uS)  &  e(uD)  &  r(uD)  &  e(pS)  &  r(pS)  &  e(pD)  &  r(pD)     ")
    println("======================================================================================================")
    for nk in 1:nkmax
        @printf("%7d & %.4f & %.2e & %.3f & %.2e & %.3f & %.2e & %.3f & %.2e & %.3f \n", nn[nk], hh[nk], euS[nk], ruS[nk], euD[nk], ruD[nk], epS[nk], rpS[nk], epD[nk], rpD[nk]);
    end
    println("======================================================================================================")
  end
  convergence_test(4)
end
