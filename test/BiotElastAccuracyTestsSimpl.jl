module BiotElastAccuracyTestsSimpl
  using Gridap
  using GridapStokesDarcyNitsche
  import Gridap: ∇
  using Printf

  const λ_P = 1.
  const λ_E = 1.
  const μ_P = 1.
  const μ_E = 1.
  const c_0 = 1.
  const α = 1.
  const κ = 1.
  const η = 1.

  Id = TensorValue(1,0,0,1)
  
  # manufactured solutions
  u_ex(x)  = VectorValue(sin(π*(x[1]+x[2])), cos(π*(x[1]^2+x[2]^2)))
  p_ex(x)  = sin(π*x[1]+x[2])*sin(π*x[2])

  #model=generate_model_unit_square_biot_stokes(4;
  #                                                       simplexify_model=true,
  #                                                       bcs_type=:full_dirichlet)

  function solve_biot_elasticity_taylor_hood(model)
    
    wallEDtag="wallED"
    wallPDtag="wallPD"
 
    labels = get_face_labeling(model)
    tags   = Gridap.Geometry.get_face_tag(labels,num_dims(model))

    φE_ex(x) = -λ_E*(∇⋅u_ex)(x)
    φP_ex(x) = α*p_ex(x)-λ_P*(∇⋅u_ex)(x)

    # Method of manufactured solutions (RHS)

    σE_ex(x) = 2*μ_E*ε(u_ex)(x)- φE_ex(x)*Id
    σP_ex(x) = 2*μ_P*ε(u_ex)(x)- φP_ex(x)*Id

    σpP(x) = κ/η*∇(p_ex)(x)

    buP(x) = -(∇⋅σP_ex)(x)
    buE(x) = -(∇⋅σE_ex)(x)
    bpP(x) = (c_0 + α^2/λ_P)*p_ex(x) - α/λ_P*φP_ex(x) - (∇⋅σpP)(x)

    # Reference FEs
    order = 2
    reffe_u   = ReferenceFE(lagrangian,VectorValue{2,Float64},order)
    reffe_p   = ReferenceFE(lagrangian,Float64,order-1)
    reffe_phi = ReferenceFE(lagrangian,Float64,order-2)

    Ω   = Interior(model)
    Ω_E = Interior(model,tags="fluid")
    Ω_P = Interior(model,tags="porous")

    # Boundary triangulations and outer unit normals

    Γ_P = BoundaryTriangulation(model,tags=wallPDtag)
    n_ΓP = get_normal_vector(Γ_P)
    Γ_E = BoundaryTriangulation(model,tags=wallEDtag)
    n_ΓE = get_normal_vector(Γ_E)

    # Triangulation and normal on interface
    # (pointing outwards with respect to Ω_P)
    Σ = InterfaceTriangulation(Ω_P,Ω_E)
    n_Σ  = get_normal_vector(Σ)

    # Numerical integration
    degree = 2*order
    dΩ = Measure(Ω,degree)
    dΩ_P = Measure(Ω_P,degree)
    dΩ_E = Measure(Ω_E,degree)

    bdegree = 2*order
    dΓ_E = Measure(Γ_E,bdegree)
    dΓ_P = Measure(Γ_P,bdegree)

    idegree = 2*order
    dΣ = Measure(Σ, idegree)

    # FE spaces
    V  = TestFESpace(Ω, reffe_u  , conformity=:H1, dirichlet_tags=[wallEDtag,wallPDtag])
    Q  = TestFESpace(Ω_P, reffe_p  , conformity=:H1, dirichlet_tags=[wallPDtag])
    Z = TestFESpace(Ω, reffe_phi, conformity=:L2)
   
    U = TrialFESpace(V,[u_ex, u_ex])
    W = TrialFESpace(Z)
    P = TrialFESpace(Q, p_ex)

    Y = MultiFieldFESpace([V,Q,Z])
    X = MultiFieldFESpace([U,P,W])

    # Galerkin terms (no stabilisation)
    avPuP(v,u) = ∫(2.0*ε(v)⊙(μ_P*ε(u)))dΩ_P
    avPφP(v,φ) = ∫(-(∇⋅v)*φ)dΩ_P

    aqPpP(q,p) = ∫( (c_0+α^2/λ_P)*(p*q)+ κ/η*(∇(p)⋅∇(q)))dΩ_P
    aqPφP(q,φ) = ∫(-α/λ_P*(φ*q))dΩ_P

    aψPφP(ψ,φ) = ∫(1.0/λ_P*φ*ψ)dΩ_P
    aψPpP(ψ,p) = ∫(-α/λ_P*(ψ*p))dΩ_P
    aψPuP(ψ,u) = ∫(ψ*(∇⋅u))dΩ_P

    avEuE(v,u) = ∫(2.0*ε(v)⊙(μ_E*ε(u)))dΩ_E
    avEφE(v,φ)  = ∫(-(∇⋅v)*φ)dΩ_E

    aψEφE(ψ,φ) = ∫(1.0/λ_E*φ*ψ)dΩ_E
    aψEuE(ψ,u) = ∫(ψ*(∇⋅u))dΩ_E

    # Method of manufactured solutions (RHS)
    lvP(v) = ∫(v⋅buP)dΩ_P + ∫(v.⁺⋅((n_Σ.⁺)⋅σP_ex))dΣ
    lvE(v) = ∫(v⋅buE)dΩ_E + ∫(v.⁻⋅((n_Σ.⁻)⋅σE_ex))dΣ
    lqP(q) = ∫(q*bpP)dΩ_P + ∫(q.⁺*((n_Σ.⁺)⋅σpP))dΣ #+ ∫(q*(n_ΓPN⋅σpP))dΓ_PN
   
    # Forms
    lhs((u,p,φ),(v,q,ψ)) =  avPuP(v,u)+avPφP(v,φ)-
                            aqPpP(q,p)-aqPφP(q,φ)-
                            aψPφP(ψ,φ)-aψPpP(ψ,p)-aψPuP(ψ,u)+
                            avEuE(v,u)+avEφE(v,φ)-
                            aψEφE(ψ,φ)-aψEuE(ψ,u)

    rhs((v,q,ψ)) = lvP(v)+lvE(v)-lqP(q) #+
                   #∫(v⋅(n_ΓEN⋅σE_ex))dΓ_EN +
                   #∫(v⋅(n_ΓPN⋅σP_ex))dΓ_PN

    # Build affine FE operator
    op = AffineFEOperator(lhs,rhs,X,Y) 
    
    uh, ph, phih = solve(op)

    eu=u_ex-uh
    ep=p_ex-ph
    eφP=φP_ex-phih
    eφE=φE_ex-phih

    error_u = sqrt(sum(∫(2*ε(eu)⊙ε(eu))*dΩ))
    error_p = sqrt(sum(∫((ep)*(ep) + κ/η* ∇(ep)⋅∇(ep))*dΩ))
    error_φ = sqrt(sum(∫((eφP)*(eφP))*dΩ_P+
                       ∫((eφE)*(eφE))*dΩ_E))

    
    error_u, error_p, error_φ, Gridap.FESpaces.num_free_dofs(X) 
  end

  function convergence_test(nkmax)


    eu   = Float64[]
    ep   = Float64[]
    ephi = Float64[]
    ru   = Float64[]
    rp   = Float64[]
    rphi = Float64[]
    hh   = Float64[]
    nn   = Int[]

    push!(ru,0.)
    push!(rp,0.)
    push!(rphi,0.)

    
    for nk in 1:nkmax
        println("******** Refinement step: $nk")
        model=generate_model_unit_square_biot_stokes(nk;
                                                         simplexify_model=true,
                                                         bcs_type=:full_dirichlet)
        push!(hh,sqrt(2)/2^nk)
        error_u,error_p,error_φ,ndofs=solve_biot_elasticity_taylor_hood(model)
        push!(nn,ndofs)
        println("******** Total DoFs: ", nn[nk])
        push!(eu,error_u)
        push!(ep,error_p)
        push!(ephi,error_φ)
        if nk>1
            push!(ru, log(eu[nk]/eu[nk-1])/log(hh[nk]/hh[nk-1]))
            push!(rp, log(ep[nk]/ep[nk-1])/log(hh[nk]/hh[nk-1]))
            push!(rphi, log(ephi[nk]/ephi[nk-1])/log(hh[nk]/hh[nk-1]))
        end
    end
    println("================================================================================")
    println("   DoF  &    h   &   e(u)   &  r(u)  &  e(p)  &  r(p)  &  e(phi)  &  r(phi)     ")
    println("================================================================================")
    for nk in 1:nkmax
        @printf("%7d & %.4f & %.2e & %.3f & %.2e & %.3f & %.2e & %.3f \n", nn[nk], hh[nk], eu[nk], ru[nk], ep[nk], rp[nk], ephi[nk], rphi[nk]);
    end
    println("================================================================================")
  end
  convergence_test(4)
end
