module StokesDarcyNitscheAccuracyTestsTHRT

  using Gridap
  using GridapStokesDarcyNitsche
  import Gridap: ∇
  using Printf
  
  (Gridap.CellData.cross)(a::Gridap.CellData.CellField,
  b::Gridap.CellData.SkeletonPair{<:Gridap.CellData.CellField}) = Operation(cross)(a,b)
  
  Id  = TensorValue(1,0,0,1)
  const μ   = 10
  const α   = 1.0
  const κ   = 1.0
  const δ   = 40
  const μP  = 10
  const λP  = 10
  const K   = 1.0
  const φ   = 0.1
  const ρf  = 1.0
  const ρp  = 1.0
  const θ_sink = 0
  
  uf_ex(x)  = VectorValue(-cos(π*x[1])*sin(π*x[2]), sin(π*x[1])*cos(π*x[2]))
  pf_ex(x)  = cos(π*x[1])*cos(π*x[2])*cos(π*x[1])*cos(π*x[2])
  ur_ex(x)  = VectorValue(x[1]^3*x[2]+x[1]^2, -x[2]^2*x[1])
  ph_ex(x)  = sin(π*x[1])*sin(π*x[2])
  ys_ex(x) = VectorValue(-cos(π*x[1])*sin(π*x[2]), sin(π*x[1])*cos(π*x[2]))
  
  function solve_stokes_darcy_TH_RT(model; generate_output=false)

    labels = get_face_labeling(model)
    tags   = Gridap.Geometry.get_face_tag(labels,num_dims(model))
    
    # Method of manufactured solutions (RHS)
    σS_ex(x)  =  2*μ*(ε(uf_ex))(x) - pf_ex(x)*Id  
    σP_ex(x)  =  2*μ*φ*(ε(ur_ex))(x) + 2*μ*φ*(ε(ys_ex))(x) - φ*ph_ex(x)*Id 
    σEP_ex(x) =  2*μP*(ε(ys_ex))(x) + λP*(∇⋅ys_ex)(x) - (1-φ)*ph_ex(x)*Id
    
    FS(x)     =  ρf*uf_ex(x) -(∇⋅σS_ex)(x)
    QS(x)     =  (∇⋅uf_ex)(x)
    F1(x)     =  ρf*φ*ur_ex(x)+ρf*φ*ys_ex(x)-(∇⋅σP_ex)(x)+φ^2*κ^(-1)*ur_ex(x)-θ_sink*(ur_ex(x)+ys_ex(x))
    F2(x)     =  (1-φ)^2*K^(-1)*ph_ex(x)+(∇⋅ys_ex)(x)+φ*(∇⋅ur_ex)(x)
    F3(x)     =  ρf*φ*ur_ex(x)+ρp*ys_ex(x) -(∇⋅σP_ex)(x) -(∇⋅σEP_ex)(x)-θ_sink*(ur_ex(x)+ys_ex(x))
    
    # Reference FEs
    order = 2 
    reffe_uf = ReferenceFE(lagrangian,VectorValue{2,Float64},2)
    reffe_pf = ReferenceFE(lagrangian,Float64,1)
    reffe_ur = ReferenceFE(lagrangian,VectorValue{2,Float64},2)
    reffe_ph = ReferenceFE(lagrangian,Float64,1)
    reffe_ys = ReferenceFE(lagrangian,VectorValue{2,Float64},2)
 
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
    
    idegree = 2*order
    dΣ = Measure(Σ, idegree)
    h_e_Σ = CellField(get_array(∫(1)dΣ),Σ)
    
    # FE spaces
    Vf  = TestFESpace(Ω_S, reffe_uf, conformity=:H1, dirichlet_tags=["wallED"])
    Vr  = TestFESpace(Ω_D, reffe_ur, conformity=:H1, dirichlet_tags=["wallPD"])
    Vs  = TestFESpace(Ω_D, reffe_ys, conformity=:H1, dirichlet_tags = ["wallPD"])
    Qf = TestFESpace(Ω_S, reffe_pf , conformity=:C0)
    Qh = TestFESpace(Ω_D, reffe_ph, conformity=:C0)

    L_= ConstantFESpace(Ω)
    L = TrialFESpace(L_) 

    Uf = TrialFESpace(Vf,[uf_ex])
    Ur = TrialFESpace(Vr,[ur_ex])
    Us = TrialFESpace(Vs,[ys_ex])
    Pf = TrialFESpace(Qf)
    Ph = TrialFESpace(Qh)
    
    Y = MultiFieldFESpace([Vf,Vr,Vs,Qf,Qh,L_])
    X = MultiFieldFESpace([Uf,Ur,Us,Pf,Ph,L])
    

    a1(uf,vf)   = ∫((2*μ)*(ε(uf)⊙ε(vf)))dΩ_S + 
                  ∫((μ*α/sqrt(κ))*(uf.⁺×n_Σ.⁺)*(vf.⁺×n_Σ.⁺))dΣ + 
                  ∫((δ*μ/h_e_Σ)*(uf.⁺⋅n_Σ.⁺)*(vf.⁺⋅n_Σ.⁺))dΣ - 
                  ∫((2*μ)*((ε(uf).⁺⋅n_Σ.⁺)⋅n_Σ.⁺)*(vf.⁺⋅n_Σ.⁺))dΣ - 
                  ∫((2*μ)*((ε(vf).⁺⋅n_Σ.⁺)⋅n_Σ.⁺)*(uf.⁺⋅n_Σ.⁺))dΣ
    
    a2(ur,vf)   = ∫((δ*μ/h_e_Σ)*(ur.⁻⋅n_Σ.⁻)*(vf.⁺⋅n_Σ.⁺))dΣ -
                  ∫((2*μ)*((ε(vf).⁺⋅n_Σ.⁺)⋅n_Σ.⁺)*(ur.⁻⋅n_Σ.⁻))dΣ
                  
    a3(pf,vf)   = ∫(-(∇⋅vf)*pf)dΩ_S + ∫(pf.⁺*(vf.⁺⋅n_Σ.⁺))dΣ 
    a4(uf,qf)   = ∫((∇⋅uf)*qf)dΩ_S - ∫(qf.⁺*(uf.⁺⋅n_Σ.⁺))dΣ 
    a5(ys,ws)   = ∫((2*μP)*(ε(ys)⊙ε(ws)))dΩ_D + ∫(λP*(∇⋅ys)*(∇⋅ws))dΩ_D
    a6(ph,ws)   = ∫(-(∇⋅ws)*ph)dΩ_D 
    
    a7(uf,ws)   = ∫(-(μ*α/sqrt(κ))*(uf.⁺×n_Σ.⁺)*(ws.⁻×n_Σ.⁺))dΣ + 
                  ∫((δ*μ/h_e_Σ)*(uf.⁺⋅n_Σ.⁺)*(ws.⁻⋅n_Σ.⁻))dΣ -
                  ∫((2*μ)*((ε(uf).⁺⋅n_Σ.⁺)⋅n_Σ.⁺)*(ws.⁻⋅n_Σ.⁻))dΣ
                  
    a8(ur,ws)   = ∫((δ*μ/h_e_Σ)*(ur.⁻⋅n_Σ.⁻)*(ws.⁻⋅n_Σ.⁻))dΣ +
                  ∫((2*μ*φ)*(ε(ur)⊙ε(ws)))dΩ_D -
                  ∫(θ_sink*(ur⋅ws))dΩ_D 
    
    a9(ur,vr)   =  ∫((δ*μ/h_e_Σ)*(ur.⁻⋅n_Σ.⁻)*(vr.⁻⋅n_Σ.⁻))dΣ + 
                   ∫((μ*α/sqrt(κ))*(ur.⁻×n_Σ.⁻)*(vr.⁻×n_Σ.⁻))dΣ + 
                   ∫((2*μ*φ)*(ε(ur)⊙ ε(vr)))dΩ_D + ∫(φ^2*κ^(-1)*(ur⋅vr))dΩ_D -
                   ∫(θ_sink*(ur⋅vr))dΩ_D 
                   
    a10(ph,vr)  = ∫(-φ*(∇⋅vr)*ph)dΩ_D
    a11(uf,vr)  = ∫((δ*μ/h_e_Σ)*(uf.⁺⋅n_Σ.⁺)*(vr.⁻⋅n_Σ.⁻))dΣ - ∫((2*μ)*((ε(uf).⁺⋅n_Σ.⁺)⋅n_Σ.⁺)*(vr.⁻⋅n_Σ.⁻))dΣ
    a12(ur,qh)  = ∫(φ*(∇⋅ur)*qh)dΩ_D 
    a13(pf,vr)  = ∫(pf.⁺*(vr.⁻⋅n_Σ.⁻))dΣ 
    a14(pf,ws)  = ∫(pf.⁺*(ws.⁻⋅n_Σ.⁻))dΣ
    a15(ur,qf)  = ∫(-qf.⁺*(ur.⁻⋅n_Σ.⁻))dΣ 
    a16(ys,vr)  = ∫(-θ_sink*(ys⋅vr))dΩ_D 
    a17(ys,ws)  = ∫(-θ_sink*(ys⋅ws))dΩ_D 
    a18(ph,ψ)   = ∫(ψ*ph)dΩ_D
    a19(φ1,qh)   = ∫(φ1*qh)dΩ_D
  
    b1(ys,vf)   = ∫(-(μ*α/sqrt(κ))*(ys.⁻×n_Σ.⁺)*(vf.⁺×n_Σ.⁺))dΣ +
                  ∫((δ*μ/h_e_Σ)*(ys.⁻⋅n_Σ.⁻)*(vf.⁺⋅n_Σ.⁺))dΣ - 
                  ∫((2*μ)*((ε(vf).⁺⋅n_Σ.⁺)⋅n_Σ.⁺)*(ys.⁻⋅n_Σ.⁻))dΣ
                  
    b2(ys,ws)   = ∫((μ*α/sqrt(κ))*(ys.⁻×n_Σ.⁺)*(ws.⁻×n_Σ.⁺))dΣ + 
                  ∫((δ*μ/h_e_Σ)*(ys.⁻⋅n_Σ.⁻)*(ws.⁻⋅n_Σ.⁻))dΣ + 
                  ∫((2*μ*φ)*(ε(ys)⊙ ε(ws)))dΩ_D 
                  
    b3(ys,vr)   = ∫((δ*μ/h_e_Σ)*(ys.⁻⋅n_Σ.⁻)*(vr.⁻⋅n_Σ.⁻))dΣ +  ∫((2*μ*φ)*(ε(ys)⊙ε(vr)))dΩ_D 
    b4(ys,qh)   = ∫((∇⋅ys)*qh)dΩ_D
    b5(ph,qh)   = ∫(((1-φ)^2*K^(-1))*ph*qh)dΩ_D
    b6(ys,qf)   = ∫(-qf.⁺*(ys.⁻⋅n_Σ.⁻))dΣ
    b7(uf,vf)   = ∫(ρf*(uf⋅vf))dΩ_S
    b8(ur,vr)   = ∫((ρf*φ)*(ur⋅vr))dΩ_D
    b9(ys,vr)   = ∫((ρf*φ)*(ys⋅vr))dΩ_D
    b10(ur,ws)   = ∫((ρf*φ)*(ur⋅ws))dΩ_D
    b11(ys,ws)   = ∫((ρp)*(ys⋅ws))dΩ_D
    
    lvf(vf)     = ∫(vf⋅FS)dΩ_S + 
                  ∫((δ*μ/h_e_Σ)*(uf_ex⋅n_Σ.⁺+ur_ex⋅n_Σ.⁻+ ys_ex⋅n_Σ.⁻)*(vf.⁺⋅n_Σ.⁺))dΣ  - 
                  ∫((-(σS_ex⋅n_Σ.⁺)×n_Σ.⁺ -(μ*α/sqrt(κ))*(uf_ex×n_Σ.⁺)+(μ*α/sqrt(κ))*(ys_ex ×n_Σ.⁺))*(vf.⁺×n_Σ.⁺))dΣ  - 
                  ∫((2*μ)*((ε(vf).⁺⋅n_Σ.⁺)⋅n_Σ.⁺)*(uf_ex⋅n_Σ.⁺+ur_ex⋅n_Σ.⁻+ ys_ex⋅n_Σ.⁻))dΣ
                    
    lvr(vr)     = ∫(vr⋅F1)dΩ_D +
                  ∫((δ*μ/h_e_Σ)*(uf_ex⋅n_Σ.⁺+ur_ex⋅n_Σ.⁻+ ys_ex⋅n_Σ.⁻)*(vr.⁻⋅n_Σ.⁻))dΣ  + 
                  ∫((-(σS_ex⋅n_Σ.⁺)⋅n_Σ.⁺ + (σP_ex⋅n_Σ.⁻)⋅n_Σ.⁻)*(vr.⁻⋅n_Σ.⁻))dΣ -
                  ∫((-(σP_ex⋅n_Σ.⁻)×n_Σ.⁻ -(μ*α/sqrt(κ))*(ur_ex×n_Σ.⁻))*(vr.⁻×n_Σ.⁻))dΣ 
                    
    lys(ws)     = ∫(F3⋅ws)dΩ_D + 
                  ∫((δ*μ/h_e_Σ)*(uf_ex⋅n_Σ.⁺+ur_ex⋅n_Σ.⁻+ ys_ex⋅n_Σ.⁻)*(ws.⁻⋅n_Σ.⁻))dΣ + 
                  ∫((-(σS_ex⋅n_Σ.⁺)×n_Σ.⁺ -(μ*α/sqrt(κ))*(uf_ex×n_Σ.⁺)+(μ*α/sqrt(κ))*(ys_ex ×n_Σ.⁺))*(ws.⁻×n_Σ.⁺))dΣ + 
                  ∫((σS_ex⋅n_Σ.⁺ + σP_ex⋅n_Σ.⁻ + σEP_ex⋅n_Σ.⁻)⋅ws.⁻)dΣ   
    
    lqf(qf)     = ∫(qf*QS)dΩ_S - ∫(qf.⁺*(uf_ex⋅n_Σ.⁺+ur_ex⋅n_Σ.⁻+ ys_ex⋅n_Σ.⁻))dΣ
    lqh(qh)     = ∫(qh*F2)dΩ_D 
    lψ(ψ)       = ∫(ψ*ph_ex)dΩ_D
    
    lhs((uf,ur,ys,pf,ph,φ1),(vf,vr,ws,qf,qh,ψ)) =
               a1(uf,vf) + a2(ur,vf) + a3(pf,vf) + a4(uf,qf) + a5(ys,ws) +
               a6(ph,ws) + a7(uf,ws) + a8(ur,ws) + a9(ur,vr) + a10(ph,vr) +
               a11(uf,vr) + a12(ur,qh) + a13(pf,vr) + a14(pf,ws) + a15(ur,qf) +
               a16(ys,vr) + a17(ys,ws) + a18(ph,ψ) + a19(φ1,qh) +
               b1(ys,vf) + b2(ys,ws) + b3(ys,vr) + b4(ys,qh) + b5(ph,qh) +
               b6(ys,qf) + b7(uf,vf) + b8(ur,vr) + b9(ys,vr) + b10(ur,ws) + b11(ys,ws)

    
    rhs((vf,vr,ws,qf,qh,ψ)) = lvf(vf) + lvr(vr) + lys(ws) + lqf(qf) +lqh(qh) + lψ(ψ) 
    
    # Build affine FE operator
    op = AffineFEOperator(lhs,rhs,X,Y) 
    ufh, urh, ysh, pfh, ph, ψh = solve(op) 

    euf  = uf_ex-ufh
    epf  = pf_ex-pfh
    eys  = ys_ex-ysh
    eur  = ur_ex-urh
    eph  = ph_ex-ph

    error_uf = sqrt(sum(∫(euf⋅euf + ε(euf)⊙ε(euf))*dΩ_S))
    error_pf = sqrt(sum(∫(epf*epf)*dΩ_S))
    error_ys = sqrt(sum(∫(eys⋅eys + ε(eys)⊙ε(eys))*dΩ_D))
    error_ur = sqrt(sum(∫(eur⋅eur)*dΩ_D))
    error_ph = sqrt(sum(∫(eph*eph)*dΩ_D))

    error_uf, error_ur, error_ys, error_pf, error_ph,  Gridap.FESpaces.num_free_dofs(X)  
  end

  function convergence_test(nkmax)

    euf   = Float64[]
    eur   = Float64[]
    eys   = Float64[]
    epf   = Float64[]
    eph   = Float64[]
    
    ruf   = Float64[]
    rur   = Float64[]
    rys   = Float64[]
    rpf   = Float64[]
    rph   = Float64[]
    hh    = Float64[]
    nn    = Int[]

    push!(ruf,0.)
    push!(rur,0.)
    push!(rys,0.)
    push!(rpf,0.)
    push!(rph,0.)
    
    for nk in 1:nkmax
        println("******** Refinement step: $nk")
        model=generate_model_unit_square_biot_stokes(nk;simplexify_model=true,bcs_type=:full_dirichlet)
        push!(hh,1/2^nk)

        error_uf,error_ur,error_ys,error_pf,error_ph,ndofs = solve_stokes_darcy_TH_RT(model; generate_output=true)
        push!(nn,ndofs)
        println("******** Total DoFs: ", nn[nk])
        push!(euf,error_uf)
        push!(eur,error_ur)
        push!(eys,error_ys)
        push!(epf,error_pf)
        push!(eph,error_ph)
        
        if nk>1
            push!(ruf, log(euf[nk]/euf[nk-1])/log(hh[nk]/hh[nk-1]))
            push!(rur, log(eur[nk]/eur[nk-1])/log(hh[nk]/hh[nk-1]))
            push!(rys, log(eys[nk]/eys[nk-1])/log(hh[nk]/hh[nk-1]))
            push!(rpf, log(epf[nk]/epf[nk-1])/log(hh[nk]/hh[nk-1]))
            push!(rph, log(eph[nk]/eph[nk-1])/log(hh[nk]/hh[nk-1]))
        end
    end
    println("====================================================================================================================================")
    println("   DoF  &    h   &   e(uf)  & r(uf) &   e(ur)  & r(ur) &   e(pf)  & r(pf) &   e(ph)  & r(ph) &   e(ys)  & r(ys)  ")
    println("====================================================================================================================================")
    for nk in 1:nkmax
        @printf("%7d & %.4f & %.2e & %.3f & %.2e & %.3f & %.2e & %.3f & %.2e & %.3f & %.2e & %.3f   \n", nn[nk], hh[nk], euf[nk], ruf[nk], eur[nk], rur[nk], epf[nk], rpf[nk], eph[nk], rph[nk],eys[nk], rys[nk]);
    end
    println("=====================================================================================================================================")
  end
  convergence_test(4)
end
