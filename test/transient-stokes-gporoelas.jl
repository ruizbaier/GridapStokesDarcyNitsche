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
  
  # manufactured solutions 
  uf_ex(t)  = x -> VectorValue(t*sin(pi*x[1])*sin(pi*x[1])*sin(pi*x[2])*cos(pi*x[2]),-t*sin(pi*x[2])*sin(pi*x[2])*sin(pi*x[1])*cos(pi*x[1]))
  pf_ex(t)  = x -> t*t*sin(pi*x[1])*sin(pi*x[2])
  ur_ex(t)  = x -> VectorValue(t*t*sin(π*x[2])*sin(π*x[2])-t*x[1]*x[1]*x[1]*cos(π*x[2]),t*t*sin(π*x[2])*sin(π*x[2])+2*t*x[1]*x[1]*x[1]*sin(π*x[2]))
  ph_ex(t)  = x -> t*t*(1-sin(π*x[1])*sin(π*x[2]))
  ys_ex(t)  = x -> VectorValue(0.5*t*t*x[1]*x[1]*x[1]*cos(π*x[2]),-t*t*x[1]*x[1]*x[1]*sin(π*x[2]))
  us_ex(t)  = x -> VectorValue(t*x[1]*x[1]*x[1]*cos(π*x[2]),-2*t*x[1]*x[1]*x[1]*sin(π*x[2]))
  
  duf_ex(t)  = x -> VectorValue(sin(pi*x[1])*sin(pi*x[1])*sin(pi*x[2])*cos(pi*x[2]), -sin(pi*x[2])*sin(pi*x[2])*sin(pi*x[1])*cos(pi*x[1]))
  dur_ex(t)  = x -> VectorValue(2*t*sin(π*x[2])*sin(π*x[2])-x[1]*x[1]*x[1]*cos(π*x[2]),2*t*sin(π*x[2])*sin(π*x[2])+2*x[1]*x[1]*x[1]*sin(π*x[2]))
  dys_ex(t)  = x -> VectorValue(t*x[1]*x[1]*x[1]*cos(π*x[2]),-2*t*x[1]*x[1]*x[1]*sin(π*x[2]))
  dph_ex(t)  = x -> 2*t*(1-sin(π*x[1])*sin(π*x[2]))
  dus_ex(t)  = x -> VectorValue(x[1]*x[1]*x[1]*cos(π*x[2]),-2*x[1]*x[1]*x[1]*sin(π*x[2]))
 

  function solve_stokes_darcy_TH_RT(model,Δt; generate_output=false)

    labels = get_face_labeling(model)
    tags   = Gridap.Geometry.get_face_tag(labels,num_dims(model))
    
    # Method of manufactured solutions (RHS)
    σS_ex(t)  = x -> 2*μ*(ε(uf_ex(t)))(x) - pf_ex(t)(x)*Id  
    σP_ex(t)  = x -> 2*μ*φ*(ε(ur_ex(t)))(x) + 2*μ*φ*(ε(dys_ex(t)))(x) - φ*ph_ex(t)(x)*Id 
    σEP_ex(t) = x -> 2*μP*(ε(ys_ex(t)))(x) + λP*(∇⋅ys_ex(t))(x) - (1-φ)*ph_ex(t)(x)*Id
    
    FS(t)     = x -> ρf*duf_ex(t)(x) -(∇⋅σS_ex(t))(x)
    QS(t)     = x -> (∇⋅uf_ex(t))(x)
    F1(t)     = x -> ρf*φ*dur_ex(t)(x)+ρf*φ*dus_ex(t)(x)-(∇⋅σP_ex(t))(x)+φ^2*κ^(-1)*ur_ex(t)(x)-θ_sink*(ur_ex(t)(x)+us_ex(t)(x))
    F2(t)     = x -> (1-φ)^2*K^(-1)*dph_ex(t)(x)+(∇⋅dys_ex(t))(x)+φ*(∇⋅ur_ex(t))(x)
    F3(t)     = x -> ρf*φ*dur_ex(t)(x)+ρp*dus_ex(t)(x) -(∇⋅σP_ex(t))(x) -(∇⋅σEP_ex(t))(x)-θ_sink*(ur_ex(t)(x)+us_ex(t)(x))
    F4(t)     = x -> ρp*us_ex(t)(x)-ρp*dys_ex(t)(x)
    
    # Reference FEs
    order = 2 
    reffe_uf = ReferenceFE(lagrangian,VectorValue{2,Float64},2)
    reffe_pf = ReferenceFE(lagrangian,Float64,1)
    reffe_ur = ReferenceFE(lagrangian,VectorValue{2,Float64},2)
    reffe_ph = ReferenceFE(lagrangian,Float64,1)
    reffe_ys = ReferenceFE(lagrangian,VectorValue{2,Float64},2)
    reffe_us = ReferenceFE(lagrangian,VectorValue{2,Float64},1)
 
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
    Ws = TestFESpace(Ω_D, reffe_us, conformity=:C0)
    Qf = TestFESpace(Ω_S, reffe_pf , conformity=:C0)
    Qh = TestFESpace(Ω_D, reffe_ph, conformity=:C0)

    L_= ConstantFESpace(Ω)
    L = TrialFESpace(L_) 

    Uf = TransientTrialFESpace(Vf,[uf_ex])
    Ur = TransientTrialFESpace(Vr,[ur_ex])
    Us = TransientTrialFESpace(Vs,[ys_ex])
    Ys = TransientTrialFESpace(Ws)
    Pf = TransientTrialFESpace(Qf)
    Ph = TransientTrialFESpace(Qh)
    
    Y = MultiFieldFESpace([Vf,Vr,Vs,Ws,Qf,Qh,L_])
    X = MultiFieldFESpace([Uf,Ur,Us,Ys,Pf,Ph,L])
    
    # terms without time derivative
    a1(t,uf,vf)   = ∫((2*μ)*(ε(uf)⊙ε(vf)))dΩ_S + 
                    ∫((μ*α/sqrt(κ))*(uf.⁺×n_Σ.⁺)*(vf.⁺×n_Σ.⁺))dΣ + 
                    ∫((δ*μ/h_e_Σ)*(uf.⁺⋅n_Σ.⁺)*(vf.⁺⋅n_Σ.⁺))dΣ - 
                    ∫((2*μ)*((ε(uf).⁺⋅n_Σ.⁺)⋅n_Σ.⁺)*(vf.⁺⋅n_Σ.⁺))dΣ - 
                    ∫((2*μ)*((ε(vf).⁺⋅n_Σ.⁺)⋅n_Σ.⁺)*(uf.⁺⋅n_Σ.⁺))dΣ
                    
    a2(t,ur,vf)   = ∫((δ*μ/h_e_Σ)*(ur.⁻⋅n_Σ.⁻)*(vf.⁺⋅n_Σ.⁺))dΣ - 
                    ∫((2*μ)*((ε(vf).⁺⋅n_Σ.⁺)⋅n_Σ.⁺)*(ur.⁻⋅n_Σ.⁻))dΣ
                    
    a3(t,pf,vf)   = ∫(-(∇⋅vf)*pf)dΩ_S + ∫(pf.⁺*(vf.⁺⋅n_Σ.⁺))dΣ 
    a4(t,uf,qf)   = ∫((∇⋅uf)*qf)dΩ_S - ∫(qf.⁺*(uf.⁺⋅n_Σ.⁺))dΣ 
    a5(t,ys,ws)   = ∫((2*μP)*(ε(ys)⊙ε(ws)))dΩ_D + ∫(λP*(∇⋅ys)*(∇⋅ws))dΩ_D
    a6(t,ph,ws)   = ∫(-(∇⋅ws)*ph)dΩ_D 
    
    a7(t,uf,ws)   = ∫(-(μ*α/sqrt(κ))*(uf.⁺×n_Σ.⁺)*(ws.⁻×n_Σ.⁺))dΣ +
                    ∫((δ*μ/h_e_Σ)*(uf.⁺⋅n_Σ.⁺)*(ws.⁻⋅n_Σ.⁻))dΣ -
                    ∫((2*μ)*((ε(uf).⁺⋅n_Σ.⁺)⋅n_Σ.⁺)*(ws.⁻⋅n_Σ.⁻))dΣ
                    
    a8(t,ur,ws)   = ∫((δ*μ/h_e_Σ)*(ur.⁻⋅n_Σ.⁻)*(ws.⁻⋅n_Σ.⁻))dΣ +
                    ∫((2*μ*φ)*(ε(ur)⊙ε(ws)))dΩ_D -∫(θ_sink*(ur⋅ws))dΩ_D 
                    
    a9(t,ur,vr)   =  ∫((δ*μ/h_e_Σ)*(ur.⁻⋅n_Σ.⁻)*(vr.⁻⋅n_Σ.⁻))dΣ + 
                     ∫((μ*α/sqrt(κ))*(ur.⁻×n_Σ.⁻)*(vr.⁻×n_Σ.⁻))dΣ + 
                     ∫((2*μ*φ)*(ε(ur)⊙ ε(vr)))dΩ_D + 
                     ∫(φ^2*κ^(-1)*(ur⋅vr))dΩ_D -∫(θ_sink*(ur⋅vr))dΩ_D 
                     
    a10(t,ph,vr)  = ∫(-φ*(∇⋅vr)*ph)dΩ_D
    a11(t,uf,vr)  = ∫((δ*μ/h_e_Σ)*(uf.⁺⋅n_Σ.⁺)*(vr.⁻⋅n_Σ.⁻))dΣ - 
                    ∫((2*μ)*((ε(uf).⁺⋅n_Σ.⁺)⋅n_Σ.⁺)*(vr.⁻⋅n_Σ.⁻))dΣ
                    
    a12(t,ur,qh)  = ∫(φ*(∇⋅ur)*qh)dΩ_D 
    a13(t,pf,vr)  = ∫(pf.⁺*(vr.⁻⋅n_Σ.⁻))dΣ 
    a14(t,pf,ws)  = ∫(pf.⁺*(ws.⁻⋅n_Σ.⁻))dΣ
    a15(t,ur,qf)  = ∫(-qf.⁺*(ur.⁻⋅n_Σ.⁻))dΣ 
    a16(t,us,vr)  = ∫(-θ_sink*(us⋅vr))dΩ_D 
    a17(t,us,ws)  = ∫(-θ_sink*(us⋅ws))dΩ_D 
    a18(t,us,vs)  = ∫(ρp*(us⋅vs))dΩ_D  
    a19(t,ph,ψ)   = ∫(ψ*ph)dΩ_D
    a20(t,φ1,qh)   = ∫(φ1*qh)dΩ_D
  
    b1(t,dys,vf)   = ∫(-(μ*α/sqrt(κ))*(dys.⁻×n_Σ.⁺)*(vf.⁺×n_Σ.⁺))dΣ + 
                     ∫((δ*μ/h_e_Σ)*(dys.⁻⋅n_Σ.⁻)*(vf.⁺⋅n_Σ.⁺))dΣ -
                     ∫((2*μ)*((ε(vf).⁺⋅n_Σ.⁺)⋅n_Σ.⁺)*(dys.⁻⋅n_Σ.⁻))dΣ
                     
    b2(t,dys,ws)   = ∫((μ*α/sqrt(κ))*(dys.⁻×n_Σ.⁺)*(ws.⁻×n_Σ.⁺))dΣ + 
                     ∫((δ*μ/h_e_Σ)*(dys.⁻⋅n_Σ.⁻)*(ws.⁻⋅n_Σ.⁻))dΣ + 
                     ∫((2*μ*φ)*(ε(dys)⊙ ε(ws)))dΩ_D 
                     
    b3(t,dys,vr)   = ∫((δ*μ/h_e_Σ)*(dys.⁻⋅n_Σ.⁻)*(vr.⁻⋅n_Σ.⁻))dΣ +  ∫((2*μ*φ)*(ε(dys)⊙ε(vr)))dΩ_D 
    b4(t,dys,qh)   = ∫((∇⋅dys)*qh)dΩ_D
    b5(t,dph,qh)   = ∫(((1-φ)^2*K^(-1))*dph*qh)dΩ_D
    b6(t,dys,qf)   = ∫(-qf.⁺*(dys.⁻⋅n_Σ.⁻))dΣ
    b7(t,duf,vf)   = ∫(ρf*(duf⋅vf))dΩ_S
    b8(t,dur,vr)   = ∫((ρf*φ)*(dur⋅vr))dΩ_D
    b9(t,dus,vr)   = ∫((ρf*φ)*(dus⋅vr))dΩ_D
    b10(t,dur,ws)   = ∫((ρf*φ)*(dur⋅ws))dΩ_D
    b11(t,dus,ws)   = ∫((ρp)*(dus⋅ws))dΩ_D
    b12(t,dys,vs)   = ∫(-(ρp)*(dys⋅vs))dΩ_D
    
    lvf(t,vf)     = ∫(vf⋅FS(t))dΩ_S +  
                    ∫((δ*μ/h_e_Σ)*(uf_ex(t)⋅n_Σ.⁺+ur_ex(t)⋅n_Σ.⁻+ dys_ex(t)⋅n_Σ.⁻)*(vf.⁺⋅n_Σ.⁺))dΣ  -
                    ∫((-(σS_ex(t)⋅n_Σ.⁺)×n_Σ.⁺ -(μ*α/sqrt(κ))*(uf_ex(t)×n_Σ.⁺)+(μ*α/sqrt(κ))*(dys_ex(t) ×n_Σ.⁺))*(vf.⁺×n_Σ.⁺))dΣ -
                    ∫((2*μ)*((ε(vf).⁺⋅n_Σ.⁺)⋅n_Σ.⁺)*(uf_ex(t)⋅n_Σ.⁺+ur_ex(t)⋅n_Σ.⁻+ dys_ex(t)⋅n_Σ.⁻))dΣ
                    
    lvr(t,vr)     = ∫(vr⋅F1(t))dΩ_D +
                    ∫((δ*μ/h_e_Σ)*(uf_ex(t)⋅n_Σ.⁺+ur_ex(t)⋅n_Σ.⁻+ dys_ex(t)⋅n_Σ.⁻)*(vr.⁻⋅n_Σ.⁻))dΣ  + 
                    ∫((-(σS_ex(t)⋅n_Σ.⁺)⋅n_Σ.⁺ + (σP_ex(t)⋅n_Σ.⁻)⋅n_Σ.⁻)*(vr.⁻⋅n_Σ.⁻))dΣ -
                    ∫((-(σP_ex(t)⋅n_Σ.⁻)×n_Σ.⁻ -(μ*α/sqrt(κ))*(ur_ex(t)×n_Σ.⁻))*(vr.⁻×n_Σ.⁻))dΣ 
                    
    lys(t,ws)     = ∫(F3(t)⋅ws)dΩ_D +
                    ∫((δ*μ/h_e_Σ)*(uf_ex(t)⋅n_Σ.⁺+ur_ex(t)⋅n_Σ.⁻+ dys_ex(t)⋅n_Σ.⁻)*(ws.⁻⋅n_Σ.⁻))dΣ +
                    ∫((-(σS_ex(t)⋅n_Σ.⁺)×n_Σ.⁺ -(μ*α/sqrt(κ))*(uf_ex(t)×n_Σ.⁺)+(μ*α/sqrt(κ))*(dys_ex(t) ×n_Σ.⁺))*(ws.⁻×n_Σ.⁺))dΣ +
                    ∫((σS_ex(t)⋅n_Σ.⁺ + σP_ex(t)⋅n_Σ.⁻ + σEP_ex(t)⋅n_Σ.⁻)⋅ws.⁻)dΣ      
                    
    lqf(t,qf)     = ∫(qf*QS(t))dΩ_S - ∫(qf.⁺*(uf_ex(t)⋅n_Σ.⁺+ur_ex(t)⋅n_Σ.⁻+ dys_ex(t)⋅n_Σ.⁻))dΣ
    lqh(t,qh)     = ∫(qh*F2(t))dΩ_D 
    lvs(t,vs)     = ∫(vs⋅F4(t))dΩ_D
    lψ(t,ψ)       = ∫(ψ*ph_ex(t))dΩ_D
    
    lhs(t,(uf,ur,ys,us,pf,ph,φ1),(vf,vr,ws,vs,qf,qh,ψ)) = 
                                     a1(t,uf,vf)+a2(t,ur,vf)+a3(t,pf,vf)+a4(t,uf,qf)+a5(t,ys,ws)+
                                     a6(t,ph,ws)+a7(t,uf,ws)+a8(t,ur,ws)+a9(t,ur,vr)+a10(t,ph,vr)+
                                     a11(t,uf,vr)+a12(t,ur,qh)+a13(t,pf,vr)+a14(t,pf,ws)+a15(t,ur,qf)+
                                     a16(t,us,vr)+a17(t,us,ws)+a18(t,us,vs)+a19(t,ph,ψ)+a20(t,φ1,qh)
    
    mass(t,(duf,dur,dys,dus,pf,dph,φ1),(vf,vr,ws,vs,qf,qh,ψ)) = 
                                     b1(t,dys,vf)+b2(t,dys,ws)+b3(t,dys,vr)+b4(t,dys,qh)+b5(t,dph,qh)+
                                     b6(t,dys,qf)+b7(t,duf,vf)+b8(t,dur,vr)+b9(t,dus,vr)+b10(t,dur,ws)+
                                     b11(t,dus,ws)+b12(t,dys,vs)  
    
    rhs(t,(vf,vr,ws,vs,qf,qh,ψ)) = lvf(t,vf) + lvr(t,vr) + lys(t,ws) + lvs(t,vs) + lqf(t,qf) +lqh(t,qh) + lψ(t,ψ) 
    
    op = TransientLinearFEOperator((lhs,mass),rhs,X,Y)
    
    ls = LUSolver()
    θ = 1
    solver = ThetaMethod(ls, Δt, θ)
    t0, tF = 0.0, 0.01
    
    xh0 = interpolate_everywhere([uf_ex(t0),ur_ex(t0),ys_ex(t0),us_ex(t0),pf_ex(t0),ph_ex(t0),1.0],X(t0))
    fesltn = solve(solver, op, t0, tF, xh0)
    
    el2uf_accum  = 0.0
    el2ur_accum  = 0.0
    el2ys_accum  = 0.0
    el2us_accum  = 0.0
    el2pf_accum  = 0.0
    el2ph_accum  = 0.0
    
    last_xhs_n = nothing  # to hold the last solution

    for (t_n, xhs_n) in fesltn
        last_xhs_n = xhs_n  # save last time step solution
        
        eh_u1 = uf_ex(t_n)  - xhs_n[1]
        eh_u2 = ur_ex(t_n)  - xhs_n[2]
        eh_u3 = ys_ex(t_n)  - xhs_n[3]
        eh_u4 = us_ex(t_n)  - xhs_n[4]
        eh_u5 = pf_ex(t_n)  - xhs_n[5]
        eh_u6 = ph_ex(t_n)  - xhs_n[6] 
        
        spatial_err2_uf = sum(∫(eh_u1⋅eh_u1 + ε(eh_u1)⊙ε(eh_u1))*dΩ_S)
        spatial_err2_ur = sum(∫(eh_u2⋅eh_u2  + ε(eh_u2)⊙ε(eh_u2))*dΩ_D)                         
        spatial_err2_ys = sum(∫(eh_u3⋅eh_u3 + ε(eh_u3)⊙ε(eh_u3))*dΩ_D)
        spatial_err2_us  = sum(∫(eh_u4⋅eh_u4)*dΩ_D)
        spatial_err2_pf  = sum(∫(eh_u5*eh_u5)*dΩ_S)
        spatial_err2_ph  = sum(∫(eh_u6*eh_u6)*dΩ_D)
        
        el2uf_accum  += spatial_err2_uf  * Δt
        el2ur_accum  += spatial_err2_ur  * Δt
        el2ys_accum  += spatial_err2_ys  * Δt
        el2us_accum  += spatial_err2_us  * Δt
        el2pf_accum  += spatial_err2_pf  * Δt
        el2ph_accum  += spatial_err2_ph  * Δt
    end
    
    error_uf  = sqrt(el2uf_accum)
    error_ur  = sqrt(el2ur_accum)
    error_ys  = sqrt(el2ys_accum)
    error_us  = sqrt(el2us_accum)
    error_pf  = sqrt(el2pf_accum)
    error_ph  = sqrt(el2ph_accum)

    error_uf, error_ur, error_ys, error_us, error_pf, error_ph, Gridap.FESpaces.num_free_dofs(X)  
  end
  
  function convergence_test(nkmax)

    euf   = Float64[]
    eur   = Float64[]
    eys   = Float64[]
    eus   = Float64[]
    epf   = Float64[]
    eph   = Float64[]
    
    ruf   = Float64[]
    rur   = Float64[]
    rys   = Float64[]
    rus   = Float64[]
    rpf   = Float64[]
    rph   = Float64[]
    hh    = Float64[]
    nn    = Int[]

    push!(ruf,0.)
    push!(rur,0.)
    push!(rys,0.)
    push!(rus,0.)
    push!(rpf,0.)
    push!(rph,0.)
    
    for nk in 1:nkmax
        println("******** Refinement step: $nk")
        model=generate_model_unit_square_biot_stokes(nk;simplexify_model=true,bcs_type=:full_dirichlet)
        push!(hh,1/2^nk)
        Δt = 1e-5

        error_uf,error_ur,error_ys,error_us,error_pf,error_ph,ndofs = solve_stokes_darcy_TH_RT(model,Δt; generate_output=true)
        push!(nn,ndofs)
        println("******** Total DoFs: ", nn[nk])
        push!(euf,error_uf)
        push!(eur,error_ur)
        push!(eys,error_ys)
        push!(eus,error_us)
        push!(epf,error_pf)
        push!(eph,error_ph)

        if nk>1
            push!(ruf, log(euf[nk]/euf[nk-1])/log(hh[nk]/hh[nk-1]))
            push!(rur, log(eur[nk]/eur[nk-1])/log(hh[nk]/hh[nk-1]))
            push!(rys, log(eys[nk]/eys[nk-1])/log(hh[nk]/hh[nk-1]))
            push!(rus, log(eus[nk]/eus[nk-1])/log(hh[nk]/hh[nk-1]))
            push!(rpf, log(epf[nk]/epf[nk-1])/log(hh[nk]/hh[nk-1]))
            push!(rph, log(eph[nk]/eph[nk-1])/log(hh[nk]/hh[nk-1]))
        end
    end
    println("====================================================================================================================================")
    println("   DoF  &    h   &   e(uf)  & r(uf) &   e(ur)  & r(ur) &   e(pf)  & r(pf) &   e(ph)  & r(ph) &   e(ys)  & r(ys)  &   e(us)  & r(us) ")
    println("====================================================================================================================================")
    for nk in 1:nkmax
        @printf("%7d & %.4f & %.2e & %.3f & %.2e & %.3f & %.2e & %.3f & %.2e & %.3f & %.2e & %.3f & %.2e & %.3f  \n", nn[nk], hh[nk], euf[nk], ruf[nk], eur[nk], rur[nk], epf[nk], rpf[nk], eph[nk], rph[nk],eys[nk], rys[nk],eus[nk], rus[nk]);
    end
    println("=====================================================================================================================================")
  end
  convergence_test(6)
end





