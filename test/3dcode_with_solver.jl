using Gridap
using GridapGmsh
import Gridap: ∇
using Printf
using GridapPETSc
using LineSearches: BackTracking 

# needed for the 3D cross product 
(Gridap.CellData.cross)(a::Gridap.CellData.CellField,
  b::Gridap.CellData.SkeletonPair{<:Gridap.CellData.CellField}) = Operation(cross)(a,b)

options = ["-ksp_type", "gmres",
           "-pc_type", "ilu",
           "-pc_factor_levels", "3",
           "-ksp_monitor"]


GridapPETSc.with(args=options) do

    Id  = TensorValue(1,0,0,0,1,0,0,0,1)
    μ     = 0.01 
    ρf    = 1.0     
    α     = 1.0       
    δ     = 10.0      
    κ     = 1e-3      
    φ     = 0.30    
    λP    = 4.9364e4 
    μP    = 1.0336e3
    K     = 0.49e3         
    ρp    = 1.2
    θ_sink = 0.0

    uin(t) = x -> VectorValue(20*(x[2]-0.5)*(1.5-x[2])*x[3]*(0.8-x[3])/0.64, 0.0, 0.0)

    ap = 0
    ff = VectorValue(0.0,0.0,0.0)
    vzero(t) = x -> VectorValue(0.0,0.0,0.0)
    model = GmshDiscreteModel("meshes/mesh3d.msh")
    labels = get_face_labeling(model)
    
    add_tag_from_tags!(labels,"FluidIn",["Inlet"])
    add_tag_from_tags!(labels,"FluidOut",["Outlet"])
    add_tag_from_tags!(labels,"FluidWall",["LateralWalls"])
    add_tag_from_tags!(labels,"ChannelT",["ChannelTop"])
    add_tag_from_tags!(labels,"ChannelB",["ChannelBottom"])

    add_tag_from_tags!(labels,"CylinderT",["CylinderTop"])
    add_tag_from_tags!(labels,"CylinderB",["CylinderBottom"])

    add_tag_from_tags!(labels,"Interface",["Interf"])
    add_tag_from_tags!(labels,"Fluid",["FluidDomain"])
    add_tag_from_tags!(labels,"Porous",["PorousDomain"])

    # Reference FEs
    order = 2 
    reffe_uf = ReferenceFE(lagrangian,VectorValue{3,Float64},2)
    reffe_pf = ReferenceFE(lagrangian,Float64,1)
    reffe_ur = ReferenceFE(lagrangian,VectorValue{3,Float64},2)
    reffe_ph = ReferenceFE(lagrangian,Float64,1)
    reffe_ys = ReferenceFE(lagrangian,VectorValue{3,Float64},2)
    reffe_us = ReferenceFE(lagrangian,VectorValue{3,Float64},1)

    Ω   = Interior(model)
    Ω_S = Interior(model,tags="Fluid")
    Ω_D = Interior(model,tags="Porous")

    # Triangulation and normal on interface (pointing outwards with respect to Ω_S)
    Σ = InterfaceTriangulation(Ω_S,Ω_D)
    n_Σ  = get_normal_vector(Σ)

    # Numerical integration
    degree = 3*order
    dΩ = Measure(Ω,degree)
    dΩ_S = Measure(Ω_S,degree)
    dΩ_D = Measure(Ω_D,degree)

    idegree = 3*order
    dΣ = Measure(Σ, idegree)
    h_e_Σ = CellField(get_array(∫(1)dΣ),Σ)

    # FE spaces
    Vf  = TestFESpace(Ω_S, reffe_uf, conformity=:H1, dirichlet_tags=["FluidIn","FluidWall","ChannelT", "ChannelB"])
    Vr  = TestFESpace(Ω_D, reffe_ur, conformity=:H1, dirichlet_tags=["CylinderT","CylinderB"])
    Vs  = TestFESpace(Ω_D, reffe_ys, conformity=:H1, dirichlet_tags=["CylinderT","CylinderB"])
    Ws = TestFESpace(Ω_D, reffe_us, conformity=:C0)
    Qf = TestFESpace(Ω_S, reffe_pf, conformity=:C0)
    Qh = TestFESpace(Ω_D, reffe_ph, conformity=:C0)

    Uf = TransientTrialFESpace(Vf,[uin,vzero,vzero,vzero])
    Ur = TransientTrialFESpace(Vr,[vzero,vzero])
    Us = TransientTrialFESpace(Vs,[vzero,vzero])
    Ys = TransientTrialFESpace(Ws)
    Pf = TransientTrialFESpace(Qf)
    Ph = TransientTrialFESpace(Qh)

    Y = MultiFieldFESpace([Vf,Vr,Vs,Ws,Qf,Qh])
    X = MultiFieldFESpace([Uf,Ur,Us,Ys,Pf,Ph])
    
    # Galerkin and Nitsche terms 
    a1(t,uf,vf)   = ∫((2*μ)*(ε(uf)⊙ε(vf)))dΩ_S + 
                    ∫((μ*α/sqrt(κ))*((uf.⁺×n_Σ.⁺)⋅(vf.⁺×n_Σ.⁺)))dΣ + 
                    ∫((δ*μ/h_e_Σ)*(uf.⁺⋅n_Σ.⁺)*(vf.⁺⋅n_Σ.⁺))dΣ - 
                    ∫((2*μ)*((ε(uf).⁺⋅n_Σ.⁺)⋅n_Σ.⁺)*(vf.⁺⋅n_Σ.⁺))dΣ - 
                    ∫((2*μ)*((ε(vf).⁺⋅n_Σ.⁺)⋅n_Σ.⁺)*(uf.⁺⋅n_Σ.⁺))dΣ
                
    a2(t,ur,vf)   = ∫((δ*μ/h_e_Σ)*(ur.⁻⋅n_Σ.⁻)*(vf.⁺⋅n_Σ.⁺))dΣ - 
                    ∫((2*μ)*((ε(vf).⁺⋅n_Σ.⁺)⋅n_Σ.⁺)*(ur.⁻⋅n_Σ.⁻))dΣ
                
    a3(t,pf,vf)   = ∫(-(∇⋅vf)*pf)dΩ_S + ∫(pf.⁺*(vf.⁺⋅n_Σ.⁺))dΣ 
    a4(t,uf,qf)   = ∫((∇⋅uf)*qf)dΩ_S - ∫(qf.⁺*(uf.⁺⋅n_Σ.⁺))dΣ 
    a5(t,ys,ws)   = ∫((2*μP)*(ε(ys)⊙ε(ws)))dΩ_D + ∫(λP*(∇⋅ys)*(∇⋅ws))dΩ_D 

    a6(t,ph,ws)   = ∫(-(∇⋅ws)*ph)dΩ_D 

    a7(t,uf,ws)   = ∫(-(μ*α/sqrt(κ))*((uf.⁺×n_Σ.⁺)⋅(ws.⁻×n_Σ.⁺)))dΣ +
                    ∫((δ*μ/h_e_Σ)*(uf.⁺⋅n_Σ.⁺)*(ws.⁻⋅n_Σ.⁻))dΣ -
                    ∫((2*μ)*((ε(uf).⁺⋅n_Σ.⁺)⋅n_Σ.⁺)*(ws.⁻⋅n_Σ.⁻))dΣ
        
    a8(t,ur,ws)   = ∫((δ*μ/h_e_Σ)*(ur.⁻⋅n_Σ.⁻)*(ws.⁻⋅n_Σ.⁻))dΣ +
                    ∫((2*μ*φ)*(ε(ur)⊙ε(ws)))dΩ_D -∫(θ_sink*(ur⋅ws))dΩ_D 
            
    a9(t,ur,vr)   =  ∫((δ*μ/h_e_Σ)*(ur.⁻⋅n_Σ.⁻)*(vr.⁻⋅n_Σ.⁻))dΣ + 
                     ∫((μ*α/sqrt(κ))*((ur.⁻×n_Σ.⁻)⋅(vr.⁻×n_Σ.⁻)))dΣ + 
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
    a18(t,us,vs)  = ∫(-(ρp)*(us⋅vs))dΩ_D  
                  
    b1(t,dys,vf)   = ∫(-(μ*α/sqrt(κ))*((dys.⁻×n_Σ.⁺)⋅(vf.⁺×n_Σ.⁺)))dΣ + 
                     ∫((δ*μ/h_e_Σ)*(dys.⁻⋅n_Σ.⁻)*(vf.⁺⋅n_Σ.⁺))dΣ -
                     ∫((2*μ)*((ε(vf).⁺⋅n_Σ.⁺)⋅n_Σ.⁺)*(dys.⁻⋅n_Σ.⁻))dΣ
         
    b2(t,dys,ws)   = ∫((μ*α/sqrt(κ))*((dys.⁻×n_Σ.⁺)⋅(ws.⁻×n_Σ.⁺)))dΣ + 
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
    b12(t,dys,vs)   = ∫(ρp*(dys⋅vs))dΩ_D
        
    lvf(t,vf)     = ∫(vf⋅ff)dΩ_S 
    lvr(t,vr)     = ∫(vr⋅ff)dΩ_D 
    lys(t,ws)     = ∫(ff⋅ws)dΩ_D               
    lqf(t,qf)     = ∫(qf*ap)dΩ_S 
    lqh(t,qh)     = ∫(qh*ap)dΩ_D 
    lvs(t,vs)     = ∫(vs⋅ff)dΩ_D
     
    lhs(t,(uf,ur,ys,us,pf,ph),(vf,vr,ws,vs,qf,qh)) =  a1(t,uf,vf)+a2(t,ur,vf)+a3(t,pf,vf)+a4(t,uf,qf)+a5(t,ys,ws)+
                                    a6(t,ph,ws)+a7(t,uf,ws)+a8(t,ur,ws)+a9(t,ur,vr)+a10(t,ph,vr)+
                                    a11(t,uf,vr)+a12(t,ur,qh)+a13(t,pf,vr)+a14(t,pf,ws)+a15(t,ur,qf)+
                                    a16(t,us,vr)+a17(t,us,ws)+a18(t,us,vs)
    
    mass(t,(duf,dur,dys,dus,pf,dph),(vf,vr,ws,vs,qf,qh)) = b1(t,dys,vf)+b2(t,dys,ws)+b3(t,dys,vr)+b4(t,dys,qh)+b5(t,dph,qh)+
                                    b6(t,dys,qf)+b7(t,duf,vf)+b8(t,dur,vr)+b9(t,dus,vr)+b10(t,dur,ws)+
                                    b11(t,dus,ws)+b12(t,dys,vs) 
    
    
    rhs(t,(vf,vr,ws,vs,qf,qh)) = lvf(t,vf) + lvr(t,vr) + lys(t,ws) + lvs(t,vs) + lqf(t,qf) +lqh(t,qh) 
    
    c(t,u,w,v) =  ∫((∇(u)⋅w)⋅v)dΩ_S 
    
    jac_m(t,(uf,ur,ys,us,pf,ph),(duf,dur,dys,dus,dpf,dph),(vf,vr,ws,vs,qf,qh)) = mass(t,(duf,dur,dys,dus,dpf,dph),(vf,vr,ws,vs,qf,qh))
    jac_s(t,(uf,ur,ys,us,pf,ph),(duf,dur,dys,dus,dpf,dph),(vf,vr,ws,vs,qf,qh)) = lhs(t,(duf,dur,dys,dus,dpf,dph),(vf,vr,ws,vs,qf,qh)) + c(t,duf,uf,vf) + c(t,uf,duf,vf) 
    res(t,(uf,ur,ys,us,pf,ph),(vf,vr,ws,vs,qf,qh)) = mass(t,(uf,ur,ys,us,pf,ph),(vf,vr,ws,vs,qf,qh)) + lhs(t,(uf,ur,ys,us,pf,ph),(vf,vr,ws,vs,qf,qh)) + c(t,uf,uf,vf) -rhs(t,(vf,vr,ws,vs,qf,qh))

    jac(t,(uf,ur,ys,us,pf,ph),(duf,dur,dys,dus,dpf,dph),(vf,vr,ws,vs,qf,qh))   = jac_s(t,(uf,ur,ys,us,pf,ph),(duf,dur,dys,dus,dpf,dph),(vf,vr,ws,vs,qf,qh))
    jac_t(t,(uf,ur,ys,us,pf,ph),(duf,dur,dys,dus,dpf,dph),(vf,vr,ws,vs,qf,qh)) =  jac_m(t,(uf,ur,ys,us,pf,ph),(duf,dur,dys,dus,dpf,dph),(vf,vr,ws,vs,qf,qh))

    # Build FE operator
    op = TransientSemilinearFEOperator(mass, res, (jac, jac_t), X, Y)
    ksp = PETScLinearSolver()
    nl_solver = NLSolver(ksp, method=:newton, iterations=50, show_trace=false, linesearch=BackTracking())
    θ = 1.0
    Δt = 5e-3
    solver = ThetaMethod(nl_solver, Δt, θ)
    t0, tF = 0.0, 1

    zero_vec = VectorValue(0.0, 0.0, 0.0)
    zero_scalar = 0.0

    xh0 = interpolate_everywhere([zero_vec, zero_vec, zero_vec, zero_vec, zero_scalar, zero_scalar], X(t0))
    fesltn = solve(solver, op, t0, tF, xh0)

    allsol = collect(fesltn)  # collect all time steps
    tfinal, xfinal_state = allsol[end]  # get final time and solution

    ufh = xfinal_state[1] 
    urh = xfinal_state[2]  
    ysh = xfinal_state[3] 
    ush = xfinal_state[4]  
    pfh = xfinal_state[5] 
    pph = xfinal_state[6]

    writevtk(Ω_S, "NS_GP_Channel_fluid", order=1, cellfields=["uf"=>ufh, "pf"=>pfh])
    writevtk(Ω_D, "NS_GP_Channel_porous", order=1, cellfields=["ur"=>urh, "ph"=>pph, "ys"=>ysh, "us"=>ush])

end # GridapPESTc
