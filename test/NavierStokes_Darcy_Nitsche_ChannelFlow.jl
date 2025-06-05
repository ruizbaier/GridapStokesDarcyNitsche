using Gridap
using GridapGmsh
using LineSearches: BackTracking 
import Gridap: ∇
using Printf

# needed for the 2D cross product 
(Gridap.CellData.cross)(a::Gridap.CellData.CellField,
  b::Gridap.CellData.SkeletonPair{<:Gridap.CellData.CellField}) = Operation(cross)(a,b)

Id  = TensorValue(1,0,0,1)
const μ   = 0.01
const α   = 1.0
const κ   = 1.0E-4
const δ   = 1.0

# inlet Stokes and Darcy velocities
uin(x) = VectorValue(0.01*tanh(40*(x[2]+0.2)*(1-x[2])),0)

ff = VectorValue(0,0)
fp = VectorValue(0,0)
v0 = VectorValue(0,0)

model = GmshDiscreteModel("meshes/ComplexChannelWithInterface.msh")
writevtk(model,"meshes/ComplexChannel")
labels = get_face_labeling(model)

add_tag_from_tags!(labels,"FluidIn",["inF"])
add_tag_from_tags!(labels,"FluidOut",["outF"])
add_tag_from_tags!(labels,"FluidWall",["wallF"])
add_tag_from_tags!(labels,"FluidObstacles",["obsF"])

add_tag_from_tags!(labels,"PorousIn",["inP"])
add_tag_from_tags!(labels,"PorousOut",["outP"])
add_tag_from_tags!(labels,"PorousWall",["wallP"])

add_tag_from_tags!(labels,"Interface",["interf"])
add_tag_from_tags!(labels,"Fluid",["fluid"])
add_tag_from_tags!(labels,"Porous",["porous"])

# Reference FEs
order = 2 
reffe_uf = ReferenceFE(lagrangian,VectorValue{2,Float64},order)
reffe_pf = ReferenceFE(lagrangian,Float64,order-1)
reffe_up = ReferenceFE(raviart_thomas,Float64,order-2)
reffe_pp = ReferenceFE(lagrangian,Float64,order-2)

Ω   = Interior(model)
Ω_f = Interior(model,tags="Fluid")
Ω_p = Interior(model,tags="Porous")

#Boundary triangulations and outer unit normals are needed if imposing weak BCs there
#Γp_out = BoundaryTriangulation(model,tags="PorousOut")
#n_Γp_out = get_normal_vector(Γp_out)
#... and others if needed 

# Triangulation and normal on interface (pointing outwards with respect to Ω_S)
Σ = InterfaceTriangulation(Ω_f,Ω_p)
n_Σ  = get_normal_vector(Σ)

# Numerical integration
degree = 2*order
dΩ = Measure(Ω,degree)
dΩ_f = Measure(Ω_f,degree)
dΩ_p = Measure(Ω_p,degree)

bdegree = 2*order
dΓp_out = Measure(Γp_out,bdegree)
# ...

idegree = 2*order
dΣ = Measure(Σ, idegree)
h_e_Σ = CellField(get_array(∫(1)dΣ), Σ)

# FE spaces
Vf  = TestFESpace(Ω_f, reffe_uf, conformity=:H1, dirichlet_tags=["FluidIn","FluidWall","FluidObstacles"])
Vp  = TestFESpace(Ω_p, reffe_up, conformity=:HDiv, dirichlet_tags=["PorousIn","PorousWall"])

Qf = TestFESpace(Ω_f, reffe_pf, conformity=:C0)
Qp = TestFESpace(Ω_p, reffe_pp, conformity=:L2)

Uf = TrialFESpace(Vf,[uin,v0,v0])
Up = TrialFESpace(Vp,[uin,v0])
Pf = TrialFESpace(Qf)
Pp = TrialFESpace(Qp)

Y = MultiFieldFESpace([Vf,Vp,Qf,Qp])
X = MultiFieldFESpace([Uf,Up,Pf,Pp])

# Galerkin and Nitsche terms 
a11(u,v) = ∫((2*μ)*(ε(u)⊙ε(v)))dΩ_f + 
           ∫((μ*α/sqrt(κ))*(u.⁺×n_Σ.⁺)*(v.⁺×n_Σ.⁺))dΣ +
           ∫((δ/h_e_Σ)*(u.⁺⋅n_Σ.⁺)*(v.⁺⋅n_Σ.⁺))dΣ 
                     
a12(v,u) = ∫((δ/h_e_Σ)*(u.⁻⋅n_Σ.⁻)*(v.⁺⋅n_Σ.⁺))dΣ

b13(v,q) = ∫(-(∇⋅v)*q)dΩ_f #+ ∫(q.⁺*(v.⁺⋅n_Σ.⁺))dΣ

b14(v,q) = ∫(q.⁻*(v.⁺⋅n_Σ.⁺))dΣ
    
b24(v,q) = ∫(-(∇⋅v)*q)dΩ_p + ∫(q.⁻*(v.⁻⋅n_Σ.⁻))dΣ 
    
a22(u,v) = ∫((μ/κ)*(u⋅v))dΩ_p + ∫((δ/h_e_Σ)*(u.⁻⋅n_Σ.⁻)*(v.⁻⋅n_Σ.⁻))dΣ

b45(ψ,q) = ∫(ψ*q)dΩ_p           

lhs((uf,up,pf,pp),(vf,vp,qf,qp)) = 
    #   uf           up           pf           pp  
    # ---------------------------------------------------|---
    a11(uf,vf) + a12(vf,up) + b13(vf,pf) + b14(vf,pp) + #|vf
    a12(uf,vp) + a22(up,vp)              + b24(vp,pp) + #|vp
    b13(uf,qf)                                        + #|qf
    b14(uf,qp) + b24(up,qp)                             #|qp

c(u,w,v) =  ∫((∇(u)⋅w)⋅v)dΩ_f 

resid((uf,up,pf,pp),(vf,vp,qf,qp)) = lhs((uf,up,pf,pp),(vf,vp,qf,qp)) + c(uf,uf,vf) - ∫(vf⋅ff)dΩ_f - ∫(vp⋅fp)dΩ_p 

jacob((uf,up,pf,pp),(duf,dup,dpf,dpp),(vf,vp,qf,qp)) = lhs((duf,dup,dpf,dpp),(vf,vp,qf,qp)) + c(duf,uf,vf) + c(uf,duf,vf) 
        
# Build FE operator
oper = FEOperator(resid,jacob,X,Y) 

nls = NLSolver(show_trace=true, method=:newton, linesearch=BackTracking())
solver = FESolver(nls)
                     
ufh, uph, pfh, pph = solve(solver,oper)

writevtk(Ω_f,"NavierStokesDarcy_Channel_fluid",order=1,cellfields=["uf"=>ufh,"pf"=>pfh])
writevtk(Ω_p,"NavierStokesDarcy_Channel_porous",order=1,cellfields=["up"=>uph,"pp"=>pph])
