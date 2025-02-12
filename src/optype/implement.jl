
"""
This file is the actual implementation of the information in types.jl & baseops.jl
"""

using DelimitedFiles

include("@__DIR__/../type.jl")
include("@__DIR__/../baseops.jl")
include("@__DIR__/../file_in.jl")

const csl = 29979.2458

function westereng(molnam::String)
   ctrl = ctrlinp(molnam)
   val, errs, ℋ, stgs = secordinp(molnam,ctrl)
   ℋ, stgs, errs = opreader(molnam,ctrl,ℋ,stgs,errs)
   ℋ = stgvalset(ℋ,stgs)
   #initialize vals
   #initialize vecs
   #initialize tvecs
   #initialize qns
end

function jvdest2(j,s,vtm)
"""
This returns the first and final indices for a certain J value for a given S.
   This is used to place the eigenvalues & vectors in the final large arrays
"""
   snd = convert(Int, (vtm+1)*(2*s+1)*sum(2 .*collect((0.5*isodd(2*s)):(j-1)) .+1))+1
   fnd = convert(Int, (vtm+1)*(2*s+1)*sum(2 .*collect((0.5*isodd(2*s)):j) .+1))
   return snd:fnd
end

#tsrcalc2(prm,stg,cdo,nf,ctrl,jlist)
function tsrcalc!(vals,vecs,qns,jlist,ctrl,prm,stg,ℋ)
   nfold = ctrl["NFOLD"]
   mc = ctrl["mcalc"]
   s = ctrl["S"]
   msd = convert(Int,(2s+1)*prod(2 .*mc.+1))
   σs = σgen_indef(nfold)
   σcnt = size(σs,2)
for j in jlist
   jd = Int(2.0*j) + 1
   dest = jvdest2(j,s,vtm) 
   ψ = Psi(J=j,S=ctrl["S"])
   Hrot = hrot2(prm[1:4],ψ)
   if S≥1.0
      Hrot += Hsr(prm[5:9],ψ.J,ψ.S,ψ) + Hqua(prm[10:12],ψ.J,ψ.S,ψ)
   elseif S==0.5
      Hrot += Hsr(prm[5:9],ψ.J,ψ.S,ψ)
   end
   Hrot = Symmetric(Hrot)
   for sc in 1:σcnt
      ψ = Psi(Psi(J=j,S=s),nf=nfold,σ=σs[σd],mc=mc)
      vals[dest,sc],vecs[1:jd*msd,dest,sc] = tsrdiag_1(Hrot,ℋ,ψ)
      qns[dest,sc] = qnlabv(j,s,nf,vtm,σ)
   end#σs
end#j
end#f

#tsr_diag variants: 0 generic universal, 1 preloads Hrs, 2 two-stage preloads Hrs
function tsrdiag_0(ℋ::Vector{Op},ψ::Psi)
   H = enact(ℋ,ψ)
   U = sparse(ones(1))
   for i in 1:length(ψ.nf)
#having U at the end of kron maintians the tight block structure of the lower index tops
   if ψ.σ[i] == 0
      U = kron(ur(mcalc),U)
   else
      U = kron(sparse(1.0I,2mcalc+1,2mcalc+1),U)
   end;end
   U = kron(u,ur(ψ.J,ψ.S))
   H = droptol!(U*H*U,2*eps())
   vals,vecs = eigen!(Symmetric(Matrix(H)))
   if (ctrl["assign"]=="ram36")||(ctrl["assign"]=="RAM36")
      perm = ramassign(vecs,j,s,mcalc,vtm)
      vals = vals[perm]
      vecs = vecs[:,perm]
   elseif ctrl["assign"]=="expectk"
      vals, vecs = expectkassign!(vals,vecs,j,s,nf,mcalc,σ)      
   elseif ctrl["assign"]=="eeo"
      vals, vecs = eeoassign!(vals,vecs,j,s,nf,mcalc,σ)
   else
      vals, vecs = expectassign!(vals,vecs,j,s,nf,mcalc,σ)
   end
   return vals, vecs
end

function tsrdiag_1(Hr::SparseMatrixCSC{Float64,Int},ℋ::Vector{Op},ψ::Psi)
   H = Hr + enact(ℋ,ψ)
   U = sparse(ones(1))
   for i in 1:length(ψ.nf)
   if ψ.σ[i] == 0
      U = kron(u,ur(mcalc))
   else
      U = kron(u,sparse(1.0I,2mcalc+1,2mcalc+1))
   end;end
   U = kron(u,ur(ψ.J,ψ.S))
   H = U*H*U
   vals,vecs = eigen!(Symmetric(Matrix(H)))
   if (ctrl["assign"]=="ram36")||(ctrl["assign"]=="RAM36")
      perm = ramassign(vecs,j,s,mcalc,vtm)
      vals = vals[perm]
      vecs = vecs[:,perm]
   elseif ctrl["assign"]=="expectk"
      vals, vecs = expectkassign!(vals,vecs,j,s,nf,mcalc,σ)      
   elseif ctrl["assign"]=="eeo"
      vals, vecs = eeoassign!(vals,vecs,j,s,nf,mcalc,σ)
   else
      vals, vecs = expectassign!(vals,vecs,j,s,nf,mcalc,σ)
   end
   return vals, vecs
end