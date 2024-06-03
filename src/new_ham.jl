
#using LinearAlgebra, SparseArrays, WIGXJPFjl, BenchmarkTools
"""
eh(N) returns the square root of the N²|NK⟩ matrix element, √(N(N+1)). See eh2
for the non-square rooted version.
"""
eh(x::Number)::Float64 = √(x*(x+1))
"""
eh2(N) returns the the N²|NK⟩ matrix element, N(N+1). See eh for the
automatically square rooted version.
"""
eh2(x::Number)::Float64 = x*(x+1)
□rt(x::Number)::Float64 =√(x*(x>zero(x)))
fh(x::Number,y::Number)::Float64 = □rt((x-y)*(x+y+1))
jnred(j::Number,n::Number)::Float64 = √((2*j+1)*(2*n+1))
nred(n::Number)::Float64 = √(n*(n+1)*(2*n+1))
"""
powneg1(x) takes a number and returns (-1)^x. I realize this is a stupid looking
function to have but it evalutes every so slightly faster than just (-1)^x
"""
powneg1(k::Number)::Int = isodd(k) ? -1 : 1
"""
δ(x,y) takes two number and returns the Kronecker delta as a float. See δi for
the integer version
"""
δ(x::Number,y::Number)::Float64 = x==y
"""
δi(x,y) takes two number and returns the Kronecker delta as an integer. See δ
for the float version
"""
δi(x::Number,y::Number)::Int = x==y
T(l::Int,q::Int)::Int = l*(l+1) + q + 1
Tq(q::Int)::Int = 3 + q #quadrupole variant (only has 2nd rank components)
Tsr(l::Int,q::Int)::Int = δi(l,2) + abs(q) + 1 #Cs sr version, no 1st rk, & symm

"""
tplus!(a) replaces the matrix a with the sum of it and it's transpose. A dense
and sparse variant are available
"""
tplus!(a::Array{Float64,2})::Array{Float64,2} = a .+= permutedims(a)
function tplus!(a::SparseMatrixCSC{Float64, Int64})::SparseMatrixCSC{Float64, Int64}
   a .+= permutedims(a)
end

"""
qngen(j,s) generates the quantum numbers that the J dependent parts of the 
Hamiltonian processes. J is the total angular moment and S is the spin.
Returns a 2D array with Ns in the first column and Ks in the second 
"""
function qngen(j,s)
   ns, nd, ni, jsd = srprep(j,s)
   out = zeros(Int,jsd,2)
   for i in 1:length(ns)
      out[ni[i,1]:ni[i,2],1] .= ns[i]
      out[ni[i,1]:ni[i,2],2] = collect(Int,-ns[i]:ns[i])
   end
   #[n k]
   return out
end
function qnlabv(j,s,nf,vtm,σ)
   σt = σtype(nf,σ)
   nlist = Δlist(j,s)
   jsd = Int((2*j+1)*(2*s+1))
   vd = Int(vtm+1)
   out = zeros(Int,0,3)
   for n in nlist
      nd = Int(2*n+1)
      part = zeros(Int,nd,3)
      part[:,1] = fill(n,nd)
      part[:,2] = collect(Int,-n:n)
      part[:,3] = k2kc.(part[:,1],part[:,2])
      out = vcat(out,part)
   end
   out[:,2] = abs.(out[:,2])
   out = kron(ones(Int,vd),out)
   vtrray = kron(nf .* vtcoll(vtm,σt) .+ σ,ones(Int,jsd))
   out = hcat(fill(Int(2*j),size(out,1)),out,vtrray,fill(σ,jsd*vd))
   return out
end
####   new 2nd order   ####
hr2on(ns,ks,bk,bn) = @. bn*eh2(ns) + bk*ks^2 
hr2of1(ns,ks,dab) = @. dab*(ks-0.5)*fh(ns,ks-1)
hr2of2(ns,ks,bpm) = @. bpm*fh(ns,ks-1)*fh(ns,ks-2)
function hrot2(pr,qns)::SparseMatrixCSC{Float64, Int64}
   ns = view(qns,:,1)
   ks = view(qns,:,2)
   out = spzeros(size(ns,1),size(ns,1))
   #p0 = hr2on(ns,ks,pr[1],pr[2])
   out[diagind(out)] .= hr2on(ns,ks,pr[1],pr[2])
   #p1 = hr2of1(ns[2:end],ks[2:end], pr[4])
   out[diagind(out,1)] .= hr2of1(ns[2:end],ks[2:end], pr[4])
   #p2 = hr2of2(ns[3:end],ks[3:end], pr[3])
   out[diagind(out,2)] .= hr2of2(ns[3:end],ks[3:end], pr[3])
   #out = spdiagm(0=>p0,1=>p1,2=>p2)
   return Symmetric(dropzeros!(out))
end
"""
hrotest(pr,j,s) generates the 2nd order rotational Hamiltonian for the given J
and S pair. pr is an array of length 4 with values of BK, BN, B±, and Dab 
respectively
"""
function hrotest(pr,j,s)
   qns = qngen(j,s)
   out = hrot2(pr,qns)
   return out
end
function hrotest(pr,n)
   qns = qngen(n,0)
   out = hrot2(pr,qns)
   return out
end

function nsred(l::Int,nb,nk)
   return 0.5*( 
   √(nk*(nk+1)*(2*nk+1))*
      wig6j( 1, 1, l,
            nb,nk,nk)*powneg1(l) + 
   √(nb*(nb+1)*(2*nb+1))*
      wig6j( 1, 1, l,
            nk,nb,nb))
end
function jsred(j,s,nb,nk)
   return wig6j(nk, s, j,
                 s,nb, 1)*√((2*nb+1)*(2*nk+1))
end
function srelem(pr::Float64,l::Int,q::Int,j,s,nb,kb,nk,kk)#::Array{Float64,2}
   return pr*wig3j( nb,l,nk,
                   -kb,q,kk)*√(2*l+1)*
       nsred(l,nb,nk)*jsred(j,s,nb,nk)*powneg1(nb-nk-kb)
end
function hsr(pr,j,s,qns)::SparseMatrixCSC{Float64, Int64}
   ns = view(qns,:,1)
   ks = view(qns,:,2)
   le = size(ns,1)
   out = spzeros(le,le)
   #awkward special array of rank, component, prm ind, & sign
   #each col is a different parameter
   Ts = SA[0 1 1 2 2 2 2 2; 0 -1 1 0 -1 1 -2 2;
           1 2 2 3 4 4 5 5; 1  1 1 1 -1 1  1 1]
   for i in 1:8
      tv = view(Ts,:,i)
      prm = pr[tv[3]]*tv[4]
      if prm ≠ 0.0
      for a in 1:le, b in 1:le
         nb = ns[b]
         kb = ks[b]
         nk = ns[a]
         kk = ks[a]
         if abs(nb-nk)≤1 && (tv[2]+kk-kb)==0
            out[b,a] += srelem(prm,tv[1],tv[2], j,s,nb,kb,nk,kk)
         end#selection rule if
      end#sr ind for loop
      end#prm chck if
   end#sr term for loop
   dropzeros!(out)
   out .*= nred(s)*powneg1(j+s)
   return out
end

function wiginv(s::Number)::Float64
   if s<one(s)
      return 0.0
   else
      return inv(wig3j( s,2,s,
                       -s,0,s))
   end
end
function qured(j,s,nb,nk)
   return 0.25*jnred(nb,nk)*
               wig6j(j, s,nb,
                     2,nk, s)
end
function qulm(pr,q,j,s,nb,kb,nk,kk)#::Array{Float64,2}
   return pr*qured(j,s,nb,nk)*
#             δ(nb,nk)* #This line can be used to emulate the perturbative 
             wig3j( nb, 2,nk,
                   -kb, q,kk)*powneg1(nb+nk-kb+s+j)
end
function hqu(pr,j,s,qns)::SparseMatrixCSC{Float64, Int64}
   ns = view(qns,:,1)
   ks = view(qns,:,2)
   le = size(ns,1)
   out = spzeros(le,le)
   #awkward special array of rank, component, prm ind, & sign
   #each col is a different parameter
   Tq = SA[0 -1 1 -2 2; 
           1  2 2  3 3;
           1  1 1  1 1]
   for i in 1:5
      tv = view(Tq,:,i)
      prm = pr[tv[2]]*tv[3]
      if prm ≠ 0.0
      for a in 1:le, b in 1:le
         nb = ns[b]
         kb = ks[b]
         nk = ns[a]
         kk = ks[a]
         if abs(nb-nk)≤2 && (tv[2]+kk-kb)==0
            out[b,a] += qulm(prm,tv[1], j,s,nb,kb,nk,kk)
         end#selection rule if
      end#qu ind for loop
      end#prm chck if
   end#qu term for loop
   dropzeros!(out)
   out .*= nred(s)*wiginv(s)*powneg1(j+s)
   return out
end

function htor2(sof::Array{Float64},ms::Array{Int})::SparseMatrixCSC{Float64, Int64}
   out = sof[1]*pa_op(ms,2)
   out += sof[4].*(I(size(out,1)) .- cos_op(ms,1))
   return out
end

####   individual operators   ####

function nnss_check(a,b)::Int
   a = a*iseven(a) + (a-1)*isodd(a)
   b = b*iseven(b) + (b-1)*isodd(b)
   return min(a,b)
end
ns_el(j,s,p,n)::Float64 = (0.5*eh2(j) - eh2(n) - eh2(s))^p
function nnss_op(j,s,qns,a,b)::Diagonal{Float64, Vector{Float64}}
   c = nnss_check(a,b)
   a -= c
   b -= c
   @views out = eh.(qns[:,1]).^a .* ns_el.(j,s,c,qns[:,1]) .* eh(s)^b
   return Diagonal(out)
end

nz_op(qns,p)::Diagonal{Float64, Vector{Float64}} = @views out = Diagonal(qns[:,2].^p)

function np_op(qns,p::Int)::SparseMatrixCSC{Float64, Int64}
   ns = qns[1+p:end,1]
   part = ones(length(ns))
   if p ≤ length(ns)
      ks = qns[1+p:end,2]
      part = ones(length(ks))
      for o in 1:p
         part .*= fh.(ns,ks.-o)
      end
   end
   out = spzeros(size(qns,1),size(qns,1))
   out[diagind(out,-p)] = part
   return out
end
npm_op(qns::Matrix{Int64},p::Int) = Symmetric(np_op(qns,p),:L)
function iny_op(qns::Matrix{Int64},p::Int)
   out = np_op(qns,1-δi(p,0))
   if p≠0
      out .-= permutedims(out)
   end
   return dropzeros!(out)
end

function sqpart(j,s,q,bqn,kqn)::Float64
   nb = bqn[1]
   kb = bqn[2]
   nk = kqn[1]
   kk = kqn[2]
   return wig3j(nb,1,nk,-kb,q,kk)*wig6j(s,nb,j,nk,s,1)*jnred(nb,nk)*powneg1(-kb)
end
function sz_op(j,s,qns,p)#::SparseMatrixCSC{Float64, Int64}
   l = size(qns,1)
   out = spzeros(l,l)
   if s != zero(s)
      for a ∈ 1:l, b ∈ a:l
         if abs(qns[a,1]-qns[b,1])≤1 && qns[a,2]==qns[b,2]
            @views out[b,a] = sqpart(j,s,0,qns[b,:],qns[a,:])
         end
      end
      dropzeros!(out)
      out .*= nred(s)*powneg1(s+j+1)
   else
      out[diagind(out)] .+= 1.0
   end
   return Symmetric(out,:L)^p
end

function sq_op(j,s,q,qns)::SparseMatrixCSC{Float64, Int64}
   l = size(qns,1)
   out = spzeros(l,l)
   if s != zero(s)
      for a ∈ 1:l, b ∈ 1:l
         if abs(qns[a,1]-qns[b,1])≤1 && (q+qns[a,2]-qns[b,2])==0
            @views out[b,a] = sqpart(j,s,q,qns[b,:],qns[a,:])
         end
      end
      dropzeros!(out)
      out .*= nred(s)*powneg1(s+j+1+δ(1,q))*√2
      #the √2 is included to convert from spherical to cylinderical 
   else
      out[diagind(out)] .+= 1.0
   end
   return out
end

sp_op(j,s,qns,p)::SparseMatrixCSC{Float64, Int64} = sq_op(j,s,1,qns)^p
sm_op(j,s,qns,p)::SparseMatrixCSC{Float64, Int64} = sq_op(j,s,-1,qns)^p
spm_op(j,s,qns,p)::SparseMatrixCSC{Float64, Int64} = sp_op(j,s,qns,p) + sm_op(j,s,qns,p)

pa_op(ms,p)::Diagonal{Float64, Vector{Float64}} = Diagonal(ms.^p)
function cos_op(ms,p)::SparseMatrixCSC{Float64, Int64}
   if p==0
   out = I(size(ms,1))
   else
   out = fill(0.5, length(ms)-p)
   out = spdiagm(-p=>out, p=>out)
   end
   return out
end
function sin_op(ms,p)::SparseMatrixCSC{Float64, Int64}
   #this is actually sin/2i as we are moving immediately multiplying it by
   #the i from Ny
   out = fill(0.25, length(ms)-p)
   out = spdiagm(-p=>out, p=>out)
   return out
end

####   collected operators   ####

function rsr_op(j::Number,s::Number,qns::Array{Int,2},a::Int,b::Int,
         c::Int,d::Int,e::Int,f::Int,jp::Int)::SparseMatrixCSC{Float64, Int64}
   out = np_op(qns,e)
   out *= sp_op(j,s,qns,f)
   tplus!(out)
   out = sz_op(j,s,qns,d)*out*iny_op(qns,jp)
   out = nnss_op(j,s,qns,a,b)*nz_op(qns,c)*out
   return dropzeros!(out)
end

tor_op(ms,g,h,j)::SparseMatrixCSC{Float64, Int64} = pa_op(ms,g)*
                                                      cos_op(ms,h)*sin_op(ms,j)

function tsr_op(prm::Float64,j::Number,s::Number,qns::Array{Int,2},ms::Array{Int},
                  a::Int,b::Int,c::Int,d::Int,e::Int,f::Int,g::Int,h::Int,
                  jp::Int)::SparseMatrixCSC{Float64, Int64}
   out = 0.25*prm*rsr_op(j,s,qns,a,b,c,d,e,f,jp)
   out = kron(tor_op(ms,g,h,jp),out)
   return dropzeros!(tplus!(out))
end
function tsr_op(prm::Float64,j::Number,s::Number,qns::Array{Int,2},
            ms::Array{Int},plist::Array{Int})::SparseMatrixCSC{Float64, Int64}
   #1/2 from a+a', 1/2 from np^0 + nm^0
   out = 0.25*prm*rsr_op(j,s,qns,plist[1],plist[2],plist[3],
                           plist[4],plist[5],plist[6],plist[9])
   out = kron(tor_op(ms,plist[7],plist[8],plist[9]),out)
   return dropzeros!(tplus!(out))
end

####   final construction and collection functions   ####

function hjbuild(sof,cdf::Array,cdo::Array,j,s,nf,mc,σ)
   qns = qngen(j,s)
   ms = msgen(mc,nf,σ)
   ℋ = hrot2(sof[1:4],qns) 
   if s==0.5
      ℋ .+= hsr(sof[5:9],j,s,qns)
   elseif s≥1
      ℋ .+= hsr(sof[5:9],j,s,qns)
      ℋ .+= hqu(sof[10:12],j,s,qns)
   end
   ℋ = kron(I(length(ms)), ℋ)
   ℋ .+= kron(htor2(sof[13:16], ms), I(size(qns,1)))
   #if s≥0.5
   #ℋ .+= kron(pa_op(1,ms), sof[14]*nz_op(qns,1)) #+ sof[15]*npm_op(qns,1) + 
   #            sof[17]*sz_op(j,s,qns,1) + sof[18]*spm_op(j,s,qns,1))
   #else
   ℋ += kron(pa_op(ms,1), sof[14]*nz_op(qns,1))
   ℋ += kron(pa_op(ms,1), sof[15]*npm_op(qns,1))
   ℋ += kron(pa_op(ms,1), sof[17]*sz_op(j,s,qns,1))
   ℋ += kron(pa_op(ms,1), sof[17]*spm_op(j,s,qns,1))
   #end
   for i in 1:length(cdf)
      ℋ .+= tsr_op(cdf[i],j,s,qns,ms,cdo[:,i] )
   end
   return dropzeros!(ℋ)
end

function tsrdiag(ctrl,sof,cdf,cdo,nf,mcalc,j,s,σ,vtm)
   H = hjbuild(sof,cdf,cdo,j,s,nf,mcalc,σ)
   if true ∈ isnan.(H)
      @warn "FUCK!!! j=$j, σ=$σ, NaN in H"
   end
   #if σtype(nf,σ) != 1 #A & B states have more symmetry
   #   U = ur(j,s,mcalc,σtype(nf,σ))*ut(mcalc,σtype(nf,σ),j,s)
   #else
   #   U = ur(j,s,mcalc,σtype(nf,σ))
   #end
   U = kron(ut(mcalc,σtype(nf,σ)),ur(j,s))
   H = (U*H*U)
   ### All lines commented with ### are for the Jacobi routine
   ###perm = kperm(j,s,mcalc)
   ###H = permute!(H,perm,perm)
   ###H, rvecs = jacobisparse(H, 3)#Int(j+s)+mcalc)
   ###rvecs = U*rvecs
   vals, vecs = eigen!(Symmetric(Matrix(H)))
   #@show vals[1:2*Int(2j+1)] ./csl
   ###perm = assignperm(vecs)
   if ctrl["assign"]=="RAM36"
      perm = ramassign(vecs,j,s,mcalc,σtype(nf,σ),vtm)
      vals = vals[perm]
      vecs = vecs[:,perm]
   elseif ctrl["assign"]=="expectk"
      vals, vecs = expectkassign!(vals,vecs,j,s,nf,mcalc,σ)      
   elseif ctrl["assign"]=="eeo"
      vals, vecs = eeoassign!(vals,vecs,j,s,nf,mcalc,σ)
   else
      vals, vecs = expectassign!(vals,vecs,j,s,nf,mcalc,σ)
   end
   ###vecs = rvecs*vecs 
   vecs = U*vecs
   pasz = zero(vals)
   if s != zero(s)#add η
      pasz = diag(vecs' * tsrop(1.0,0,0,0,0,1,1,0,0,
         j,s,permutedims(ngen(j,s)),
         permutedims(kgen(j,s)),mb,ngen(j,s),kgen(j,s),mk) * vecs)
   end
   return vals, vecs, pasz
end

function tsrcalc(ctrl,prm,stg,cdo,nf,vtm,mcalc,jlist,s,sd,σ)
   sof = prm[1:18]
   cdf = prmsetter(prm[19:end],stg)
   mcd = Int(2*mcalc+(σtype(nf,σ)==2)+1)
   vtd = Int(vtm+1)
   σt = σtype(nf,σ)
   jmin = 0.5*iseven(sd)
   jmax = jlist[end]
   jfd = sd*Int(sum(2.0 .* collect(Float64,jmin:jmax) .+ 1.0))
   msd = sd*mcd
   #mstrt, mstop = mslimit(nf,mcalc,σ)
   outvals = zeros(Float64,jfd*vtd)
   outpasz = zeros(Float64,jfd*vtd)
   outquns = zeros(Int,jfd*vtd,6)
   outvecs = zeros(Float64,Int(sd*(2*jmax+1)*mcd),jfd*vtd)
   for j in jlist #thread removed for troubleshooting purposes
#   @threads for j in jlist
      jd = Int(2.0*j) + 1
      sind, find = jvdest(j,s,vtm) 
      tvals, tvecs, tpasz = tsrdiag(ctrl,sof,cdf,cdo,nf,mcalc,j,s,σ,vtm)
      outvals[sind:find] = tvals
      outpasz[sind:find] = tpasz
      outquns[sind:find,:] = qnlabv(j,s,nf,vtm,σ)
      outvecs[1:jd*msd,sind:find] = tvecs###[:,pull]
   end
   return outvals, outvecs, outquns, outpasz
end

function tsrcalc2(prm,stg,cdo,nf,ctrl,jlist)
   s = ctrl["S"]
   mcalc = ctrl["mcalc"]
   vtm = ctrl["vtmax"]
   #sd = Int(2*s + 1)
   sof = prm[1:18]
   cdf = prmsetter(prm[19:end],stg)
   #println(cdf)
   vtd = Int(vtm+1)
   jmin = 0.5*iseven((2*s + 1))
   jmax = jlist[end,1]
   jfd = Int((2s+1)*sum(2.0 .* collect(Float64,jmin:jmax) .+ 1.0))
   σcnt = σcount(nf)
   mcd = Int(2*mcalc+(iseven(nf))+1)
   fvls = zeros(Float64,jfd*vtd,σcnt)
   fqns = zeros(Int,jfd*vtd,6,σcnt)
   fvcs = zeros(Float64,Int((2*s + 1)*(2*jmax+2)*mcd),jfd*vtd,σcnt)
   @time for sc in 1:σcnt
      σ = sc - 1
      #mcd = Int(2*mcalc+(σtype(nf,σ)==2)+1)
      mcd = Int(2*mcalc+(σtype(nf,σ)==2)+1)
      msd = Int(2*s + 1)*mcd
      mstrt, mstop = mslimit(nf,mcalc,σ)
      jmsd = Int(mcd*(2s+1)*(2*jmax+1))
      jsvd = Int(jfd*vtd)
      jsublist = jlist[isequal.(jlist[:,2],σ), 1] .* 0.5
      @threads for j in jsublist
         jd = Int(2.0*j) + 1
         #pull = indpuller(vtm,mcalc,σt,Int(jd*sd))
         sind, find = jvdest(j,s,vtm) 
         fvls[sind:find,sc], fvcs[1:jd*msd,sind:find,sc], = tsrdiag(ctrl,
            sof,cdf,cdo,nf,mcalc,j,s,σ,vtm)
         #fvls[sind:find,sc] = tvals#[pull]
         #fvcs[1:jd*msd,sind:find,sc] = tvecs#[:,pull]
         fqns[sind:find,:,sc] = qnlabv(j,s,nf,vtm,σ)
      end
   end
   return fvls, fvcs, fqns
end
