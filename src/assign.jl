"""
This is where the assignment routines are located for westerfit.
   The is a simple assigner based on RAM36. There are also the start of a 
   Jacobi eigenvalue routine based version. It is not done yet
"""


### EXPECTATION
#july 28 24: fun discovery! this breaks for NFOLD=2
kperm(n::Int)::Array{Int} = sortperm(Int.(cospi.(collect(-n:n).+isodd(n))) .* collect(-n:n))
keperm(n::Int)::Array{Int} = sortperm(sortperm(collect(-n:n), by=abs))[kperm(n)]

function mexpect(vecs,jsd,nf,mc,σ)
   m2 = kron(pa_op(msgen(nf,mc,σ),2) ^2, I(jsd))
   mlist = diag(vecs' * m2 * vecs)
   #println(mlist)
   mlist = sort(sortperm(mlist, by=abs)[1:jsd])
   return mlist
end

function nexpect(vecs,mcd,j,s,jsd,ns,nd,ni)
   n2 = diag(vecs' * kron(I(mcd),nt2_op(qngen(j,s),2)) * vecs)
   nlist = zeros(Int,jsd)
   for i in 1:length(ns)
      n = ns[i]
      part = sort(sortperm(n2 .- eh(n)^2, by=abs)[1:nd[i]])
      #part = part[kexpect(vecs[:,part], mcd,n,nd[i],j,s)]
      nlist[ni[i,1]:ni[i,2]] = part[keperm(n)]
      n2[part] .= 0.0
   end
   return nlist
end
function expectassign!(vals,vecs,j,s,nf,mc,σ)
   ns, nd, ni, jsd = srprep(j,s)
   list = mexpect(vecs,jsd,nf,mc,σ)
   #println(list)
   md = 2*mc + 1
   vals = vals[list]
   vecs = vecs[:,list]
   list = nexpect(vecs,md,j,s,jsd,ns,nd,ni)
   vals = vals[list]
   vecs = vecs[:,list]   
   return vals, vecs
end

### EXPECTATION WITH SYMMETRY ON K
function ksymer(vecs,npart,n)
   vecs = vecs[kperm,:] # symmetrerizes vectors into A2 & A1 blocks
   A1 = zeros(Int,n+iseven(n))
   A2 = zeros(Int,n+isodd(n))
   for i in 1:length(npart)

   end
end

### EXPECATION WITH K
kperm2(n::Int)::Array{Int} = sortperm(Int.(cospi.(collect(-n:n).+iseven(n))) .* collect(-n:n))
kaperm(n::Int)::Array{Int} = sortperm(sortperm(collect(-n:n), by=abs))[kperm2(n)]

function kexpect(vecs,mcd,n,nd,j,s)
   nz = diag(vecs' * kron(I(mcd),nz_op(qngen(j,s),2)) * vecs)
   perm = sortperm(nz, by=abs)[1:nd]
   perm = perm[kaperm(n)]
   return perm
end
function nkexpect(vecs,mcd,j,s,jsd,ns,nd,ni)
   n2 = diag(vecs' * kron(I(mcd),nt2_op(qngen(j,s),1)) * vecs)
   nlist = zeros(Int,jsd)
   for i in 1:length(ns)
      n = ns[i]
      part = sort(sortperm(n2 .- eh2(n), by=abs)[1:nd[i]])
      part = part[kexpect(vecs[:,part], mcd,n,nd[i],j,s)]
      nlist[ni[i,1]:ni[i,2]] = part[keperm(n)]
      n2[part] .= 0.0
   end
   return nlist
end

function expectkassign!(vals,vecs,j,s,nf,mc,σ)
   ns, nd, ni, jsd = srprep(j,s)
   list = mexpect(vecs,jsd,nf,mc,σ)
   #println(list)
   md = 2*mc + 1 
   vals = vals[list]
   vecs = vecs[:,list]
   list = nkexpect(vecs,md,j,s,jsd,ns,nd,ni)
   vals = vals[list]
   vecs = vecs[:,list]   
   return vals, vecs
end

"""
the new k assignment
avecs = abs.(vecs)
perm = zeros(Int,2N+1)
for i in 1:(2N+1)
   a,b = Tuple(argmax(avecs))
   perm[b] = a
   arv5[:,b] *= 0.0
   arv5[a,:] *= 0.0
end
"""
### Expect-Expect-Overlap
function eeoassign!(vals,vecs,j,s,nf,mc,σ)
   ns, nd, ni, jsd = srprep(j,s)
   list = mexpect(vecs,jsd,nf,mc,σ)
   #println(list)
   md = 2*mc + 1
   vals = vals[list]
   vecs = vecs[:,list] #simplifies down to just ground tor-state
   list = neko(vecs,md,j,s,jsd,ns,nd,ni)
   vals = vals[list]
   vecs = vecs[:,list]   
   return vals, vecs
end
function neko(vecs,mcd,j,s,jsd,ns,nd,ni)
   n2 = diag(vecs' * kron(I(mcd),nt2_op(qngen(j,s),1)) * vecs)
   nlist = zeros(Int,jsd)
   for i in 1:length(ns)
      n = ns[i]
      part = sort(sortperm(n2 .- eh2(n), by=abs)[1:nd[i]])
      part = part[koverlap(vecs,part,ni[i,1],ni[i,2],nd[i],n)]
      nlist[ni[i,1]:ni[i,2]] = part
      n2[part] .= 0.0
   end
   return nlist
end
function koverlap(vecs,part,ns,nf,nd,n)
   avecs = abs.(vecs[ns:nf,part])[kperm(n),:]
   piece = zeros(Int,nd)
   for j in 1:nd 
      a,b = Tuple(argmax(avecs))
      piece[a] = b
      avecs[:,b] *= 0.0
      avecs[a,:] *= 0.0
   end
   piece[1:(n+isodd(n))] = sort(piece[1:(n+isodd(n))],rev=true)
   piece[(n+isodd(n)+1):end] = sort(piece[(n+isodd(n)+1):end])
   return piece
end

### RAM36
#sum across each m to find dominate torsional state
function mfinder(svcs,jsd::Int,md::Int,mcalc,vtmax)
   ovrlp = zeros(md,size(svcs,2))
   @simd for i in 1:md 
      ovrlp[i,:] = sum(svcs[(jsd*(i-1)+1):(jsd*i), :], dims=1)
   end
   #ovrlp = argmax(ovrlp,dims=1)
   mind = zeros(Int,size(svcs,1))
   #@simd for i in 1:length(mind)
   #   mind[i] = ovrlp[i][1]
   #end
   cap = min(vtmax+4,md)
   for v in 0:vtmax#cap #THIS HAS TO BE SERIAL DON'T SIMD THIS ONE FUTURE WES
      mg = mcalc + vt2m(v) + 1
      perm = sort(sortperm(ovrlp[mg,:], rev=true)[1:jsd])
      ovrlp[:,perm] .= 0.0
      ovrlp[mg,:] .= 0.0
      mind[perm] .= mg
   end
   #println(mind)
   return mind
end
#4/9 11:30pm Something is fucked up for spin+low barrier. FUCK QNS
function mfinderv2(svcs,nind,ns,jsd,md,mcalc,vtmax)
   mind = zeros(Int,size(svcs,1))
   for ng in 1:length(ns)
      n = ns[ng]
      mind = mfinderforagivenn(svcs,mind,nind,ng,n,jsd,md,mcalc,vtmax)
   end
   return mind
end
function mfinderforagivenn(svcs,mind,nind,ng,n,jsd,md,mcalc,vtmax)
   list = collect(1:size(svcs,1))[nind .== ng]
   tvcs = svcs[:,list]
   ovrlp = zeros(md,size(tvcs,2))
   @simd for i in 1:md
      ovrlp[i,:] = sum(tvcs[(jsd*(i-1)+1):(jsd*i), :], dims=1)
   end
   part = zeros(Int,length(list))
   nd = 2*n+1
   for v in 0:vtmax
      mg = mcalc + vt2m(v) + 1
      perm = sort(sortperm(ovrlp[mg,1:nd*(vtmax+2)], rev=true)[1:nd])
      ovrlp[mg,:] .= 0.0 #prevents reassigning
      ovrlp[:,perm] .= 0.0 #prevents reassigning
      #println(perm)
      part[perm] .= mg
   end
   #println(size(part))
   #println(size(mind[list]))
   mind[list] .= part
   return mind
end
function nfinder(svcs,vtmax,md,jd,sd,ns,ni)
#this needs to be fully reworked it is only looking at the top part of the vector
   jsd = jd*sd
   ovrlp = zeros(length(ns),size(svcs,2))
   for i in 1:length(ns) 
      nd = 2*ns[i]+1
      for m in 1:md
         frst = ni[i,1] + jsd*(m-1)
         last = ni[i,2] + jsd*(m-1)
         ovrlp[i,:] += transpose(sum(svcs[frst:last, :], dims=1))
      end
   end
   nind = zeros(Int,size(svcs,1))
   for i in 1:length(ns)
      nd = (2*ns[i]+1)
      count = nd*min((vtmax+4),md)
      #vlimit = min(vtmax+3,md) 
      #for v in 0:vlimit
         perm = sort(sortperm(ovrlp[i,:], rev=true)[1:count])#[1:nd]
         nind[perm] .= i
         ovrlp[:,perm] .= 0.0 
         ovrlp[i,:] .= 0.0
      #end
   end
   #println(nind)
   return nind
end
function nfinderform(svcs,mind,mg,vtmax,md,jd,sd,ns,ni)
#this needs to be fully reworked it is only looking at the top part of the vector
#as of sept24, I'm very sure this is the problem. No! the issue is the lack of update & expansion
   jsd = jd*sd
   list = collect(1:size(svcs,1))[mind .== mg]
   tvcs = svcs[:,list]
   ovrlp = zeros(length(ns),size(tvcs,2))
   for i in 1:length(ns) 
      nd = 2*ns[i]+1
      for m in 1:md
         frst = ni[i,1] + jsd*(m-1)
         last = ni[i,2] + jsd*(m-1)
         ovrlp[i,:] += transpose(sum(tvcs[frst:last, :], dims=1))
      end
   end
   nind = zeros(Int,size(svcs,1))
   part = zeros(Int,length(list))
   for i in 1:length(ns)
      nd = (2*ns[i]+1)#*(vtmax+1)
      #count = min(nd*(vtmax+1),nd*(md))
      #vlimit = min(vtmax+3,md) 
      #for v in 0:vlimit
      perm = sort(sortperm(ovrlp[i,:], rev=true)[1:nd])#[1:nd]
      nind[perm] .= i
      ovrlp[:,perm] .= 0.0 
      #end
   end
   #println(nind)
   return nind
end
function nfinderv2(svcs,mind,md,mc,vtmax,jd,sd,ns,ni)
   nind = zeros(Int,size(svcs,1))
   for v in 0:vtmax
      mg = mc + vt2m(v) + 1
      nind = nfinderform(svcs,mind,mg,vtmax,md,jd,sd,ns,ni)
   end
   return nind
end
function nfinderv3(svcs,mind,md,mc,vtmax,jd,sd,ns,ni)
   jsd = jd*sd
   vlist = (mc+1) .+ vt2m.(collect(0:vtmax))
   list = collect(1:size(svcs,1))[Bool.(sum(mind .∈ vlist',dims=2))[:]]
   tvcs = svcs[:,list]
   ovrlp = zeros(length(ns),size(tvcs,2))
   for i in 1:length(ns) 
      nd = 2*ns[i]+1
      for m in 1:md
         frst = ni[i,1] + jsd*(m-1)
         last = ni[i,2] + jsd*(m-1)
         ovrlp[i,:] += transpose(sum(tvcs[frst:last, :], dims=1))
      end
   end
   nind = zeros(Int,size(svcs,1))
   part = zeros(Int,length(list))
   for i in 1:length(ns)
      nd = (2*ns[i]+1)*min((vtmax+1),md)
      #count = min(nd*(vtmax+1),nd*(md))
      #vlimit = min(vtmax+3,md) 
      #for v in 0:vlimit
      perm = sort(sortperm(ovrlp[i,:], rev=true)[1:nd])#[1:nd]
      nind[list[perm]] .= i
      ovrlp[:,perm] .= 0.0 
      #end
   end
   #println(nind)
   return nind
end

#kperm(n::Int)::Array{Int} = sortperm(Int.(cospi.(collect(-n:n).+isodd(n))) .* collect(-n:n))
#keperm(n::Int)::Array{Int} = sortperm(sortperm(collect(-n:n), by=abs))[kperm(n)]

function ramassign(vecs,j::Float64,s::Float64,mcalc::Int,vtmax)
   jd = Int(2.0*j) + 1
   sd = Int(2.0*s) + 1
   ns, nd, ni, jsd = srprep(j,s)
   #println(ns)
   #println(ni)
   md = 2*mcalc + 1 
   count = min(vtmax+4,md)
   svcs = abs.(vecs[:,1:jsd*count]).^2

   mind = mfinder(svcs,jsd,md,mcalc,vtmax)
   nind = nfinderv3(svcs,mind,md,mcalc,vtmax,jd,sd,ns,ni) 
   #nfinder(svcs,vtmax,md,jd,sd,ns,ni)
   #println(mind)
   col = collect(1:size(vecs,1))
   perm = zeros(Int,size(vecs,1)) #initalize big because easier
   for ng in 1:length(ns)
      nfilter = (nind .== ng)
      for v in 0:vtmax
         mg = mcalc + vt2m(v) + 1
         frst = jsd*(mg-1) + ni[ng,1]
         last = jsd*(mg-1) + ni[ng,2]
         part = col[nfilter .* (mind .== mg)]
         part = part[keperm(ns[ng])]
         perm[frst:last] = part#col[filter][keperm(ns[ng])]
      end
   end
   perm = perm[perm .!= 0]
   return perm
end
#Attempting to generalize ramassign to indefinite fold problem
#I can simply use a much larger mcalc value and essentially search
#for the combined m states

### JACOBI
#Reorganize matrix to better match my conception of the quantum numbers
function kperm(j,s)::Array{Int}
   perm = zeros(Int,Int((2*j+1)*(2*s+1)))
   shift = 0
   for n in Δlist(j,s)
      nd = 2*n+1
      perm[(1+shift):(nd+shift)] = kperm(n) .+ shift
      shift += nd
   end
   return perm
end
function kperm(j,s,shift::Int,jsd::Int,Δl::Array,perm)::Array{Int}
   #perm = zeros(Int,jsd)
   for n in Δl
      nd = 2*n+1
      perm[(1+shift):(nd+shift)] = kperm(n) .+ shift
      shift += nd
   end
   return perm#, shift
end
function kperm(j,s,m)
   jsd = Int((2*j+1)*(2*s+1))
   Δlst = Δlist(j,s)
   shift = 0
   perm = zeros(Int,jsd*(2*m+1))
   for i in 0:(2*m)
      perm = kperm(j,s,shift,jsd,Δlst,perm)
      shift += jsd
   end
   return perm
end

#This takes the eigenvectors & builds the permutation
assignperm(vec) = sortperm([iamax(vec[:,i]) for i in 1:size(vec,2)])

