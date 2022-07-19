
function build_jcbn!(jcbn,inds,vecs,params,perm)
"""
This builds the Jacobian based on the Hellmann–Feynman theorem.
"""
   Threads.@threads for a in 1:size(inds)[1]
      ju = 0.5*inds[a,1]
      jl = 0.5*inds[a,4]
      σu = inds[a,2]
      σl = inds[a,5]
      vecu = vecs[1:Int((2*S+1)*(2*ju+1)*(2*mcalc+1)),inds[a,3],σu+1]
      vecl = vecs[1:Int((2*S+1)*(2*jl+1)*(2*mcalc+1)),inds[a,6],σl+1]
      for i in 1:length(perm)
         b = perm[i]
         jcbn[a,i] = anaderiv(jl,S,σl,vecl,params,b) - anaderiv(ju,S,σu,vecu,params,b)
      end
   end
   return jcbn
end
function build_hess!(hssn,dk,jcbn,β)
   hssn = transpose(jcbn)*jcbn
   for i in size(hssn)[1]
      dk[i,i] = norm(hssn[:,i])
   end
   return hssn, dk
end

function lbmq_acc(g,jcb,dk,β,omc)
   println(size(g))
   println(size(transpose(jcb)))
   println(size(dk))
   a = (dk^2)*omc
   a = -0.5*inv(g)*transpose(jcb)
   check = 2*norm(a)/norm(β)
   if check > 0.75
      a = zeros(Float64,size(a))
   end
   return a
end
function lbmq_step!(β,jcbn,dk, weights, omc, λ)
   jtw = transpose(jcbn)*weights
   jtj = jtw*jcbn
   A = jtj + λ*Diagonal(jtj)#transpose(dk)*dk
   B = factorize(Symmetric(A))
   X = jtw*omc
   β = ldiv!(β, B, -X)
   return β,X
end

function lbmq_opttr(nlist,ofreqs,uncs,inds,params,scales,λ)
   vals,vecs = limeigcalc(nlist, inds, params)
   rms, omc = rmscalc(vals, inds, ofreqs)
   perm,n = findnz(sparse(scales))
   println("Initial RMS = $rms")
   goal = sum(uncs)/length(uncs)
   newparams = copy(params)
   W = diagm(0=>(uncs .^ -1))
   converged = false
   THRESHOLD = 1.0E-8
   RHOTHRES = -1.0E-6
   ϵ0 = 0.1E-2
   LIMIT = 50
   λlm = 0.0
   Δlm = 100.0
   counter = 0
   rms, omc = rmscalc(vals,inds,ofreqs)
   nparams = copy(params)
   β = zeros(Float64,size(perm))
   J = zeros(Float64,size(inds)[1],length(perm))
   A = zeros(Float64,size(J))
   H = zeros(Float64,length(perm),length(perm))
   D = zeros(Float64,size(H))
   J = build_jcbn!(J,inds,vecs,params,perm)
   H, D = build_hess!(H,D,J,β)
   converged=false
   while converged==false
      β,g = lbmq_step!(β,J,D,W,omc,λlm)
      #β += lbmq_acc(g,J,D,β,omc)
      if norm(D*β)>Δlm
         β *= Δlm/norm(D*β)
      end
      nparams[perm] = params[perm] + β
      vals, nvecs = limeigcalc(nlist, inds, nparams)
      nrms, nomc = rmscalc(vals,inds,ofreqs)
      check = abs(nrms-rms)/rms
      if nrms < rms
         rms = nrms
         omc = nomc
         params = nparams
         vecs = nvecs
         J = build_jcbn!(J,inds,vecs,params,perm)
         H, D = build_hess!(H,D,J,β)
         λlm *= 1.0/3.0
      else
         λlm *= 2.0
         λlm = max(λlm,1.0E-13)
      end
      if rms < ϵ0
         converged = true
      end
      ρlm = lbmq_gain(β,λlm,g,rms,nrms)
      if ρlm ≥ 0.75
         Δlm *= 2.0
      else
         Δlm *= 0.8
      end
      srms = (@sprintf("%0.4f", rms))
      slλ = (@sprintf("%0.4f", log10(λlm)))
      sΔ = (@sprintf("%0.6f", Δlm))
      counter += 1
      scounter = lpad(counter,3)
      println("After $scounter interations, RMS = $srms, log₁₀(λ) = $slλ, Δₖ = $sΔ")
      if (check < THRESHOLD)||(rms ≤ goal)#&&(counter > 1)
         println("A miracle has come to pass. The fit has converged")
         break
      elseif counter ≥ LIMIT
         println("Alas, the iteration count has exceeded the limit")
         println(omc)
         break
      else
         #write update to file
      end #check if
   end#while
   return params, vals
end
