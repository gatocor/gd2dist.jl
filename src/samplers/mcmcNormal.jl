mutable struct MCMCNormal <: MCMCsample

    k::Array{Int}
    samples::Dict

end

function mcmcSampleNormalStep!( data, datac, k, w, wc, μ, μc, σ2, σ2c, steps = 1)
    
    #Declare auxiliar variables
    wAux = zeros(k[1])    
    wAuxc = zeros(k[1]*k[2])
    nAuxSum = zeros(k[1])
    nAuxSumc = zeros(k[2])
    nAux = zeros(k[1],k[2]+1)
    xAux = zeros(k[1],k[2]+1)
    x2Aux = zeros(k[1],k[2]+1)
    
    for i in 1:steps
        #Clean auxiliar
        nAux .= 0. 
        xAux .= 0.
        x2Aux .= 0.

        #SAMPLE IDENTITIES DATA
        for (j,d) in enumerate(data)
            #Weight normals
            for i in 1:k[1]
                wAux[i] = log(w[i])+lognorm(d,μ[i],σ2[i])
            end
            wMax = maximum(wAux)
            wAux = exp.(wAux.-wMax)
            wAux /= sum(wAux)
            #Sample
            mult = Multinomial(1,wAux)
            chosen = findfirst(rand(mult).==1)
            #Add to auxiliar
            nAux[chosen,end] += 1
            nAuxSum[chosen] += 1
            xAux[chosen,end] += d
            x2Aux[chosen,end] += (d-μ[chosen])^2
        end

        #SAMPLE IDENTITIES DATA CONVOLVED
        for (j,d) in enumerate(datac)
            #Weight normals
            for i in 1:k[1]
                for ic in 1:k[2]
                    wAuxc[k[1]*(ic-1)+i] = log(w[i])+log(wc[ic])+lognorm(d,μ[i]+μc[ic],σ2[i]+σ2c[ic])
                end
            end
            wMax = maximum(wAuxc)
            wAuxc = exp.(wAuxc.-wMax)
            wAuxc /= sum(wAuxc)
            #Sample
            mult = Multinomial(1,wAuxc)
            chosenc = findfirst(rand(mult).==1)
            chosen = mod(chosenc-1,k[1])+1
            chosenc = Int((chosenc-chosen)/k[1])+1
            #Add to auxiliars
            nAuxSum[chosen] += 1
            nAuxSumc[chosenc] += 1
            nAux[chosen,chosenc] += 1
            xAux[chosen,chosenc] += d
            x2Aux[chosen,chosenc] += (d-μ[chosen]-μc[chosenc])^2            
        end

        #SAMPLE WEIGHTS
        nAuxSum .+= 0.00000001
        dirich = Dirichlet(nAuxSum)
        w = rand(dirich)
        nAuxSumc .+= 0.00000001
        dirich = Dirichlet(nAuxSumc)
        wc = rand(dirich)

        #SAMPLE VARIANCES
        for i in 1:k[1]
            step = sum(nAux[i,:])
            if !(step ≈ 0.)
                σ2[i]=sliceSamplingLog(logσ2sum,σ2[i],[x2Aux[i,:],nAux[i,:],σ2c])#,step=1/sqrt(step))
            end
        end
        for i in 1:k[2]
            step = sum(nAux[:,i])
            if !(step ≈ 0.)
                step = sum(nAux[:,i])
                σ2c[i]=sliceSamplingLog(logσ2sum,σ2c[i],[x2Aux[:,i],nAux[:,i],σ2])#,step=1/sqrt(step))
            end
        end

        #SAMPLE MEANS
        effVar = 0
        effMean = 0
        for i in 1:k[1]
            effVar = 0.
            effMean = 0.
            for j in 1:k[2]
                effMean += (xAux[i,j]-nAux[i,j]*μc[j])/(σ2[i]+σ2c[j])
                effVar += nAux[i,j]/(σ2[i]+σ2c[j])
            end
            effMean += xAux[i,end]/(σ2[i])
            effVar += nAux[i,end]/(σ2[i])
            if !(effVar ≈ 0.)
                effMean /= effVar
                effVar = 1/effVar

                norm = Normal(effMean,sqrt(effVar))
                μ[i] = rand(norm)
            end
        end
        for i in 1:k[2]
            effVar = 0.
            effMean = 0.
            for j in 1:k[1]
                effMean += (xAux[i,j]-nAux[i,j]*μ[j])/(σ2c[i]+σ2[j])
                effVar += nAux[i,j]/(σ2c[i]+σ2[j])
            end
            if !(effVar ≈ 0.)
                effMean /= effVar
                effVar = 1/effVar

                norm = Normal(effMean,sqrt(effVar))
                μc[i] = rand(norm)
            end
        end
             
    end
        
    return
end

function mcmcSampleNormal( data, datac, k; chains::Int = 1, stopCriterion = :steps, steps::Int = 1000, ignore::Int=1000, gap::Int=1, basis = :Normal, initialCondition = Dict())

    #Declare parameters
    w = zeros(k[1])
    wc = zeros(k[2])
    μ = zeros(k[1])
    μc = zeros(k[2])
    σ = zeros(k[1])
    σc = zeros(k[2])
        
    #Declare storing dictionary
    out = Dict()
    out["w"] = zeros(steps*chains,k[1])
    out["μ"] = zeros(steps*chains,k[1])
    out["σ"] = zeros(steps*chains,k[1])
    out["wc"] = zeros(steps*chains,k[2])
    out["μc"] = zeros(steps*chains,k[2])
    out["σc"] = zeros(steps*chains,k[2])
    
    #Initialise model
    if initialCondition != Dict()
        w = initialCondition["w"]
        wc = initialCondition["wc"]
        μ = initialCondition["μ"]
        μc = initialCondition["μc"]
        σ2 = initialCondition["σ"]
        σ2c = initialCondition["σc"]    
    end

    pos = 0:steps:chains*steps

    #Sample
    Threads.@threads for chainId in 1:chains
                
        #Make ignorable steps
        mcmcSampleNormalStep!( data, datac, k, w, wc, μ, μc, σ2, σ2c, ignore)
        
        #Make steps in gap
        for i in 1:steps
            mcmcSampleNormalStep!( data, datac, k, w, wc, μ, μc, σ2, σ2c, gap)
            #Save
            out["w"][pos[chainId]+i,:] .= w
            out["μ"][pos[chainId]+i,:] = μ
            out["σ"][pos[chainId]+i,:] = σ2
            out["wc"][pos[chainId]+i,:] = wc
            out["μc"][pos[chainId]+i,:] = μc
            out["σc"][pos[chainId]+i,:] = σ2c           
        end 
    end
    
    return out

end
