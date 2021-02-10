mutable struct MCMCGamma <: MCMCsample

    k::Array{Int}
    samples::Dict

end


function mcmcSampleGamma( data, datac, k; chains::Int = 1, stopCriterion = :steps, steps::Int = 1000, ignore::Int=1000, gap::Int=1, basis = :Normal, initialCondition = Dict())

    #Declare parameters
    w = zeros(k[1])
    wc = zeros(k[2])
    θ = zeros(k[1])
    θc = zeros(k[2])
    κ = zeros(k[1])
    κc = zeros(k[2])
        
    #Declare storing dictionary
    out = Dict()
    out["w"] = zeros(steps*chains,k[1])
    out["θ"] = zeros(steps*chains,k[1])
    out["κ"] = zeros(steps*chains,k[1])
    out["wc"] = zeros(steps*chains,k[2])
    out["θc"] = zeros(steps*chains,k[2])
    out["κc"] = zeros(steps*chains,k[2])
    
    #Initialise model
    if initialCondition != Dict()
        w = initialCondition["w"]
        wc = initialCondition["wc"]
        θ = initialCondition["θ"]
        θc = initialCondition["θc"]
        κ = initialCondition["κ"]
        κc = initialCondition["κc"]    
    end

    pos = 0:steps:chains*steps

    #Sample
    Threads.@threads for chainId in 1:chains
                
        #Make ignorable steps
        mcmcSampleGammaStep!( data, datac, k, w, wc, θ, θc, κ, κc, ignore)
        
        #Make steps in gap
        for i in 1:steps
            mcmcSampleGammaStep!( data, datac, k, w, wc, θ, θc, κ, κc, gap)
            #Save
            out["w"][pos[chainId]+i,:] .= w
            out["θ"][pos[chainId]+i,:] = θ
            out["κ"][pos[chainId]+i,:] = κ
            out["wc"][pos[chainId]+i,:] = wc
            out["θc"][pos[chainId]+i,:] = θc
            out["κc"][pos[chainId]+i,:] = κc           
        end 
    end
    
    return out

end
