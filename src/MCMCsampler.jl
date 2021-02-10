abstract type MCMCsample end

function mcmcSample( data, datac, k; stopCriterion = :steps, steps = 1000, basis = :Normal, initialCondition = Dict())

    if basis == :Normal
        return mcmcSampleNormal( data, datac, k, stopCriterion = stopCriterion, steps = steps, basis = :Normal, initialCondition = initialCondition)
    elseif basis == :Gamma
        return mcmcSampleGamma( data, datac, k, stopCriterion = stopCriterion, steps = steps, basis = :Normal, initialCondition = initialCondition)
    end

end
