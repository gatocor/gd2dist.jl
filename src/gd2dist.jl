module gd2dist

    using Random
    using Distributions
    import SpecialFunctions

    export norm, lognorm, logσ2, logσ2sum, gamma, loggamma, loggammaθ, loggammaκ, loggammaθsum, loggammaκsum
    export sliceSamplingLog
    export mcmcSample

    include("./distributions.jl")
    include("./sliceSampling.jl")
    include("./MCMCsampler.jl")
    include("./samplers/mcmcNormal.jl")
    include("./samplers/mcmcGamma.jl")

end