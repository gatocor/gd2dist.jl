function sliceSamplingLog(f,x0,args;it=10,step=0.1)
    
    xOld = x0
    dist = Uniform(0,1)
        
    for i in 1:it
    
        xTry = xOld
        xMin = xOld
        xMax = xOld

        p = f(xOld,args...) + log(rand(dist))
        #Expand
        xMin = xOld - step
        while f(xMin,args...) > p
            xMin -= step
        end
        counter = 0
        xMax = xOld + step
        while f(xMax,args...) > p && counter < 2000
            xMax += step
            counter += 1
        end
        #Sample
        xTry = (xMax-xMin)*rand(dist)+xMin
        while f(xTry,args...) < p
            xTry = (xMax-xMin)*rand(dist)+xMin
            if xTry < xOld
                xMin = xTry
            else
               xMax = xTry 
            end
        end
        
        xOld = xTry
    end
    
    return xOld
end