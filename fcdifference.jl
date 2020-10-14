using FileIO, JLD, DataFrames

function calculateFDRUpDown(samples,fcup=1.05,fcdown=0.95)
  up = 0
  down = 0
  for i in 1:size(samples)[1]
      samples[i]<log(fcdown) ? down += 1 : down = down
      samples[i]>log(fcup) ? up += 1 : up = up
  end

  upfdr = 1 - up/(size(samples)[1])
  downfdr = 1 - down/(size(samples)[1])

  upfdr,downfdr
end

function calculateFDRSame(samples,fcup=1.5,fcdown=0.666666)
  same = 0
  for i in 1:size(samples)[1]
      samples[i]<log(fcup) && samples[i]>log(fcdown) ? same += 1 : same = same
  end

  samefdr = 1 - same/(size(samples)[1])

  samefdr
end


ctlsameresults = readtable("Ctl-Same.csv")
ctlsameresults[:ctlmeanlogfc] = ctlsameresults[:meanlogfc]
ctlsameresults[:ctltotalmiss] = [ctlsameresults[i,:nmissC1] + ctlsameresults[i,:nmissC2] for i in 1:size(ctlsameresults)[1]]

kosameresults = readtable("KO-Same.csv")
kosameresults[:komeanlogfc] = kosameresults[:meanlogfc]
kosameresults[:kototalmiss] = [kosameresults[i,:nmissC1] + kosameresults[i,:nmissC2] for i in 1:size(kosameresults)[1]]

combinedResults = join(ctlsameresults[:,[:Index,:ctlmeanlogfc,:ctltotalmiss]],kosameresults[:,[:Index,:komeanlogfc,:kototalmiss]],on=[:Index])

datoriginal = readtable("Progenesis.csv",skipstart=3)[:,1:2]
datoriginal[:Index] = [i for i in 1:size(datoriginal)[1]]

#origCols = [names(datoriginal)[1:11];names(datoriginal)[28]]
fcdiffresults = join(datoriginal[:,[1;3]],combinedResults, on = [:Index])
fcdiffFDR = Array{Float64,1}(size(fcdiffresults)[1])
fcdiffDirection = Array{String,1}(size(fcdiffresults)[1])
fcdiffmeans = Array{Float64,1}(size(fcdiffresults)[1])
fcdifflq = Array{Float64,1}(size(fcdiffresults)[1])
fcdiffuq = Array{Float64,1}(size(fcdiffresults)[1])

#Read MCMC samples to calculate differences in log-fold-change
ctlsamps = jldopen("samples/Ctl-samps.jld","r")
kosamps = jldopen("samples/KO-samps.jld","r")
for i in 1:size(fcdiffresults)[1]
  s = fcdiffresults[i,:Index]

  try
    ctlfc = read(ctlsamps,"sim-"*string(s))
    kofc = read(kosamps,"sim-"*string(s))

    nSamples = 10000

    fcdiff = Array{Float64,1}(nSamples)

    for j in 1:nSamples
      ctlind = rand(1:length(ctlfc))
      koind = rand(1:length(kofc))

      fcdiff[j] = kofc[koind] - ctlfc[ctlind]
    end

    up,down = calculateFDRUpDown(fcdiff,1.25,1/1.25)

    fcdifflq[i] = quantile(fcdiff,0.05)
    fcdiffmeans[i] = mean(fcdiff)
    fcdiffuq[i] = quantile(fcdiff,0.95)

    fcdiffDirection[i] = (up<down ? "UP" : "DOWN")
    fcdiffFDR[i] = minimum([up;down])
    println(fcdiffDirection[i]*": " * string(fcdiffFDR[i]))

  catch e
    fcdifflq[i] = -1e9
    fcdiffmeans[i] = 0.0
    fcdiffuq[i] = 1e9

    fcdiffDirection[i] = "NoData"
    fcdiffFDR[i] = 1.0
    println(fcdiffDirection[i]*": " * string(fcdiffFDR[i]))


  end
end

close(ctlsamps)
close(kosamps)

fcdiffresults[:Q5Difference] = fcdifflq
fcdiffresults[:MeanDifference] = fcdiffmeans
fcdiffresults[:Q95Difference] = fcdiffuq
fcdiffresults[:Direction] = fcdiffDirection
fcdiffresults[:localFDR] = fcdiffFDR

#Calculate cumulative mean localFDR/PEP
sort!(fcdiffresults,cols=[:localFDR])
cumSumlocalFDR = cumsum(fcdiffresults[:localFDR])
fcdiffresults[:globalFDR] = [cumSumlocalFDR[i]/i for i in 1:size(fcdiffresults)[1]]

writetable("fcdifference-1.25.csv",fcdiffresults)
