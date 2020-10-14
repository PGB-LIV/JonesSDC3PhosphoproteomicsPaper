using DataFrames, Stan, Mamba, Distributions, FileIO, JLD

#Reshape MCMC samples from multiple chains into one array
function unfold(sim,variable)
  a = Array{Float64,1}(size(sim.value)[1] * size(sim.value)[3])

  for k in 1:size(sim.value)[3]
    a[1+(k-1)*size(sim.value)[1]:k*size(sim.value)[1]] = sim[:,variable,k].value
  end

  return a
end

function unfold(sim)
  a = Array{Float64,1}(size(sim.value)[1] * size(sim.value)[3])

  for k in 1:size(sim.value)[3]
    a[1+(k-1)*size(sim.value)[1]:k*size(sim.value)[1]] = sim[:,sim.names[1],k].value
  end

  return a
end

#Calculate quantitative localFDR/PEP from MCMC samples
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

#Calculate quantitative localFDR/PEP from MCMC samples
function calculateFDRSame(samples,fcup=1.5,fcdown=0.666666)
  same = 0
  for i in 1:size(samples)[1]
      samples[i]<log(fcup) && samples[i]>log(fcdown) ? same += 1 : same = same
  end

  samefdr = 1 - same/(size(samples)[1])

  samefdr
end


#Run model for each comparison
for comparison in ["Ten", "Zero", "Ctl", "KO"]

  dat = readtable(comparison*"Dat.csv")

  conditions = levels(dat[:Condition])
  nConditions = length(conditions)

  dat = dat[:,2:size(dat)[2]]

  indices = levels(dat[:Index])


  newdf = unique(dat[dat[:Index] .== indices[1],names(dat)[2:3]])
  newdf = vcat(newdf,newdf)
  newdf = hcat(newdf,unstack(dat[dat[:Index] .== indices[1],:], :Condition,:Rep,:value))

  for i in 2:length(indices)
    temp = unique(dat[dat[:Index] .== indices[i],names(dat)[2:3]])
    temp = vcat(temp,temp)
    temp = hcat(temp,unstack(dat[dat[:Index].==indices[i],:],:Condition,:Rep,:value))

    newdf = vcat(newdf,temp)
  end

  nmiss = Array{Int64,1}(size(newdf)[1])
  for i in 1:size(newdf)[1]
    nmiss[i] = 0
    for j in size(newdf)[2]-3:size(newdf)[2]
      if isna(newdf[i,names(newdf)[j]])
        nmiss[i] += 1
      end
    end
  end

  newdf[:nmiss] = nmiss

  meanabundance = Array{Float64,1}(size(newdf)[1])
  for i in 1:size(newdf)[1]
    total = 0
    for j in 4:7
      if !isna(newdf[i,names(newdf)[j]])
        total += newdf[i,names(newdf)[j]]
      end
    end
    meanabundance[i] = total/(4-newdf[i,:nmiss])
  end

  newdf[:mean] = meanabundance

  #Fitting lognormals to mean abundances for each condition/nmiss combination
  priors = Array{Distributions.Normal,2}(5,2)
  for i in 0:3
    for j in 1:2
      t = newdf[newdf[:nmiss].==i,:]
      t = t[t[:Condition].==conditions[j],:]

      priors[i+1,j] = Distributions.fit(Distributions.Normal,log(t[:mean]))
    end
  end

  means = Array{Float64,2}(4,2)
  for i in 1:4
    for j in 1:2
      means[i,j] = priors[i,j].μ
    end
  end

  sds = Array{Float64,2}(4,2)
  for i in 1:4
    for j in 1:2
      sds[i,j] = priors[i,j].σ
    end
  end

  #linear regression on means of lognormals to infer prior for n = 4 missing
  means = reshape(means,(8,1))
  interceptSelect = [1 0; 1 0; 1 0; 1 0; 0 1; 0 1; 0 1; 0 1]
  conditionMatrix = [4 0; 3 0; 2 0; 1 0; 0 4; 0 3; 0 2; 0 1]

  linregData = [Dict{String,Any}(
    "y" => means[:,1],
    "N" => length(means),
    "C" => nConditions,
    "X" => interceptSelect,
    "conditionMatrix" => conditionMatrix
  )]

  #Stan model for linear regression
  linregmodel = """
    data {
      int<lower=1> N; //num data points
      vector<lower=0>[N] y;
      int<lower=2> C;
      matrix[N,C] X;
      matrix[N,C] conditionMatrix;
    }

    parameters {
      vector<lower=0>[C] intercept;
      vector[C] beta;
      vector<lower=0>[C] sigma;
    }

    transformed parameters {
    }

    model {
      intercept ~ normal(0,10);
      beta ~ normal(0,10);
      sigma ~ student_t(3,0,5);
      y ~ normal(X*intercept + conditionMatrix*beta, X*sigma);
    }

    generated quantities {
    }
  """

  initialIters = 10000
  nChains = 4
  warmup = 0.5
  stanmodel = Stan.Stanmodel(Stan.Sample(save_warmup = false),
                    name=comparison*"linregmodel",
                    nchains=nChains,
                    model=linregmodel,
                    adapt = round(Int, warmup*initialIters),
                    update = round(Int, initialIters*(1-warmup)),
                    thin = 1)

  sim = Stan.stan(stanmodel, linregData, pwd(), CmdStanDir=Stan.CMDSTAN_HOME)

  #Use maximum sd across all subsets to be conservative
  priors[5,1] = Normal(mean(unfold(sim[:,"intercept.1",:])),maximum(sds[:,1]))
  priors[5,2] = Normal(mean(unfold(sim[:,"intercept.2",:])),maximum(sds[:,2]))

  #Stan model for Poisson regression
  missmodel = """
    data{
      int<lower=0> N_obs;
      vector[N_obs] y_obs;
      vector [N_obs] X_obs;
      real priorMu;
      real priorSigma;
    }
    parameters{
      real<lower=0> peptideSigma;
      real intercept;
      real conditionEffect;
      vector[N_obs] lambda_obs;
    }
    transformed parameters{
    }
    model{
      intercept ~ normal(priorMu,priorSigma);
      conditionEffect ~ normal(0,50);
      peptideSigma ~ student_t(3,0,5);

      lambda_obs ~ normal(intercept+X_obs*conditionEffect,peptideSigma);

      //y_obs ~ poisson_log(lambda_obs);
      for (i in 1:N_obs){
        target += lmultiply(y_obs[i], exp(lambda_obs[i])) - exp(lambda_obs[i]);
      }
    }
  """

  #Stan model for Poisson regression (for case where there is only single observation)
  missmodelsinglepoint = """
    data{
      int<lower=0> N_obs;
      real y_obs;
      real X_obs;
      real priorMu;
      real<lower=0> priorSigma;
    }
    parameters{
      //vector[N_obs] pepR;
      real<lower=0> peptideSigma;
      real intercept;
      real conditionEffect;
      real lambda_obs;
    }
    transformed parameters{

    }
    model{
      intercept ~ normal(priorMu,priorSigma);
      conditionEffect ~ normal(0,50);
      peptideSigma ~ student_t(3,0,5);

      lambda_obs ~ normal(intercept+X_obs*conditionEffect,peptideSigma);

      //y_obs ~ poisson_log(lambda_obs);
      target += lmultiply(y_obs, exp(lambda_obs)) - exp(lambda_obs);
    }
  """

  #Create tables for results of up/down tests and same test
  udresults = unique(dat[:,[:Index,:Unique_Id]])
  udtestDirections = Array{String,1}(length(indices))
  udlocalFDRs = Array{Float64,1}(length(indices))

  sameresults = unique(dat[:,[:Index,:Unique_Id]])
  samelocalFDRs = Array{Float64,1}(length(indices))

  meanlogfc = Array{Float64,1}(length(indices))
  absmeanlogfc = Array{Float64,1}(length(indices))
  uq = Array{Float64,1}(length(indices))
  lq = Array{Float64,1}(length(indices))

  missC1 = Array{Int64,1}(length(indices))
  missC2 = Array{Int64,1}(length(indices))

  for s in 1:length(indices)
    d = dat[dat[:Index].==s,:]
    tmp = newdf[newdf[:Index].==s,:]
    t=tmp[tmp[:Condition].==conditions[1],:]
    missC1[s] = t[:nmiss][1]
    t=tmp[tmp[:Condition].==conditions[2],:]
    missC2[s] = t[:nmiss][1]
  end

  totalMiss = [missC1[i] + missC2[i] for i in 1:length(indices)]

  Rhat = Array{Float64,1}(length(indices))

  Eff = Array{Float64,1}(length(indices))

  #Open JLD file to save MCMC samples
  samps = jldopen("samples/"*comparison*"-samps.jld","w")

  #Run model for each peptide
  @time for s in 1:length(indices)
    d = dat[dat[:Index].==s,:]

    println("INDEX: " * string(s) * "   missC1: " * string(missC1[s]) * "   missC2: " * string(missC2[s]))

    dObs = d[!isna(d[:value]),:]
    dMis = d[isna(d[:value]),:]

    N_obs = size(dObs)[1]
    if N_obs < 1
      #No observed data
      udtestDirections[s] = "NoData"
      udlocalFDRs[s] = 1.0
      samelocalFDRs[s] = 1.0

      meanlogfc[s] = 0.0
      absmeanlogfc[s] = 0.0
      uq[s] = 1e9
      lq[s] = -1e9

      Rhat[s] = 1e9
      Eff[s] = 0

    else
      y_obs = [dObs[i,:value]+0.001 for i in 1:N_obs]

      j = (missC2[s] == 4 ? 1 : 2)

      priorMu = 0
      priorSigma = 100

      if missC1[s] == 4
        priorMu = priors[5,1].μ
        priorSigma = priors[5,1].σ
      elseif missC2[s] == 4
        priorMu = priors[5,2].μ
        priorSigma = priors[5,2].σ
      end

      X_obs = Array{UInt8,}(N_obs)
      for i in 1:N_obs
        X_obs[i] = (dObs[i,:Condition] == conditions[j] ? 1 : 0)
      end


      standata = [Dict{String,Any}(
        "N_obs" => N_obs,
        "y_obs" => y_obs,
        "X_obs" => X_obs,
        "priorMu" => priorMu,
        "priorSigma" => priorSigma
      )]

      initialIters = 10000
      numChains = 8
      warmup = 0.5
      thinning =  1

      stanmodel = Stan.Stanmodel(Stan.Sample(save_warmup = false),
                        name=(N_obs > 1 ? comparison*"missmodel" : comparison*"missmodelsinglepoint"),
                        nchains=numChains,
                        model=(N_obs > 1 ? missmodel : missmodelsinglepoint),
                        adapt = round(Int, warmup*initialIters),
                        update = round(Int, initialIters*(1-warmup)),
                        thin = thinning)

      sim = Stan.stan(stanmodel, standata, pwd(), CmdStanDir=Stan.CMDSTAN_HOME)

      Rhat[s] = gelmandiag(sim[:,"conditionEffect",:]).value[1,2,1]
      Eff[s] = summarystats(sim[:,"conditionEffect",:]).value[1,5,1]

      println("Rhat: " * string(Rhat[s]) * "   ESS: "*string(Eff[s]))

      fc = unfold(sim[:,"conditionEffect",:])

      if j == 1
        fc = -fc
      end
    
    #Write MCMC samples to file
      write(samps,"sim-"*string(s),fc)

    end
  end

  #Write Rhat diagnostics and effective sample sizes to file
  write(samps,"Rhat",Rhat)
  write(samps,"Eff",Eff)

  close(samps)

  #Open file with samples and calculate summary stats
  f = jldopen("samples/"*comparison*"-samps.jld","r")
  for s in 1:length(indices)
    if missC1[s]+missC2[s] < 8
      fc = read(f,"sim-"*string(s))
    
    #Calculate localFDR/PEP
      up,down = calculateFDRUpDown(fc,1.5,1/1.5)
      same = calculateFDRSame(fc,1.5,1/1.5)

      #println("UP: " * string(up) * "DOWN: " * string(down) * "SAME: " * string(same))

      udtestDirections[s] = (up < down ? "UP" : "DOWN")
      udlocalFDRs[s] = minimum([up,down])
      samelocalFDRs[s] = same

      meanlogfc[s] = mean(fc)
      absmeanlogfc[s] = abs(meanlogfc[s])
      uq[s] = quantile(fc,0.95)
      lq[s] = quantile(fc,0.05)
    else
      udtestDirections[s] = "NoData"
      udlocalFDRs[s] = 1.0
      samelocalFDRs[s] = 1.0

      meanlogfc[s] = 0.0
      absmeanlogfc[s] = 0.0
      uq[s] = 1e9
      lq[s] = -1e9

      Rhat[s] = 1e9
      Eff[s] = 0
    end
  end
  close(f)

  #Add columns to tables
  udresults[:Q5logfc] = lq
  udresults[:meanlogfc] = meanlogfc
  udresults[:Q95logfc] = uq
  udresults[:absmeanlogfc] = absmeanlogfc
  udresults[:nmissC1] = missC1
  udresults[:nmissC2] = missC2

  udresults[:Rhat] = Rhat
  udresults[:ESS] = Eff

  sameresults[:Q5logfc] = lq
  sameresults[:meanlogfc] = meanlogfc
  sameresults[:Q95logfc] = uq
  sameresults[:absmeanlogfc] = absmeanlogfc
  sameresults[:nmissC1] = missC1
  sameresults[:nmissC2] = missC2

  sameresults[:Rhat] = Rhat
  sameresults[:ESS] = Eff

  udresults[:Test] = udtestDirections
  udresults[:localFDR] = udlocalFDRs

  sameresults[:Test] = ["SAME" for i in 1:length(indices)]
  sameresults[:localFDR] = samelocalFDRs

  datoriginal = readtable("Progenesis.csv",skipstart=3)[:,[1;6;7;8;9;14;15;16;17]]

  datoriginal[:Index] = [i for i in 1:size(datoriginal)[1]]

  udresults2 = join(datoriginal,udresults,on=[:Index])
  sameresults2 = join(datoriginal,sameresults,on=[:Index])

  for i in 1:size(sameresults2)[1]
    if sameresults2[i,:nmissC1] == 4 && sameresults2[i,:nmissC2] == 4
      sameresults2[i,:Test] = "NoData"
    end
  end

  udresults2 = udresults2[udresults2[:Test].!="NoData",:]
  sameresults2 = sameresults2[sameresults2[:Test].!="NoData",:]

  #Calculate cumulative mean localFDR/PEP
  udresults3 = sort(udresults2,cols=[:localFDR,:absmeanlogfc],rev=[false,true])
  c = cumsum(udresults3[:localFDR])
  globalFDRs = [c[i]/i for i in 1:length(c)]
  udresults3[:globalFDR] = globalFDRs

  sameresults3 = sort(sameresults2,cols=[:localFDR,:absmeanlogfc],rev=[false,false])
  c = cumsum(sameresults3[:localFDR])
  globalFDRs = [c[i]/i for i in 1:length(c)]
  sameresults3[:globalFDR] = globalFDRs

  delete!(udresults3,:Unique_Id_1)
  delete!(sameresults3,:Unique_Id_1)

  delete!(udresults3,[:Rhat,:ESS])
  delete!(sameresults3,[:Rhat,:ESS])

  #Write output tables
  writetable(comparison*"-UpDown.csv",udresults3)
  writetable(comparison*"-Same.csv",sameresults3)

end