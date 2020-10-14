using DataFrames, Stan, Mamba, Distributions

d = readtable("Progenesis.csv",skipstart=3)

d[:Index] = [i for i in 1:size(levels(d[:Unique_Id]))[1]]

#Replace zeros with NA
for j in 2:17
  for i in 1:size(d)[1]
    if d[i,j] == 0.0
      d[i,j] = NA
    end
  end
end

d = melt(d,[:Unique_Id,:Index],names(d)[2:17])

#Add columns for replicate and condition
replicate = Dict{Symbol,String}(
 :Sum_260416_Sample5 => "R1",
 :Sum_260416_Sample9 => "R2",
 :Sum_260416_Sample13 => "R3",
 :Sum_260416_Sample17 => "R4",
 :Sum_260416_Sample6 => "R1",
 :Sum_260416_sample10 => "R2",
 :Sum_260416_Sample14 => "R3",
 :Sum_260416_Sample18 => "R4",
 :Sum_260416_Sample7 => "R1",
 :Sum_260416_Sample11 => "R2",
 :Sum_260416_Sample15 => "R3",
 :Sum_260416_Sample19 => "R4",
 :Sum_260416_Sample8 => "R1",
 :Sum_260416_Sample12 => "R2",
 :Sum_260416_Sample16 => "R3",
 :Sum_260416_Sample20 => "R4"
)

condition = Dict{Symbol,String}(
 :Sum_260416_Sample5 => "Ctl0",
 :Sum_260416_Sample9 => "Ctl0",
 :Sum_260416_Sample13 => "Ctl0",
 :Sum_260416_Sample17 => "Ctl0",
 :Sum_260416_Sample6 => "Ctl10",
 :Sum_260416_sample10 => "Ctl10",
 :Sum_260416_Sample14 => "Ctl10",
 :Sum_260416_Sample18 => "Ctl10",
 :Sum_260416_Sample7 => "KO0",
 :Sum_260416_Sample11 => "KO0",
 :Sum_260416_Sample15 => "KO0",
 :Sum_260416_Sample19 => "KO0",
 :Sum_260416_Sample8 => "KO10",
 :Sum_260416_Sample12 => "KO10",
 :Sum_260416_Sample16 => "KO10",
 :Sum_260416_Sample20 => "KO10"
)

d[:Condition] = [condition[d[i,:variable]] for i in 1:size(d)[1]]
d[:Rep] = [replicate[d[i,:variable]] for i in 1:size(d)[1]]

#Write out csv for each comparison

writetable("CtlDat.csv",vcat(d[d[:Condition].=="Ctl0",:],d[d[:Condition].=="Ctl10",:]))

writetable("KODat.csv",vcat(d[d[:Condition].=="KO0",:],d[d[:Condition].=="KO10",:]))

writetable("0Dat.csv",vcat(d[d[:Condition].=="Ctl0",:],d[d[:Condition].=="KO0",:]))

writetable("10Dat.csv",vcat(d[d[:Condition].=="Ctl10",:],d[d[:Condition].=="KO10",:]))
