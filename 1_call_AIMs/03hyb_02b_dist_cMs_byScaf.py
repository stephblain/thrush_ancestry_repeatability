
#!/usr/bin/env python3


############
#estimates the distance in cM between bp based on a recomb map
#inputs:
    #Ne = effective population size
    #mapFile is a recombination map with window size of 100,000 and columns named  
    ##"mean_rho" (mean rho for window in 4Ner/kb), "midpos" (midpoint of window, starting from 50000)
    ##"chr" (chromosome)
############

print("start py script")

Ne1=570338
#mapFile="C:/Users/Steph/OneDrive - Texas A&M University/Thrushes/generate_AIMs/data/rhomap_res_win_100kb.csv"
#aimsFile="C:/Users/Steph/OneDrive - Texas A&M University/Thrushes/generate_AIMs/data/current_aims_test.txt"
mapFile="/scratch/user/sblain/thrush_hybrids/rhomap_res_win_100kb_noUnderscores.csv"
aimsFile="/scratch/user/sblain/get_AIMs/hmm_AIMs/lowCov/get_f50_cMs/current_aims_f50"

#set range from chrom list to run through
j1=0
j2=200

import pandas as pd

recomb=pd.read_csv(mapFile)
aims=pd.read_csv(aimsFile,header=None,sep="\t")
aims.columns=["chr","pos","Ref","Alt"]

#function to calculate Morgans per bp given Ne and rho
#assumes rho is equal to 4Ner/kb
def calc_M_per_bp(Ne,rho):
  return rho/(4*Ne*1000)

recomb["startpos"]=recomb["midpos"]-50000
aims["cM"]=float()

print("subset aims by recomb")

aims=aims[aims["chr"].isin(recomb.chr.unique())]
print("outputting scaffolds:",aims.chr.unique()[j1:j2])

for chr1 in aims.chr.unique()[j1:j2]:
    
    aimsChr=aims[aims['chr']==chr1]
    recombChr=recomb[recomb['chr']==chr1]
    aimsChr=aimsChr.reset_index(drop=True)   
    print("start scaffold:",chr1)
 
    for i in range(1,len(aimsChr)): #start at 1 to skip first position
        
        pos1=aimsChr.iloc[i]["pos"] #get position
        pos2=round(pos1,-5) #round to nearest 100000 to get correct start position
                
        if pos1 < pos2: #if rounded up, subtract 100,000
            pos2=pos2-100000       
        
        #if position is past maximum in recomb map, take the last estimate
        if pos2 > max(recombChr["startpos"]):
            focalRecomb=recombChr.iloc[[len(recombChr)-1]]
        else:    
            # print("chr:",chr1,"position:",pos1,pos2)
            #find in recombination map
            focalRecomb=recombChr[recombChr["startpos"]==pos2]
            while len(focalRecomb)==0: #if recomb rate is missing, move to next chunk
                pos2=pos2+100000
                focalRecomb=recombChr[recombChr["startpos"]==pos2]
                print("used recombination rate from",chr1,pos2,"for",chr1,pos1)

        dist_bp=aimsChr.iloc[i]["pos"]-aimsChr.iloc[i-1]["pos"] #get distance in bp from previous site
        recomb_M_bp=calc_M_per_bp(Ne1,focalRecomb["mean_rho"]) #get recomb rate in morgans per bp
        # print("recomb rate:",recomb_M_bp)
        dist_cM=dist_bp*recomb_M_bp*100 #get distance in cM
        dist_cM=dist_cM.to_frame(name="col0").reset_index()
        # print("distance:",dist_cM)
        aimsChr.loc[i,"cM"]=dist_cM.loc[0,"col0"]
        # print("end pos:",pos1)
    
    outFile=aimsFile+"_"+chr1+".cMs"
    aimsChr.to_csv(outFile,index=False,header=False,sep="\t",float_format='%.15f')

print("outputted scaffolds:",aims.chr.unique()[j1:j2])
