# Use it : give in argument one txt file with the accession ID / or a csv file + path where put output + analyse name(without space or special charactère)

#BALISE # PATHDEF where a define path
#==> see them when ai write this script to be used by nextflow or on cluster


# ---- Import

import os
import sys
import subprocess
import csv
import statistics
from math import log10

# ---- Set path

#inputFile=sys.argv[1] #Used at the final version
inputdir="/home/smile/Bureau/TestautoPip/ManonData-PRJNA371" # PATHDEF
#outdir=sys.argv[2] #Used at the final version
outdir="/home/smile/Bureau/TestautoPip/ManonData-PRJNA371"

#To suppres to didn't keep file generated that user don't want
#subprocess.run(["pwd"], output=True) #Todo
#currentdir= 'Command permit to have pwd result' #TODO
currentdir='/home/smile/Bureau/autoPip/'
tmpdir=outdir+'/working'
#TODO suppress tmpdir content (generated in the previous run) VOIR COMMENT ENV INTERACTIF DE CODIUM GERE

# ---- Other variables

#Construct files path here

#Calcul percentage value : RULES THAT WE SET

#.. Minimal quality required for ech base on the reads
#Set the maximal error percentage authorized. 45 bp or less : 0.39811 (Q4), 200 bp or more : 0.1 (Q10)
minLim=45
maxLim=200
minFpercent=0.39811
maxFpercent=0.1
limDifNt=maxLim-minLim
limDifPercent=minFpercent-maxFpercent
if limDifNt<=0 or limDifPercent<=0 :
    print('\n&&&&&&&&&& Warning calculation error 1\n') #TODO delete after debug
FpercentByNt=limDifPercent/limDifNt

#.. Percentage of sequence length we want to trim according to reads length
LminLim=50
LmaxLim=150
minPercentL=0.4
maxPercentL=0.6
limDifNtL=LmaxLim-LminLim
limDifPercentL=maxPercentL-minPercentL
if limDifNtL<=0 or limDifPercentL<=0 :
    print('\n&&&&&&&&&& Warning calculation error 2\n') #TODO delete after debug
LpercentByNt=limDifPercentL/limDifNtL

############################
# FOUNCTIONS
############################


############################
# CHECK USER INPUT
############################

# arg number
# arg file : fomat, not empty, readable with ',' si csv, exist
# arg repertory : exist
# arg name : not exist output directory on output path with this name, if it's case asking to user if he want erase this output analyse
# Sample to analyse : more than one (else bad list reading)

############################
# Create symbolic links with Acc ID file + analyse if there are missing datas
############################

#TEMP
multiqcPath=inputdir+"/multiqc/multiqc_data/multiqc_fastqc.txt" # PATHDEF
sampleList=inputdir+'/Acc.txt' #TODO make the file with alla accession, put the good path here and define givedAcc wi
file=open(sampleList,"r")
givedAcc=(file.read()).split("\n") # PATHDEF
file.close()
#givedAcc.sort()
totAccNumber=len(givedAcc)

#print("\n\n ----- SAMPLE ANALYZED --------\n",givedAcc[0]) #TEMP

############################
# Perform the first QC
############################

############################
# Analysing the multiQC
############################

# ---- Paired=true ?

# LINK LIST ON LINKED_FILES REPERTORY
#TODO suppr when work : now files in fastq.zip directory
#TODO
dataFiles=['SRR5229900.fastqc.gz']

# Files names construct with the accession ID
paired=[]
single=""

# TODO comparer len acc id in file + len fastqc_zip listing

for acc in range(totAccNumber-1) : #-1 because a \n symbol at the end of the list file

    #TODO for linked files repertory adapt the extension name
    
    # .. Test the paired reads condition
    pairedFile=givedAcc[acc]+"_1"+".fastqc.gz"  # TODO : keep [0] ?
    paired=[pairedFile]
    pairedFile=givedAcc[acc]+"_2"+".fastqc.gz"
    paired.append(pairedFile)
    if (paired[0] in dataFiles) and (paired[1] in dataFiles):
        single=givedAcc[acc]+".fastqc.gz"
        # .... If it's case you must have only two file for the accession ID
        if single in dataFiles:
            print("The database seem have three file for this accession ["+givedAcc[acc]+"]. Please verify them.")
            # TODO exit program
       # else:
          #  added=[givedAcc[acc]+'_1']
           # added.append(givedAcc[acc]+'_2')
           # totPairedSample.append(added)

        # .... Other tested cases files contains : 0 paired (p) + 0 single (s), 0p+1s, 1p+0s, 1p+1s
    else :
        single=givedAcc[acc]+".fastqc.gz"
        if single in dataFiles:
            if (paired[0] in dataFiles) or (paired[1] in dataFiles):
                print("The database seems contain one fastq file for paired reads data AND one fastq file for single reads [ACCESSION : "+givedAcc[acc]+"]. Please manually check if paired or single reads and put the right data files.")
                # TODO exit program
          #  else :
              #  totSingleSample.append(givedAcc[acc][0])
        else:
            if (paired[0] in dataFiles) or (paired[1] in dataFiles):
                print("The sample seem to be paired reads but database contain only one of the two fastq file for this accession ["+givedAcc[acc]+"]. Please manually add the second.")
                # TODO exit program
            else:
                print("The database seem don't contain the data for this accession ["+givedAcc[acc]+"]. Please manually add it.")
                # TODO exit program



# ---- Read multiQC data

#If run CQ always from the same CQ directory
# WARNING : my CQ not generate plot data text file --> multiQC -p option
# But multiqc_fastqc.txt is generated all times and contains general information for each sample

file=open(multiqcPath,"r")                               # PATHDEF
lines=(file.read()).split("\n")

# .. Check the header : read correctly data
header=(lines[0]).split("\t")
multiqcNormalHeader=['Sample', 'Filename', 'File type', 'Encoding', 'Total Sequences', 'Sequences flagged as poor quality', 'Sequence length', '%GC', 'total_deduplicated_percentage', 'avg_sequence_length', 'basic_statistics', 'per_base_sequence_quality', 'per_sequence_quality_scores', 'per_base_sequence_content', 'per_sequence_gc_content', 'per_base_n_content', 'sequence_length_distribution', 'sequence_duplication_levels', 'overrepresented_sequences', 'adapter_content']

#TODO to delete
# TAB EXPLAINED : multiQC content
#print('\n\n********* multiqcHeader ***********\n\n')
#for element in range(0,(len(multiqcNormalHeader))):
#    print(element,'\t',multiqcNormalHeader[element])
#print('**************************************\n\n')

if header!=multiqcNormalHeader :
    print("MultiQC report file doesn't seem like expected !\nExpected header : "+multiqcNormalHeader)
    print("Header on "+filePath+" : "+header)
    print("\nPlease fix the multiQC output or the analyse script.")
    # TODO exit program


# .. Make tab with sample important informations
# Tab with multiQC output analyse 
totSample=[]



for line in range(1,len(lines)): 
    infos=(lines[line]).split("\t")

    if infos != [''] :

        newEntry=[]

        # ....newEntry[0] Read sample name
        newEntry.append(infos[0])

        # ....newEntry[1] Technology used for sequencing
        newEntry.append(infos[3])

        # ....newEntry[2] Total sequences count
        newEntry.append(float(infos[4]))

        # ....newEntry[3] Sequences average Length
        newEntry.append(float(infos[9]))

        # ....newEntry[4] There is adaptator ?
        info=infos[19]
        if (info =="pass") or (info =="Pass") or (info =="PASS") :
            newEntry.append(False)
        else :
            newEntry.append(True)

        # ....newEntry[5] There is overrepresented sequence ? [T/F, PolyX seq overrepresented]
        info=infos[18]
        added=[]
        if (info !="pass") and (info !="Pass") and (info !="PASS") :  #TODO TODO Improve my script : always lower case ? il it's case not need all tests in if condition and all other in the script
            added.append(True)
        else :
            added.append(False)
        newEntry.append(added)

  

        info =infos[13]
        # ....newEntry[6] There is WARN in Per base sequence Content ?
        # ....newEntry[7] There is FAIL in Per base sequence Content ?
        newEntryFdata=[]
        newEntryWdata=[]
        if (info =="pass") or (info =="Pass") or (info =="PASS") :
            newEntryWdata.append(False)
            newEntryFdata.append(False)
        else :
            if (info =="warn") or (info =="Warn") or (info =="WARN") :
                newEntryWdata.append(True)
                newEntryFdata.append(False)
                #Analyse only Warn to add list of warned pos in read
            else : 
                newEntryFdata.append(True)
                newEntryWdata.append("Not evaluated")
                #Analyse Warn and Fail 
                # Add list of failed pos in read
                #Set newEntryWdata T or F with pos len (0 is F)
                # add list of warn pos in read
        newEntry.append(newEntryWdata)
        newEntry.append(newEntryFdata)

        # .... Put the new entry on the tab
        totSample.append(newEntry)


file.close()

############################
# Analysing the fastQC
############################

# .. Open others files to complete the analyse 
for sample in range (0,(len(totSample))) :

    fastqc_zip=inputdir+'/fastqc_zip' #Todo path to modify when automatize thz pipeline # PATHDEF
    fastqc_zip_dir=fastqc_zip+'/'+totSample[sample][0]+'_fastqc.zip'
    subprocess.run(["unzip","-q",fastqc_zip_dir,"-d",tmpdir])

    # .. OPEN 
    fastqc_zip_file = tmpdir+'/'+totSample[sample][0]+'_fastqc'+'/fastqc_data.txt'
    file=open(fastqc_zip_file,'r')

    fastqcInfo=(file.read()).split('#')

    # .... Adding bad quality in read - Entry[8] & Entry[9]  & Entry[10]

    #TODO or #test. if F result length = 0, else result=[T,(position,meanQualValue,minQualValue),(position,meanQualValue,minQualValue)]   
    fastqcPart=fastqcInfo[4]
    fastqcReaded=((fastqcPart.split('>'))[0].split('\n'))
    #print('----------- What read for fastqc called fastqc info ----------\n',fastqcReaded[0],'\n\n\n')

    # Position(s) with median number of read < 28 - [True , (pos, pos, pos)] OR [False , ()]
    median=[]
    # Position(s) with 10th percentile quality of read < 20 - [True , (pos, pos, pos)] OR [False , ()]
    minimal20=[]
    # Position(s) with 10th percentile quality of read < 10 - [True , (pos, pos, pos)] OR [False , ()]
    minimal10=[]
    # There is position(s) with wrong quality
    targetMed=False
    targetMin20=False
    targetMin10=False
    MedPos=[]
    Min20Pos=[]
    Min10Pos=[]

    for pos in range(1,(len(fastqcReaded)-1)) :
        infos=(fastqcReaded[pos]).split('\t')
        if (float(infos[2])<28) or (float(infos[5])<20) :
            if (float(infos[2])<28) :

                if '-' in infos[0] :
                    suppPos=range(int((infos[0].split('-'))[0]),int((infos[0].split('-'))[1]))
                    MedPos.extend(suppPos)
                else:
                    MedPos.append(float(infos[0]))

                if targetMed==False :
                    targetMed=True

            if (float(infos[5])<20) :
                if targetMin20==False:
                    targetMin20=True
                if (float(infos[5])<10):
                    if targetMin10==False:
                        targetMin10=True

                    if '-' in infos[0] :  
                        suppPos=range(int((infos[0].split('-'))[0]),int((infos[0].split('-'))[1]))
                        Min10Pos.extend(suppPos)
                        Min20Pos.extend(suppPos)
                    else:
                        Min10Pos.append(float(infos[0]))
                        Min20Pos.append(float(infos[0]))

                else:

                    if '-' in infos[0] : 
                        suppPos=range(int((infos[0].split('-'))[0]),int((infos[0].split('-'))[1]))
                        Min20Pos.extend(suppPos)
                    else:
                        Min20Pos.append(float(infos[0]))
    
    median.append(targetMed)
    median.append(MedPos)
    minimal20.append(targetMin20)
    minimal20.append(Min20Pos)
    minimal10.append(targetMin10)
    minimal10.append(Min10Pos)


    totSample[sample].append(median)
    totSample[sample].append(minimal20)
    totSample[sample].append(minimal10)

    # .... Adding bad content in read - Entry[5] & Entry[6]

# ++++++ REMAINDER +++++++  #TODO delete after debug
# newEntry[5] There is WARN in Per base sequence Content ? newEntryWdata
# newEntry[6] There is FAIL in Per base sequence Content ? newEntryFdata

#if (info =="pass") or (info =="Pass") or (info =="PASS") :
    #newEntryWdata.append(False)    -     newEntryFdata.append(False)
#else :
    #if (info =="warn") or (info =="Warn") or (info =="WARN") :
        #newEntryWdata.append(True)    -     newEntryFdata.append(False)
    #else : 
        #newEntryFdata.append(True)    -     newEntryWdata.append("Not evaluated")
                #Analyse Warn and Fail 
                # Add list of failed pos in read
                #Set newEntryWdata T or F with pos len (0 is F)
                # add list of warn pos in read
# +++++++++++++++++++++

    Wcontent=totSample[sample][6]


    Fcontent=totSample[sample][7]


    posW=[]
    posF=[]
    W=[False] #The evaluated level during analyses
    F=[False]

    if (Wcontent[0]==True) or (Wcontent[0]=="Not evaluated") :   

        fastqcPart=(fastqcInfo[6].split('>'))[0]
        fastqcReaded=fastqcPart.split('\n')

        for pos in range(1,(len(fastqcReaded)-1)) :
            
            line=fastqcReaded[pos].split('\t')
            values=sorted([float(line[2]),float(line[3])])    
            diffAT=values[1]-values[0]
            values=sorted([float(line[1]),float(line[4])])
            diffGC=values[1]-values[0]

            if (diffAT>10) or (diffGC>10) :
        
                if (diffAT>20) or (diffGC>20) : 
                    F=[True]
                    if '-' in line[0] :
                        suppPos=range(int((line[0].split('-'))[0]),int((line[0].split('-'))[1]))
                        posF.extend(suppPos)
                    else :
                        posF.append(int(line[0]))

                else:
                    W=[True]
                    if '-' in line[0] :
                        suppPos=range(int((line[0].split('-'))[0]),int((line[0].split('-'))[1]))
                        posW.extend(suppPos)
                    else :
                        posW.append(int(line[0]))
           
                        
                        
    totSample[sample][6]=W
    totSample[sample][6].append(posW)  

    totSample[sample][7]=F
    totSample[sample][7].append(posF)

    # .... Adding polyX seq in overrepresented - Entry[5][1] 

    polyX=False
    if totSample[sample][5][0]==True:
        fastqcPart=(fastqcInfo[12].split('>'))[0]
        fastqcReaded=fastqcPart.split('\n')

        for line in range(1,len(fastqcReaded)):
            seq=((fastqcReaded[line]).split('\t'))[0]
            print(seq)
            nucleotide=set((list(seq)))
            print(nucleotide)
            if len(nucleotide)==0:
                print('\nPROBLEM ON OVERREPRESENTED SEQUENCE READING !!!!!\n') #TODO treat them
            if len(nucleotide)==1 :
                polyX=True
            #if len(nucleotide)==2 :  #TODO : -y in this case ?
                #polyX='Maybe'
    totSample[sample][5].append(polyX)



    file.close()

    #dirToSup=fastqc_zip+'/'+totSample[sample][0]+'_fastqc' # PATH
    dirToSup=tmpdir+'/'+totSample[sample][0]+'_fastqc'
    subprocess.run(["rm","-r",dirToSup])

# ....Entry[11] RNA Type

# List of csv file the contain the AccID for each type of RNA
listType=['transcriptomic_circRNAAccession.csv','transcriptomic_lncRNAAccession.csv','transcriptomic_RandomRNAAccession.csv','transcriptomic_SMRTAccession.csv','transcriptomic_sRNAAccession.csv','transcriptomic_totalAccession.csv']
accType=[]

for listedCSV in listType :
    createdList=[]
    #.. Define the first element of this list
    split1=listedCSV.split('_')
    split2=(split1[1]).split('Accession')
    type=split2[0]
    createdList.append(type)

    #.. Fill the list
    path=currentdir+listedCSV
    with open(path, mode='r') as csvFile :
        file=csv.reader(csvFile, delimiter=',')
        for row in file :
            createdList.append(row)
    accType.append(createdList)


for sample in range(len(totSample)) :

    totSample[sample].append('Not set')
    for type in range(len(accType)) :
        if totSample[sample][0] in accType[type] :
            totSample[sample][11]=accType[type][0]
    
    if totSample[sample][11]=='Not set':
        totSample[sample][11]='mRNA'

    #TODO how evaluated if random RNA is mRNA or sRNA ???? Or other
    if totSample[sample][11]=='RandomRNA':
        print('\n\nSome RANDOM RNA have been found !!!!! \n\n')

############################
# Define the fastP options                    #TODO TODO TODO base op add -n ==> see how fastP do by default
############################


# ....Entry[13] fastP options
# (or put after file close and tab make)
# --> in that case
# + if overrepresented sequence in data where minlenght >100 bp : see if polyX nuceotidique sequence --> open fastqc zipped dir and read the info OR read on multiqc meta data
# + if true qual at pos < 28 : see if some read on sample with bas qual < 28 OR red limit --> open fastqc zipped dir and read the info OR read on multiqc meta data
#           : at the end juste -3 -5 ok, if in read peak ou si après trimming put quality treshord to red limit and % of base not on trheshold to 0
dfopt='-5 -3 -M 28 -u 0 -D -n 0'
trashed=[]

for sample in range(len(totSample)) :

    newEntry=[]
    # If the sample can't be treated by fastP and it's analysis ended
    sampleTrash=False

    name=totSample[sample][0]
    if '_' in name :
        totSample[sample].append('paired')
    else :
        totSample[sample].append('single')

    opt=dfopt
    

    #.... Did I put the --adapter-for-pe option ?  #TODO TODO TODO si single read et adapter ok ==> desactiver le trimming des adpter pour fastp
    if totSample[sample][4]==True :
        if totSample[sample][12]=='paired' :
            opt=opt+' --detect_adapter_for_pe'

    length=round(totSample[sample][3])
        
    #.... Did I cut reads front ? #TODO define for the other type of RNA
    cutInOpt=0
    if (totSample[sample][11]=='mRNA') and ( (totSample[sample][6][0]==True) or (totSample[sample][7][0]==True) ):
        if (totSample[sample][6][0]==True) and (totSample[sample][7][0]==True):
            #pos=totSample[sample][5][1].extend(totSample[sample][6][1]) #Notwork
            pos=totSample[sample][6][1]
            for value in totSample[sample][7][1]:
                pos.append(value)
        else :
            if totSample[sample][6][0]==True :
                pos=totSample[sample][6][1]
            if totSample[sample][7][0]==True :
                pos=totSample[sample][7][1]
            # Keep greater value < 26 (read begin)
        position=0
        for value in pos :
            if value < 26: 
                if value >position:
                    position=value
        cutInOpt=position+1 # We must indicate to fastp the first nucleotide to didn't trimm
        add=' -f '+str(cutInOpt)
        opt=opt+add
    # if totSample[sample][11]=='sRNA':, do nothing for that opt


    #.... set the minimal (and maximal for sRNA) lenght of trimmed reads (default : 15) #TODO define for the other type of RNA
    if totSample[sample][11]=='mRNA':
        
        if length!=LminLim:
            diffLength=max(length,LminLim)-min(length,LminLim)
            if length<LminLim:
                percentage=minPercentL-(diffLength*LpercentByNt)
            else :
                percentage=minPercentL+(diffLength*LpercentByNt)
        else:
            percentage=minPercentL

        finalLength=int(percentage*(length-position))
        if finalLength<15: #15 is default setting
            finalLength=15
        added=' -l '+str(finalLength)
        opt=opt+added

    if totSample[sample][11]=='sRNA':
        if totSample[sample][7][0]==True:
            pos=totSample[sample][7][1]
            # All position > 17 (default : -l 18)
            trimmedPos=[]
            for target in pos:
                # Value I expected don't have fail because normal RNA length
                if 17<target<=25 :
                    sampleTrash=True
                    print('\n\n Problem with sample ',totSample[sample][0],' fail on pos 17-25 on per base sequence content')
                    trashed.append(totSample[sample][0])
                if 25<target :
                    trimmedPos.append(target)
            if sampleTrash==False:
                trimmedPos=trimmedPos.sort()
                cutUpperPos=trimmedPos[0]-1
                added=' -l 18 -b '+str(cutUpperPos)
                opt=opt+added

    #.. If pass the sRNA sequence content analyse (18-32 nt long un sRNA)
    if sampleTrash==False :  

        #.... Set the quality threshold (if different for paired read ==> take the max)
        if length>=maxLim :
            maxFalsePercent=maxFpercent # maxFalsePercent is the maximum probability of nt is not the real (error during sequencing)
        else:
            if length<=minLim :
                maxFalsePercent=minFpercent
            else:
                diffLength=maxLim-length
                maxFalsePercent=maxFpercent+(diffLength*FpercentByNt)
        
        #Cacul the according quality Threshord
        minQTot=-10*(log10(maxFalsePercent))
        minQ=round(minQTot)
        if minQ>10 :
            minQ=10
        #Put result on options
        opt=opt+' -q '+str(minQ)



        #.... Did i put the -y option ? #Don't do for now, see later
        #if length>100 :
            #opt=opt+' -y -Y 1'



            #diversityOption=True  #TODO TODO for length<100, change the % of minimal diversity according with length / change maximal authorized number of position with quality bad
        #else :
            #print('/////////////// ',totSample[sample][4])
            #if totSample[sample][5][1]==True : #TOOD if i keep maybe, change here with != false
                #opt=opt+' -y'
        #if diversityOption==True:
            #opt=opt+' -y'
      
        #.... Now just print sample with 90th lower percentile < 10 for quality, after treat them with position (change percentage of base authorized with quality < treshold (warning threshold))
        # And after see if I can refine the window size for read quality # TODO TODO
    
        if totSample[sample][10][0]==True :
            print('\n$$$$$$$$$$$$$$$$$ 90th lower percentile < 10 quality for sample '+totSample[sample][0]+' $$$$$$$$$$$$$$$$$')
            print('If result = [], this value is upper for all reads positions')
            print('Position with Median quality <28 : ', totSample[sample][8][1])
            print('Position with 90th Lower Percentile quality <20 : ', totSample[sample][9][1])
            print('Position with 90th Lower Percentile quality <10 : ', totSample[sample][10][1])
            print('$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n\n')

 


    totSample[sample].append(opt)



print('\n-----------  SAMPLE(s) TRASHED FOR ANY REASON  ------')
print(trashed)
print('\nSample I cant deddup :')

############################
# Output file fill
############################


#TODO ++++++++
# to delete aftre debug

outputFile=outdir+'/LastAnalyseRes.txt'              #PATH
output=open(outputFile,'w')
output.write('')
#Header
outputHeader='Sample\tSequencing technology\tTotal Seq\tMean Length\tAdaptor\tOverrepresentedSeq/polyXSeq\tPerBaseSeqContent WARN (T/F,pos)\tPerBaseSeqContent FAIL\tMedQualLow28\tMinQualLow20\tMinQualLow10\tDupplicated(level,percentageNotDuplicated)\tRNAtype\tSingleOrPairedRead\tfastPoptions\n'
output.write(outputHeader)

for element in totSample :
    output.write(str(element))
    output.write('\t')
    output.write('\n')

output.close



outputFile=outdir+'/FastoOpt.txt'              #PATH
output=open(outputFile,'w')
output.write('')

for element in totSample :
    entry=element[0]+'\t'+element[-1]
    output.write(entry)
    output.write('\n')

output.close

#++++++++++


#TAB EXPLAINED : totSample
print('\n\n********* LastAnalyseRes Header ***********\n\n')
resHeader=outputHeader.split('\t')
for element in range(0,(len(resHeader))):
    print(element,'\t',resHeader[element])
print('**************************************\n\n')


############################
# Make the trimming sub data set
############################

## .. Do paired reads opt the same for pair of reads are on the same dataset at the end

#TODO TODO TODO make a founction to define between thze two value the one keeped

for sample in range(len(totSample)-1):

    if '_1' in totSample[sample][0]:
        if totSample[sample][13] != totSample[sample+1][13]:
            optR1=(totSample[sample][13]).split(" ")
            optR2=(totSample[sample+1][13]).split(" ")

            print('~~~~~~~~~~ ONE AJUSTMENT MAKE : ',totSample[sample][0],totSample[sample][-2],totSample[sample][-1],' & ',totSample[sample+1][0],totSample[sample+1][-2],totSample[sample+1][-1])

            finalopt=dfopt
            maxInd=max(len(optR1),len(opt2))
            R1f=''
            R1l=''
            R1b=''
            R1y=''
            R1q=''
            R2f=''
            R2l=''
            R2b=''
            R2y=''
            R2q=''
        # Compare values can be différent and decide the one to keep from each
            for i in range(10,(maxInd-1),2):

                if i<=(len(optR1)-1):
                #.. put values on the good variable
                    if optR1[i]=='-f':
                        R1f=int(optR1[i+1])
                    else :
                        if optR1[i]=='-l':
                            R1l=int(optR1[i+1])
                        else:
                            if optR1[i]=='-b':
                                R1b=int(optR1[i+1])
                            else:
                                if optR1[i]=='-y':
                                    R1y=True
                                else:
                                    if optR1[i]=='-q':
                                        R1q=int(optR1[i+1])
                                    else:
                                        print('\n\n @@@@@@@@@@@@@@@@@@ BE AWARE : unknow option or fastp otion indice not set correctly') #TODO delete after debug

                if i<=(len(optR2)-1):
                    if optR2[i]=='-f':
                        R2f=int(optR2[i+1])
                    else :
                        if optR2[i]=='-l':
                            R2l=int(optR2[i+1])
                        else:
                            if optR2[i]=='-b':
                                R2b=int(optR2[i+1])
                            else:
                                if optR2[i]=='-y':
                                    R2y=True
                                else:
                                    if optR2[i]=='-q':
                                        R2q=int(optR2[i+1])
                                    else:
                                        print('\n\n @@@@@@@@@@@@@@@@@@ BE AWARE : unknow option or fastp option indice not set correctly') #TODO delete after debug
        
            #.. compare each different options
            nbDiff=0
            optYadd=False

            if R1f!='' or R2f!='':
                if R1f==R2f:
                    finalopt=finalopt+' -f '+str(R1f)
                if R1f!=R2f and R1f!='' and R2f!='':
                    nbDiff=nbDiff+1
                    finalopt=finalopt+' -f '+str(max(R1f,R2f))
                else:
                    print('\n\n @@@@@@@@@@@@@@@@@@ BE AWARE : This case not must happen for paired read and -f OPT')
                    if R1f=='' :
                        finalopt=finalopt+' -f '+str(R2f)
                    if R2f=='' :   # So case R1f==R2f with '' is not treat
                        finalopt=finalopt+' -f '+str(R1f)


            if R1l!='' or R2l!='':
                if R1l==R2l:
                    finalopt=finalopt+' -l '+str(R1l)
                if R1l!=R2l and R1l!='' and R2l!='':
                    nbDiff=nbDiff+1
                    finalopt=finalopt+' -l '+str(mean(R1l,R2l))
                else:
                    print('\n\n @@@@@@@@@@@@@@@@@@ BE AWARE : This case not must happen for paired read and -l OPT')
                    if R1l=='' :
                        finalopt=finalopt+' -l '+str(R2l)
                    if R2l=='' :  
                        finalopt=finalopt+' -l '+str(R1l)


            if R1b!='' or R2b!='':
                if R1b==R2b:
                    finalopt=finalopt+' -b '+str(R1b)
                if R1b!=R2b and R1b!='' and R2b!='':
                    nbDiff=nbDiff+1
                    finalopt=finalopt+' -b '+str(min(R1b,R2b))
                else:
                    print('\n\n @@@@@@@@@@@@@@@@@@ BE AWARE : This case not must happen for paired read and -b OPT')
                    if R1b=='' :
                        finalopt=finalopt+' -b '+str(R2b)
                    if R2b=='' :  
                        finalopt=finalopt+' -b '+str(R1b)

            if R1y!='' or R1y!='':
                if R1y==R2y:
                    finalopt=finalopt+' -y -Y 1'
                else:
                    finalopt=finalopt+' -y -Y 1'
                    optYadd=True
                    print('\n\n @@@@@@@@@@@@@@@@@@ BE AWARE : -y option différent btwenn this two sample (median length one <100 other >100)',totSample[sample][0],totSample[sample+1][0])


            if R1q!='' or R2q!='':
                if R1q==R2q:
                    finalopt=finalopt+' -q '+str(R1q)
                if R1q!=R2q and R1q!='' and R2q!='':
                    nbDiff=nbDiff+1
                    finalopt=finalopt+' -q '+str(max(R1q,R2q))
                else:
                    print('\n\n @@@@@@@@@@@@@@@@@@ BE AWARE : This case not must happen for paired read and -q OPT')
                    if R1q=='' :
                        finalopt=finalopt+' -q '+str(R2q)
                    if R2q=='' :  
                        finalopt=finalopt+' -q '+str(R1q)


            if nbDiff==0 and optYadd==False:
                print('\n\n @@@@@@@@@@@@@@@@@@ BE AWARE : I find différent fastP option bus after parse them not diff !!!', optR1,optR2)

            print('final gived option : ', finalopt)





# TAB EXPLAINED : rst tab
print('\n\n********* dataset table index ***********\n\n')
index='DatasetX (with X the ID of dataset)\tPipeline paired option (if \'paired : 1\', reads are paired)\tfastp options\tList of sample in this dataset'
resIndex=index.split('\t')
for element in range(0,(len(resIndex))):
    print(element,'\t',resIndex[element])
print('**************************************\n\n')


datasetID=1
rst=[] #See above
sampledst=''
# Variable to link fastp option already found with corresponding dataset ID
# [ [dataset1, opt] , [dataset2, opt] ]  --> Have matched indice
# OR
# **** [dataset1, dataset2]  [opt1, opt2] with dataset and opt adding in the same  order --> Have matched indice
# OR
# [opt1, opt2]   [ [dataset1, opt] , [dataset2, opt] ]
opt=[]
dtset=[]
dstP=[]

for sample in range(0,len(totSample)):

    # .. Give the dataset name to the sample

    if len(opt)==0:
        # Update existing dst
        opt.append(totSample[sample][13])
        dtset.append('dataset'+str(datasetID))
        if '_' in totSample[sample][0]:
            dstP.append(1)
        else:
            dstP.append(0)

        #Save sample dst info
        sampledst='dataset'+str(datasetID)
        datasetID=datasetID+1

        #Init rst tab
        Entry=[]
        spList=[]
        Entry.append(sampledst)
        Entry.append(dstP[0])
        Entry.append(opt[0])
        spList.append(totSample[sample][0])
        Entry.append(spList)
        rst.append(Entry)

        print(opt,dtset,dstP)

    else:
        #Search if sample opyt in opt tab with opt.index(sampleOPT) --> treat the result, see what founction return if not find (if error, use if opsample in opt before)
        if totSample[sample][13] in opt :  #TODO sauv - info : sauv juste les options qui ne sont pas par défautin opt seulemnt
            finded=False
            indexFinded=''
            occ=[]

            for i in range(len(dtset)):
                if totSample[sample][13]==opt[i]:
                    occ.append(int(i))
            if len(occ)>2: #Delete after debug TODO
                print('\n\nêêêêêêêêêêêêêê - Resultat pas possible !!!! 1 ')

            print('--- Finded occurences : ',occ)

            if len(occ)==2 : #Two dataset with the same options : can occur when read is paired
                if '_' in totSample[sample][0]:
                    spP=1
                else:
                    spP=0
                
                for occurence in occ:
                    if dstP[occurence]==spP:
                        indexFinded==occurence
                        datasetID=dtset[occurence]
                        finded=True
            else:
                indexFinded=occ[0]
                datasetID=dtset[occ[0]]
                finded=True

            if finded==True:
                print(totSample[sample][13],' --> ',datasetID)
            else:
                print('\n\nêêêêêêêêêêêêêê - Resultat pas possible !!!! 2 ')

            #Put sample in rst tab
            rst[indexFinded][3].append(totSample[sample][0])

        else :
            opt.append(totSample[sample][13])
            dtset.append('dataset'+str(datasetID))
            sampledst='dataset'+str(datasetID)
            datasetID=datasetID+1

            #Add on rst tab
            Entry=[]
            spList=[]
            Entry.append(sampledst)
            Entry.append(dstP[-1]) # The last I put on tabs
            Entry.append(opt[-1])
            spList.append(totSample[sample][0])
            Entry.append(spList)
            rst.append(Entry)

  
print(rst)


############################
# Performe the trimming
############################

############################
# Analyse each quality of trimmed data : one not OK message to user
############################

# ---- Coverage analysis

# ---- Multiqc analyse

############################
# Performe mapping
############################


