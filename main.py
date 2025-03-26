# Use it : give in argument one txt file with the accession ID / or a csv file + path where put output + analyse name(without space or special charactère)

#BALISE # PATHDEF where a define path
#==> see them when ai write this script to be used by nextflow or on cluster


# ---- Import

import os
import sys
import subprocess
import csv

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

        # ....newEntry[1] Total sequences count
        newEntry.append(float(infos[4]))

        # ....newEntry[2] Sequences average Length
        newEntry.append(float(infos[9]))

        # ....newEntry[3] There is adaptator ?
        info=infos[19]
        if (info =="pass") or (info =="Pass") or (info =="PASS") :
            newEntry.append(False)
        else :
            newEntry.append(True)

        # ....newEntry[4] There is overrepresented sequence ? [T/F, PolyX seq overrepresented]
        info=infos[18]
        added=[]
        if (info =="fail") or (info =="Fail") or (info =="FAIL") :
            added.append(True)
        else :
            added.append(False)
        newEntry.append(added)

  

        info =infos[13]
        # ....newEntry[5] There is WARN in Per base sequence Content ?
        # ....newEntry[6] There is FAIL in Per base sequence Content ?
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

    # .... Adding bad quality in read - Entry[7] & Entry[8]  & Entry[9]

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

    Wcontent=totSample[sample][5]


    Fcontent=totSample[sample][6]


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
           
                        
                        
    totSample[sample][5]=W
    totSample[sample][5].append(posW)  

    totSample[sample][6]=F
    totSample[sample][6].append(posF)

    # .... Adding dupplication analyse - Entry[10]

    fastqcPart=(fastqcInfo[9]).split('>')
    fastqcReaded=(fastqcPart[-1]).split('\t')
    # Pass, Warn or Fail level
    duplLevel=(fastqcReaded[-1].split('\n'))[0]

    fastqcPart=fastqcInfo[10]
    fastqcReaded=(fastqcPart.split('\t'))
    # Percentage of sequence of total sequence if we dedupplicated them
    percentageDedup=float((fastqcReaded[-1].split('\n'))[0])

    newEntry=[duplLevel,percentageDedup]
    totSample[sample].append(newEntry)

    file.close()

    #dirToSup=fastqc_zip+'/'+totSample[sample][0]+'_fastqc' # PATH
    dirToSup=tmpdir+'/'+totSample[sample][0]+'_fastqc'
    subprocess.run(["rm","-r",dirToSup])

    # .... Adding polyX seq in overrepresented - Entry[4][1] 

    polyX=False
    if totSample[sample][4][0]==True:
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
    totSample[sample][4].append(polyX)



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
# Define the fastP options
############################


# ....Entry[13] fastP options
# (or put after file close and tab make)
# --> in that case
# + if overrepresented sequence in data where minlenght >100 bp : see if polyX nuceotidique sequence --> open fastqc zipped dir and read the info OR read on multiqc meta data
# + if true qual at pos < 28 : see if some read on sample with bas qual < 28 OR red limit --> open fastqc zipped dir and read the info OR read on multiqc meta data
#           : at the end juste -3 -5 ok, if in read peak ou si après trimming put quality treshord to red limit and % of base not on trheshold to 0

fastp=[]
undedupp=[]
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

    dfopt='-5 -3 -M 28'
    opt=dfopt
    
    #.. Dedupplication and list all sample that can't be dedupplicated without poor sample number in analysis
    undedupp=[]

    if (totSample[sample][10][0]=='warn') or (totSample[sample][10][0]=='fail') :
        afTrimRead=(float(totSample[sample][10][1])/100)*totSample[sample][1]
        if totSample[sample][11]=='mRNA' :
            if (afTrimRead<8000000) or (30000000<afTrimRead<45000000):
                undedupp.append(totSample[sample][0])
                undedupp.append(afTrimRead)
                sampleTrash=True
                trashed.append(totSample[sample][0])
        if totSample[sample][11]=='sRNA' :
            if (afTrimRead<900000) or (6000000<afTrimRead):
                undedupp.append(totSample[sample][0])
                undedupp.append(afTrimRead)
                sampleTrash=True
                trashed.append(totSample[sample][0])
        #TODO define for the other type of RNA
    
    #.. If pass the dedupplication analyse
    if sampleTrash==False :

        #.... Did I put the --adapter-for-pe option ?
        if totSample[sample][3]==True :
            if totSample[12]=='paired' :
                opt=opt+' --detect_adapter_for_pe'

        length=round(totSample[sample][2])

        #.... Did i put the -y option ? #TODO define for the other type of RNA
        diversityOption=False
        if length<=100 :
            diversityOption=True
        else :
            if totSample[sample][4][1]==True : #TOOD if i keep maybe, change here with != false
                opt=opt+' -y'
        if diversityOption==True:
            opt=opt+' -y'
        
        #.... Did I cut reads front ? #TODO define for the other type of RNA
        cutInOpt=0
        print(totSample[sample][5][1])
        print(totSample[sample][6][1])
        if (totSample[sample][11]=='mRNA') and ( (totSample[sample][5][0]==True) or (totSample[sample][6][0]==True) ):
            print('jentre in boucle 1')
            if (totSample[sample][5][0]==True) and (totSample[sample][6][0]==True):
                print('suis in boucle2')
                #pos=totSample[sample][5][1].extend(totSample[sample][6][1])
                pos=totSample[sample][5][1]
                for value in totSample[sample][6][1]:
                    pos.append(value)
                print('######### ',pos)
            else :
                if totSample[sample][5][0]==True :
                    pos=totSample[sample][5][1]
                if totSample[sample][6][0]==True :
                    pos=totSample[sample][6][1]
                # Trim pos to keep them < 25 (read begin)
            begPos=[]
            for value in pos :
                print('------ ',value)
                if value < 26:
                    begPos.append(int(value))
            print(begPos[-1])
            begPos=begPos.sort()
            print(begPos)
            cutInOpt=int(begPos[-1])+1 # We must indicate to fastp the first nucleotide to didn't trimm
            add=' -f '+str(cutInOpt)
            opt=opt+add
        # if totSample[sample][11]=='sRNA':, do nothing for that opt

        #.... Dedupplicate
        if totSample[sample][10][0]!='pass':
            opt=opt+' -D'


        #.... set the minimal (and maximal for sRNA) lenght of trimmed reads (default : 15) #TODO define for the other type of RNA
        if totSample[sample][11]=='mRNA':
            percentage=((length-int(begPos[-1]))*0.6)/150
            finalLength=int(percentage*(length-int(begPos[-1])))
            print('Final length : ',percentage,' * (',length,' - ',int(BegPos[-1]),')') #TODO delete after debugg

        if totSample[sample][11]=='sRNA':
            if totSample[sample][6][0]==True:
                pos=totSample[sample][6][1]
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

    if sampleTrash==False :        


        #.... Now just print sample with 90th lower percentile < 10 for quality, after treat them with position (change percentage of base authorized with quality < treshold (warning threshold))
        # And after see if I can refine the window size for read quality # TODO TODO
    
        if totSample[sample][9][0]==True :
            print('\n$$$$$$$$$$$$$$$$$ 90th lower percentile < 10 quality for sample '+totSample[sample][0]+' $$$$$$$$$$$$$$$$$')
            print('If result = [], this value is upper for all reads positions')
            print('Position with Median quality <28 : ', totSample[sample][7][1])
            print('Position with 90th Lower Percentile quality <20 : ', totSample[sample][8][1])
            print('Position with 90th Lower Percentile quality <10 : ', totSample[sample][9][1])
            print('$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n\n')

        newEntry.append(totSample[sample][0])
        newEntry.append(opt)
        fastp.append(newEntry)


print('\n-----------  SAMPLE(s) TRASHED FOR ANY REASON  ------')
print(trashed)
print('\nSample I cant deddup :')
print(undedupp) #TOOD delete after debug

############################
# Output file fill
############################


#TODO ++++++++
# to delete aftre debug

outputFile=outdir+'/LastAnalyseRes.txt'              #PATH
output=open(outputFile,'w')
output.write('')
#Header
outputHeader='Sample\tTotal Seq\tMean Length\tAdaptor\tOverrepresented Seq\tPerBaseSeqContent WARN (multiqc,fastqc,pos)\tPerBaseSeqContent FAIL(multiqc,fastqc,pos)\tMedQualLow28\tMinQualLow20\tMinQualLow10\tDupplicated\tRNAtype\tSingleOrPairedRead\tfastPoptions\n'
output.write(outputHeader)

for element in totSample :
    output.write(str(element))
    output.write('\t')
    output.write('\n')

output.close



outputFile=outdir+'/FastoOpt.txt'              #PATH
output=open(outputFile,'w')
output.write('')

for element in fastp :
    entry=element[0]+'\t'+element[1]
    output.write(entry)
    output.write('\n')

output.close

#++++++++++



print('\n\n********* LastAnalyseRes Header ***********\n\n')
resHeader=outputHeader.split('\t')
for element in range(0,(len(resHeader))):
    print(element,'\t',resHeader[element])
print('**************************************\n\n')


############################
# Make the trimming sub data set
############################

# Read the csv if I do one

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


