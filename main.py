# Use it : give in argument one txt file with the accession ID / or a csv file + path where put output + analyse name(without space or special charactère)

#BALISE # PATHDEF where a define path
#==> see them when ai write this script to be used by nextflow or on cluster


# ---- Import

import os
import sys

# ---- Set path

#inputFile=sys.argv[1] #Used at the final version
inputdir="/home/smile/Bureau/TestautoPip/ManonData" # PATHDEF
outdir=sys.argv[2]

#To suppres to didn't keep file generated that user don't want
currentdir= 'Command permit to have pwd result' #TODO
tmpdir=currentdir+'/work'
#TODO suppress tmpdir content (generated in the previous run) VOIR COMMENT ENV INTERACTIF DE CODIUM GERE

# ---- Other variables

#Construct files path here


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
multiqcPath=inputdir+"/multiqc_fastqc.txt" # PATHDEF
sampleList='' #TODO make the file with alla accession, put the good path here and define givedAcc wi
file=open(sampleList,"r")
givedAcc=(file.read()).split("\n") # PATHDEF
file.close()
givedAcc.sort()
totAccNumber=len(givedAcc)

print("\n\n ----- SAMPLE ANALYZED --------\n",givedAcc[0]) #TEMP

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
dataFiles=[]

# Files names construct with the accession ID
paired=[]
single=""

# TODO comparer len acc id in file + len fastqc_zip listing

for acc in range(totAccNumber) :
    
    # .. Test the paired reads condition
    pairedFile=givedAcc[acc][0]+"_1"+".fastq.gz"  # TODO : keep [0] ?
    paired=[pairedFile]
    pairedFile=givedAcc[acc][0]+"_2"+".fastq.gz"
    paired.append(pairedFile)
    print(paired) #TEMP
    if (paired[0] in dataFiles) and (paired[1] in dataFiles):
        single=givedAcc[acc][0]+".fastq.gz"
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
        single=givedAcc[acc][0]+".fastq.gz"
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
if header!=multiqcNormalHeader :
    print("MultiQC report file doesn't seem like expected !\nExpected header : "+multiqcNormalHeader)
    print("Header on "+filePath+" : "+header)
    print("\nPlease fix the multiQC output or the analyse script.")
    # TODO exit program


# .. Make tab with sample important informations
# Tab with multiQC output analyse 
totSample=[]

#TODO ++++++++
# to delete aftre debug
subprocess.run(["touch","\'LastAnalyseRes.txt\'"])
outputFile=currentdir+'/LastAnalyseRes.txt'              #PATH
output=open(outputFile,'w')
#++++++++++

for line in range(1,len(lines)): 
    infos=(lines[line]).split("\t")

    newEntry=[]

    # ....newEntry[0] Read sample name
    newEntry.append(infos[0])

    # ....newEntry[1] Total sequences count
    newEntry.append(infos[4])

    # ....newEntry[2] Sequences average Length
    newEntry.append(infos[9])  #TODO or take the median : chiffre rond

    # ....newEntry[3] There is adaptator ?
    info=infos[19]
    if (info =="pass") or (info =="Pass") or (info =="PASS") :
        newEntry.append(False)
    else :
        newEntry.append(True)

    # ....newEntry[4] There is overrepresented sequence ?
    info=infos[18]
    if (info =="fail") or (info =="Fail") or (info =="FAIL") :
        newEntry.append(True)
    else :
        newEntry.append(False)

  

    info =infos[14]
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
    totSample.append[newEntry]

    #Todo sup ap debug
    sauv=''
    for element in newEntry :
        sauv=sauv+'\t'+element

    .write()

file.close()
output.close()


# --------------------------------------------------------------------------

# .. Open others files to complete the analyse 
for sample in range (len(totSample)) :

    fastqc_zip=inputdir+'/fastqc_zip' #Todo path to modify when automatize thz pipeline # PATHDEF
    fastqc_zip_dir=fastqc_zip+'/'+totSample[sample][0]+'_fastqc.zip'
    subprocess.run(["unzip",fastqc_zip_dir])

    # .. OPEN 
    fastqc_zip_file = fastqc_zip+'/'+totSample[sample][0]+'_fastqc'+'/fastqc_data.txt'
    file=open(fastqc_zip_file,'r')



    # .... Adding bad quality in read - Entry[7] & Entry[8]

    # Position(s) with median number of read < 28 - [True , (pos, pos, pos)] OR [False , ()]
    #TODO or #test. if F result length = 0, else result=[T,(position,meanQualValue,minQualValue),(position,meanQualValue,minQualValue)]   
    median=[]
    # Position(s) with minimal quality of read < 10 - [True , (pos, pos, pos)] OR [False , ()]
    minimal=[]


    file.close()

# ....Entry[9] RNA Type
# ....Entry[10] fastP options
# (or put after file close and tab make)
# --> in that case
# + if overrepresented sequence in data where minlenght >100 bp : see if polyX nuceotidique sequence --> open fastqc zipped dir and read the info OR read on multiqc meta data
# + if true qual at pos < 28 : see if some read on sample with bas qual < 28 OR red limit --> open fastqc zipped dir and read the info OR read on multiqc meta data
#           : at the end juste -3 -5 ok, if in read peak ou si après trimming put quality treshord to red limit and % of base not on trheshold to 0


    dirToSup=fastqc_zip+'/'+totSample[sample][0]+'_fastqc' # PATH
    subprocess.run(["rm","-r","-i",dirToSup])


# To gain space idea : write tabs on a csv/txt file (For user permit to keep trimming option if is useful)
# Can do after two tabs completed or after each column 



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