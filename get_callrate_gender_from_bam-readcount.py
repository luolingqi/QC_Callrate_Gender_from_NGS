#!/usr/bin/python
####################################################################################################################
# This script reads in the BAM-READCOUNT output for the 96 SNPs of Biomarks (93 if Y is absent for Female) and do the following two things:
# 1) calculate the call rate based on depth cutoff
# 2) Call gender using biomarks gender calling rule and compare with the true gender in samplesheet
#
#
####################################################################################################################

import sys

in_file=sys.argv[1] # bam-readcount output
samplesheet_in = sys.argv[2] # read in the sample sheet

# Declare a list to record the information for the 96 SNPs
aList = []
out_file=in_file.replace("bam_readcount.txt","get_callrate_gender.out")
output = open(out_file,"w")

out_log=in_file.replace("bam_readcount.txt","get_callrate_gender.log")
output_log = open(out_log,"w")


def extract_true_gender():
    true_gender = "Unknown"
    gender_col = -1
    
    with open(samplesheet_in) as ss:
        for line in ss:
            line = line.strip()
            arr = line.split(",")
            # determine the column No. for true gender
            if line.startswith("Sample_ID"):
                gender_col = arr.index("Patient Gender")
            #print line.split(",")[1]
            if arr[1] in in_file: # sample name in samplesheet matches partially with bam-readcount file
                #print line
                #true_gender = line.split(",")[14]
                if gender_col != -1:
                    true_gender = arr[gender_col]
#    print "gender column No. is: ", gender_col
#    print "true gender is: ", true_gender
    return true_gender
                
def determine_genotype(info): #info is the count information for A,C,G,T
    geno = []
    base = []
    base_count = []
    arr = info.split("\t")
    for i in range(0,4):
        base.append(arr[i].split(":")[0])
        base_count.append(int(arr[i].split(":")[1]))

    # determine genotype as the combination of bases with 'Non-zero' count
    [ geno.append(base[i]) for i,count in enumerate(base_count) if count > 5 ]
    if len(geno) == 1:
        return geno[0]+":"+geno[0]
    else:
        return ":".join(geno) 


def determine_gender(d_sex):
    sex = "NA"
    # determine female
    x = set()
    y = set()
    for key in d_sex:
        if key.endswith("X"):
            x.add(d_sex[key])
    for key in d_sex:
        if key.endswith("Y"):
            y.add(d_sex[key])
    
    # determine female
    if len(y) == 1 and list(y)[0] == "No Call:No Call":
        sex = "Female"
    # determine male
    if "No Call:No Call" not in y and all([len(set(i.split(":")))==1 for i in list(x)]): # check all 3 X-SNP are homozygous
        sex = "Male"
    if "No Call:No Call" not in y and any([len(set(i.split(":")))!=1 for i in list(x)]): # check any one of the 3 S-SNPs is heterozygous
        sex = "Klinefelter"
    return sex


# All the 6 gender SNPs (chromosome & hg38 position)
aDict = {}
aDict['chrY-7000077'] = ['hu98Y','No Call','No Call']
aDict['chrX-91341160'] = ['hu103X','No Call','No Call']
aDict['chrX-117469443'] = ['hu107X','No Call','No Call']
aDict['chrX-149485920'] = ['hu109X','No Call','No Call']
aDict['chrY-12914512'] = ['hu111Y','No Call','No Call']
aDict['chrY-15174113'] = ['hu209Y','No Call','No Call']


sex_dict = {}
sex_dict['hu98Y'] = "No Call"+":"+ "No Call"
sex_dict['hu103X'] = "No Call"+":"+ "No Call"
sex_dict['hu107X'] = "No Call"+":"+ "No Call"
sex_dict['hu109X'] = "No Call"+":"+ "No Call"
sex_dict['hu111Y'] = "No Call"+":"+ "No Call"
sex_dict['hu209Y'] = "No Call"+":"+ "No Call"
nocall = 0

sex_snps = ['chrY-7000077','chrX-91341160','chrX-117469443','chrX-149485920','chrY-12914512','chrY-15174113']
x_snps = ['chrX-91341160','chrX-117469443','chrX-149485920']
y_snps = ['chrY-7000077','chrY-12914512','chrY-15174113']
print >>output_log, "Chromesome\tPosition\tReference\tRead Count"
print >>output, "Call rate,Pass callrate,Gender,Pass gender"

with open(in_file) as fin:
    for line in fin:
        line = line.strip()
        arr = line.split("\t")
        if arr[3] < 5 and '-'.join(arr[0],arr[1]) not in y_snps: # depth < 5 means no call
            nocall += 1
        aKey = "-".join(arr[0:2])
        print>>output_log, arr[0]+"\t"+arr[1]+"\t"+arr[2]+"\t"+arr[3]
        if aKey in sex_snps:
            sex_dict[aDict[aKey][0]] = determine_genotype("\t".join(arr[5:9]))



gender = determine_gender(sex_dict)
true_gender = extract_true_gender()

#print true_gender
pass_call = "FALSE"
pass_gender = "TRUE"
call_rate = (93-nocall)/93.0*100
if call_rate > 0.95:
    pass_call = "TRUE"
if true_gender.lower() == gender.lower():
    pass_gender = "TRUE"
              
print>>output, "{0:.2f}%".format((93-nocall)/93.0*100) + "," + pass_call + "," + gender + "," + pass_gender
print>>output_log, "True gender is: "+true_gender

        #print arr[0]+"\t"+arr[1]+"\t"+determine_genotype("\t".join(arr[5:9]))

output.close()
output_log.close()
        
