#!/usr/bin/env python
import sys
import os
from os import listdir
from os.path import isfile, join
import gzip
import json
import tqdm
import pandas as pd
import argparse


parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--indir", required=True, help="Path to input folder - all input files should be in this folder (no subfolders will be scanned)")
parser.add_argument("--outdir", required=False, default="./outdir", help="Output directory")
parser.add_argument("--mode", required=False, default="single", help="'single': create one parsed table pr input nirvana json.gz. 'joint': create one joint table from individual tables created in 'single' mode ")

args = parser.parse_args()
indir= args.indir
outdir = args.outdir
mode = args.mode

# Nirvana dir check
if os.path.exists(indir):
    print("- Nirvana input path: %s" % indir )
else:
    sys.exit("> Error! Nirvana input path (%s) does not exist!")

# Nirvana dir check
if os.path.exists(outdir):
    print("- Outdir: %s" % outdir)
else:
    os.mkdir(outdir)
    print("- Outdir: %s (created)" % outdir)

# Function to replace duplicate strings
# by alphanumeric strings to make all
# strings in the array unique
def replaceDuplicates(names):

    # Store the frequency of strings
    hash = {}

    # Iterate over the array
    for i in range(0, len(names)):

        # For the first occurrence,
        # update the frequency count
        if names[i] not in hash:
            hash[names[i]] = 1

        # Otherwise
        else:
            count = hash[names[i]]
            hash[names[i]] += 1

            # Append frequency count
            # to end of the string
            names[i] = names[i] + "_" + str(count)

    return(names)

print("-Mode: %s" % mode )

if mode == "single":
    # Read all vcf nirvana files
    nirvs = [indir + "/" + f for f in listdir(indir) if isfile(join(indir, f))]

    if len(nirvs) == 0:
        sys.exit("> Error: No files found in nirvanapath (%s)" % indir)
    else:
        print("- Nirvana files detected:")
        print(nirvs)


    # Function to take unique values from list, and conserve order
    def f7(seq):
        seen = set()
        seen_add = seen.add
        return [x for x in seq if not (x in seen or seen_add(x))]

    # Set joint table structure
    pdjoint=pd.DataFrame(columns=["var_id","chr","position","ref_allele","alt_allele",
                                   "globalMinorAllele","globalMinorAlleleFreq",
                                   "gene_symbol","gene_consequence","variantType","aminoAcids",
                                   "non_coding","clinVar_significance","clinVar_phenotype",
                                   "hgvsc","hgvsp","dbsnp","phylopScore","primateAI",
                                   "regulatory_region","revel","polyPhenPreds",
                                   "polyPhenScores","siftPreds","siftScores","codons","exons"])

    # scan all nirvana files and parse
    for nirv in nirvs[0:3]:
        with gzip.open(nirv, 'r') as fin:        # 4. gzip
            json_bytes = fin.read()
        json_str = json_bytes.decode('utf-8')            # 2. string (i.e. JSON)
        data = json.loads(json_str)                      # 1. data
        sample=data['header']['samples'][0]
        print("- Sample: %s" % sample)
        outfile=outdir + "/%s.snv.annotation.nirvana.table.csv" % sample
        print("- Outfile: %s " % outfile)
        # anontation per transcript
        annot=["source","bioType","consequence","hgnc",
               "transcript","Ensembl_transcript","RefSeq_transcript",
               "aminoAcids","hgvsc","hgvsp",
               "polyPhenPred","polyPhenScore",
               "siftPred","siftScore",
               "codons","exons"]
        pd.set_option('display.max_columns', None)
        pd.set_option('display.max_rows', None)
        vcf_dict = pd.DataFrame(columns=["var_id","chr","position","ref_allele","alt_allele","variantType",
                                                 "genotype","AF","total_depth","alt_allele_depth","ref_allele_depth",
                                                 "hgnc","Ensembl_transcript","RefSeq_transcript",
                                                 "bioType",
                                                 "globalMinorAllele","globalMinorAlleleFreq",
                                                 "filterVCdragen","SomaticQuality","clinVar_significance","clinVar_phenotype",
                                                 "regulatory_region","revel",
                                                 "dbsnp","phylopScore","primateAI",
                                                 "transcript_consequence","aminoAcids",
                                                 "hgvsc","hgvsp",
                                                 "polyPhenScore","polyPhenPred",
                                                 "siftPred","siftScore",
                                                 "codons","exons","introns"])
        pd.set_option('display.max_columns', None)
        pd.set_option('display.max_rows', None)
        for pos in tqdm.tqdm(data['positions'][0:50]):
            verb=0
            vid=""
            chrom=""
            position=""
            varType=""
            hgvsc=""
            conseqPr=""
            btypes=[]
            conseqs=[]
            bt_cs=[]
            protein_codings=[]
            hgvscs=[]
            hgvsps=[]
            bt_cs=""
            pr_coding=""
            gt=""
            filt=""
            varfreq=""
            sq=""
            protein_coding=""
            totalDepth=""
            alleleDepths=[]
            clinset=""
            phylop=""
            dbsnp=""
            regulatories=""
            revel=""
            clphenset=""
            primateAI=""
            polyPhenPreds=[]
            polyPhenScores=[]
            siftPreds=[]
            siftScores=[]
            codons=[]
            exons=[]
            refallele=""
            altallele=""
            filt = pos['filters']
            samples=pos["samples"][0]
            altallele=pos["altAlleles"]
            refallele=pos["refAllele"]
            chrom=pos["chromosome"]
            position=pos["position"]
            if "variantFrequencies" in samples:
                varfreq=samples["variantFrequencies"][0]
            if "genotype" in samples:
                gt=samples["genotype"]
            if "somaticQuality" in samples:
                sq=samples["somaticQuality"]
            if "totalDepth" in samples:
                totalDepth=samples["totalDepth"]
            if "alleleDepths" in samples:
                alleleDepths=samples["alleleDepths"]
            if "variants" in pos:
                for var in pos['variants']:
                    #print("---- new variant -----")
                    if "vid" in var:
                        vid=var['vid']
                    if "dbsnp" in var:
                        dbsnp=var['dbsnp']
                    if "globalAllele" in var:
                        globalMin=var["globalAllele"]["globalMinorAllele"]
                        globalFreq=var["globalAllele"]["globalMinorAlleleFrequency"]
                    if "phylopScore" in var:
                        phylop=var["phylopScore"]
                        #print(phylop)
                    if "primateAI" in var:
                        primateAI=var["primateAI"][0]["scorePercentile"]
                        #print("primateIAscore: " + str(primateAI))
                    if "regulatoryRegions" in var:
                        regulatories=[]
                        for reg in var["regulatoryRegions"]:
                            regulatories.append(reg["type"])
                    if "revel" in var:
                        revel=var["revel"]["score"]
                        #print(revel)
                    if "variantType" in var:
                        varType=str(var["variantType"])
                    # CLINVAR
                    if "clinvar" in var:
                        clinvars=[]
                        clinvarphenos=[]
                        for clin in var["clinvar"]:
                            clid=clin["id"]
                            clsign=clin["significance"]
                            #print(clsign)
                            if "patho" in clsign[0] or "risk" in clsign or "assoc" in clsign or "drug" in clsign or "protect" in clsign or "affect" in clsign:
                                clinvars.insert(0,str(clsign[0]))
                                verb=1
                                #print("verb:" + str(verb))
                            else:
                                clinvars.append(clsign[0])
                            if "phenotypes" in clin:
                                clphen=clin["phenotypes"]
                                for phen in clphen:
                                    clinvarphenos.append(phen)
                        clinset=f7(clinvars)
                        clphenset=f7(clinvarphenos)
                        if "pathogenic" in clinset:
                            clinset.insert(0, clinset.pop(clinset.index("pathogenic")))
                    # Transcript
                    canons=pd.DataFrame(columns=annot)
                    main_canons=pd.DataFrame(columns=annot)
                    canot={}
                    if "transcripts" in var:
                        for tr in var["transcripts"]:
                            # Get canonicals
                            if "isCanonical" in tr:
                                canot={}
                                if tr["isCanonical"]:
                                    if "source" in tr:
                                        canot["source"]=tr["source"]
                                    if "bioType" in tr:
                                        btype=tr["bioType"]
                                        canot["bioType"]=btype
                                    if "consequence" in tr:
                                        conseq=tr["consequence"]
                                        canot["consequence"]=conseq
                                    if "hgnc" in tr:
                                        canot["hgnc"]=tr["hgnc"]
                                    if "transcript" in tr:
                                        canot["transcript"]=tr["transcript"]
                                    if "aminoAcids" in tr:
                                        canot["aminoAcids"]=tr["aminoAcids"]
                                    if "hgvsc" in tr:
                                        canot["hgvsc"]=tr["hgvsc"]
                                    if "hgvsp" in tr:
                                        canot["hgvsp"]=tr["hgvsp"]
                                    if "polyPhenPrediction" in tr:
                                        canot["polyPhenPred"]=tr["polyPhenPrediction"]
                                    if "polyPhenScore" in tr:
                                        canot["polyPhenScore"]=tr["polyPhenScore"]
                                    if "siftPrediction" in tr:
                                        canot["siftPred"]=tr["siftPrediction"]
                                    if "siftScore" in tr:
                                        canot["siftScore"]=tr["siftScore"]
                                    if "codons" in tr:
                                        canot["codons"]=tr["codons"]
                                    if "exons" in tr:
                                        canot["exons"]=tr["exons"]
                                    if "introns" in tr:
                                        canot["introns"]=tr["introns"]
                                    else:
                                        canot["introns"]=""
                                    canons=canons.append(canot,ignore_index=True)
                        # Go through all canonical transcripts
                        # remove Ensembl/RefSeq redundancies (Ensembl have priority if hgnc is identical between the two)
                        # write one row in output file pr canonical transcript (main_canons)
                        # 1. get all unique hgnc ids
                        hgnc=set(list(canons["hgnc"]))
                        # 2. iterate through all hgnc ids
                        for h in hgnc:
                            main=pd.DataFrame(columns=annot)
                            # 3. Get all canonical transcripts with that id
                            curr=canons[canons["hgnc"]==h]
                            # if current is Ensembl, set as main
                            if "Ensembl" in list(curr["source"]):
                                # If also there is a RefSeq transcript with same hgnc ID
                                # add refseq_transcript id to variable
                                if "RefSeq" in list(curr["source"]):
                                    refseq_transcript=str(curr[curr["source"]=="RefSeq"]["transcript"][curr[curr["source"]=="RefSeq"]["transcript"].index[0]])
                                else:
                                    refseq_transcript=""
                                main=curr[curr["source"]=="Ensembl"]
                                idx=main.index[0]
                                main.loc[idx,"RefSeq_transcript"]=refseq_transcript
                                main.loc[idx,"Ensembl_transcript"]=str(curr[curr["source"]=="Ensembl"]["transcript"][curr[curr["source"]=="Ensembl"]["transcript"].index[0]])
                            elif "RefSeq" in list(curr["source"]):
                                main=curr[curr["source"]=="RefSeq"]
                                idx=main.index[0]
                                refseq_transcript=str(curr[curr["source"]=="RefSeq"]["transcript"][curr[curr["source"]=="RefSeq"]["transcript"].index[0]])
                                main.loc[idx,"RefSeq_transcript"]=refseq_transcript
                                main.loc[idx,"Ensembl_transcript"]=""
                            if len(main)>0:
                                main_canons=main_canons.append(main,ignore_index=True)
                        for idx,row in main_canons.iterrows():
                            hgnc_curr=main_canons.loc[idx,"hgnc"]
                            if len(alleleDepths)<2:
                                alleleDepths=["na","na"]
                            vcf_dict = vcf_dict.append({"var_id":vid,
                                                 "chr": chrom,
                                                 "position": position,
                                                 "variantType": varType,
                                                 "genotype": gt,
                                                 "AF": varfreq,

                                                 "Ensembl_transcript": row["Ensembl_transcript"],
                                                 "RefSeq_transcript": row["RefSeq_transcript"],
                                                 "bioType": row["bioType"],
                                                 "hgnc": row["hgnc"],
                                                 "transcript_consequence": str(row["consequence"]).replace("[","").replace("]","").replace("'",""),
                                                 "aminoAcids": row["aminoAcids"],
                                                 "hgvsp": row["hgvsp"],
                                                 "hgvsc": row["hgvsc"],
                                                 "polyPhenPred": row["polyPhenPred"],
                                                 "polyPhenScore":row["polyPhenScore"],
                                                 "siftPred":row["siftPred"],
                                                 "siftScore":row["siftScore"],
                                                 "codons":row["codons"],
                                                 "exons":row["exons"],
                                                 "introns": row["introns"],

                                                 "ref_allele": refallele,
                                                 "alt_allele": altallele[0],
                                                 "primateAI":primateAI,
                                                 "total_depth": totalDepth,
                                                 "alt_allele_depth": alleleDepths[1],
                                                 "ref_allele_depth": alleleDepths[0],
                                                 "filterVCdragen": str(filt).replace("[","").replace("]","").replace("'",""),
                                                 "SomaticQuality": sq,
                                                 "globalMinorAllele": globalMin,
                                                 "globalMinorAlleleFreq": globalFreq,
                                                 "clinVar_significance": str(clinset).replace("[","").replace("]","").replace("'",""),
                                                 "clinVar_phenotype": str(clphenset).replace("[","").replace("]","").replace("'",""),
                                                 "phylopScore": phylop,
                                                 "dbsnp":str(dbsnp).replace("[","").replace("]","").replace("'","").replace("{","").replace("}",""),
                                                 "regulatory_region":str(regulatories).replace("[","").replace("]","").replace("'",""),
                                                 "revel":revel,
                                                 },ignore_index=True)
                            vcf_dict=vcf_dict.set_index(replaceDuplicates(vcf_dict["var_id"]+"_"+vcf_dict["hgnc"]))
        vcf_dict.to_csv(outfile,index=True)
        print("Table with written to %s" % outfile)

##################################################################################################################
## JOINT TABLES
##################################################################################################################



if mode == "joint":
    # Get output table from each samples
    samples=os.listdir(indir)
    # remove potential joint table files
    samples = [x for x in samples if not x.startswith("joint")]
    print("-- Fetching individual tables to create joint table: ")
    for samp in samples:
        print("- %s" % samp)
    ## PASS variants ONLY
    indxs=[]
    vids=[]
    passonly=1
    for samp in samples:
        vcf_dict=pd.read_csv(outdir+"/"+samp,index_col=0)
        for idx,row in vcf_dict.iterrows():
            if passonly==0:
                indxs.append(idx)
                vids.append(vcf_dict.loc[idx]["var_id"])
            elif passonly==1:
                if vcf_dict.loc[idx]["filterVCdragen"]=="PASS":
                    vids.append(vcf_dict.loc[idx]["var_id"])
                    indxs.append(idx)
    # Get all indexes

    collist=["var_id","chr","position","ref_allele","alt_allele","variantType",
                     "hgnc","Ensembl_transcript","RefSeq_transcript",
                     "bioType","globalMinorAllele","globalMinorAlleleFreq",
                     "clinVar_significance","clinVar_phenotype",
                     "regulatory_region","revel",
                     "dbsnp","phylopScore","primateAI",
                     "transcript_consequence","aminoAcids",
                     "hgvsc","hgvsp",
                     "polyPhenScore","polyPhenPred",
                     "siftPred","siftScore",
                     "codons","exons","introns"]

    pdjoint=pd.DataFrame(columns=collist)
    pdjoint["indx"]=list(set(indxs))
    pdjoint["var_id"]=[vids.split("_")[0] for vids in list(set(indxs))]
    pdjoint=pdjoint.set_index(pdjoint["indx"])
    pdjoint=pdjoint.drop(columns="indx")

    # Add sample AFs
    for sample in samples:
        af=sample+"_AF"
        pdjoint[af]=""
        vcf_dict=pd.read_csv(outdir+"/"+sample,index_col=0)
        for idx,row in vcf_dict.iterrows():
            if passonly==1:
                if vcf_dict.loc[idx]["filterVCdragen"]=="PASS":
                    # Make sure not to overvrite with empty rows
                    if pdjoint.loc[idx,"chr"]!=pdjoint.loc[idx,"chr"]:
                        pdjoint.loc[idx,collist]=row
                        pdjoint.loc[idx][af]=row["AF"]
                    else:
                        pdjoint.loc[idx][af]=row["AF"]
            else:
                pdjoint.loc[idx,collist]=row
                pdjoint.loc[idx][af]=row["AF"]

    if passonly==0:
        for sample in samples:
            af=sample+"_AF"
            flag=sample+"_dragenFlag"
            pdjoint[flag]=""
            vcf_dict=pd.read_csv(outdir+"/"+sample,index_col=0)
            for idx,row in vcf_dict.iterrows():
                if row["filterVCdragen"]!="PASS":
                    pdjoint.loc[idx][af]=str(row["AF"])+" (f)"
                    pdjoint.loc[idx][flag]=str(row["filterVCdragen"])
        jointout=outdir + "/joint_variant_AF_table.FILTER.csv"
    else:
        jointout=outdir + "/joint_variant_AF_table.PASS.csv"
    pdjoint=pdjoint.sort_values(by=["chr","position"])
    pdjoint.to_csv(jointout,index=False)
    print("Joint PASS Table with written to %s" % jointout)









    ## FILTER variants
    indxs=[]
    vids=[]
    passonly=0
    for samp in samples:
        vcf_dict=pd.read_csv(outdir+"/"+samp,index_col=0)
        for idx,row in vcf_dict.iterrows():
            if passonly==0:
                indxs.append(idx)
                vids.append(vcf_dict.loc[idx]["var_id"])
            elif passonly==1:
                if vcf_dict.loc[idx]["filterVCdragen"]=="PASS":
                    vids.append(vcf_dict.loc[idx]["var_id"])
                    indxs.append(idx)
    # Get all indexes

    collist=["var_id","chr","position","ref_allele","alt_allele","variantType",
                     "hgnc","Ensembl_transcript","RefSeq_transcript",
                     "bioType","globalMinorAllele","globalMinorAlleleFreq",
                     "clinVar_significance","clinVar_phenotype",
                     "regulatory_region","revel",
                     "dbsnp","phylopScore","primateAI",
                     "transcript_consequence","aminoAcids",
                     "hgvsc","hgvsp",
                     "polyPhenScore","polyPhenPred",
                     "siftPred","siftScore",
                     "codons","exons","introns"]

    pdjoint=pd.DataFrame(columns=collist)
    pdjoint["indx"]=list(set(indxs))
    pdjoint["var_id"]=[vids.split("_")[0] for vids in list(set(indxs))]
    pdjoint=pdjoint.set_index(pdjoint["indx"])
    pdjoint=pdjoint.drop(columns="indx")

    # Add sample AFs
    for sample in samples:
        af=sample+"_AF"
        pdjoint[af]=""
        vcf_dict=pd.read_csv(outdir+"/"+sample,index_col=0)
        for idx,row in vcf_dict.iterrows():
            if passonly==1:
                if vcf_dict.loc[idx]["filterVCdragen"]=="PASS":
                    # Make sure not to overvrite with empty rows
                    if pdjoint.loc[idx,"chr"]!=pdjoint.loc[idx,"chr"]:
                        pdjoint.loc[idx,collist]=row
                        pdjoint.loc[idx][af]=row["AF"]
                    else:
                        pdjoint.loc[idx][af]=row["AF"]
            else:
                pdjoint.loc[idx,collist]=row
                pdjoint.loc[idx][af]=row["AF"]

    if passonly==0:
        for sample in samples:
            af=sample+"_AF"
            flag=sample+"_dragenFlag"
            pdjoint[flag]=""
            vcf_dict=pd.read_csv(outdir+"/"+sample,index_col=0)
            for idx,row in vcf_dict.iterrows():
                if row["filterVCdragen"]!="PASS":
                    pdjoint.loc[idx][af]=str(row["AF"])+" (f)"
                    pdjoint.loc[idx][flag]=str(row["filterVCdragen"])
        jointout=outdir + "/joint_variant_AF_table.FILTER.csv"
    else:
        jointout=outdir + "/joint_variant_AF_table.PASS.csv"
    pdjoint=pdjoint.sort_values(by=["chr","position"])
    pdjoint.to_csv(jointout,index=False)
    print("Joint FILTER Table with written to %s" % jointout)
