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


parser = argparse.ArgumentParser()
parser.add_argument("--nirvpath", required=True, help="Path to nirvana files - all nirvana files should be in this folder (no subfolders will be scanned)")
parser.add_argument("--passonly", required=False, default="1", help="If include non-PASS variants in joint table (with all sample AFs). Set to 1 if only PASS, 0 if all variants (with filters)")
parser.add_argument("--outdir", required=False, default="./outdir", help="Output directory")

args = parser.parse_args()
nirvpath= args.nirvpath
passonly = args.passonly
outdir = args.outdir

# Nirvana dir check
if os.path.exists(nirvpath):
    print("- Nirvana input path: %s" % nirvpath )
else:
    sys.exit("> Error! Nirvana input path (%s) does not exist!")

# Nirvana dir check
if os.path.exists(outdir):
    print("- Outdir: %s" % outdir)
else:
    os.mkdir(outdir)
    print("- Outdir: %s (created)" % outdir)

# Read all vcf nirvana files
nirvs = [nirvpath + "/" + f for f in listdir(nirvpath) if isfile(join(nirvpath, f))]

if len(nirvs) == 0:
    sys.exit("> Error: No files found in nirvanapath (%s)" % nirvpath)
else:
    print("- Nirvana files detected:")
    print(nirvs)

print("PASS variants only in joint table: ")
if passonly=="1":
    print("- Yes")
else:
    print("- No")

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

# Dictionary to store all pd.dataframes (vcf tables)
vcf_dict={}
# scan all nirvana files and parse
for nirv in nirvs:
    with gzip.open(nirv, 'r') as fin:        # 4. gzip
        json_bytes = fin.read()
    json_str = json_bytes.decode('utf-8')            # 2. string (i.e. JSON)
    data = json.loads(json_str)                      # 1. data
    sample=data['header']['samples'][0]
    print("Sample: %s" % sample)
    outfile=outdir + "/%s.snv.annotation.nirvana.table.csv" % sample
    print(outfile)
    vcf_dict[sample] = pd.DataFrame(columns=["var_id","chr","position","ref_allele","alt_allele",
                                  "globalMinorAllele","globalMinorAlleleFreq",
                                  "gene_symbol","gene_consequence","variantType","aminoAcids",
                                  "non_coding","genotype","AF","total_depth",
                                  "alt_allele_depth","ref_allele_depth","SomaticQuality",
                                  "filterVCdragen","clinVar_significance","clinVar_phenotype",
                                  "hgvsc","hgvsp","dbsnp","phylopScore","primateAI",
                                  "regulatory_region","revel","polyPhenPreds",
                                  "polyPhenScores","siftPreds","siftScores","codons","exons"])
    pd.set_option('display.max_columns', None)
    pd.set_option('display.max_rows', None)
    for pos in tqdm.tqdm(data['positions'][0:500]):
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
                    #print(globallMin)
                    #print(globallFreq)
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
                    varType=var["variantType"]
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
                if "transcripts" in var:
                    #print("---new var---")
                    btypes=[]
                    conseqs=[]
                    bt_cs=[]
                    protein_codings=[]
                    hgvscs=[]
                    hgvsps=[]
                    for tr in var["transcripts"]:
                        if "bioType" in tr:
                            btype=tr["bioType"]
                            btypes.append(btype)
                            if "consequence" in tr:
                                conseq=tr["consequence"]
                                conseqs.append(conseq)
                            if "protein" in btype:
                                if "hgnc" in tr:
                                    hgnc=tr["hgnc"]
                                if "transcript" in tr:
                                    trid=tr["transcript"]
                                if "aminoAcids" in tr:
                                    protein_coding=str(conseq[0]) + ":" + hgnc + ":" + tr["aminoAcids"]
                                else:
                                    protein_coding=str(conseq[0]) + ":" + hgnc + ":" + ""
                                #print(str(conseq[0]) + ":" + hgnc)
                                protein_codings.append(protein_coding)
                                if "hgvsc" in tr:
                                    hgvsc=tr["hgvsc"].split(":")[1]
                                    hgvscs.append(hgvsc)
                                if "hgvsp" in tr:
                                    hgvsp=tr["hgvsp"].split(":")[1]
                                    hgvsps.append(hgvsp)
                                if "polyPhenPrediction" in tr:
                                    polyPhenPred=tr["polyPhenPrediction"]
                                    polyPhenPreds.append(polyPhenPred)
                                if "polyPhenScore" in tr:
                                    polyPhenScore=tr["polyPhenScore"]
                                    polyPhenScores.append(polyPhenScore)
                                if "siftPrediction" in tr:
                                    sift=tr["siftPrediction"]
                                    siftPreds.append(sift)
                                if "siftScore" in tr:
                                    sift=tr["siftScore"]
                                    siftScores.append(sift)
                                if "codons" in tr:
                                    codon=tr["codons"]
                                    codons.append(codon)
                                if "exons" in tr:
                                    exon=tr["exons"]
                                    exons.append(exon)
                            else:
                                bt_cs.append(str(btype) + ":" + str(conseq[0]) + "")
                if len(bt_cs)>0:
                    bt_cs=f7(bt_cs)
                else:
                    bt_cs=""
                if len(protein_codings)>0:
                    pr_coding=f7(protein_codings)
                else:
                    pr_coding=""
                if len(polyPhenScores)>0:
                    polyPhenScores=max(polyPhenScores)
                if len(polyPhenPreds)>0:
                    polyPhenPreds=f7(polyPhenPreds)
                if len(siftScores)>0:
                    siftScores=max(siftScores)
                    if siftScores==1:
                        siftScores=1.0
                    if siftScores==0:
                        siftScores=0.0
                if len(siftPreds)>0:
                    siftPreds=f7(siftPreds)
                if len(codons)>0:
                    codons=f7(codons)
                if len(exons)>0:
                    exons=f7(exons)
        gene=""
        conseqPr=""
        aminoacid=""
        for pc in pr_coding:
            gene=pc.split(":")[1]
            conseqPr=pc.split(":")[0]
            aminoacid=pc.split(":")[2]
        if len(alleleDepths)<2:
            alleleDepths=["na","na"]

        vcf_dict[sample] = vcf_dict[sample].append({"var_id":vid,
                                 "chr": chrom,
                                 "position": position,
                                 "variantType": varType,
                                 "gene_symbol": gene,
                                 "gene_consequence": conseqPr,
                                 "aminoAcids": aminoacid,
                                 "genotype": gt,
                                 "AF": varfreq,
                                  "ref_allele": refallele,
                                  "alt_allele": altallele[0],
                                  "primateAI":primateAI,
                                 "total_depth": totalDepth,
                                 "alt_allele_depth": alleleDepths[1],
                                 "ref_allele_depth": alleleDepths[0],
                                 "filterVCdragen": str(filt).replace("[","").replace("]","").replace("'",""),
                                 "SomaticQuality": sq,
                                 "non_coding":str(bt_cs).replace("[","").replace("]","").replace("'",""),
                                 "globalMinorAllele": globalMin,
                                 "globalMinorAlleleFreq": globalFreq,
                                 "hgvsp": str(f7(hgvsps)).replace("[","").replace("]","").replace("'",""),
                                 "hgvsc": str(f7(hgvscs)).replace("[","").replace("]","").replace("'",""),
                                 "clinVar_significance": str(clinset).replace("[","").replace("]","").replace("'",""),
                                 "clinVar_phenotype": str(clphenset).replace("[","").replace("]","").replace("'",""),
                                 "phylopScore": phylop,
                                 "dbsnp":str(dbsnp).replace("[","").replace("]","").replace("'","").replace("{","").replace("}",""),
                                 "regulatory_region":str(regulatories).replace("[","").replace("]","").replace("'",""),
                                 "revel":revel,
                                 "polyPhenPreds":str(polyPhenPreds).replace("[","").replace("]","").replace("'",""),
                                  "polyPhenScores":str(polyPhenScores).replace("[","").replace("]","").replace("'",""),
                                  "siftPreds":str(siftPreds).replace("[","").replace("]","").replace("'",""),
                                  "siftScores":str(siftScores).replace("[","").replace("]","").replace("'",""),
                                  "codons":str(codons).replace("[","").replace("]","").replace("'",""),
                                  "exons":str(exons).replace("[","").replace("]","").replace("'","")
                                 },
                      ignore_index=True)
        vcf_dict[sample]=vcf_dict[sample].set_index(vcf_dict[sample].var_id)
    vcf_dict[sample].to_csv(outfile,index=False)
    print("Table written to %s" % outfile)


# In[39]:


vids=[]
for samp in vcf_dict.keys():
    for idx,row in vcf_dict[samp].iterrows():
        if passonly==0:
            vids.append(idx)
        elif passonly==1:
            if vcf_dict[samp].loc[idx]["filterVCdragen"]=="PASS":
                vids.append(idx)

collist=["var_id","chr","position","ref_allele","alt_allele","globalMinorAllele","globalMinorAlleleFreq",
                               "gene_symbol","gene_consequence","variantType","aminoAcids",
                               "non_coding","clinVar_significance","clinVar_phenotype",
                               "hgvsc","hgvsp","dbsnp","phylopScore","primateAI",
                               "regulatory_region","revel","polyPhenPreds",
                                 "polyPhenScores","siftPreds","siftScores","codons","exons"]
pdjoint=pd.DataFrame(columns=collist)
pdjoint["var_id"]=list(set(vids))
pdjoint=pdjoint.set_index(pdjoint["var_id"])
pdjoint.index.name="index"

# Add sample AFs
for sample in vcf_dict.keys():
    af=sample+"_AF"
    pdjoint[af]=""
    for idx,row in vcf_dict[sample].iterrows():
        if passonly==1:
            if vcf_dict[sample].loc[idx]["filterVCdragen"]=="PASS":
                pdjoint.loc[idx,collist]=row
                pdjoint.loc[idx][af]=row["AF"]
        else:
            pdjoint.loc[idx,collist]=row
            pdjoint.loc[idx][af]=row["AF"]

if passonly==0:
    for sample in vcf_dict.keys():
        af=sample+"_AF"
        flag=sample+"_dragenFlag"
        pdjoint[flag]=""
        for idx,row in vcf_dict[sample].iterrows():
            if row["filterVCdragen"]!="PASS":
                pdjoint.loc[idx][af]=str(row["AF"])+" (f)"
                pdjoint.loc[idx][flag]=str(row["filterVCdragen"])
    jointout=outdir + "/joint_variant_AF_table.nonFilter.csv"
else:
    jointout=outdir + "/joint_variant_AF_table.PASS.csv"
pdjoint.to_csv(jointout,index=False)


# In[ ]:
