import vcf
import pandas as pd
import sys
import os

# SAMPLES = ["MPNST_0001_Pb_R_Control_1001", 
# "MPNST_0002_Pb_R_Control_1508", 
# "MPNST_0003_Pb_R_Control_1543", 
# "MPNST_0004_Pb_R_Control_1547", 
# "MPNST_0005_Pb_R_Control_1782", 
# "MPNST_0006_Pb_R_Control_1832", 
# "MPNST_0007_Pb_R_Control_2130", 
# "MPNST_0008_Pb_R_Control_2138", 
# "MPNST_0009_Pb_R_Control_2178", 
# "MPNST_0010_Pb_R_Control_2179", 
# "MPNST_0011_Pb_R_Control_2215", 
# "MPNST_0012_Pb_R_Control_2314", 
# "MPNST_0013_Pb_R_Control_2414", 
# "MPNST_0014_Pb_R_Control_2474", 
# "MPNST_0015_Pb_R_Control_2491", 
# "MPNST_0016_Pb_R_Control_2797", 
# "MPNST_0017_Pb_R_Control_3048", 
# "MPNST_0018_Pb_R_Control_3175", 
# "MPNST_0019_Pb_R_Control_882", 
# "MPNST_0020_Pb_R_Control_M803"]

apprisISO = "/Volumes/TGL/gsi/databases/appris/appris_data.principal_20190121.txt"
apprisInfo = pd.read_csv(apprisISO, sep = "\t")
transcripts = apprisInfo['TRANSCRIPT'].str.split('.', n = 0, expand = True)[0].tolist()
genes = apprisInfo['HUGO'].tolist()
txDict = {}
for g in list(set(genes)):
	if g not in txDict.keys():
		txDict[g] = apprisInfo.loc[apprisInfo['HUGO'] == g]['TRANSCRIPT'].str.split('.', n = 0, expand = True)[0].tolist()

# ## allowed variant types
Allowed_VarTypes = ["frameshift_variant", "missense_variant", "splice_acceptor", "splice_donor", "start_lost", "stop_gained", "synonymous_variant",
"inframe_insertion", "inframe_deletion", "frameshift_insertion", "frameshift_deletion"]


vcf_file = "/Volumes/gsiprojects/external/zadehglab/MPNST/data/MPNST/EXOME/Seqware_GATKHaplotypeCaller/JointCalling/MPNST_jointGT.annotated.vcf"

vcf_reader = vcf.Reader(open(vcf_file, 'r'))
for record in vcf_reader:
	SAMPLES = [str(s).split(',')[0].split("=")[1] for s in record.samples]
	break

CSQ = "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|VARIANT_CLASS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|SIFT|PolyPhen|HGVS_OFFSET|CLIN_SIG|SOMATIC|PHENO|1000Gp3_AF|CADD_phred|ExAC_AF|FATHMM_pred|GERP++_RS|MutationTaster_pred|Polyphen2_HDIV_pred|REVEL_score|SIFT_pred|clinvar_clnsig|clinvar_golden_stars|genename|gnomAD_exomes_AF|gnomAD_genomes_AF|LoFtool".split("|")
ENTRY = CSQ
writeREC = "\t".join(str(k) for k in ['chrom', 'pos', 'ref', 'alt', 'qual'] + ENTRY + SAMPLES)
outfile = "/Volumes/gsiprojects/external/zadehglab/MPNST/data/MPNST/EXOME/Seqware_GATKHaplotypeCaller/JointCalling/MPNST_jointGT.annotated2.tsv"
outStream = open(outfile, "w")
outStream.write(writeREC + "\n")

# vcf_reader = vcf.Reader(open(vcf_file, 'r'))
for record in vcf_reader:
	# print record
	GT = list()
	DP = list()
	GQ = list()
	AD = list()
	# sample = "NA"
	for s in record.samples:
		# print s
		# print s['sample']
		GT.append(str(s['GT']))
		DP.append(str(s['DP']))
		GQ.append(str(s['GQ']))
		AD.append(str(s['AD']))
	# break
		# sample = s
	# print GT
	sample_deets = [':'.join(i) for i in zip(GT,DP,GQ,AD)]
	INFO = record.INFO['CSQ']
	# print INFO
	# break
	varType = []
	for info in INFO:
		# print info.split("|")[1]
		varType.append(info.split("|")[1])
	availVarTypes = list(set(varType))
	printInf = []
	for av in availVarTypes: # variant type filtering
		if av in Allowed_VarTypes:
			printInf = INFO
	# if record.QUAL < 10: # variant quality < 10 is discarded
	# 	continue
	# QUAL = record.QUAL
		# print printInf
	if not printInf:
		continue
	for inf in printInf:
		# print inf
		tx = inf.split("|")[6]
		# print tx
		genename = inf.split("|")[3]
		# print genename
		if genename not in txDict.keys():
			continue
		if tx in txDict[genename]:
			chrom = record.CHROM
			pos = record.POS
			ref = record.REF
			alt = ",".join([str(r) for r in record.ALT])
			qual = record.QUAL
			infoDeets = inf
			# gt = GT
			# # dp = DP
			# gq = GQ
			# # ad = AD
			infoMap = dict(zip(CSQ,infoDeets.split("|")))
			# print infoMap
			# print infoMap['Consequence']
			if infoMap['Consequence'] not in Allowed_VarTypes:
				continue # variant type filtering
			infoLog = []
			for csq in ENTRY:
				infoLog.append(infoMap[csq])
				# print csq,"=",infoMap[csq]
			# print "\t".join(str(k) for k in [chrom, pos, ref, alt, qual, gt, gq, dp] + infoLog + [sample])

			writeREC = "\t".join(str(k) for k in [chrom, pos, ref, alt, qual] + infoLog + sample_deets)
			outStream.write(writeREC + "\n")
	# break
outStream.close()
