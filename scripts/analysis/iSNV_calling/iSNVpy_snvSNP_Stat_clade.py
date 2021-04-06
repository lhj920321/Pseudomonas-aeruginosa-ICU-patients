
#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys,os
from optparse import OptionParser
import matplotlib.pyplot as plt
import numpy as np 



def iSNV_SNP_tableRead(iSNV_SNP_table):
	iSNVSNPtablels = open(iSNV_SNP_table,'r').readlines()
	HeaderLs = iSNVSNPtablels[0].split("\n")[0].split("\t")
	HeaderL = []
	for HeaderID in HeaderLs:
		if HeaderID != "#Posi" and HeaderID !=  "snv/SNP":
			
			HeaderL.append(HeaderID.split(".sort")[0])
		else:
			HeaderL.append(HeaderID)

	lie = 0
	lieD = {}
	for Header in HeaderL:
		lieD[Header] = lie
		lie += 1
	lsD = {}
	for iSNVSNPtablel in iSNVSNPtablels:
		if "#" not in iSNVSNPtablel and iSNVSNPtablel != "\n" :
			lL = iSNVSNPtablel.split("\n")[0].split("\t")
			Pos = int(lL[0])
			lsD[Pos] = lL
	return lieD,lsD



def samp_Freq_Dic(samples,lieD,lsD):
	samp_FreqDic = {}
	sampNA_PosiDic = {}
	print(lieD)
	for sample in samples:
		Lie = lieD[sample]
		samp_FreqDic[sample] = {}
		sampNA_PosiDic[sample] = []
		for Posi in lsD:
			Freq = lsD[Posi][Lie]
			if Freq != "NA" :
				if Freq != "NO":				
					Freq = float(lsD[Posi][Lie])
					if Freq >= 0.05:
						samp_FreqDic[sample][Posi] = Freq
				else:
					Freq = 0
			else:
				sampNA_PosiDic[sample].append(Posi)
	#print(samp_FreqDic.keys())
	return samp_FreqDic,sampNA_PosiDic



def Person_NA_Posi(samples,sampNA_PosiDic):
	personNAPosi = {}
	for sample in samples:
		sampleNAPosi = sampNA_PosiDic[sample]
		for posi in sampleNAPosi:
			if posi not in personNAPosi:
				personNAPosi[posi] = ''
	return personNAPosi


def Person_commonSNP_Posi(samples,samp_FreqDic,lsD):
	personCommonPosi = {}
	
	for posi in lsD.keys():
		count = 0
		for sample in samples:
			sampleFreqDic = samp_FreqDic[sample]
			if posi in sampleFreqDic:
				FREQ = sampleFreqDic[posi]
				if FREQ >= snv_MaxFreq:
					count += 1
		if count == len(samples):
			personCommonPosi[posi] = ''
			
	return personCommonPosi



def SNPPosi_withsnv_Posi(samples,samp_FreqDic,lsD,personCommonPosi):
	SNPPosi_withsnv = {}
	for posi in lsD.keys():
		SNPcount = 0
		snvCount = 0
		if posi not in personCommonPosi:	
			for sample in samples:
				sampleFreqDic = samp_FreqDic[sample]
				if posi in sampleFreqDic:
					FREQ = sampleFreqDic[posi]
					if FREQ >= snv_MaxFreq:
						SNPcount += 1
					if FREQ >= snv_MinFreq and  FREQ < snv_MaxFreq:
						snvCount += 1
		if snvCount >= 1  and SNPcount >= 1 and SNPcount +  snvCount == len(samples):
			SNPPosi_withsnv[posi] = ''
			
	return SNPPosi_withsnv




def readInfo(InfoF,person):
	print(person)
	cladeSampDic = {}
	for InfoFl in open(InfoF).readlines():
		if "#" not in InfoFl:
			InfoFl_Tags = InfoFl.split("\n")[0].split("\t")
			seqID = InfoFl_Tags[18]
			clade = InfoFl_Tags[20]
			#if seqID not in ["P1-S1","P1-S2","P1-S3","P2-S1","P5-S17","P7-S3"]:	
			if person in clade:
				if clade  not in cladeSampDic:
					cladeSampDic[clade] = []
				cladeSampDic[clade].append(seqID)			
		
	return cladeSampDic








def tongjiFREQ(sample,samp_FreqDic,personNAPosi,personCommonPosi,snvSNPcount_Dic,personSNP_Lst,personsnv_Lst,SNPPosi_withsnv,EcoliFiltedPosiDic,repeatRegPosi,effDic):

	FreqDic = samp_FreqDic[sample]
	tongji_Dic = {}
	#print(len(personCommonPosi))
	#print("len common")
	#for XX in range(0,101,step):
		#tongji_Dic[XX] = 0

	snvCount = 0
	SNPCount = 0
	SNPPosi_withsnv_Num = 0
	EcoliHomo_SNPCount = 0
	EcoliHomo_SNPPosi_withsnv_Num = 0
	EcoliHomo_snvCount = 0
	SNPSynNum,SNPMissNum,SNPStopNum,SNP_interNum = 0,0,0,0
	RepeatSNV,RepeatSNP,RepeatSNPPosi_withsnv = 0,0,0

	for Posi in FreqDic:
		if Posi not in personNAPosi and Posi not in personCommonPosi :
			if Posi not in EcoliFiltedPosiDic and Posi not in repeatRegPosi:
				Freq = FreqDic[Posi]
				if Freq >= snv_MinFreq and Freq < snv_MaxFreq:
					snvCount += 1
					if Posi not in personsnv_Lst:
						personsnv_Lst.append(Posi)

				elif Freq >= snv_MaxFreq and Posi not in SNPPosi_withsnv:
					SNPCount += 1
					if "synonymous_variant" in effDic[Posi]:
						SNPSynNum += 1
					if "missense_variant" in effDic[Posi]:
						SNPMissNum += 1
					if "stop_gained" in effDic[Posi]:
						SNPStopNum += 1
					if"upstream_gene_variant" in effDic[Posi] or "downstream_gene_variant" in effDic[Posi]:
						SNP_interNum += 1

					if Posi not in personSNP_Lst:
						personSNP_Lst.append(Posi)
				if Freq >= snv_MaxFreq and Posi in SNPPosi_withsnv:
					SNPPosi_withsnv_Num += 1
			elif Posi in EcoliFiltedPosiDic:
				Freq = FreqDic[Posi]
				if Freq >= snv_MinFreq and Freq < snv_MaxFreq:
					EcoliHomo_snvCount += 1
				elif Freq >= snv_MaxFreq and Posi not in SNPPosi_withsnv:
					EcoliHomo_SNPCount += 1
				if Freq >= snv_MaxFreq and Posi in SNPPosi_withsnv:
					EcoliHomo_SNPPosi_withsnv_Num += 1
			elif Posi in repeatRegPosi :

				Freq = FreqDic[Posi]
				if Freq >= snv_MinFreq and Freq < snv_MaxFreq:
					RepeatSNV += 1
				elif Freq >= snv_MaxFreq and Posi not in SNPPosi_withsnv:
					RepeatSNP += 1
				if Freq >= snv_MaxFreq and Posi in SNPPosi_withsnv:
					RepeatSNPPosi_withsnv += 1




	
	snvSNPcount_Dic[sample] = "\t".join([sample,str(snvCount),str(SNPCount),str(SNPSynNum),str(SNPMissNum),str(SNPStopNum),str(SNP_interNum),str(len(personCommonPosi)),str(SNPPosi_withsnv_Num),\
		str(EcoliHomo_snvCount),str(EcoliHomo_SNPCount),str(EcoliHomo_SNPPosi_withsnv_Num),str(len(EcoliFiltedPosiDic)),str(RepeatSNV),str(RepeatSNP),str(RepeatSNPPosi_withsnv),str(len(repeatRegPosi))])

	return snvSNPcount_Dic,personSNP_Lst,personsnv_Lst





	

def Posi_lines(cladeSamps,lieD,lsD,Cite_Lst,effDic):
	SampLst = []
	for cladeSamp in cladeSamps:
		SampLst.append(cladeSamp)
	outlines = ["\t".join(["chromosome","position in chromosome","Reference nucleotide",\
		"Query nucleotide","Type of Mutation","Level of Mutation effect",\
		"Mutation nucleotide level","Mutation AA level","Gene locus tag(PAO1)","Gene name(PAO1)","no.SNP"] + SampLst)]
	
	for POSI in Cite_Lst:
		lineLst = [effDic[POSI],'']
		SNPno = 0
		for cladeSamp in cladeSamps:
			Samp_col = lieD[cladeSamp]
			Freq = lsD[POSI][Samp_col]
			lineLst.append(Freq)
			if Freq != "NO":
				SNPno += 1
		lineLst[1] = str(SNPno)
		line = "\t".join(lineLst)
		outlines.append(line)
	lines = "\n".join(outlines)
	return lines




def FiltEcoliPosi(FiltEcoliPosiF):
	EcoliFiltedPosiDic = {}
	for FiltEcoliPosiFl in open(FiltEcoliPosiF).readlines():
		if FiltEcoliPosiFl != "\n":
			EcoliFiltedPosiDic[FiltEcoliPosiFl.split("\n")[0].split("\t")[0]] = ''
	return EcoliFiltedPosiDic


def FiltRefRepeatRegions1(RepeatMaskedGnm):
	for RepeatMaskedGnml in open(RepeatMaskedGnm).readlines():
		if ">" in RepeatMaskedGnml:
			gnm = ''
		else:
			gnm += RepeatMaskedGnml.split("\n")[0]
	repeatRegPosi={}
	for gnmPosiIndex in range(0,len(gnm)):
		if gnm[gnmPosiIndex] == "N":
			repeatRegPosi[gnmPosiIndex+1] = ''
	return repeatRegPosi



def FiltRefRepeatRegions(RepeatPosiF):
	repeatRegPosi = {}
	for RepeatPosiFl in open(RepeatPosiF).readlines():
		repeatRegPosi[int(RepeatPosiFl.split("\t")[0])] = ''
	return repeatRegPosi



def readSnpEff(snpEffPath):
	effDic = {}
	fileList = os.listdir(snpEffPath)
	files = []
	for file in fileList:
		if os.path.splitext(file)[1] == ".vcf":
			files.append(snpEffPath + "/" + file)
	for sampleSnpAnnF in files:
		for sampleSnpAnnFl in open(sampleSnpAnnF).readlines():
			if "#" not in sampleSnpAnnFl:
				outl = []
				Tags = sampleSnpAnnFl.split("\n")[0].split("\t")
				RefID = Tags[0]
				Posi = int(Tags[1])
				ref = Tags[3]
				alt = Tags[4]
				infos = Tags[7]
				outl = [RefID,str(Posi),ref,alt]
				
				Eff,chengdu,aaChange,geneTag,gene = "","","","",""
				if Posi not in effDic:
					if "ANN=" in infos:
						InfoTags = infos.split("ANN=")[1].split("|")
						Eff = InfoTags[1]
						chengdu = InfoTags[2]
						gene = InfoTags[3]
						geneTag = InfoTags[4]
						ntChange = InfoTags[9].split("c.")[1]
						aaChange = InfoTags[10]
						outl.append(Eff)
						outl.append(chengdu)
						outl.append(ntChange)						
						outl.append(aaChange)
						outl.append(geneTag)
						outl.append(gene)
					line = "\t".join(outl)
					effDic[Posi] = line

	return effDic











def main():
	person = iSNV_SNP_table.split("/")[-2]
	if person == "P5_withP2":
		person = "P5"
	lieD,lsD = iSNV_SNP_tableRead(iSNV_SNP_table)
	print(lieD)
	EcoliFiltedPosiDic = FiltEcoliPosi(FiltEcoliPosiF)
	repeatRegPosi = FiltRefRepeatRegions(RepeatMaskedGnm)
	#iSNVtable_lst = list(lieD.keys())
	cladeSampDic = readInfo(InfoF,person)
	print("clade samps :")
	print(cladeSampDic)
	effDic = readSnpEff(snpEffPath)
	for clade in cladeSampDic:
		cladeSamps = cladeSampDic[clade]
		print(cladeSamps)


		samp_FreqDic,sampNA_PosiDic = samp_Freq_Dic(cladeSamps,lieD,lsD)
		personNAPosi = Person_NA_Posi(cladeSamps,sampNA_PosiDic)
		personCommonPosi = Person_commonSNP_Posi(cladeSamps,samp_FreqDic,lsD)
		SNPPosi_withsnv = SNPPosi_withsnv_Posi(cladeSamps,samp_FreqDic,lsD,personCommonPosi)
		
		print(SNPPosi_withsnv)
		
		snvSNPcount_Dic = {}
		personSNP_Lst = []	
		personsnv_Lst = []	
		for sample in cladeSamps:
			#print(sample)
			snvSNPcount_Dic,personSNP_Lst,personsnv_Lst = tongjiFREQ(sample,samp_FreqDic,personNAPosi,personCommonPosi,snvSNPcount_Dic,personSNP_Lst,personsnv_Lst,SNPPosi_withsnv,EcoliFiltedPosiDic,repeatRegPosi,effDic)		

		#person = iSNV_SNP_table.split("/")[-2]
		out_snvSNPCountF = out_P + "/" + clade + ".FiltSNPwithsnv.snvSNPcount.txt"
		out_snvSNPCountF_O=open(out_snvSNPCountF,'a')
		out_snvSNPCountF_O.write("\t".join(["sample","snvCount","SNPCount_noSNV","SNPSynNum","SNPMissNum","SNPStopNum","SNP_interNum","personCommonCitesNum","SNPPosi_withsnv_Num","OtherBacteriaHomo_snvCount","OtherBacteriaHomo_SNPCount","OtherBacteriaHomo_SNPPosi_withsnv_Num","OtherBacteria_FiltedPosiNum","RepeatSNV","RepeatSNP","RepeatSNPPosi_withsnv","RepeatPosis"]) + "\n")
		for sampleID in snvSNPcount_Dic:	
			out_snvSNPCountF_O.write(snvSNPcount_Dic[sampleID] + "\n")
		out_snvSNPCountF_O.close()	


		personSNP_Lst.sort()
		SNPCiteslines = Posi_lines(cladeSamps,lieD,lsD,personSNP_Lst,effDic)
		out_SNP_F = out_P + "/" + clade + ".FiltSNPwithsnv.SNP_Posi_list.txt"
		if ( os.path.exists(out_SNP_F)):
			os.remove(out_SNP_F)
		out_SNP_F_O=open(out_SNP_F,'a')
		out_SNP_F_O.write(SNPCiteslines + "\n")
		out_SNP_F_O.close()	

		print(len(personSNP_Lst))
		print(len(personsnv_Lst))
		personsnv_Lst.sort()
		snvCiteslines = Posi_lines(cladeSamps,lieD,lsD,personsnv_Lst,effDic)
		out_snv_F = out_P + "/" + clade + ".FiltSNPwithsnv.snv_Posi_list.txt"
		if ( os.path.exists(out_snv_F)):
			os.remove(out_snv_F)
		out_snv_F_O=open(out_snv_F,'a')
		out_snv_F_O.write(snvCiteslines + "\n")
		out_snv_F_O.close()	








if __name__ == "__main__":
	usage = "usage:python  %prog [option]"
	parser = OptionParser(usage=usage)

	parser.add_option("-i","--iSNV_SNP_table",
					  dest = "iSNV_SNP_table",
					  default = "",
					  metavar = "file",
					  help = "iSNV table.  [required]")
	parser.add_option("-o","--out_P",
					  dest = "out_P",
					  default = "",
					  metavar = "path",
					  help = "output stat file Path of snv Freq distribution.  [required]")

	parser.add_option("-f","--InfoF",
					  dest = "InfoF",
					  default = "",
					  metavar = "file",
					  help = "Info File.  [required]")

	parser.add_option("-e","--snpEffPath",
					  dest = "snpEffPath",
					  default = "",
					  metavar = "path",
					  help = "snpEff file Path.  [required]")

	parser.add_option("-R","--RepeatMaskedGnm",
					  dest = "RepeatMaskedGnm",
					  default = "",
					  metavar = "file",
					  help = "RepeatMasker marked Gnm.  [required]")




	parser.add_option("-m","--snv_MinFreq",
					  dest = "snv_MinFreq",
					  default = "0.05",
					  metavar = "float",
					  help = "min Freq of snv to calculate (0-100%).  [required]")
	parser.add_option("-M","--snv_MaxFreq",
					  dest = "snv_MaxFreq",
					  default = "0.95",
					  metavar = "float",
					  help = "max Freq of snv to calculate (0-100%).  [required]")
	parser.add_option("-E","--FiltEcoliPosiF",
					  dest = "FiltEcoliPosiF",
					  default = "",
					  metavar = "file",
					  help = "Filt Ecoli homo Posi file.  [required]")


	(options,args) = parser.parse_args()
	iSNV_SNP_table   = os.path.abspath(options.iSNV_SNP_table)
	out_P		  = os.path.abspath(options.out_P)
	InfoF   = os.path.abspath(options.InfoF)
	snpEffPath = os.path.abspath(options.snpEffPath)
	RepeatMaskedGnm = os.path.abspath(options.RepeatMaskedGnm)
	snv_MinFreq		= float(options.snv_MinFreq)
	snv_MaxFreq		= float(options.snv_MaxFreq)
	FiltEcoliPosiF   = os.path.abspath(options.FiltEcoliPosiF)

	if ( not os.path.exists(out_P)):
		os.mkdir(out_P)


	main()




