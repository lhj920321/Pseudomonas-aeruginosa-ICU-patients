
#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys,os
from optparse import OptionParser
import matplotlib.pyplot as plt
import numpy as np 



def readInfo(InfoF):
	#cladeSampDic = {}
	allcladeSampDic = {}
	for InfoFl in open(InfoF).readlines():
		if "#" not in InfoFl:
			InfoFl_Tags = InfoFl.split("\n")[0].split("\t")
			seqID = InfoFl_Tags[0]
			clade = InfoFl_Tags[20]
			if clade  not in allcladeSampDic:
				allcladeSampDic[clade] = []
			allcladeSampDic[clade].append(seqID)

	return allcladeSampDic




def readPindelSV(pindelSVFs):
	allSamp_pindelDic = {}	
	for samppindelSVF in pindelSVFs:
		ID = samppindelSVF.split("/")[-1].split("_BDM")[0].split("BJ13-")[1]
		allSamp_pindelDic[ID] = {}

		for samppindelSVFl in open(samppindelSVF).readlines():
			if "#" not in samppindelSVFl and samppindelSVFl != "\n":
				Tags = samppindelSVFl.split("\t")
				infoTags = Tags[7].split(";")
				Posi = int(Tags[1])
				SVType = infoTags[3]

				for infTag in infoTags:
					if "SVTYPE" in infTag:
						SVType = infTag.split("=")[-1]
					if "SVLEN" in infTag:
						SVLen = abs(int(infTag.split("=")[-1]))

				if SVType == "DEL" and SVLen <= maxIndelLen:
					allSamp_pindelDic[ID][Posi] = SVLen
					
	return allSamp_pindelDic



def readVarscanSV(varscanFs):
	allSamp_dellyDic = {}
	for varscanF in varscanFs:
		ID = varscanF.split("/")[-1].split("_BDM")[0].split("BJ13-")[1]
		varscanF = varscanP + "/"+ varscanF	
		allSamp_dellyDic[ID] = {}
		for varscanFl in open(varscanF).readlines():
			if "#" not in varscanFl and varscanFl != "\n":
				Tags = varscanFl.split("\t")
				infoTags = Tags[7].split(";")
				PASS = Tags[6]
				Posi = int(Tags[1])
				for infTag in infoTags:
					if "SVTYPE" in infTag:
						SVType = infTag.split("=")[-1]
					if infTag.split("=")[0] ==  "END":
						SV_end =  int(infTag.split("=")[-1])
						SVLen = SV_end - Posi + 1
				if SVType == "DEL" and SVLen >= maxIndelLen:		
					allSamp_dellyDic[ID][Posi] = SVLen

	return allSamp_dellyDic




def read_file_name(file_dir):
	Files = []
	for root,dirs,files in os.walk(file_dir):
		for file in files:
			if "_D.vcf" in file:
				Files.append(os.path.join(root,file))
	return Files






def pindel_dellyCommonSV(CladeSamps,allSamp_pindelDic,allSamp_dellyDic,sampCommonSVdic,allSV_over50notCommom):
	#sampCommonSVdic = {}
	AllcommSVPosiLst = []
	print(CladeSamps)
	print(allSamp_dellyDic["0016"][3841810])
	print(allSamp_pindelDic["0016"][3841810])
	for sampleID in CladeSamps:
		sampCommonSVdic[sampleID] = {}
		
		for SV_pindel_Posi in allSamp_pindelDic[sampleID]:			
			if SV_pindel_Posi in allSamp_dellyDic[sampleID] and  abs(allSamp_pindelDic[sampleID][SV_pindel_Posi] -  allSamp_dellyDic[sampleID][SV_pindel_Posi] ) <= 50: 
				Flag = "Same"
				sampCommonSVdic[sampleID][SV_pindel_Posi] = allSamp_pindelDic[sampleID][SV_pindel_Posi]
				if SV_pindel_Posi not in AllcommSVPosiLst:
					AllcommSVPosiLst.append(SV_pindel_Posi)
					#print(sampleID)
				if "0016" in sampleID and SV_pindel_Posi == 3841810:
					print(SV_pindel_Posi)

	
	for commSVPosi in AllcommSVPosiLst:
		shareCount = 0	 
		for sampID in CladeSamps:
			Len = ""
			if commSVPosi in sampCommonSVdic[sampID]:
				Len = str(sampCommonSVdic[sampID][commSVPosi])
				shareCount  += 1
		#print(commSVPosi)
		#print(shareCount)
		#print(CladeSamps)
		#print()
		if shareCount != len(CladeSamps):
			if commSVPosi not in allSV_over50notCommom:
				allSV_over50notCommom.append(commSVPosi)


	return allSV_over50notCommom,sampCommonSVdic

			
			




def main():

	pindelSVFs = read_file_name(pindelP)
	varscanFs = os.listdir(varscanP)

	allSamp_pindelDic = readPindelSV(pindelSVFs)
	allSamp_dellyDic = readDellySV(varscanFs)


	allcladeSampDic = readInfo(InfoF)
	print(allcladeSampDic)

	allSamps = []	
	for Samp_delly in pindelSVFs:
		sampID = Samp_delly.split("/")[-1].split("_BDM")[0].split("BJ13-")[1]
		allSamps.append(sampID)
	print(allSamps)


	sampCommonSVdic = {}
	allSV_over50notCommom = []

	for clade in allcladeSampDic:
		cladeSamps = allcladeSampDic[clade]	
		allSV_over50notCommom,sampCommonSVdic = pindel_dellyCommonSV(cladeSamps,allSamp_pindelDic,allSamp_dellyDic,sampCommonSVdic,allSV_over50notCommom)


	print(sampCommonSVdic)
	header = ["Posi"]
	allSamps.sort()
	for samp in allSamps:
		header.append(samp)
	Head = "\t".join(header)	

	outL = []
	#allSV_over50notCommom.sort()
	print(allSV_over50notCommom)
	for SV_over50notCommom in allSV_over50notCommom:
		outl = [str(SV_over50notCommom)]
		for samp in allSamps:
			Len = ''
			if SV_over50notCommom in sampCommonSVdic[samp]:
				Len = str(sampCommonSVdic[samp][SV_over50notCommom])
			outl.append(Len)
		print(outl)
		outL.append("\t".join(outl))

	out_F = out_P + "/allPerson.Deletion.pindel_delly_common.Stat.txt"
	if ( os.path.exists(out_F)):
		os.remove(out_F)
	out_SNP_F_O=open(out_F,'a')
	out_SNP_F_O.write(Head + "\n")
	out_SNP_F_O.write("\n".join(outL) + "\n")
	out_SNP_F_O.close()	





if __name__ == "__main__":
	usage = "usage:python  %prog [option]"
	parser = OptionParser(usage=usage)

	parser.add_option("-f","--InfoF",
					  dest = "InfoF",
					  default = "",
					  metavar = "file",
					  help = "Info File.  [required]")

	parser.add_option("-p","--pindelP",
					  dest = "pindelP",
					  default = "",
					  metavar = "path",
					  help = "pindel Path.  [required]")

	parser.add_option("-d","--varscanP",
					  dest = "varscanP",
					  default = "",
					  metavar = "path",
					  help = "varscanP.  [required]")

	parser.add_option("-o","--out_P",
					  dest = "out_P",
					  default = "",
					  metavar = "path",
					  help = "out_P. (.common-) [required]")

	parser.add_option("-m","--maxIndelLen",
					  dest = "maxIndelLen",
					  default = "50",
					  metavar = "int",
					  help = "max indels length.  [optional]")


	(options,args) = parser.parse_args()
	InfoF   = os.path.abspath(options.InfoF)
	pindelP = os.path.abspath(options.pindelP)
	varscanP  = os.path.abspath(options.varscanP)
	out_P  = os.path.abspath(options.out_P)
	maxIndelLen =int(options.maxIndelLen)



	if ( not os.path.exists(out_P)):
		os.mkdir(out_P)


	main()





