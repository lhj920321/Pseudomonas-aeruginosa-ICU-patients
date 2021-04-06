
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
			
			HeaderL.append(HeaderID.split("_BDM")[0].split("BJ13-")[1])
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
			Pos = lL[0]
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
					if Freq >= snv_MinFreq:
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





def tongjiFREQ(sample,samp_FreqDic,personNAPosi,out_F,allSampsOverMinFreqDic,allPosi):
	FreqDic = samp_FreqDic[sample]
	tongji_Dic = {}
	for XX in range(10,91,step):
		tongji_Dic[XX] = 0


	for Posi in FreqDic:
		if Posi not in personNAPosi:
			Freq = FreqDic[Posi]
			if Freq >= snv_MinFreq and Freq < snv_MaxFreq:
				for XX in range(10,91,step):
					if Freq *100 >= XX and Freq *100 < XX + step:
						tongji_Dic[XX] += 1
			
			if Freq >= snv_MinFreq:
				allSampsOverMinFreqDic[sample][Posi] = Freq
			if Posi not in allPosi:
				allPosi.append(Posi)


	freq_flagLst=[]
	snvCount_list=[]
	outlst = []
	for XX in range(10,91,step):
		outlst.append(sample +"\t" + str(XX) + "\t" + str(tongji_Dic[XX]))
		freq_flagLst.append(XX)
		snvCount_list.append(tongji_Dic[XX])



	outlines = "\n".join(outlst) 
	out_F_O=open(out_F,'a')
	out_F_O.write(outlines)
	out_F_O.close()	
	
	out_plotF = out_F.split(".txt")[0] + ".pdf"
	plt.figure(dpi=300,figsize=(8,3))
	plt.bar(range(len(snvCount_list)), snvCount_list,color='red',tick_label=freq_flagLst)  


	plt.title(sample.split(".sort")[0])
	plt.xlabel('Freq', fontsize = 13)
	plt.ylabel('count', fontsize = 13)
	plt.xticks(rotation=90)
	plt.savefig(out_plotF)
	plt.savefig(out_plotF.split(".pdf")[0] + ".png")

	return allSampsOverMinFreqDic,allPosi



def main():

	lieD,lsD = iSNV_SNP_tableRead(iSNV_SNP_table)

	allSamps = []
	for Samp in lieD:
		if Samp != "#Posi" and Samp != "snv/SNP":
			allSamps.append(Samp)

	samp_FreqDic,sampNA_PosiDic = samp_Freq_Dic(allSamps,lieD,lsD)
	personNAPosi = Person_NA_Posi(allSamps,sampNA_PosiDic)	
	allSampsOverMinFreqDic = {}
	allPosi = []
	for sample in allSamps:
		out_F = out_P + "/" + sample + ".snv" + str(snv_MinFreq) + "-" + str(snv_MaxFreq) + ".txt"

		allSampsOverMinFreqDic[sample] = {}
		allSampsOverMinFreqDic,allPosi = tongjiFREQ(sample,samp_FreqDic,personNAPosi,out_F,allSampsOverMinFreqDic,allPosi)		


	samp1 = allSamps[0]
	samp2 = allSamps[1]
	samp1_FreqD = allSampsOverMinFreqDic[samp1]
	samp2_FreqD = allSampsOverMinFreqDic[samp2]

	OverMinFreqPosi_FlagDic = {}
	for OverMinFreqPosi in allPosi:
		if OverMinFreqPosi in  samp1_FreqD and OverMinFreqPosi not in  samp2_FreqD :
			Flag = "0028_single"
		elif OverMinFreqPosi in  samp1_FreqD and OverMinFreqPosi in  samp2_FreqD :
			Flag = "share"
		elif OverMinFreqPosi not in  samp1_FreqD and OverMinFreqPosi in  samp2_FreqD :
			Flag = "0134_single"
		OverMinFreqPosi_FlagDic[OverMinFreqPosi] = Flag

	SampFlags = {allSamps[0]:["share","0028_single"],allSamps[1]:["share","0134_single"]}



	for sample in allSamps:
		tongji_Dic = {}
		for XX in range(10,101,step):
			tongji_Dic[XX] = {}

		for Position in allSampsOverMinFreqDic[sample]:
			Freq = allSampsOverMinFreqDic[sample][Position]
			if Freq < snv_MaxFreq:
				Flag = OverMinFreqPosi_FlagDic[Position]

				for XX in range(10,101,step):
					if Freq *100 >= XX and Freq *100 < XX + step:
						if Flag not in tongji_Dic[XX]:
							tongji_Dic[XX][Flag] = 0 
						
						tongji_Dic[XX][Flag] += 1

		print(tongji_Dic)

		freq_flagLst = []
		for XX in range(10,101,step):
			freq_flagLst.append(XX)

		outLst = []
		FlagNumDic = {}
		for Flag in SampFlags[sample]:
			print(Flag)
			FlagNumDic[Flag] = []
			for XX in range(10,101,step):
				if Flag in tongji_Dic[XX]:
					FlagNumDic[Flag].append(tongji_Dic[XX][Flag])
				else:
					FlagNumDic[Flag].append(0)
				outline = "\t".join( [sample,Flag,str(XX),str(FlagNumDic[Flag][-1])])
				outLst.append(outline)

		out = "\n".join(outLst)
		out_putF = out_P + "/" + sample + ".0028-0134.duiji.txt"
		out_F_O=open(out_putF,'a')
		out_F_O.write(out)
		out_F_O.close()	
	
		
		out_plotF = out_P + "/" + sample + ".0028-0134.duiji.pdf"
		plt.figure(dpi=300,figsize=(8,3))

		plt.bar(freq_flagLst, FlagNumDic[SampFlags[sample][0]],color='red', alpha=0.3,tick_label=freq_flagLst)
		plt.bar(freq_flagLst, FlagNumDic[SampFlags[sample][1]],color='blue', bottom=FlagNumDic[SampFlags[sample][0]])

	
		
		#plt.bar(range(len(snvCount_list)), snvCount_list,color='red',tick_label=freq_flagLst)  


		plt.title(sample)
		plt.xlabel('Freq', fontsize = 13)
		plt.ylabel('count', fontsize = 13)
		plt.xticks(rotation=90)
		plt.savefig(out_plotF)
		plt.savefig(out_plotF.split(".pdf")[0] + ".png")




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


	parser.add_option("-m","--snv_MinFreq",
					  dest = "snv_MinFreq",
					  default = "0.1",
					  metavar = "float",
					  help = "min Freq of snv to calculate (0-100%).  [required]")
	parser.add_option("-M","--snv_MaxFreq",
					  dest = "snv_MaxFreq",
					  default = "0.9",
					  metavar = "float",
					  help = "max Freq of snv to calculate (0-100%).  [required]")

	parser.add_option("-s","--step",
					  dest = "step",
					  default = "2 ",
					  metavar = "int",
					  help = "step of Freq Stat.  [required]")

	(options,args) = parser.parse_args()
	iSNV_SNP_table   = os.path.abspath(options.iSNV_SNP_table)
	out_P		  = os.path.abspath(options.out_P)
	InfoF   = os.path.abspath(options.InfoF)
	snv_MinFreq		= float(options.snv_MinFreq)
	snv_MaxFreq		= float(options.snv_MaxFreq)
	step = int(options.step)
	if ( not os.path.exists(out_P)):
		os.mkdir(out_P)


	main()




