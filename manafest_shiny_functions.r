# Functions for analysis of Manafest data in shiny app
# 
#
geti = function(x,i){return(x[i])}

# modified version of the parse.immunoseq function from the tcR package to read new Adaptive format
parse.immunoseq4 = function (.filename, reads, aa.seq) 
{
  require(tcR)
  filename <- .filename
  nuc.seq <- "nucleotide"
  barcodes <- "vIndex"
  vgenes <- "vGeneName"
  jgenes <- "jGeneName"
  dgenes <- "dFamilyName"
  vend <- "n1Index"
  jstart <- "jIndex"
  dalignments <- c("dIndex", "n2Index")
  vd.insertions <- "n1Insertion"
  dj.insertions <- "n2Insertion"
  total.insertions <- NA
  .skip = 0
  .sep = "\t"
  df <- parse.cloneset(.filename = filename, .nuc.seq = nuc.seq, 
                       .aa.seq = aa.seq, .reads = reads, .barcodes = barcodes, 
                       .vgenes = vgenes, .jgenes = jgenes, .dgenes = dgenes, 
                       .vend = vend, .jstart = jstart, .dalignments = dalignments, 
                       .vd.insertions = vd.insertions, .dj.insertions = dj.insertions, 
                       .total.insertions = total.insertions, .skip = .skip, 
                       .sep = .sep)
  df$CDR3.nucleotide.sequence <- substr(df$CDR3.nucleotide.sequence, 
                                        df$Umi.count + 1, nchar(df$CDR3.nucleotide.sequence) - 6)
  df$CDR3.amino.acid.sequence <- bunch.translate(df$CDR3.nucleotide.sequence)
  .fix.genes <- function(.col) {
    .col <- gsub(",", ", ", .col, fixed = T, useBytes = T)
    .col <- gsub("([0])([0-9])", "\\2", .col, useBytes = T)
    .col <- gsub("TCR", "TR", .col, fixed = T, useBytes = T)
    .col
  }
  df$V.gene <- .fix.genes(df$V.gene)
  df$D.gene <- .fix.genes(df$D.gene)
  df$J.gene <- .fix.genes(df$J.gene)
  .fix.poses <- function(.col) {
    df[df[[.col]] != -1, .col] <- df[df[[.col]] != -1, .col] - 
      df$Umi.count[df[[.col]] != -1]
    df
  }
  df <- .fix.poses("V.end")
  df <- .fix.poses("D3.end")
  df <- .fix.poses("D5.end")
  df <- .fix.poses("J.start")
  df$Umi.count <- NA
  df$Umi.proportion <- NA
  return(df)
}

# reads files, removes non-productive sequencies, extracts counts, 
# creates all necessary objects, and saves them
# input: a path to folder with input files
# output:
# mergedData is AA level counts (merged by AA sequences)
# ntData is nucleotide level counts
# productiveReadCounts is the total number of reads of productive sequencies
readMergeSave = function(files, filenames = NULL)
{
	if (length(files) == 0)
	{
		print('There is no files to read')
			return(NULL);
	}
		require(tools)
		require(immunarch)
		
  # output objects
		mergedData = ntData = list()
		
		# read all files with immunarch functionality
		repertoire = repLoad(.path = files)
		
		# the names of files that were actually read in
		readFiles = c()
		for(i in names(repertoire$data))
		{
			print(i)
		  dat = repertoire$data[[i]]
				#count reads of productive sequences only
				mergedData[[i]] = tapply(dat[,'Clones'], dat[,'CDR3.aa'], sum, na.rm = T)
				# nucleotide level data
				ntData[[i]] = tapply(dat[,'Clones'], dat[,'CDR3.nt'], sum, na.rm = T)
				readFiles = c(readFiles, i)
		}
		if (length(mergedData) == 0)
		{
			print(paste('There are no data to read'))
			return(NULL)
		}
		# if file names are not supplied, use internal shiny server file names (0,1,2,..)
		if (is.null(filenames)) 
		{
		  filenames = sapply(unlist(lapply(readFiles,basename)),file_path_sans_ext)
		} else {
		  # it files names are supplied (actual file names that were loaded)
		  filenames = sapply(unlist(filenames),file_path_sans_ext)
		}	
		# assign file names as names to objects
		names(mergedData) = names(ntData) = filenames
		return(list(mergedData = mergedData,ntData = ntData))
}

# fisher's test for a pair of samples
# NB: reference sample is the second sample in pair if dir = T (suggesting expansion compare to baseline and corresponding OR > 1)
# Use filter for the number of reads (nReadFilter) if a set of clones to analyze is not provided
runFisher = function(pair, mergedData,totalReadCounts, nReadFilter = c(10,0), dir = T, clones = NULL)
{
	sampID1 = pair[1]
	sampID2 = pair[2]

	if (dir) {clonesToTest = setdiff(names(mergedData[[sampID1]])[which(mergedData[[sampID1]] >= nReadFilter[1])],"")
	}else{
		clonesToTest = union(names(mergedData[[sampID1]])[which(mergedData[[sampID1]] >= nReadFilter[1])],
			names(mergedData[[sampID2]])[which(mergedData[[sampID2]] >= nReadFilter[2])])
		clonesToTest = setdiff(clonesToTest,"")
	}
	
#print(length(clonesToTest))
	
# if set of clones is specified take intersection
	if(!is.null(clones)) clonesToTest = intersect(clonesToTest, clones)
	
	if (length(clonesToTest) == 0)
	{
		print('There are no clones to test')
		return(NULL)
	}
	
	dat = matrix(0,ncol = 2, nrow = length(clonesToTest))
	colnames(dat) = c(sampID1,sampID2)
	rownames(dat) = clonesToTest
	s = intersect(clonesToTest,names(mergedData[[sampID1]]))
	dat[s,sampID1] = mergedData[[sampID1]][s]
	s = intersect(clonesToTest,names(mergedData[[sampID2]]))
	dat[s,sampID2] = mergedData[[sampID2]][s]
#	print(head(dat))
#	print(dim(dat))
	
	res = matrix(nrow = length(clonesToTest), ncol = 4)
	rownames(res) = clonesToTest
#print(sampID1)
	for (clone in clonesToTest)
	{
#print(clone)
		counts1 = dat[clone,sampID1]
		counts2 = dat[clone,sampID2]
		tab = rbind(c(counts1,counts2),c(totalReadCounts[sampID1]-counts1, totalReadCounts[sampID2]-counts2))
		fres = fisher.test(tab)
#print(fres)
		res[clone,] = c(fres$p.value,fres$estimate,fres$conf.int)	
	}
#print(dim(res))
	if(length(clonesToTest)>1) 
	{
		res = cbind(res[,1],p.adjust(res[,1], method = 'BH'),res[,c(2:4)])
	}else{
		res = data.frame(matrix(c(res[,1],p.adjust(res[,1], method = 'BH'),res[,c(2:4)]), nrow = 1))
		rownames(res) = clonesToTest
	}
#print(dim(res))
	colnames(res) = c('p.val','FDR','odds.ratio','CI_low','CI_up')
	res = res[order(res[,'FDR']),]
	return(res)
}

# write output
# p-value, OR from fisher test + abundance + frequency
createResTable = function(res,mergedData,totalReadCountPerSample, orThr = 1, FDR_threshold = 0.05, 
	saveCI = T, significanceTable = F)
{
	clones = c()
	notNullComp = names(res)[!sapply(res,is.null)]
	# get all clones for output
	for (i in notNullComp){
#print(i)
		clones = union(clones,rownames(res[[i]])[which(as.numeric(res[[i]][,'odds.ratio']) >= orThr & 
			as.numeric(res[[i]][,'CI_low']) > 1 & as.numeric(res[[i]][,'FDR']) < FDR_threshold)])
	}
	if(length(clones)==0){print('There is no significant clones'); return(NULL)}
	# create output matrices
	output_counts = output_freq = matrix(0,nrow = length(clones), ncol = length(mergedData)) 
	output_fdr = output_OR = output_CI = matrix(nrow = length(clones), ncol = length(res))
	rownames(output_counts) = rownames(output_freq) = rownames(output_fdr) = rownames(output_OR) = rownames(output_CI) = clones
	colnames(output_fdr) = colnames(output_OR) = colnames(output_CI) = names(res)
	colnames(output_counts) = colnames(output_freq) = names(mergedData)
	for (i in notNullComp)
	{
		rows = intersect(rownames(res[[i]]),rownames(output_fdr))
		output_fdr[rows,i] = res[[i]][rows,'FDR']
		output_OR[rows,i] = round(res[[i]][rows,'odds.ratio'],3)
		if (saveCI) {output_CI[rows,i] = paste(round(res[[i]][rows,'CI_low'],3),'-',round(res[[i]][rows,'CI_up'],3), sep = '')}
	}
	for (i in names(mergedData))
	{
		rows = intersect(names(mergedData[[i]]),rownames(output_counts))
		output_counts[rows,i] = mergedData[[i]][rows]
		output_freq[rows,i] = round((mergedData[[i]][rows]/totalReadCountPerSample[i])*100,3)
	}	
	colnames(output_fdr) = paste('FDR:',names(res))
	colnames(output_OR) = paste('OR:',names(res))
	colnames(output_CI) = paste('CI95%:',names(res))
	colnames(output_counts) = paste(names(mergedData),'abundance', sep = '_')
	colnames(output_freq) = paste(names(mergedData),'percent', sep = '_')
	if(saveCI) {
		output_OR_CI = c()
		for(i in 1:ncol(output_OR))
		{
		  output_OR_CI = cbind(output_OR_CI, output_OR[,i],output_CI[,i])
		}
		colnames(output_OR_CI) = paste(rep(c('OR:','CI95%:'),length(res)),rep(names(res),each = 2))
	}else {output_OR_CI = output_OR}
	tab = data.frame(sequence = clones,output_fdr,significant_comparisons = apply((as.numeric(output_fdr) < FDR_threshold & output_OR > orThr),1,sum,na.rm = T),
		output_OR_CI,output_counts,output_freq, check.names = F)
	tab = tab[which(tab[,'significant_comparisons'] > 0),]
	# if significanceTable return a binary table clones vs conditions specifying which clone is significant in what condition
	if (significanceTable)
	{
		signTab = matrix(as.numeric(output_fdr[rownames(tab),]) < FDR_threshold & output_OR[rownames(tab),] > orThr, 
			nrow = nrow(tab), ncol = length(res), dimnames = list(rownames(tab),names(res)))
		return(signTab)
	} else return(tab)
}
# write results
# if there are only comparisons to the reference sample, saves the results to one text file
# if there are all pairwise comparisons, create a text file per peptide and create an archive with them
writeResToFile = function(res,fileName = 'res.txt')
{ 
		if (is.null(dim(res)))
		{
			write.table(t(res),file =fileName, sep = '\t',row.names = F,quote = F, eol = '\r')
		}else{
			write.table(res,file = fileName, sep = '\t',row.names = F,quote = F, eol = '\r')
		}
}


# returns read counts and frequencies for selected clones and samples
getCountsFreq = function(clones, mergedData,totalReadCountPerSample, samp = names(mergedData))
{
	output_counts = output_freq = matrix(0,nrow = length(clones), ncol = length(samp)) 
	rownames(output_counts) = rownames(output_freq) = clones
	colnames(output_counts) = colnames(output_freq) = samp
	for (i in samp)
	{
		rows = intersect(names(mergedData[[i]]),clones)
		output_counts[rows,i] = mergedData[[i]][rows]
		output_freq[rows,i] = (mergedData[[i]][rows]/totalReadCountPerSample[i])*100
	}	
	colnames(output_counts) = paste(samp,'abundance', sep = '_')
	colnames(output_freq) = paste(samp,'percent', sep = '_')
	tab = cbind(output_counts,output_freq)
	return(tab)
}

# returns difference in percent for selected clones and samples
getDiff = function(clones, mergedData, totalReadCountPerSample, samp = names(mergedData), refSamp)
{
	output_freq = matrix(0,nrow = length(clones), ncol = length(samp)+1, dimnames = list(clones, c(samp,refSamp)))
	# get frequincies for all samples
	for (i in c(samp,refSamp))
	{
		rows = intersect(names(mergedData[[i]]),clones)
		output_freq[rows,i] = (mergedData[[i]][rows]/totalReadCountPerSample[i])*100
	}
	if (ncol(output_freq) == 2)
	{
		output_diff = data.frame(output_freq[,samp] - output_freq[,refSamp])
		colnames(output_diff) = paste0('DIFF:',samp)
	}else{
		output_diff = apply(output_freq[,samp],2,function(x) x-output_freq[,refSamp])
		colnames(output_diff) = paste0('DIFF:',colnames(output_diff))	
	}
	return(output_diff)
}

# return frequencies in percent for selected clones and samples
getFreq = function(clones, mergedData,totalReadCountPerSample, samp = names(mergedData), colSuf = 'percent', minRead = 0)
{
	output_freq = matrix(minRead,nrow = length(clones), ncol = length(samp)) 
	rownames(output_freq) = clones
	colnames(output_freq) = samp
	for (i in samp)
	{
		rows = intersect(names(mergedData[[i]]),clones)
		output_freq[rows,i] = mergedData[[i]][rows]
		output_freq[,i] = (output_freq[,i]/totalReadCountPerSample[i])*100		
	}	
	if (colSuf != '') colnames(output_freq) = paste(samp,colSuf, sep = '_')
	else colnames(output_freq) = samp
	return(data.frame(output_freq,check.names = F))
}
# return frequencies or abundancies in percent for selected clones and samples
getFreqOrCount = function(clones, mergedData,totalReadCountPerSample, samp = names(mergedData), colSuf = 'percent', minRead = 0, returnFreq = T)
{
	output_freq = matrix(minRead,nrow = length(clones), ncol = length(samp)) 
	rownames(output_freq) = clones
	colnames(output_freq) = samp
	for (i in samp)
	{
		rows = intersect(names(mergedData[[i]]),clones)
		output_freq[rows,i] = mergedData[[i]][rows]
		if (returnFreq) output_freq[,i] = (output_freq[,i]/totalReadCountPerSample[i])*100		
	}	
	if (colSuf != '') colnames(output_freq) = paste(samp,colSuf, sep = '_')
	else colnames(output_freq) = samp
	return(data.frame(output_freq,check.names = F))
}

# return fold change of frequencies for selected clones and samples
getFC = function(clones, mergedData, totalReadCountPerSample, refSamp, samp = names(mergedData), colSuf = 'FC:')
{
	samp = setdiff(samp,refSamp)
	freq = getFreq(clones, mergedData, totalReadCountPerSample, union(refSamp,samp), colSuf = '', minRead = 1)
	# reference samples frequencies
	ref = freq[,refSamp]
	# divide to reference samples frequencies
	fc = round(sweep(data.frame(freq[,samp]),1, ref, "/"))
	colnames(fc) = paste0(colSuf,samp)
	rownames(fc) = clones
	return(fc)
}

# returns clones that unique for a selected condition
getUniqueClones = function(samp, mergedData, readCountThr = 0)
{
	res = names(mergedData[[samp]])[which(mergedData[[samp]]>= readCountThr)]
	for(i in setdiff(names(mergedData),samp)) { res = setdiff(res,names(mergedData[[i]])[which(mergedData[[i]]> readCountThr)])}
	return(res)
}

# make heatmaps of frequencies (if refSamp is not specified) of each element of input list of clones and samples and plot FC if refSamp is specified
makeHeatmaps = function(listOfclones, mergedData, totalReadCountPerSample, samp = names(mergedData), refSamp = NULL, 
	fileName = 'heatmap.pdf',size = 7)
{
	if (length(listOfclones)==0)
	{
		print('There are no clones to plots')
		return;
	}
	require(gplots)
	pdf(fileName, width = size, height = size)
	samp = setdiff(samp, refSamp)
	for(i in 1:length(listOfclones))
	{
		clones = listOfclones[[i]]
		if(length(clones)< 2) next;
		plotData = getFreq(clones, mergedData,totalReadCountPerSample, samp)
		if (!is.null(refSamp))
		{
			freqRef = getFreq(clones, mergedData,totalReadCountPerSample, refSamp)
			freqRef[freqRef == 0] = min(plotData[plotData > 0]) - 1e-10
			plotData = sweep(plotData,1,unlist(freqRef),'/')
		}
#		par(oma = c(.1,.1,.1,3))
		heatmap(as.matrix(plotData), col = bluered(100),labCol = lapply(strsplit(colnames(plotData),'_percent'),geti,1),
			scale = 'row', cexCol = .5, cexRow = .3)
		title(main = names(listOfclones)[[i]], cex=0.2)
	}
	dev.off()
}

runSingleFisher = function(clone, pair, mergedData,totalReadCounts)
{
	samp1 = pair[1]
	samp2 = pair[2]
	counts1 = mergedData[[samp1]][clone]
	counts2 = mergedData[[samp2]][clone]
	if(is.na(counts2)) counts2 = 0
	tab = rbind(c(counts1,counts2),c(totalReadCounts[samp1]-counts1, totalReadCounts[samp2]-counts2))
	fres = fisher.test(tab)
	return(c(fres$p.value,fres$estimate,fres$conf.int))
}

getUniqSignClones = function(tab, FDR_threshold = 0.05)
{
	clones = rownames(tab)[which(tab[,'significant_comparisons'] == 1)]
	n = grep('FDR:',colnames(tab), fixed = T, value= T)
	if(length(n) == 1)
	{
		return(rownames(tab)[which(as.numeric(tab[,n]) < FDR_threshold)])
	}else{
		fdrs = t(apply(tab[clones,n],1,function(x){as.numeric(x)< FDR_threshold}))
		colnames(fdrs) = n
		return(apply(fdrs,2,function(x){res = names(x)[x];return(res[!is.na(res)])}))
	}
}

# returns a vector with positive clones as names and conditions, in which a clone is significant, as values
getPositiveClones = function(analysisRes, mergedData, totalReadCounts, samp = names(mergedData), orThr = 1, fdrThr = 0.05, nReads = 10)
{
	resTable = createResTable(analysisRes, mergedData,totalReadCounts, orThr = orThr,
				FDR_threshold=fdrThr, saveCI =F) 
	if(is.null(resTable)) return(NULL)
	
	# find significant expansions
	clones = rownames(resTable)[which(resTable[,'significant_comparisons'] == 1)]
	#================
	# find condition in which it's expanded
	signTable = createResTable(analysisRes, mergedData,totalReadCounts, orThr = orThr,
				FDR_threshold=fdrThr, saveCI =F, significanceTable = T)
#	signTable = data.frame(signTable[clones,])
	signMatrix = matrix(signTable[clones,], 
			nrow = length(clones), ncol = ncol(signTable), dimnames = list(clones,sapply(strsplit(colnames(signTable),'_vs_'), function(x)x[1])))

	# fix column names
#	colnames(signTable) = sapply(strsplit(colnames(signTable),'_vs_'), function(x)x[1])
	# returns condition in which a clone is significant

	signCond = apply(signMatrix,1,function(x) colnames(signMatrix)[which(x)])
	#=====================
	# if we have only one comparison
	if (length(analysisRes) == 1) 
	{	
		# get name of the condition
		fdrCol = grep('FDR: ', colnames(resTable), fixed = T, value = T)
		n = unlist(strsplit(fdrCol,'FDR: '))[2]
		n = unlist(strsplit(n,'_vs_'))[1]
		res = rep(n, length(clones))
		names(res) = clones
		return(res)
	}
	#=====================
	# check if they are unique by checking top two conditions with the highest number of reads
#	countMatrix = getFreqOrCount(clones,mergedData,totalReadCounts,samp, colSuf = '', returnFreq = F)
	freqMatrix = getFreqOrCount(clones,mergedData,totalReadCounts,samp, colSuf = '', returnFreq = T)
#print(countMatrix)
	# remove conditions with less than nReads reads from analysis 
#	freqMatrix[countMatrix<nReads] = 0 # removed this filter in v12
	# compare with the second highest
	fishRes1 = getFisherForNclone(freqMatrix, rownames(freqMatrix),2,mergedData,totalReadCounts)
	# compare with the third highest
	fishRes2 = getFisherForNclone(freqMatrix, rownames(freqMatrix),3,mergedData,totalReadCounts)
	# combine results
	fishResComb = cbind(fishRes1[,'FDR'],fishRes2[,'FDR'], fishRes1[,'odds.ratio'],fishRes2[,'odds.ratio'])
#print(fishResComb)
#print(dim(fishResComb))
#	print(fishRes1)
#	print(fishRes2)
	# select clones that have significant FDRs and OR higher than threshold meaning that a clone is significant and unique expansion
	# also select clones that have NAs in FDR and OR, which means that this clone appears in only one condition and there is nothing to compare
	# check if there is a condition that is also significantly expanded
	fdrClones2 = apply(fishResComb,1,function(x) any(as.numeric(x[1:2])>fdrThr|as.numeric(x[3:4])<orThr))
#print(cbind(fishResComb,fdrClones2))
	posClones = names(fdrClones2)[which(!fdrClones2|is.na(fdrClones2))]
	if(length(posClones)>0)
	{
	# get corresponding condition
		return(signCond[posClones])
	}else{return(NULL)}
}
# runs Fisher's test for the nth clone with the highest frequency/the number of reads
getFisherForNclone = function(freq, clones, n = 2,mergedData,productiveReadCounts)
{
	fishRes = c()
	freq[freq==0] = NA
	for( i in clones)
	{
		r = freq[i,]
		r = r[order(r,decreasing = T)]
		if (!is.na(r[n])) res = runSingleFisher(i,names(r)[c(1,n)],mergedData,productiveReadCounts) else res = rep(NA,4)
		fishRes = rbind(fishRes,c(names(r)[1],res, n, names(r)[n]))
	}
	rownames(fishRes) = clones
	fishResMatrix = matrix(ncol = ncol(fishRes)+1, nrow = nrow(fishRes), 
		dimnames = list(rownames(fishRes),c('condition','p.val','FDR','odds.ratio','CI_low','CI_up','nComp','condToComp')))
	fishResMatrix[,'FDR'] = p.adjust(fishRes[,2], method = 'BH')
	fishResMatrix[,setdiff(colnames(fishResMatrix),'FDR')] = fishRes
	return(fishResMatrix)
}


# create output tables to be saved in Excel
createPosClonesOutput = function(posClones,  mergedData,totalReadCounts, refSamp = NULL, baselineSamp = NULL, addDiff = T)
{
	output = vector(mode = 'list')
	clones = names(posClones)
	# write peptide summary of positive clones	
	freqMatrix = getFreq(clones,mergedData,totalReadCounts,names(mergedData), colSuf = '')
	peptLevelList = tapply(clones,posClones, FUN = function(x)return(x))
	peptideTab = matrix(nrow = length(peptLevelList), ncol = 2, 
		dimnames = list(names(peptLevelList),c('positive_clones','sum_freq')))
	peptideTab[,'positive_clones'] = sapply(peptLevelList,length)
	peptideTab[,'sum_freq'] = sapply(names(peptLevelList),function(x) sum(freqMatrix[peptLevelList[[x]],x]))
	output$condition_summary = data.frame(peptideTab)
		
	# write clone-level summary 
	if(!is.null(baselineSamp) && !is.null(refSamp))
	{
		fc_ref = getFC(clones,mergedData,totalReadCounts,refSamp, unique(posClones), "")
		fc_bl = getFC(clones,mergedData,totalReadCounts,baselineSamp, unique(posClones), "")
		tab = data.frame(condition = posClones, 
			getFreq(clones,mergedData,totalReadCounts,baselineSamp),
			sapply(clones,function(x) fc_bl[x,posClones[x]]),
			getFreq(clones,mergedData,totalReadCounts,refSamp),
			sapply(clones,function(x) fc_ref[x,posClones[x]]),check.names = F)
		colnames(tab)[c(3,5)] = paste0('FC:', c(baselineSamp,refSamp))
	} else {
		if(!is.null(refSamp))
		{
			fc_ref = getFC(clones,mergedData,totalReadCounts,refSamp, unique(posClones), "")
			tab = data.frame(condition = posClones, 
				getFreq(clones,mergedData,totalReadCounts,refSamp),
				sapply(clones,function(x) fc_ref[x,posClones[x]]),check.names = F)
			colnames(tab)[3] = paste0('FC:', refSamp)
		}else{
			tab = data.frame(condition = posClones)
		}
	}
	output$positive_clones_summary = tab
	
	# extended information for all samples
	tab = getCountsFreq(clones, mergedData,totalReadCounts, samp = names(mergedData))
	if (addDiff) 
	{
		tab = cbind(tab, getDiff(clones, mergedData,totalReadCounts, samp = setdiff(names(mergedData),c(refSamp, baselineSamp)), refSamp))		
	}
#browser()
	output$positive_clones_all_data = data.frame(condition = posClones,tab,check.names = F)
	return(output)
}

# calculate frequency threshold using the number of cells and probability
# (1-((1-p)^(1/n)))
# where n=number of cells per well and p=selected probability
getFreqThreshold = function(n, p)
{
	return (1-((1-p)^(1/n)))
}


# selects clones that pass the nReads threshold and compare top conditions with the second and third top conditions
# and returns a table with FDRs and ORs for those comparisons that will be used to select positive clones using FDR and OR threshold later
compareWithOtherTopConditions = function(mergedData, productiveReadCounts, sampForAnalysis, nReads = 10, clones = NULL)
{
	# combine clones from all conditions that have the number of clones more than nReads
	allClones = sapply(mergedData[sampForAnalysis], function(x) setdiff(names(x)[which(x >= nReads)],""))
	clonesToTest = c()
	for(i in names(allClones))
	{
		clonesToTest = union(clonesToTest, allClones[[i]])
	}

# if set of clones is specified take intersection
	if(!is.null(clones)) clonesToTest = intersect(clonesToTest, clones)

	if (length(clonesToTest) == 0)
	{
		print('There are no clones to test')
		return(NULL)
	}

	# get frequency matrix
	freqMatrix = getFreqOrCount(clonesToTest,mergedData,productiveReadCounts,sampForAnalysis, colSuf = '', returnFreq = T)
	#print(countMatrix)
	# compare with the second highest
	fishRes1 = getFisherForNclone(freqMatrix, rownames(freqMatrix),2,mergedData,productiveReadCounts)
	# compare with the third highest
	fishRes2 = getFisherForNclone(freqMatrix, rownames(freqMatrix),3,mergedData,productiveReadCounts)
	# combine results
	fishResComb = cbind(FDR1 = fishRes1[,'FDR'], FDR2 = fishRes2[,'FDR'], OR1 = fishRes1[,'odds.ratio'], OR2 = fishRes2[,'odds.ratio'], condition = fishRes1[,'condition'])
	
	return(fishResComb)
}

# select clones that have significant FDRs and OR higher than threshold meaning that a clone is significant and unique expansion
# also select clones that have NAs in FDR and OR, which means that this clone appears in only one condition and there is nothing to compare
# check if there is a condition that is also significantly expanded
# This function takes a table with FDRs, ORs, and condition in which the clone is the most abundant
# Column 1 and 2 are FDRs, 3 and 4 - ORs, and 5 is condition
getPositiveClonesFromTopConditions = function(fisherResTable, orThr = 1, fdrThr = 0.05)
{
	# apply FDR and OR thresholds
	fdrClones2 = apply(fisherResTable,1,function(x) any(as.numeric(x[1:2])>fdrThr|as.numeric(x[3:4])<orThr))
	# find positive clones
	posClones = names(fdrClones2)[which(!fdrClones2|is.na(fdrClones2))]
	if(length(posClones)>0)
	{
	# return conditions of positive clones
		return(fisherResTable[posClones,'condition'])
	}else{return(NULL)}
}
