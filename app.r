# 2017-11-14 (v6)
# add fields for the number of cells and probability to calculate baseline frequency 
# to filter out clones with low frequency
# this version uses only clones that appear in baseline above calculated frequency
# previously, we used 10 reads as a threshold and do not restrict to baseline clones

# 2018-01-23 (v7)
# Added a threshold for the number of reads into interface
# Allowed user to ignore baseline threshold
# Make only one heatmap with all positive clones
# Output all significant clones relative to the reference even there are no positive clones

# 2018-03-07 (v7.1)
# Trying to solve out of memory issue when saving large table to Excel
# Solved by using WriteXLS and saving all results at once
# fixed problem with using AA level data when selecting clones to analyze with Use nucleotide level option selected

# 2018-04-18 (v8)
# Added a tab with user's manual to UI 
# Added baseline percent and the corresponding number of templates in baseline to output parameter sheet

# 2018-06-07 (v9)
# switch to the parser from the tcR package to eliminate sequensing error on 5' end of the nt sequence

# 2018-07-20 (v10)
# fixed reading function. It doesn't save mergedData and ntData to the disc; it returns a list with two element - mergedData and ntData
# fixed creating Excel file with no significant clones

# 2019-06-04
# fixed output parameters issue
# add the OR threshold when compare to all other conditions (previously, it was only the FDR threshold)

# 2019-06-14 (v12)
# improved an input process. Now a user can select what file format to upload: Adaptive, VDJtools, or the previously saved R object
# add excluded samples to output 
# remove nRead requirement when compare with other wells

# 2020-05-03 (v13)
# switch to immunarch: can't install on shiny server
# make comparison with reference as an option (checkbox)
# if off, then each condition will be compared with top two other conditions


server <- function(input, output,session) {
	
	# increase file upload limit to 100M
	options(shiny.maxRequestSize=100*1024^2, java.parameters = "-Xmx8000m")
	require(shiny)
	require(tools)
#	require(xlsx)
	if (!require(WriteXLS)) install.packages("WriteXLS")
	require("Matrix")
	if(!require(tcR)) install.packages("tcR")
	# if (!require(foreign)) install.package('foreign')
	# if(!require(immunarch)) {
		# install.packages("devtools") # skip this if you already installed devtools
		# devtools::install_url("https://github.com/immunomind/immunarch/raw/master/immunarch.tar.gz")
	# }

	source('manafest_shiny_functions.r')
	
# read source files
observeEvent(input$sourceFiles,{  output$message = renderUI({
   # check if there is file to analyze
	if (length(input$sourceFiles) == 0)
	{
		HTML('There are no files to analyze. Please upload files')
	}else if(length(input$sourceFiles$name) == 1){
		HTML('There should be at least two files to analyze')
	}else if (length(input$sourceFiles$name) > 1){ 
		# read loaded files in and save into inputData object
		
		# if there are previously loaded data, remove it
		if (exists('mergedData', envir = .GlobalEnv)) 
		{
			rm(list = c('mergedData','productiveReadCounts','ntData'), envir = .GlobalEnv)
			updateSelectInput(session, "refSamp", choices=list('None'))
			updateSelectInput(session, "baselineSamp", choices=list('None'))
			updateSelectInput(session, "excludeSamp", choices=list('None'))
		}
		# specify input file format
		if (input$inputFiles == 'VDJtools files') fileFormat = 'vdjtools'
		if (input$inputFiles == 'Adaptive files') fileFormat = 'adaptive'
		# read in input files
		res = readMergeSave(files = unlist(input$sourceFiles$datapath), filenames = unlist(input$sourceFiles$name), inputFileFormat = fileFormat)
		
		if(!is.null(res))
		{
			assign('mergedData', res$mergedData, envir = .GlobalEnv)
			assign('ntData', res$ntData, envir = .GlobalEnv)
			assign('productiveReadCounts', sapply(res$mergedData,sum), envir = .GlobalEnv)
			
			if(file.exists('inputData.rda')) file.remove('inputData.rda')
			updateSelectInput(session, "refSamp", choices=c('None',names(mergedData)))
			updateSelectInput(session, "baselineSamp", choices=c('None',names(mergedData)))
			updateSelectInput(session, "excludeSamp", choices=names(mergedData))
	#		HTML(c('Uploaded files:<br/>',paste(unlist(input$sourceFiles$name), collapse = '<br/>')))
			HTML(c('Uploaded files:<br/>',paste(names(mergedData), collapse = '<br/>')))
		}else{
			HTML('Error in reading files')
		}
	}
			
  })
}) 

# load RDA file with the input data
observeEvent(input$inputObj,{  output$message = renderUI({
   # check if there is a file to analyze
	if (length(input$inputObj) == 0)
	{
		HTML('There are no files to analyze. Please upload files')
	}else if(length(input$inputObj$name) > 0 )
	{
		# if there are previously loaded data, remove it
		if (exists('mergedData', envir = .GlobalEnv)) 
		{
			rm(list = c('mergedData','productiveReadCounts','ntData'), envir = .GlobalEnv)
			updateSelectInput(session, "refSamp", choices=list('None'))
			updateSelectInput(session, "baselineSamp", choices=list('None'))
			updateSelectInput(session, "excludeSamp", choices=list('None'))
		}
		if (tolower(file_ext(input$inputObj$name)) == 'rda'){ 
			# upload an existing input object
			load(input$inputObj$datapath, envir = .GlobalEnv)
			if (exists('mergedData', envir = .GlobalEnv))
			{
				updateSelectInput(session, "refSamp", choices=c('None',names(mergedData)))
				updateSelectInput(session, "baselineSamp", choices=c('None',names(mergedData)))
				updateSelectInput(session, "excludeSamp", choices=names(mergedData))
				HTML(c('Uploaded samples:<br/>',paste(names(mergedData), collapse = '<br/>')))
			}else{ HTML('There are no input data to analyze. Please check the input file')}
		}else{ 
			HTML('This is not an .rda file. Please check the input file')
		}	
	}
  })
})


# save the inputData object with uploaded files for further analysis when the Save Input Object button is clicked
output$saveInputObj <- downloadHandler(
	filename='inputData.rda',
	content=function(file){
#		load('inputData.rda')
		if (exists('mergedData', envir = .GlobalEnv)) {
			save(mergedData,productiveReadCounts,ntData,file=file)
		}else{
			print('No data to save')
			output$message = renderText('No data to save')
			emptyObj = NULL
			save(emptyObj,file=file)
		}
})

	# run analysis with the Run Analysis button is clicked
observeEvent(input$runAnalysis,{
	output$message = renderText('Analysis is running...')
	# remove results of the previous analysis
	if(exists('analysisRes', envir = .GlobalEnv)) rm('analysisRes', envir = .GlobalEnv)
	# check if there are objects to run analysis
	if(exists('mergedData', envir = .GlobalEnv) && exists('ntData', envir = .GlobalEnv))
	{
	   # run analysis on aa level data
	   # load object with input data
	   if (!input$nuctleotideFlag){
				sampNames = names(mergedData)
				obj = mergedData
		}else{ # run analysis on nucleotide level data if 'Use nucleotide level data' checkbox is selected
				sampNames = names(ntData)
				obj = ntData
		}
		sampForAnalysis = setdiff(sampNames, c(input$excludeSamp,input$refSamp, input$baselineSamp))
		# get baseline frequencies to find a threshold
		baselineSamp = input$baselineSamp
		# if baseline field is empty, use refSamp to get threshold
		if(input$baselineSamp == 'None') baselineSamp = input$refSamp	
		
		#===================
		clonesToTest = NULL
		if(!input$ignoreBaseline)
		{
			baselineFreq = getFreq(clones = names(obj[[baselineSamp]]),obj,productiveReadCounts,samp=baselineSamp)
			clonesToTest = rownames(baselineFreq)[which(baselineFreq[,1] > getFreqThreshold(as.numeric(input$nCells),input$prob)*100)] # 
		}
				
		# if the comparison to reference should be included in the analysis
		if(input$compareToRef)
		{
			if (input$refSamp == 'None') 
			{
				output$message = renderText('There is no reference. Please select a reference sample')
			}else{
			
				# create comparing pairs (to refSamp)
				compPairs = cbind(sampForAnalysis,rep(input$refSamp,length(sampForAnalysis)))
				# if the Ignore baseline flag is on
				if(input$ignoreBaseline) clonesToTest = NULL
				# run Fisher's test 
				if(nrow(compPairs) == 1)
				{
					res = list(runFisher(compPairs,obj,productiveReadCounts, clones = clonesToTest, nReadFilter = c(as.numeric(input$nReads),0)))
				}else{
					res = apply(compPairs,1,runFisher,obj,productiveReadCounts, clones = clonesToTest, nReadFilter = c(as.numeric(input$nReads),0))
				}
		#			browser()
				if (!is.null(res))
				{
					names(res) = apply(compPairs,1,paste,collapse = '_vs_')
					output$message = renderText('Analysis is done. Click Download Results to save the results')
					assign('analysisRes',res, envir = .GlobalEnv)  
				}	else{
					output$message = renderText('There are no clones to analyze. Try to reduce confidence or the number of templates 1')
				}
			}
		}else{ 
		# if there is no comparison with the reference, then compare only within conditions
			if(input$ignoreBaseline) clonesToTest = NULL
			# take only clones that have the number of reads more then nReads and compare with top second and top third conditions
			fisherRes = compareWithOtherTopConditions(obj, productiveReadCounts, sampForAnalysis, 
				nReads = as.numeric(input$nReads), clones = clonesToTest)
			if (!is.null(fisherRes))
			{
				output$message = renderText('Analysis is done. Click Download Results to save the results')
				assign('analysisRes',fisherRes, envir = .GlobalEnv)  
			}	else{
				output$message = renderText('There are no clones to analyze. Try to reduce the number of templates')
			}			
		}
	}else{	
		output$message = renderText('There are no data to analyze. Please load files')
	}

})
	

 	# save results with selected thresholds when the Download Results button is clicked
	output$saveResults <- downloadHandler(
		filename=function() {
			paste0('analysisRes_',Sys.Date(),'.xlsx')
		},
		content=function(file){
		if (exists('analysisRes', envir = .GlobalEnv))
		{
			# check if analysis was done on aa or nt level
			if (input$nuctleotideFlag) obj = ntData else obj = mergedData
			sampForAnalysis = setdiff(names(obj), c(input$excludeSamp,input$refSamp, input$baselineSamp))
			
			#======================
			# get positive clones
			posClones = NULL
			if(input$compareToRef) #if there is comparison to the ref sample
			{
			posClones = getPositiveClones(analysisRes, obj, productiveReadCounts, samp = sampForAnalysis,
				orThr = as.numeric(input$orThr), fdrThr=as.numeric(input$fdrThr), nReads = as.numeric(input$nReads))
			}else{ # if there is no comparison to ref sample			
				posClones = getPositiveClonesFromTopConditions(analysisRes, orThr = as.numeric(input$orThr), fdrThr=as.numeric(input$fdrThr))
			}
			
			#===============================
			baselineSamp = input$baselineSamp
			if(input$baselineSamp == 'None') baselineSamp = NULL
			refSamp = input$refSamp
			if(input$refSamp == 'None') refSamp = NULL
			
			
			#========================
			# create object with results to write to Excel
			#=======================
			tablesToXls = vector(mode = 'list')
			# if there is no positive clones
			if (length(posClones)==0)
			{
				output$message = renderText('There are no positive clones. Try to adjust thresholds')
				# add a sheet to the output with significant clones comparing to reference 
				#clones = rownames(resTable)[which(resTable[,'significant_comparisons'] == 1)]
				tablesToXls$summary = data.frame('There is no positive clones', row.names = NULL, check.names = F)
			}else{
			# if there are positive clones, save them in to Excel file
				# create table with results 
				tablesToXls = createPosClonesOutput(posClones, obj, productiveReadCounts, refSamp, baselineSamp, addDiff = F)
			}
			
			#===================
			# add the ref_comparison_only sheet
			#===================
			if(input$compareToRef) #if there is comparison to the ref sample
			{
				# create a table with results
				resTable = createResTable(analysisRes,obj,productiveReadCounts, orThr = as.numeric(input$orThr), 
					FDR_threshold = as.numeric(input$fdrThr), saveCI = F)
				if (!is.null(resTable))
				{
					# make a numeric table to save in an Excel spreadsheet	

					refCompRes = resTable[,2:ncol(resTable)]
					refCompRes = t(do.call('rbind',refCompRes))
					rownames(refCompRes) = rownames(resTable)			
					# find clones that are significant in one comparison only
					#clones = rownames(resTable)[which(resTable[,'significant_comparisons'] == 1)]
					
					# table with results of comparison to the reference sample only
					tablesToXls$ref_comparison_only = data.frame(refCompRes,check.names = F)
				}else{
					tablesToXls$ref_comparison_only = data.frame(res = 'There is no significant clones')
				}
			}
			#============
			# add a sheet with parameters
			#============
			s = c(refSamp, baselineSamp,sampForAnalysis)
			# add the baseline threshold percentage and the corresponding number of templates in baseline sample
			baselineThrNames = baselineThrVal = NULL
			if(!input$ignoreBaseline)
			{
				 baselineThrNames =c('confidence','nCells','baseline_threshold_percent','baseline_threshold_templates')
				 freq = getFreqThreshold(as.numeric(input$nCells),input$prob)
				if(input$baselineSamp == 'None') baselineSamp = input$refSamp
				 baselineThrVal = c(input$prob,input$nCells,freq*100,floor(round(freq*productiveReadCounts[baselineSamp])))
			}
			param = c('Reference_samp','Baseline_sample','Excluded samples','Compare to reference','nTemplates_threshold','FDR_threshold','OR_threshold','Ignore_baseline_threshold',
				'Nucleotide_level',baselineThrNames,'nAnalyzedSamples',paste(s, 'nTemplates',sep = '_'))
			value = c(toString(refSamp), toString(baselineSamp),paste(input$excludeSamp, collapse = ', '),input$compareToRef, input$nReads,input$fdrThr, input$orThr,input$ignoreBaseline,
				input$nuctleotideFlag, baselineThrVal, length(sampForAnalysis),productiveReadCounts[s])
				
			tablesToXls$parameters = data.frame(param, value)
			
			#========
			# save results to xlsx
#			tablesToXls$input = data.frame(isolate(reactiveValuesToList(input)))
			WriteXLS('tablesToXls', file, SheetNames = names(tablesToXls), row.names = T)

			output$message = renderText('The results are saved')
		}else{
		
			output$message = renderText('Please click the Run Analysis button')

		}
	})
		
	# save heatmaps with selected thresholds when the Download Heatmaps button is clicked
	output$saveHeatmaps <- downloadHandler(
		filename=function(){
		 paste0('heatmaps_',Sys.Date(),'.pdf')},
		content=function(file){
		if (exists('analysisRes', envir = .GlobalEnv))
		{
			# check if analysis was done on aa or nt level
			if (input$nuctleotideFlag) obj = ntData else obj = mergedData
			sampForAnalysis = setdiff(names(obj), c(input$excludeSamp,input$refSamp, input$baselineSamp))
			posClones = getPositiveClones(analysisRes, obj, productiveReadCounts, samp = sampForAnalysis,
				orThr = as.numeric(input$orThr), fdrThr=as.numeric(input$fdrThr))

			# create a table with results
			resTable = createResTable(analysisRes,mergedData,productiveReadCounts, orThr = as.numeric(input$orThr), 
				FDR_threshold = as.numeric(input$fdrThr), saveCI = F)
			# make heatmap with all significant clones
			posClones = list(rownames(resTable))
			names(posClones) = 'All_significant'

			makeHeatmaps(posClones, obj, productiveReadCounts, samp = sampForAnalysis, refSamp = input$refSamp, 
				fileName = file, size = 7)	
			output$message = renderText('The heat map is saved')

		}
	})
}

ui <- fluidPage(
# layout with tabs
tabsetPanel(
	# tab with FEST analysis
	tabPanel("FEST data analysis",
	# headerPanel("FEST data analysis"),
	 sidebarLayout(
	   # the left side panel 
		sidebarPanel(
		# select the format of input files
		selectInput('inputFiles', 'Select an input format', choices = c('VDJtools files','Adaptive files','an R object with data'), selected = NULL,
			multiple = FALSE,selectize = TRUE, width = NULL, size = NULL),
		
		# if a previously saved R object with the input data will be uploaded
		 conditionalPanel(
			  condition = "input.inputFiles == 'an R object with data'",
			  fileInput('inputObj', 'Upload an R object with data (.rda)', multiple = FALSE,
					  accept=c('rda', '.rda'))
		 ),
		# show the fileInput control depending on the selected values of inputFiles
		# if files with raw data will be uploaded
		 conditionalPanel(
			  condition = "input.inputFiles == 'VDJtools files' | input.inputFiles == 'Adaptive files'",
			  fileInput('sourceFiles', 'Upload FEST files (.tsv)', multiple = TRUE,
					 accept=c('text/tsv', 'text/comma-separated-values,text/plain', '.tsv')),
			 downloadButton('saveInputObj', 'Save Input Object')
		 ),

		tags$hr(),
		# check box that controls the type of analysis - if the comparison with a reference should be performed or not
		checkboxInput('compareToRef','Compare to reference', value = TRUE),
		
		selectInput('refSamp', 'Select a reference sample', choices = list('None'), selected = NULL, multiple = FALSE,
			selectize = TRUE, width = NULL, size = NULL),
		
		selectInput('baselineSamp', 'Select a baseline sample', choices = list('None'), selected = NULL, multiple = FALSE,
			selectize = TRUE, width = NULL, size = NULL),
		selectInput('excludeSamp', 'Select samples to exclude', choices = list('None'), selected = NULL, multiple = TRUE,
			selectize = TRUE, width = NULL, size = NULL),
		textInput('nReads', 'Specify the minimal number of templates (increase in case of large files)', value = "1", width = NULL, placeholder = NULL),
		checkboxInput('ignoreBaseline','Ignore baseline threshold'),
		conditionalPanel(
			condition = "input.ignoreBaseline == false",
			numericInput('prob', 'Specify clone confidence', value = 0.99, max = 1, min = 0, step= 0.01)
		),
#		numericInput('prob', 'Specify clone confidence', value = 0.99, max = 1, min = 0, step= 0.01),
		conditionalPanel(
			condition = "input.ignoreBaseline == false",
		textInput('nCells', 'Estimated number of cells per well', value = "100000", width = NULL, placeholder = NULL)
		),
		# specify if analysis should be performed on nucleotide level data
		checkboxInput('nuctleotideFlag','Use nucleotide level data'),
		actionButton('runAnalysis', 'Run Analysis', width = '60%'),
		tags$hr(),
		numericInput('fdrThr', 'Specify FDR threshold', value = 0.05, max = 1, step= 0.01),
		textInput('orThr', 'Specify odds ratio threshold', value = "5", width = NULL, placeholder = NULL),
 		downloadButton('saveResults', 'Download Results'),
		downloadButton('saveHeatmaps', 'Create heatmaps')
	#		downloadLink("downloadData", "Download Results")
		),

		
		# the main panel
		mainPanel(
		 # textOutput('contents'),
		  htmlOutput("message")
		)
	 )
	 ),# end tabPanel with FEST analysis
	 
	 #===========
	 # the tab with manual
	 tabPanel("User\'s manual", includeHTML('userManual.html'))
	 
 ) # end tabsetPanel
)

shinyApp(ui = ui, server = server)
