library(shiny)
library(bslib)
library(tidyverse)
library(plotly)
library(sangeranalyseR)


ui <- page_sidebar(
  title = "Sanger Sequence Viewer v0.1",
	sidebar=sidebar(fileInput("sanger_file", "Step 1: Upload an ab1 file",
		multiple = FALSE,
		accept = c(".ab1","application/octet-stream")),
		sliderInput("quality","Step 2: Choose a desired quality for visualization:", step=1, min=0, max=60, value=20, dragRange = F),
		textOutput("interact_text"),
		#sliderInput("zoom","Zoom Level:", step=0.5, min=1, max=50, value=1, dragRange = F),
		textOutput("n_text"),
		uiOutput("range"),
		downloadButton("download_fasta", "Step 6: Download .fasta")
	),
  p("Instructions: Use this application to open .ab1 sanger sequencing files. You can visually expore the chromatogram traces with associated sequence quality information, modify bases to change them to unknown N, and download your resulting file. Follow the steps sequentially in the side bar. Updates coming for this section..."),
	plotlyOutput(outputId = "sangerPlot")
)

server <- function(input, output, session) {
	# store  whether plot has been clicked
	data <- reactiveValues(click_data=NULL, signals=NULL, zoom_Data=NULL)
	
	#metadata <- reactiveValues(sangerRead=NULL, peak_metadata=NULL)
	
	observe({
		# if click less than 0?
		data$click_data <- event_data("plotly_click", source="A", priority="input")
		print("click detected")
	})
	# observe({
	# 	data$zoom_data <-event_data("plotly_relayout", source="A", priority = "input")
	# 	print("relay detected")
	# })
	
	# observeEvent(input$zoom, {
	# 	print("hi")
	# })

	observeEvent(input$sanger_file, {
		
		sangerRead <<- new("SangerRead",inputSource = "ABIF",readFeature= "Forward Read", readFileName = input$sanger_file$datapath,
				  M1TrimmingCutoff      = 0,
				  M2CutoffQualityScore  = NULL,
				  M2SlidingWindowSize   = NULL,
				  baseNumPerRow         = 100,
				  heightPerRow          = 200,
				  signalRatioCutoff     = 0.33,
				  showTrimmed           = TRUE)
		
		# first extract chromatogram signals
		signals <<- sangerRead@traceMatrix %>% as.data.frame() %>% 
			mutate(Time=1:nrow(.)) %>% dplyr::rename(A=V1, C=V2, G=V3, `T`=V4) %>%
			pivot_longer(!Time, names_to = "Base", values_to = "Signal") %>%
			mutate(Base=factor(Base, levels=c("A","C","G","T","N")))
		
		peak_metadata <<- data.frame(
			Base=sangerRead@primarySeq%>% as.character() %>% str_split(., pattern="") %>% unlist(),
			quality=sangerRead@QualityReport@qualityPhredScores,
			peaksA=sangerRead@peakPosMatrix[,1],
			peaksC=sangerRead@peakPosMatrix[,2],
			peaksG=sangerRead@peakPosMatrix[,3],
			peaksT=sangerRead@peakPosMatrix[,4])%>%
			mutate(Time = case_when(Base=="A" ~ peaksA,
						Base=="C" ~ peaksC,
						Base=="G" ~ peaksG,
						Base=="T" ~ peaksT,
						length(is.na(c(peaksA,peaksC,peaksG,peaksT))[is.na(c(peaksA,peaksC,peaksG,peaksT))==T])!=0 ~ mean(c(peaksA,peaksC,peaksG,peaksT),na.rm=T),
						T ~ NA))%>%
			select(-c(peaksA,peaksC,peaksG,peaksT))
		
		peak_metadata$original_base <<- peak_metadata$Base
		output$download_fasta <- downloadHandler(
			filename = function() {
				sub(".ab1$", ".fasta", input$sanger_file$name)
			},
			content = function(out_file) {
				bases <- paste(peak_metadata[input$range[1]:input$range[2],]$Base, collapse = '')
				fileConn<-file(out_file)
				writeLines(paste(paste0(">",sub(".ab1$", "", input$sanger_file$name)),sep="\n",bases), fileConn)
				close(fileConn)
			}, contentType = "text/plain"
		)
	})
	
	output$interact_text <- renderText({"Step 3: Scroll to zoom in and out, and drag to move the graph to explore the sequences."})
	output$n_text <- renderText({"Step 4: Click on bases to change them to N if you are not confident in the base call."})
	
	output$range <- renderUI({
		numericRangeInput("range", "Step 5: Select range of bases to export:",
				  min = 0, max = 0, step = 1, value=c(0,0))
	})
	
	
	output$sangerPlot <- renderPlotly({
		req(input$sanger_file)
		
		if(isTruthy(input$range[1]==0)){
			output$range <- renderUI({
				numericRangeInput("range", "Step 5: Select range of bases to export:",
				min = 1, max = nrow(peak_metadata), step = 1,
				value = c(1, nrow(peak_metadata)))
			})
		}

		if(isTruthy(input$range[1]>input$range[2])){
			output$range <- renderUI({
				numericRangeInput("range", "Step 5: Select range of bases to export:",
						  min = 1, max = nrow(peak_metadata), step = 1,
						  value = c(input$range[2], input$range[1]))
			})
		}

		if(!is.null(data$click_data$x) & isTruthy(data$click_data$y<0)){
			selected_base <- peak_metadata %>% filter(Time==data$click_data$x)
			print(selected_base)
			if(selected_base$Base == selected_base$original_base){
				peak_metadata <<- peak_metadata %>%
					mutate(Base = ifelse(Time==data$click_data$x, "N", Base))
			} else if(isTruthy(data$click_data$y==-max(signals$Signal)*0.1)){
				peak_metadata <<- peak_metadata %>%
					mutate(Base = ifelse(Time==data$click_data$x, selected_base$original_base, Base))
			}
			# change y to greater than 0 so that this if statement isn't accidentally retriggered
			data$click_data$y <- 1
		}
		# x axis is time but the base positions are essentially randomly placed on time
		# this function samples the base at 10 times to ascertain the correct base
		# get_break_labels <- function(){
		# 	# true breaks
		# 	breaks <- seq(0,max(signals$Time),max(signals$Time)/nrow(peak_metadata)*10)
		# 	break_labels <- c()
		# 	for(break_time in breaks){
		# 		break_labels <- c(break_labels, peak_metadata %>% filter(Time<=break_time) %>% summarise(num_bases=n()) %>% pull(num_bases))
		# 	}
		# 	return(break_labels)
		# }
		get_break_positions <- function(){
			breaks <- seq(10,nrow(peak_metadata),10)
			break_positions <- c()
			for(base_num in breaks){
				break_positions <- c(break_positions, peak_metadata[base_num,3])
			}
			return(break_positions)
		}
		
		plot_obj<-signals %>% 
			ggplot(aes(x=Time, y=Signal,color=Base)) +
			geom_bar(data=peak_metadata, stat="identity", position="identity", aes(x=Time, y=quality*max(signals$Signal)/max(quality)/5, fill=ifelse(quality<=input$quality,"fail","pass")), alpha=.8, width = 6, color="white")+
			geom_line() +
			geom_vline(xintercept = peak_metadata[input$range[1],3]-7, alpha=.5)+
			geom_vline(xintercept = peak_metadata[input$range[2],3]+7, alpha=.5)+
			geom_hline(yintercept = input$quality*max(signals$Signal)/max(peak_metadata$quality)/5, alpha=.5, linetype="dashed", color="grey50")+
			geom_text(data=peak_metadata, aes(label=Base, y=-max(signals$Signal)*
							  	case_when(Base=="A"~0.02,Base=="C"~0.04,Base=="G"~0.06,Base=="T"~0.08,Base=="N"~0.1), x=Time)) +
			scale_x_continuous(labels = seq(10,nrow(peak_metadata),10),
					   breaks = get_break_positions())+#seq(0,max(signals$Time),max(signals$Time)/nrow(peak_metadata)*10))+
			scale_color_manual(values=c("A"="chartreuse3","G"="black","C"="blue","T"="red","N"="grey50")) +
			scale_fill_manual(values=c("gold","grey90")) +
			xlab("Base") + ylab("Signal Intensity") +
			qiime2R::theme_q2r() + theme(legend.position="none") # !!!!!!!!!!!!!!!remove dependency on theme_q2r()
		
		plot <- ggplotly(plot_obj, dynamicTicks = F, tooltip = F, source = "A")%>%
			config(plot, scrollZoom = T, displayModeBar=F) 
		# keep same zoom if plotly click event is triggered (which redraws the plot)
		plot <- plot %>%
			layout(yaxis = list(fixedrange = T), dragmode="pan")
		if(is.null(data$click_data$x)){
			plot <- plot %>%
				layout(yaxis = list(fixedrange = T), dragmode="pan", xaxis = list(range=list(0,peak_metadata[100,3])))
		}else{
			if(is.null(data$zoom_data)){
				plot <- plot %>%
					layout(yaxis = list(fixedrange = T), dragmode="pan", xaxis = list(range=list(0,peak_metadata[100,3])))
			}else{
				plot <- plot %>%
					layout(yaxis = list(fixedrange = T), dragmode="pan", xaxis = list(range = list(
						data$zoom_data[[1]], data$zoom_data[[2]])))
			}
						
		}
		print("plotting")
		return(plot)
	})
	
	observeEvent(event_data("plotly_relayout"), {
		data$zoom_data <-event_data("plotly_relayout", source="A", priority = "input")
	})
	
	
}

shinyApp(ui = ui, server = server)





