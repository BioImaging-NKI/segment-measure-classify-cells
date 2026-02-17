/*
 * Macro to quantify the fluorescent signals in the nuclei/cells for all channels
 * 
 * Brief workflow:
 * - Nuclei/cells are detected using the StarDist/Cellpose convolutional neural network model for fluorescent nuclei.
 * - Nuclear ROIs can be filtered on size and eroded to prevent edge effects.
 * - Background signal in the measurement channel is measured (several options).
 * - N.B. Intensity values shown in the Results table are already background-subtracted.
 * 
 * Input: a folder containing 2D images with at least 2 channels. Multi-series microscopy format files are also supported.
 * 
 * Required Fiji update sites:
 * - StarDist
 * - CSBDeep
 * - TensorFlow
 * - CLIJ
 * - CLIJ2
 * - SCF MPI CBG
 * - PTBIOP
 * 
 * Author: Bram van den Broek, Netherlands Cancer Institute, 2020-2026
 * b.vd.broek@nki.nl
 *
 * N.B. This script heavily relies on (mostly) CLIJ2 (https://clij.github.io/), StarDist (https://imagej.net/StarDist) and Cellpose (https://github.com/MouseLand/cellpose).
 * If you use this script in a publication, please cite them appropriately.
 * 
 * 
 * ----------------------------------
 * Changelog, from v0.9 onwards:
 * ----------------------------------
 * v1.0:
 * - Added downsampling possibility before StarDist (better for high resolution images)
 * - Nuclei numbers visible in output images 
 * 
 * v1.1:
 * - Added metrics: nucleus area & total intensity
 * - Added an option to exclude nuclei touching image edges
 * - Perform writing of nuclei numbers in output image with CLIJ2 (different method, faster)
 * - Added some visualization options
 *
 * v1.2, March 2021:
 * - Fully omit the RoiManager by:
 *   * Creating outlines (as overlay) using CLIJ
 *   * Writing cell numbers as overlay at the center of mass of the detected cells
 * - General improvements
 * 
 * v1.3, April 2021:
 * - Support for multi-series files (opened with Bio-Formats)
 * 
 * v1.4, July 2021:
 * - Fixed a critial bug: the GPU memory was not cleared after pushCurrentSlice(), causing incorrect results in all channels >1.
 * 
 * v1.5, March 2022:
 * - don't remember :-(
 * 
 * v2.0, June 2022:
 * - Added possibility to segment cells using cellpose
 * 
 * v2.1, November 2025:
 * - bugfixes concerning downscaling for StarDist
 * - Channels and slices are flipped if only one channel is detected
 * 
 * v3.0, December 2025:
 * - Added Cell classification (positive/negative)
 * 
 * v3.3, January 2026:
 * - Upgraded to new BIOP Cellpose wrapper
 * - Added Median intensity measurements
 * - Possibility to measure background as Mean or Median of a certain lower percentile pixels
 * 
 * v3.4,January 2026:
 * - Fixed bugs for classification with multiple images
 * - Added automatic downscaling for StarDist
 * - Added a top percentile measurement  
 * - Added differential treatmentper channel for classifying
 * - Added a Classification results table with % positive cells
 * 
 */


#@ File[] 	inputFiles				(label = "Input files", style = "file") 
#@ File		outputFolder			(label = "Output folder", style = "directory")
#@ String	message1				(value="N.B. The output folder should NOT be a subfolder of the input folder.", visibility="MESSAGE", required=false, persist=false) 
#@ String	fileExtension			(label = "Process files with extension", value = ".tif")
#@ String	segmentationMethod 		(label = "Segmentation method", choices={"StarDist nuclei segmentation","Cellpose segmentation"}, style="listBox") 
#@ Integer	nucleiChannel			(label = "Nucleus marker channel (0 if not present)", value = 1, min = 0)
#@ Integer	cytoChannel				(label = "Cytoplasm marker channel (0 if not present)", value = 2, min = 0)
#@ Double	medianRadius_setting 	(label = "Pre-segmentation median filter radius (um)", description="Filtering creates a more round segmentation and can be beneficial for very noisy data.", value = 0.0, min = 0.0)

#@ String	message2				(value="<html><p style='font-size:14px; color:#3366cc; font-weight:bold'>StarDist settings</p></html>", visibility="MESSAGE", required=false, persist=false)
#@ Double	downsampleFactorStarDist(label = "Stardist nuclei rescaling factor [1-n], 0 for automatic, 1 for no rescaling", value = 0, min=0, description="Stardist is trained on medium resolution images, and generally performs well on images with pixel sizes around 0.2-0.5 µm.\nSet to 0 for automatic rescaling the nuclei to an optimal pixel size of 0.25 µm, or put any other number for manual control of the rescaling.")
#@ Double	probabilityThreshold	(label = "StarDist nucleus probability threshold (0-1, higher is more strict)", value = 0.5, description="The macro uses StarDist's pre-trained deep learning model for nuclei prediction. Higher thresholds will eliminate less likely nuclei.")

#@ String	message3				(value="<html><p style='font-size:14px; color:#3366cc; font-weight:bold'>Cellpose settings</p></html>", visibility="MESSAGE", required=false, persist=false)
#@ String	cellPoseModel			(label = "Cellpose model", choices={"nuclei","cyto", "cyto2", "cyto2_omni", "bact_omni"}, style="listBox")
#@ Integer	CellposeDiameter		(label="Estimated cell diameter (in pixels, 0 for autodetect)", value=50)
#@ Double	cellposeFlowThreshold	(label="Cell flow error threshold (Cellpose default=0.4, higher -> more cells)", style="scroll bar", value=0.4, min=0.0, max=3.0, stepSize=0.1)
#@ Double	cellposeProbability		(label="Cell probability threshold (default=0, higher -> smaller cell area)", style="scroll bar", value=0.0, min=-6.0, max=6.0, stepSize=0.25) 

#@ String	message4				(value="<html><p style='font-size:14px; color:#3366cc; font-weight:bold'>Other settings</p></html>", visibility="MESSAGE", required=false, persist=false)
#@ Boolean	excludeEdges			(label = "Exclude nuclei/cells touching the edges of the image", value = false)
#@ Integer	minNucleusSize_setting	(label = "Remove nulei/cells with diameter smaller than (um)", value = 4, min = 0)
#@ Integer	maxNucleusSize_setting	(label = "Remove nulei/cells with diameter larger than (um)", value = 100) 
#@ Float	shrinkSize_setting		(label = "Shrink segmented nuclei/cells with (units) - negative means expand", value = 0.5)
#@ Integer	percentile_setting		(label = "Measure mean of top # percentile pixels", value = 25, min=1, max=99)

#@ String	message5				(value="<html><p style='font-size:14px; color:#3366cc; font-weight:bold'>Background handling</p></html>", visibility="MESSAGE", required=false, persist=false)
#@ String	background_subtraction	(label = "Background subtraction method", choices={"Calculate value for every channel and image", "Automatic rolling ball", "Manual fixed value everywhere"}, style="listBox")
#@ Integer	rollingBallRadius_setting (label = "Rolling ball radius (units) (if applicable)", value = 100)
#@ Integer	background				(label = "Manual background value (if applicable)", value = 0)

#@ String	message6				(value="<html><p style='font-size:14px; color:#3366cc; font-weight:bold'>Cell classification</p></html>", visibility="MESSAGE", required=false, persist=false)
#@ Boolean	classifyCells_boolean	(label = "Classify cells (positive/negative)?", value = false)
//#@ String	classifyChannelsString	(label="Channels to classify (comma-separated list - empty for all channels)", value="")
//#@ String	metric					(label = "Classify using metric", choices={"Mean", "StdDev", "Total Int", "Percentile"}, style="list")
#@ String	thresholdChoice			(label = "Thresholding", choices={"automatic","manual"}, style="radioButtonHorizontal", value="automatic")
#@ String	autoThresholdMethod		(label = "Threshold method (if automatic)", choices = {"Default", "Huang", "Intermodes", "IsoData", "IJ_IsoData", "Li", "MaxEntropy", "Mean", "MinError", "Minimum", "Moments", "Otsu", "Percentile", "RenyiEntropy", "Shanbhag", "Triangle", "Yen"}, style="list", value="Otsu")
#@ String	classifyChoice			(label = "Classify", choices={"each image separately","after running all images"}, style="listbox", value="after running all images")

#@ String	message7				(value="<html><p style='font-size:14px; color:#3366cc; font-weight:bold'>Visualization options</p></html>", visibility="MESSAGE", required=false, persist=false)
#@ Boolean	saveImages				(label = "Save output images with overlayed nucleus outlines", value = false)
#@ Boolean	thickOutlines			(label = "Use thick outlines", value = false)
#@ Integer	labelOpacity			(label = "Opacity% of nucleus outlines overlay", value = 100, min=0, max=100)
#@ String	overlayChoice			(label = "Overlay display option", choices={"Cell numbers and outlines as overlay", "Cell numbers and outlines imprinted as pixels (RGB output)"}, style="listBox")
#@ Boolean	addNumbersOverlay		(label = "Add cell numbers in output image", value = false)
#@ Integer	labelFontSize			(label = "Cell numbers font size", value = 12, min = 1)
#@ ColorRGB	fontColor				(label = "Cell number font color", value="yellow")
#@ Boolean	hideImages				(label = "Hide images during processing", value = false)

version = 3.4;

//Settings outside the dialog
backgroundPercentile = 0.33;			// For automatic background value per image calculation 
maxTileSize = 4096;						// Maximum StarDist tile size
measureMedian = true;					// Additionally measure the median value per cell
measurePercentile = true;				// Additionally measure the mean of the top xx percentile pixels per cell
backgroundMedianInsteadOfMean = false;	// Determine which measurement to use
convertLabelsToROIs = true;				// Convert labels to ROIs before saving

//Default metrics for classifying  cells
defaultCheckboxes = newArray(0,1,1,1);
defaultMetrics = newArray("Mean", "Percentile", "Mean", "Mean");

var nrOfImages = 0;
var current_image_nr = 0;
var processtime = 0;
var nrCells = 0;
var nrChannels = 0;
var nrHistograms = 0;
var labelmap = "";	//labelmap of the last image
colors = newArray("Blue", "Cyan", "Green", "Red", "Magenta");
var n = 0;
if(!File.exists(outputFolder)) File.makeDirectory(outputFolder);
outputSubfolder = outputFolder;		// Initialize this variable

nrBins = 81;
minLogInt = 1;
maxLogInt = 5; 

saveSettings();

run("Bio-Formats Macro Extensions");
run("Labels...", "color=white font=" + labelFontSize + " show draw");
run("Set Measurements...", "area mean median min stack redirect=None decimal=3");
run("Input/Output...", "jpeg=85 gif=-1 file=.tsv use_file copy_row save_column save_row");
if(nImages>0) run("Close All");
print("\\Clear");
run("Clear Results");
setBatchMode(true);

resultsTable = "All Results";
if(isOpen("Results_all_files.tsv")) close("Results_all_files.tsv");
if(!isOpen(resultsTable)) Table.create(resultsTable);
else Table.reset(resultsTable);
Table.update(resultsTable);
Table.showRowIndexes(true);
Table.setLocationAndSize(100, 100, 1000, 500);

classTable = "Classification results";
if(classifyCells_boolean == true && !isOpen(classTable)) Table.create(classTable);
else if(classifyCells_boolean == true) {
	Table.reset(classTable);
	Table.update;
}

for(f=0; f<inputFiles.length; f++) {
	processFile(inputFiles[f]);
	current_image_nr++;
}

//Classify after measuring all images 
//if(classifyCells_boolean == true && classifyChoice == "after running all images") {
//	if(classifyChannelsString != "") {
//		classifyChannels = split(classifyChannelsString, ",");
//		for(i=0;i<classifyChannels.length; i++) classifyChannels[i] = parseInt(classifyChannels[i]);
//	}
//	else {
//		classifyChannels = Array.getSequence(nrChannels);
//		classifyChannels = addScalarToArray(classifyChannels, 1);
//	}
//	image = getTitle;
//	for(i=0; i<classifyChannels.length; i++) {
//		if(classifyChannels[i] != nucleiChannel && classifyChannels[i] != cytoChannel) {
//			nrPositiveCells = classifyCells(image, i+1, nrBins, minLogInt, maxLogInt, metric, thresholdChoice);
//			print("Channel "+classifyChannels[i]+": Positive cells: "+nrPositiveCells+"/"+Table.size(resultsTable));
//		}
//	}
//}

print("\\Update1:Finished processing "+inputFiles.length+" files.");
print("\\Update2:Average speed: "+d2s(current_image_nr/processtime,1)+" files per minute.");
print("\\Update3:Total run time: "+d2s(processtime,1)+" minutes.");
print("\\Update4:-------------------------------------------------------------------------");

selectWindow("Results");
run("Close");
selectWindow(resultsTable);
//Table.rename(resultsTable, "Results");
saveAs("Results", outputFolder + File.separator + "Results_all_files.tsv");

if(classifyCells_boolean == true) {
	selectWindow(classTable);
	Table.save(outputFolder + File.separator + "Classification results.tsv");
}

restoreSettings;


// The actual processing
function processFile(path) {
	run("Close All");

	starttime = getTime();
	print("\\Update1:Processing file "+current_image_nr+"/"+nrOfImages+": " + path);
	print("\\Update2:Average speed: "+d2s((current_image_nr-1)/processtime,1)+" files per minute.");
	time_to_run = (nrOfImages-(current_image_nr-1)) * processtime/(current_image_nr-1);
	if(time_to_run<5) print("\\Update3:Projected run time: "+d2s(time_to_run*60,0)+" seconds ("+d2s(time_to_run,1)+" minutes).");
	else if(time_to_run<60) print("\\Update3:Projected run time: "+d2s(time_to_run,1)+" minutes. You'd better get some coffee.");
	else if(time_to_run<480) print("\\Update3:Projected run time: "+d2s(time_to_run,1)+" minutes ("+d2s(time_to_run/60,1)+" hours). You'd better go and do something useful.");
	else if(time_to_run<1440) print("\\Update3:Projected run time: "+d2s(time_to_run,1)+" minutes. ("+d2s(time_to_run/60,1)+" hours). You'd better come back tomorrow.");
	else if(time_to_run>1440) print("\\Update3:Projected run time: "+d2s(time_to_run,1)+" minutes. This is never going to work. Give it up!");
	print("\\Update4:-------------------------------------------------------------------------");

	run("Bio-Formats Macro Extensions");	//Somehow necessary for every file(?)
	Ext.setId(path);
	Ext.getSeriesCount(nr_series);

	if(endsWith(fileExtension, "tif") || endsWith(fileExtension, "jpg")) {	//Use standard opener
		open(path);
		process_current_series(path);
	}
	else {	//Use Bio-Formats
		for(s = 0; s < nr_series; s++) {
			run("Close All");
			run("Bio-Formats Importer", "open=["+path+"] autoscale color_mode=Default view=Hyperstack stack_order=XYCZT series_"+s+1);
			seriesName = getTitle();
			seriesName = replace(seriesName,"\\/","-");	//replace slashes by dashes in the seriesName
			print(File.getParent(path) + File.separator + seriesName);
	//		outputPath = output + File.separator + substring(seriesNa
			process_current_series(seriesName);
		}
	}
}


function process_current_series(image) {
	image = getTitle();
	getDimensions(width, height, channels, slices, frames);
	if(channels == 1) {
		print("WARNING: single channel detected. Switching slices and channels.");
		run("Re-order Hyperstack ...", "channels=[Slices (z)] slices=[Channels (c)] frames=[Frames (t)]");
	}
	getDimensions(width, height, channels, slices, frames);
	nrChannels = channels;
	Stack.setDisplayMode("grayscale");
	Stack.setChannel(nucleiChannel);
	if (hideImages == false) setBatchMode("show");
	run("Enhance Contrast", "saturated=0.35");

	getPixelSize(unit, pw, ph);
	minNucleusSize = PI*pow((minNucleusSize_setting / pw / 2),2);	//Calculate the nucleus area as if it were a circle
	maxNucleusSize = PI*pow((maxNucleusSize_setting / pw / 2),2);

	medianRadius = medianRadius_setting / pw;
	rollingBallRadius = rollingBallRadius_setting / pw;
	shrinkSize = shrinkSize_setting / pw;
	
	maxSelectedChannel = maxOf(nucleiChannel, cytoChannel);
	if(maxSelectedChannel > channels && current_image_nr==1) {
		exit("The image has less channels ("+channels+") than the selected channel ("+maxSelectedChannel+").");
	}
	if(segmentationMethod == "StarDist nuclei segmentation") detect_nuclei(image, nucleiChannel);
	else if(segmentationMethod == "Cellpose segmentation") detect_cells(image, nucleiChannel, cytoChannel);
	labelmap = getLabelMaps_GPU(image, unit);
//	roiManager("reset");

	//Calculate backgrounds
	selectWindow(image);
	backgrounds = newArray(channels);
	showStatus("Determining backgrounds...");
	if(background_subtraction == "Automatic rolling ball") {
		for(c=1; c<=channels; c++) {
			Stack.setChannel(c);
			showStatus("Subtracting background channel "+c+"...");
			showProgress(c, channels+1);
			run("Subtract Background...", "rolling="+rollingBallRadius+" stack");
			backgrounds[c-1] = "auto subtracted";
		}
	}
	else if(background_subtraction == "Calculate value for every channel and image") {
		for(c=1; c<=channels; c++) {
			backgrounds[c-1] = get_background(image, c, backgroundPercentile);
		}
	}
	else if(background_subtraction == "Manual fixed value everywhere") {
		for(c=1; c<=channels; c++) {
			backgrounds[c-1] = background;
		}
	}
	selectWindow(image);
	Stack.setChannel(nucleiChannel);
	run("Enhance Contrast", "saturated=0.35");
	
	labels_to_ROI_Manager(labelmap);
	measureROIs(image, labelmap, backgrounds, percentile_setting);

//	if(convertLabelsToROIs == true) labels_to_ROI_Manager(labelmap);

	//Save the output image (+ cosmetics)
	filename = File.getNameWithoutExtension(path);
	if (saveImages == true) {
		selectWindow(image);
		if(nucleiChannel != 0) {
			setSlice(nucleiChannel);
			setMetadata("Label", "nucleus marker");
		}
		if(cytoChannel != 0) {
			setSlice(cytoChannel);
			setMetadata("Label", "cytoplasm marker");
		}
		if(overlayChoice == "Cell numbers and outlines imprinted as pixels (RGB output)") run("Flatten", "stack");

		//Optional: save with ROIs
/*
		Overlay.remove;
		Stack.setDisplayMode("composite");
		roiManager("show all without labels");
*/
		saveAs("Tiff", outputFolder + File.separator + filename+"_analyzed");

	if(classifyCells_boolean == true && classifyChoice == "each image separately") {
		setTool("hand");
		image = getTitle;
		getDimensions(width, height, channels, slices, frames);
		metricChoices = newArray("Mean", "StdDev", "Total Int", "Percentile");
		Dialog.createNonBlocking("Classify channels - select metric");
		for(i = 0; i < channels; i++) {
			Dialog.addCheckbox("channel "+i+1, defaultCheckboxes[i]);
			Dialog.addToSameRow();
			Dialog.addChoice("", metricChoices, defaultMetrics[i]);
		}
		Dialog.show();
		classifyChannels = newArray(channels);
		metric = newArray(channels);
		for(i = 0; i < channels; i++) {
			classifyChannels[i] = Dialog.getCheckbox();
			metric[i] = Dialog.getChoice();
		}
		
		selectWindow(classTable);
		nrImagesProcessed = Table.size;
	
		for(i=0; i<channels; i++) {
			if(classifyChannels[i] == true) {
				nrPositiveCells = classifyCells(image, i+1, nrBins, minLogInt, maxLogInt, metric[i], thresholdChoice);
				print("Channel "+i+1+": Positive cells "+nrPositiveCells+"/"+nrCells + " ("+d2s(nrPositiveCells/nrCells*100,1) + "%)");
		
				Table.set("file name", nrImagesProcessed, image, classTable);
				Table.set("%pos ch "+i+1, nrImagesProcessed, d2s(nrPositiveCells/nrCells*100,1), classTable);
				Table.update(classTable);
			}
		}
	}

//	run("Close All");
	endtime = getTime();
	processtime = processtime+(endtime-starttime)/60000;
}


function detect_nuclei(image, nucleiChannel) {
	selectWindow(image);
	run("Duplicate...", "duplicate channels=" + nucleiChannel + " title=nuclei");
	getPixelSize(unit, pw, ph);
	
	nuclei_downscaled = getTitle;
	if(downsampleFactorStarDist != 1) {
		if(downsampleFactorStarDist == 0 && unit == "µm" || unit == "um" || unit == "microns" || unit == "micron") {
			downsampleFactorStarDist = 0.25/pw;	//scale to 0.25 um/pixel
			if(downsampleFactorStarDist < 1 &&  downsampleFactorStarDist > 0.5) downsampleFactorStarDist = 1;	//Do not upsample unless the pixel size is > 0.5 um.
			else if(downsampleFactorStarDist < 0.25) print ("[WARNING] The pixel size is very large ("+pw+" "+unit+"). Stardist may have difficulties segmenting the nuclei.");
			else if(downsampleFactorStarDist > 0.8 && downsampleFactorStarDist < 1.2) downsampleFactorStarDist = 1;	//Too small change - skip rescaling 
		}
		else if (downsampleFactorStarDist == 0) {
			print("[WARNING] Pixel size seems incorrect ("+pw+" "+unit+")"+". Cannot determine the nuclei downsample factor for Stardist. Stardist may have difficulties segmenting the nuclei.");
			downsampleFactorStarDist = 1;
		}
		else if(downsampleFactorStarDist < 0.15) {
			print("[WARNING] Overruling StarDist downsample factor ("+downsampleFactorStarDist+"). Will be set to the minimum value of 0.15");
			downsampleFactorStarDist = 0.15;	//Allow upscaling up to 6.667x
		}
		print("Pixel size: "+pw+" µm\nStarDist downsample factor (automatic): "+downsampleFactorStarDist+" (effective pixel size: "+pw*downsampleFactorStarDist+" µm)\n");		

		run("Duplicate...", "duplicate channels=" + nucleiChannel + " title=nuclei_downscaled");		
//		run("Bin...", "x="+downsampleFactorStarDist+" y="+downsampleFactorStarDist+" bin=Average");
		run("Scale...", "x="+1/downsampleFactorStarDist+" y="+1/downsampleFactorStarDist+" interpolation=Bilinear average process create");
		close(nuclei_downscaled);
		nuclei_downscaled = getTitle();
	}
//	if(medianRadius > 0) run("Median...", "radius="+medianRadius); 
	getDimensions(width, height, channels, slices, frames);

	starDistTiles = pow(floor(maxOf(width, height)/maxTileSize)+1,2);	//Determine the nr. of tiles
	// Run StarDist and output to the ROI manager (creating a label image works only when not operating in batch mode, and that is slower and more annoying.)
	run("Command From Macro", "command=[de.csbdresden.stardist.StarDist2D], args=['input':'"+nuclei_downscaled+"', 'modelChoice':'Versatile (fluorescent nuclei)', 'normalizeInput':'true', 'percentileBottom':'1.0', 'percentileTop':'99.8', 'probThresh':'"+probabilityThreshold+"', 'nmsThresh':'0.3', 'outputType':'ROI Manager', 'nTiles':'"+starDistTiles+"', 'excludeBoundary':'2', 'roiPosition':'Stack', 'verbose':'false', 'showCsbdeepProgress':'false', 'showProbAndDist':'false'], process=[false]");

	//Scale up ROIs
	if(downsampleFactorStarDist != 1) RoiManager.scale(downsampleFactorStarDist, downsampleFactorStarDist, false);
	close("nuclei");
	close(nuclei_downscaled);
}


function detect_cells(image, nucleiChannel, cytoChannel) {
	selectWindow(image);
	if(nucleiChannel != 0) run("Duplicate...", "duplicate channels=" + nucleiChannel + " title=nuclei");
	setBatchMode("show");
	selectWindow(image);
	if(cytoChannel != 0) run("Duplicate...", "duplicate channels=" + cytoChannel + " title=cytoplasm");
	setBatchMode("show");
	if(nucleiChannel != 0 && cytoChannel != 0) run("Merge Channels...", "c1=nuclei c2=cytoplasm create");
	rename("for_Cellpose");
	if(medianRadius > 0) run("Median...", "radius="+medianRadius);
	setBatchMode("show");
	getDimensions(width, height, channels, slices, frames);

	List.setCommands;
	if (List.get("Cellpose ...")!="") {
		//Get Cellpose settings
		envPath = getPref("Packages.ch.epfl.biop.wrappers.cellpose.ij2commands.Cellpose", "env_path");
		envType = getPref("Packages.ch.epfl.biop.wrappers.cellpose.ij2commands.Cellpose", "env_type");
		if(envType == "<null>") envType = "conda";	//Default is conda, but returns <null>
		print("Cellpose environment type: "+envType);
		print("Cellpose environment path: "+envPath);
		cellposeModelPath = "path\\to\\own_cellpose_model";	//currently no custom model
	}
	else exit("[ERROR] Fiji Cellpose plugin not found!\nPlease go to Help -> Update... -> Manage Update Sites,\nmake sure PTBIOP is checked, update and restart Fiji.");

	setBatchMode("exit and display");
	selectImage("for_Cellpose");
	if(cytoChannel  >  0) run("Cellpose ...", "env_path="+envPath+" env_type="+envType+" model=["+cellPoseModel+"] model_path=["+cellposeModelPath+"] diameter="+CellposeDiameter+" ch1="+cytoChannel+" ch2="+nucleiChannel+" additional_flags=[--use_gpu, --flow_threshold, "+cellposeFlowThreshold+", --cellprob_threshold, "+cellposeProbability+"]");
	else run("Cellpose ...", "env_path="+envPath+" env_type="+envType+" model=["+cellPoseModel+"] model_path=["+cellposeModelPath+"] diameter="+CellposeDiameter+" ch1="+cytoChannel+" ch2=0 additional_flags=[--use_gpu, --flow_threshold, "+cellposeFlowThreshold+", --cellprob_threshold, "+cellposeProbability+"]");
	rename("labelmap");
	run("glasbey_on_dark");
	setBatchMode("hide");

	close("for_Cellpose");
}


function getLabelMaps_GPU(image, unit) {
	run("CLIJ2 Macro Extensions", "cl_device=");
	Ext.CLIJ2_clear();
	// In case another GPU needs to be selected:
	//Ext.CLIJ2_listAvailableGPUs();
	//availableGPUs = Table.getColumn("GPUName");
	//run("CLIJ2 Macro Extensions", "cl_device=" + availableGPUs[1]);

	//Create labelmap
	if(segmentationMethod == "StarDist nuclei segmentation") run("ROI Manager to LabelMap(2D)");	//Only for StarDist; Cellpose already produces a labelmap
	else selectWindow("labelmap");
	run("glasbey_on_dark");
	labelmap_nuclei_raw = getTitle();
	Ext.CLIJ2_push(labelmap_nuclei_raw);

	//exclude labels on edges
	if(excludeEdges) Ext.CLIJ2_excludeLabelsOnEdges(labelmap_nuclei_raw, labelmap_nuclei);
	else labelmap_nuclei = labelmap_nuclei_raw;

	//Filter on area
	Ext.CLIJ2_getMaximumOfAllPixels(labelmap_nuclei, nucleiStarDist);	//count nuclei detected by StarDist
	run("Clear Results");
	Ext.CLIJ2_statisticsOfBackgroundAndLabelledPixels(labelmap_nuclei, labelmap_nuclei); //Somehow if you put (image, labelmap) as arguments the pixel count is wrong
	Ext.CLIJ2_pushResultsTableColumn(area, "PIXEL_COUNT");

	Ext.CLIJ2_excludeLabelsWithValuesOutOfRange(area, labelmap_nuclei, labelmap_nuclei_filtered, minNucleusSize, maxNucleusSize);
	Ext.CLIJ2_release(labelmap_nuclei);

	//Shrink nuclei/cells
	if(shrinkSize > 0) {
		Ext.CLIJ2_erodeLabels(labelmap_nuclei_filtered, labelmap_final, shrinkSize, false);
		Ext.CLIJ2_release(labelmap_nuclei_filtered);
	}
	if(shrinkSize < 0) {
		Ext.CLIJ2_dilateLabels(labelmap_nuclei_filtered, labelmap_final, shrinkSize);
		Ext.CLIJ2_release(labelmap_nuclei_filtered);
	}
	else if(shrinkSize == 0) labelmap_final = labelmap_nuclei_filtered;
	
	Ext.CLIJ2_getMaximumOfAllPixels(labelmap_final, nrCells);	//get the number of nuclei after filtering
	run("Clear Results");
	Ext.CLIJ2_closeIndexGapsInLabelMap(labelmap_final, labelmap_final_ordered);	//Renumber the cells from top to bottom
	Ext.CLIJ2_release(labelmap_final);
	Ext.CLIJ2_statisticsOfLabelledPixels(labelmap_final_ordered, labelmap_final_ordered); //Somehow if you put (image, labelmap) as arguments the pixel count is wrong
	print(image + " : " +nucleiStarDist+" nuclei detected by StarDist ; "+nucleiStarDist - nrCells+" nuclei with diameter outside ["+d2s(minNucleusSize_setting,0)+" - "+d2s(maxNucleusSize_setting,0)+"] range "+unit+" were removed.");

	Ext.CLIJ2_pull(labelmap_final_ordered);
	run("glasbey_on_dark");

	if(saveImages == true) {
		//Detect label outlines, dilate them and add as overlay to the original image
		Ext.CLIJ2_detectLabelEdges(labelmap_final_ordered, labelmap_edges);
		if(thickOutlines) {
			Ext.CLIJ2_dilateBox(labelmap_edges, labelmap_edges_dilated);
		}
		else labelmap_edges_dilated = labelmap_edges;
		if(thickOutlines) {	
			Ext.CLIJ2_dilateLabels(labelmap_final_ordered, labelmap_final_ordered_extended, 1);
		}
		else labelmap_final_ordered_extended = labelmap_final_ordered;
		Ext.CLIJ2_mask(labelmap_final_ordered_extended, labelmap_edges_dilated, labelmap_outlines);
		Ext.CLIJ2_release(labelmap_edges);
		if(thickOutlines) Ext.CLIJ2_release(labelmap_final_ordered_extended);
		if(thickOutlines) Ext.CLIJ2_release(labelmap_edges_dilated);
		Ext.CLIJ2_pull(labelmap_outlines);
		run("glasbey_on_dark");
		selectWindow(image);
		run("Add Image...", "image="+labelmap_outlines+" x=0 y=0 opacity="+labelOpacity+" zero");	//Add labelmap to image as overlay
//		Overlay.setPosition(nucleiChannel, 0, 0);
		Overlay.setPosition(0, 0, 0);
		if(addNumbersOverlay) {
			setFont("SansSerif", labelFontSize, "antialiased");
			color = color_to_hex(fontColor);
			setColor(color);
			for (i = 0; i < nrCells; i++) {
				x = getResult("MASS_CENTER_X", i);
				y = getResult("MASS_CENTER_Y", i);
				Overlay.drawString(i+1, x - labelFontSize/2, y + labelFontSize/2);
			}
			Overlay.show;
		}
	}

	return labelmap_final_ordered;
}

function get_background(image, channel, percentile) {
	selectWindow(image);
	Stack.setChannel(channel);
	getRawStatistics(nPixels, mean, min, max, std, histogram);
	total = 0;
	bin=0;
	while (total < nPixels*percentile) {
		total += histogram[bin];
		bin++;
	} 
	setThreshold(0,bin-1);
	if(backgroundMedianInsteadOfMean == true) background = getValue("Median limit");
	else background = getValue("Mean limit");
	resetThreshold();
	return background;
}


function measureROIs(image, labelmap, backgrounds, percentile_setting) {
	//create data arrays
	file_name_image = newArray(nrCells);
	cell_nr_image = newArray(nrCells);
	nuc_area_image = newArray(nrCells);

	for(i=0; i<nrCells; i++) {
		file_name_image[i] = image;
		cell_nr_image[i] = i+1;
	}
	//measure intensity (only background corrected if 'Rolling ball' is selected!)
	selectWindow(image);

	//Measure, get the results and put in 'All Results' table
	if(nrCells != 0) {
		selectWindow(resultsTable);
		nRows = Table.size;
		for(i=nRows; i<nRows+nrCells; i++) {
			Table.set("file name", i, image, resultsTable);
			Table.set("cell nr", i, cell_nr_image[i-nRows], resultsTable);
		}
		for (c = 1; c <= channels; c++) {
			run("Clear Results");
			selectWindow(image);
			Stack.setChannel(c);
			run("Enhance Contrast", "saturated=0.35");
			run("Clear Results");
			Ext.CLIJ2_pushCurrentSlice(image);	//Or maybe copySlice if image is still in the GPU memory?
			Ext.CLIJ2_statisticsOfLabelledPixels(image, labelmap);
			Ext.CLIJ2_release(image);
			mean_image_channel = Table.getColumn("MEAN_INTENSITY", "Results");
			stdDev_image_channel = Table.getColumn("STANDARD_DEVIATION_INTENSITY", "Results");
			nuc_area_image = Table.getColumn("PIXEL_COUNT", "Results");
			
			//Measure median per cell (using ROIs)
			if(measureMedian == true) {
				if(c==1) {
					run("Set Measurements...", "median redirect=None decimal=3");
					selectImage(image);
				}
				roiManager("Measure");
				median_image_channel = Table.getColumn("Median", "Results");
			}

			//Measure mean of highest xx percentile per ROI
			if(measurePercentile == true) {
				percentile_mean_image_channel = newArray(nrCells);
				percentile = 1 - percentile_setting/100;
				selectImage(image);
				setBatchMode("hide");
				for(i=0; i<nrCells; i++) {
					selectImage(image);
					showStatus("Measuring mean of top "+percentile_setting+"% pixels, channel "+c+", "+i+1+"/"+nrCells);
					roiManager("Select", i);
					getRawStatistics(nPixels, mean, min, max, std, histogram);
					total = 0;
					bin = 0;
					while (total < nPixels*percentile) {
						total += histogram[bin];
						bin++;
					} 
					setThreshold(bin-1, max);
					percentile_mean_image_channel[i] = getValue("Mean limit");
				}
				selectImage(image);
				showStatus("");
				run("Select None");
				resetThreshold();
			}
			
			//Add measurements to 'Results_all_files' table
			for(i=nRows; i<nRows+nrCells; i++) {
				//Area
				if(c==1) Table.set("Nucleus area ("+unit+"^2)", i, nuc_area_image[i-nRows]*pw*pw, resultsTable);	//in units^2

				//Mean, median and percentile
				if(background_subtraction == "Automatic rolling ball") {
					Table.set("Mean ch"+c, i, mean_image_channel[i-nRows], resultsTable);
					if(measureMedian == true) Table.set("Median ch"+c, i, median_image_channel[i-nRows], resultsTable);
					if(measurePercentile == true) Table.set(percentile_setting + " percentile mean ch"+c, i, percentile_mean_image_channel[i-nRows], resultsTable);
				}
				else {
					Table.set("Mean ch"+c, i, mean_image_channel[i-nRows] - backgrounds[c-1], resultsTable);
					if(measureMedian == true) Table.set("Median ch"+c, i, median_image_channel[i-nRows] - backgrounds[c-1], resultsTable);
					if(measurePercentile == true) Table.set(percentile_setting + " percentile mean ch"+c, i, percentile_mean_image_channel[i-nRows] - backgrounds[c-1], resultsTable);
				}
				Table.set("StdDev ch"+c, i, stdDev_image_channel[i-nRows], resultsTable);
				
				if(background_subtraction == "Automatic rolling ball") Table.set("Total Int ch"+c, i, nuc_area_image[i-nRows]*pw*pw * (mean_image_channel[i-nRows]), resultsTable);
				else Table.set("Total Int ch"+c, i, nuc_area_image[i-nRows]*pw*pw * (mean_image_channel[i-nRows] - backgrounds[c-1]), resultsTable);
				
				Table.set("Background ch"+c, i, backgrounds[c-1], resultsTable);
			}
		}
	}
	Table.update;
	selectImage(image);
	setBatchMode("show");
}


function color_to_hex(color) {
	colorArray = split(color,",,");
	hexcolor = "#" + IJ.pad(toHex(colorArray[0]),2) + IJ.pad(toHex(colorArray[1]),2) + IJ.pad(toHex(colorArray[2]),2);
	return hexcolor;
}


function classifyCells(image, channel, nrBins, minLogInt, maxLogInt, metric, thresholdChoice) {
	nrPositiveCells = 0;
	selectWindow(image);
	getLocationAndSize(xpos, ypos, width, height);
	if(metric == "Percentile") metric = "" + percentile_setting + " percentile mean";
	
	//Get data from ResultsTable
	selectWindow(resultsTable);
	nRows = Table.size;
	data = newArray(nrCells);
	k=0;
	for(i=nRows-nrCells; i<nRows; i++) {
		data[k] = Table.get(metric+" ch"+channel, i, resultsTable);
		k++;
	}
	dataImageID = arrayTo1DImage(data);

	rename("data ch"+channel);
	min = getValue("Min");
//	run("Subtract...", "value="+min+" slice");
//	run("Add...", "value=1 slice");
	run("Log");								//ln(Int) - intensity distributions are usually log-normal
	run("Divide...", "value=2.302585093");	//Divide by ln10 to arrive at 10log(Int)
	setAutoThreshold(autoThresholdMethod+" dark");
	getThreshold(lowerThresholdLog, upperThresholdLog);
	resetThreshold();

	run("Histogram", "bins="+nrBins+" x_min="+minLogInt+" x_max="+maxLogInt+" y_max=Auto");
	rename("(10log)Histogram of ch"+channel);
	histogram = getImageID();

	close(dataImageID);

	selectImage(histogram);
	getDimensions(hwidth, hheight, hchannels, hslices, hframes);
	setBatchMode("show");
	setOption("DisablePopupMenu", true);
	setLocation(xpos+width, ypos+(1.35*hheight));
	histogramID = getImageID();
	nrHistograms++;

	if(thresholdChoice == "manual" && classifyChoice == "each image separately" || thresholdChoice == "manual" && classifyChoice == "after running all images" && inputFiles.length == 1) {
		selectWindow(image);
		run("Duplicate...", "title=raw_ch"+channel+" duplicate channels="+channel);
		Overlay.remove;
		run("Grays");
		setMinAndMax(1, 2*meanOfArray(data));
		setBatchMode("show");
		
		data_and_background = prependToArray(1, data);	//Add background intensity - will become 0 after taking logartihm
		log_data_and_background = logArray(data_and_background);
		Ext.CLIJ2_pushArray(dataImageGPU, log_data_and_background, data_and_background.length, 1, 1);
		Ext.CLIJ2_replaceIntensities(labelmap, dataImageGPU, labelmap_log_intensities);

		thresholded_labeledges = "thresholded_labeledges";
		setLineWidth(2);
		
		Ext.CLIJ2_reduceLabelsToLabelEdges(labelmap_log_intensities, labeledgesThin_log_intensities);
		Ext.CLIJ2_dilateLabels(labeledgesThin_log_intensities, labeledges_log_intensities, 1);
		Ext.CLIJ2_threshold(labeledges_log_intensities, thresholded_labeledges, lowerThresholdLog);
		Ext.CLIJ2_threshold(dataImageGPU, data_thresholded, lowerThresholdLog);
		Ext.CLIJ2_getSumOfAllPixels(data_thresholded, nrPositiveCells);
		Ext.CLIJ2_pull(thresholded_labeledges);
		run(colors[channel-1]);
		selectImage("raw_ch"+channel);
		Overlay.clear;
		roiManager("show all without labels");
		run("Add Image...", "image=["+thresholded_labeledges+"] x=0 y=0 opacity=100 zero");
		
		//Draw a line at the threshold location
		selectImage(histogram);
		linePosition = Math.map(lowerThresholdLog, minLogInt, maxLogInt, 0, 256);
		//Overlay.drawLine(23 + linePosition, 13, 23 + linePosition, 152);
		Overlay.drawLine(20 + linePosition, 12, 20 + linePosition, 137);
		setJustification("left");
		Overlay.drawString(nrPositiveCells, 25 + linePosition, 25);
		setJustification("right");
		Overlay.drawString(nrCells - nrPositiveCells, 12 + linePosition, 25);
		Overlay.show();
		Overlay.setStrokeColor(colors[channel-1]);

//		waitForUser("Automatic log(intensity) threshold for channel "+channel+" is set at "+lowerThresholdLog+".\nFirst press OK, then adjust the threshold if necessary, using the (log)histogram (top-right), and left-click to continue.");
		modifiers = 0;
		Ext.CLIJ2_reduceLabelsToLabelEdges(labelmap_log_intensities, labeledgesThin_log_intensities);
		Ext.CLIJ2_dilateLabels(labeledgesThin_log_intensities, labeledges_log_intensities, 1);
		while(!isActive(histogramID) || modifiers != 16) {	//left mouse button	
			close(thresholded_labeledges);
			selectImage(histogram);
			getLocationAndSize(hx, hy, hwidth, hheight);
			getCursorLoc(x, y, z, modifiers);
			
			if((x > 0 && x < hwidth) && (y > 0 && y < hheight)) {	//If the cursor is on the histogram window
				lowerThresholdLog = Math.map(x, 20, 276, minLogInt, maxLogInt);
				//Draw a line at the threshold location
				//setColor(colors[i]);
				Overlay.setStrokeColor(colors[channel-1]);
				linePosition = Math.map(lowerThresholdLog, minLogInt, maxLogInt, 0, 256);
				Overlay.clear;
				//Overlay.drawLine(23 + linePosition, 13, 23 + linePosition, 152);
				Overlay.drawLine(20 + linePosition, 12, 20 + linePosition, 137);
				setJustification("left");
				Overlay.drawString(nrPositiveCells, 25 + linePosition, 25);
				setJustification("right");
				Overlay.drawString(nrCells - nrPositiveCells, 12 + linePosition, 25);
				Overlay.setStrokeColor(colors[channel-1]);
				Overlay.show();
			}
			//Create overlay of positive cells
			//This takes some time unfortunately. Try: display ROIs and only change ROI groups? (Need to sort dataImage and ROIs on metric first to know where the cutoff between the groups is).
			if(lowerThresholdLog > 0) {
				Ext.CLIJ2_threshold(labeledges_log_intensities, thresholded_labeledges, lowerThresholdLog);
				Ext.CLIJ2_threshold(dataImageGPU, data_thresholded, lowerThresholdLog);
				Ext.CLIJ2_getSumOfAllPixels(data_thresholded, nrPositiveCells);
				Ext.CLIJ2_pull(thresholded_labeledges);
				run(colors[channel-1]);
				selectImage("raw_ch"+channel);
				Overlay.clear;
				run("Add Image...", "image=["+thresholded_labeledges+"] x=0 y=0 opacity=100 zero");
			}
			selectImage(histogram);
			wait(10);
		}
//		Ext.CLIJ2_release(labeledges_log_intensities);
//		Ext.CLIJ2_release(labeledgesThin_log_intensities);
//		Ext.CLIJ2_release(labeledges_log_intensities);
//		Ext.CLIJ2_release(thresholded_labeledges);
		Ext.CLIJ2_release(data_thresholded);
		//Add overlay to the original image
		selectImage(image);
		Stack.setChannel(channel);
		run("Add Image...", "image=["+thresholded_labeledges+"] x=0 y=0 opacity=100 zero");
		Overlay.setPosition(channel, 0, 0);
		Overlay.show;
		close(thresholded_labeledges);
	}
	else if(thresholdChoice == "manual" && classifyChoice == "after running all images") {
		waitForUser("Press OK and then click in the histogram to set the threshold.");
		modifiers = 0;
		while(!isActive(histogramID) || modifiers != 16) {	//left mouse button
			getCursorLoc(x, y, z, modifiers);
			if((x > 0 && x < hwidth) && (y > 0 && y < hheight)) {	//If the cursor is on the histogram window
				Overlay.remove;
				lowerThresholdLog = Math.map(x, 20, 276, 0, 10);
	
				//Draw a red line at the threshold location
				setColor("red");
				linePosition = Math.map(lowerThresholdLog, minLogInt, maxLogInt, 0, 256);
				//Overlay.drawLine(23 + linePosition, 13, 23 + linePosition, 152);
				Overlay.drawLine(20 + linePosition, 11, 20 + linePosition, 138);
				Overlay.show();
			}
			wait(25);
		}
		lowerThresholdLog = Math.map(x, 20, 276, 0, 10);
	}
	close(histogram);

	lowerThresholdLin = pow(10,lowerThresholdLog);		//transform the threshold back to linear space
	
	for(k=0; k<data.length; k++) {
		if(data[k] > lowerThresholdLin) Table.set("Class ch"+channel, nRows-k-1, 1, resultsTable);	//positive
		else Table.set("Class ch"+channel, nRows-k-1, 0, resultsTable);								//negative
	}
	Table.update;
	print("10log of threshold (channel "+channel+") : "+lowerThresholdLog);
	print("Threshold converted to intensities: "+lowerThresholdLin);
	return nrPositiveCells;
}


//Convert a labelmap to ROIs using CLIJ2
function labels_to_ROI_Manager(labelmap) {
	selectImage(labelmap);
	if(getValue("Max") == 0) return;
	Ext.CLIJ2_push(labelmap);
	roiManager("reset");

	run("Clear Results");
	Ext.CLIJ2_statisticsOfLabelledPixels(labelmap, labelmap);
	boundingBox_X = Table.getColumn("BOUNDING_BOX_X", "Results");
	boundingBox_Y = Table.getColumn("BOUNDING_BOX_Y", "Results");
	boundingBox_width = Table.getColumn("BOUNDING_BOX_WIDTH", "Results");
	boundingBox_height = Table.getColumn("BOUNDING_BOX_HEIGHT", "Results");
	Array.getStatistics(boundingBox_width, min, boundingBoxMax_X, mean, stdDev);
	Array.getStatistics(boundingBox_height, min, boundingBoxMax_Y, mean, stdDev);
	Ext.CLIJ2_getMaximumOfAllPixels(labelmap, nr_cells);

	if(isOpen("ROI Manager")) { selectWindow("ROI Manager"); }//run("Close"); }	//This step goes faster when the ROI manager is not visible.
	for (i = 0; i < nr_cells; i++) {
		showStatus("Converting "+nr_cells+" labels to ROIs...");
		Ext.CLIJ2_crop2D(labelmap, label_cropped, boundingBox_X[i], boundingBox_Y[i], boundingBoxMax_X, boundingBoxMax_Y);
		Ext.CLIJ2_labelToMask(label_cropped, mask_label, i+1);
		Ext.CLIJ2_pullToROIManager(mask_label);
		roiManager("Select",i);
		Roi.move(boundingBox_X[i], boundingBox_Y[i]);
		roiManager("update");
	}
	run("Select None");
	Ext.CLIJ2_release(mask_label);
	Ext.CLIJ2_release(label_cropped);
	roiManager("Deselect");
	roiManager("Remove Frame Info");
	roiManager("Remove Channel Info");
	roiManager("Remove Slice Info");
	//	roiManager("Set Color", "#555555");
	roiManager("Set Color", "cyan");
	roiManager("Set Line Width", 2);
}


//Puts the values in a 1D 32-bit image and returns the Image ID
function arrayTo1DImage(array) {
	newImage("1D Image", "32-bit", array.length, 1, 1);
	imageID = getImageID();
	for(i=0; i<array.length; i++) {
		setPixel(i, 0, array[i]);
	}
	updateDisplay();
	return imageID;
}


//Prepends the value to the array
function prependToArray(value, array) {
	temparray=newArray(lengthOf(array)+1);
	for (i=0; i<lengthOf(array); i++) {
		temparray[i+1]=array[i];
	}
	temparray[0]=value;
	array=temparray;
	return array;
}


//Returns an array containing the 10log values of all elements
function logArray(array) {
	log_array=newArray(lengthOf(array));
	for (a=0; a<lengthOf(array); a++) {
		log_array[a]= Math.log10(array[a]);
	}
	return log_array;
}


//Returns the mean of the array
function meanOfArray(array) {
	Array.getStatistics(array, min, max, mean, stdDev);
	return mean;
}


//Adds a scalar to all elements of an array
function addScalarToArray(array, scalar) {
	added_array=newArray(lengthOf(array));
	for (a=0; a<lengthOf(array); a++) {
		added_array[a]=array[a] + scalar;
	}
	return added_array;
}


//Returns the first index at which a value occurs in an array
function firstIndexOfArray(array, value) {
	for (a=0; a<lengthOf(array); a++) {
		if (array[a]==value) {
			break;
		}
	}
	return a;
}


//Get the persistent value of the script parameter 'param' in class. N.B. This returns 'null' when the parameter is set to the default value!
function getPref(class, param) {
	return eval("js",
		"var ctx = Packages.ij.IJ.runPlugIn(\"org.scijava.Context\", \"\");" +
		"var ps = ctx.service(Packages.org.scijava.prefs.PrefService.class);" +
		"var " + param + " = ps.get(" + class + ".class, \"" + param + "\", \"<null>\");" +
		param + ";"
	);
}
