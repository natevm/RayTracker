/* Function to initialize all visualizations */
function InitializeAll() {
	PixelAnalysisView = new PixelAnalysisView();
	ImageView = new ImageView();
	ImageView.setParallelBarView(PixelAnalysisView.parallelBarView);
	PixelAnalysisView.parallelBarView.setImageView(ImageView);
	PixelAnalysisView.parallelBarView.setRayTreeView(PixelAnalysisView.rayTreeView);
}

function UpdateAll() {
	ImageView.update();
	PixelAnalysisView.update();
}

document.addEventListener("DOMContentLoaded", function(event) { 
	/* Initialize all visualizations */
	InitializeAll();
	UpdateAll();
	
	//load the pixel dataset
	d3.json("./Data/128/pixeldata.json", function(error, data){
		console.log(error);
		console.log(data);
		PixelAnalysisView.setData(data.PixelData);
		ImageView.setDimensions(data.PixelData.dimensions);
		UpdateAll();
	});
});

window.onresize = function(event) {
	PixelAnalysisView.resize();
    ImageView.resize();
};

d3.select("#story-toggle").on("click", function() {
	d3.selectAll(".toggled").classed("toggled", false);
	d3.select(this).classed("toggled", true);

	ImageView.toggleStoryMode();
	PixelAnalysisView.toggleStoryMode();
});

d3.select("#explore-toggle").on("click", function() {
	d3.selectAll(".toggled").classed("toggled", false);
	d3.select(this).classed("toggled", true);
	ImageView.toggleExploreMode();
	PixelAnalysisView.toggleExploreMode();
});