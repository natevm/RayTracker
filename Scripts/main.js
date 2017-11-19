/* Function to initialize all visualizations */
function InitializeAll() {
	ParallelBarView = new ParallelBarView();
	ImageView = new ImageView(ParallelBarView);
	RayTreeView = new RayTreeView();
}

function UpdateAll() {
	ImageView.update();
	RayTreeView.update();
}

document.addEventListener("DOMContentLoaded", function(event) { 
	/* Initialize all visualizations */
	InitializeAll();
	UpdateAll();
	
	//load the pixel dataset
	d3.json("./Data/64/pixeldata.json", function(error, data){
		console.log(error);
		console.log(data);
		ParallelBarView.setData(data.PixelData);

		UpdateAll();
	});

	// d3.json("./Data/64/raydata/.json", function(error, data){
	// 	console.log(error);
	// 	console.log(data);
	// 	RayTreeView.setData(data.RayData);
	// 	RayTreeView.update();
	// });

});

window.onresize = function(event) {
	ParallelBarView.resize();
    ImageView.resize();
	RayTreeView.resize();
};