/* Function to initialize all visualizations */
function InitializeAll() {
	ImageView = new ImageView();
	RayTreeView = new RayTreeView();
	ParallelBarView = new ParallelBarView();
}

function UpdateAll() {
	ImageView.update();
	RayTreeView.update();
	ParallelBarView.update();
}

document.addEventListener("DOMContentLoaded", function(event) { 
	/* Initialize all visualizations */
	InitializeAll();
	UpdateAll();
	
	//load the initial dataset
	d3.json("./Data/Ray_tracing_test_data.json", function(error, data){
		console.log(error);
		console.log(data);
		//pass it to RayTreeView
		RayTreeView.setData(data.RayData);
		ParallelBarView.setData(data.RayData);

		//call update on RayTreeView
		RayTreeView.update();
		ParallelBarView.update();
	});
});

window.onresize = function(event) {
    ImageView.resize();
	RayTreeView.resize();
	ParallelBarView.resize();
};