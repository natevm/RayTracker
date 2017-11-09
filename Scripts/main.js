/* Function to initialize all visualizations */
function InitializeAll() {
	ImageView = new ImageView();
}

function UpdateAll() {
	ImageView.update();
}

document.addEventListener("DOMContentLoaded", function(event) { 
	/* Initialize all visualizations */
	InitializeAll();
	UpdateAll();
});

window.onresize = function(event) {
    ImageView.resize();
};