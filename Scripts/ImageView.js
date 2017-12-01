/*
	Nate Morrical. 
	An "ImageView" displays an image, and allows panning, zooming, 
	and brushing pixels for analysis.
*/
class ImageView {
	constructor(parallelBarView) {
		// Create the menu
		d3.select("#image-menu").html("");
		this.menuDiv = d3.select("#image-menu");
		this.menuSVG = this.menuDiv.append("svg");
		this.menuWidth = $("#image-menu").width();
    	this.menuHeight = $("#image-menu").height();
    	this.menuSVG.attr("width", this.menuWidth).attr("height", this.menuHeight);
		this.checkpoint = 0;
		this.maxCheckpoint = 5;

		// Create the viewport SVG
		d3.select("#image-viewport").html("");
		this.parallelBarView = parallelBarView;
		this.parallelBarView.setImageView(this);
		this.viewportDiv = d3.select("#image-viewport");
		this.svg = this.viewportDiv.append("svg");
		this.width = $("#image-viewport").width();
    	this.height = $("#image-viewport").height();
    	this.imgWidth = 256; 
    	this.imgHeight = 256; 
		this.svg.attr("width", this.width).attr("height", this.height);
		this.lastTransform = {"k":1,"x":0,"y":0};
	};

	setDimensions(dimensions) {
		this.dimensions = dimensions;
	}

	update() 
	{
		this.updateMenu();
		this.updateViewport();
	}

	updateMenu() 
	{
		let buttonVertOffset = 10;
		let buttonWidth = 5;


		let self = this;
		let svg = this.menuSVG;
		let scaleHeight = d3.scaleLinear().domain([0,100]).range([0, this.menuHeight]);
		let scaleWidth = d3.scaleLinear().domain([0,100]).range([0, this.menuWidth]);

		/* Clear the menu */
		svg.selectAll("*").remove();

		/* Add text for the current checkpoint */
		let textGroup = svg.append("g").classed("textGroup", true);

		textGroup.append("text")
		.text(this.getCheckpointHeader(this.checkpoint))
		.classed("h", true).classed("unselectable", true)
		.attr("text-anchor", "middle")
		.attr("x", scaleWidth(50)).attr("y", scaleHeight(25));


		textGroup.append("text")
		.text(this.getCheckpointLine1(this.checkpoint))
		.classed("p", true).classed("unselectable", true)
		.attr("text-anchor", "middle")
		.attr("x", scaleWidth(50)).attr("y", scaleHeight(50));
		textGroup.append("text")
		.text(this.getCheckpointLine2(this.checkpoint))
		.classed("p", true).classed("unselectable", true)
		.attr("text-anchor", "middle")
		.attr("x", scaleWidth(50)).attr("y", scaleHeight(65));
		textGroup.append("text")
		.text(this.getCheckpointLine3(this.checkpoint))
		.classed("p", true).classed("unselectable", true)
		.attr("text-anchor", "middle")
		.attr("x", scaleWidth(50)).attr("y", scaleHeight(80));


		/* Add next/previous buttons */
		svg.append("rect")
		.classed("button", true)
		.attr("id", "imageMenuBackButton")
		.attr("x", scaleWidth(5)).attr("y", scaleHeight(30+buttonVertOffset))
		.attr("width", scaleWidth(buttonWidth)).attr("height", scaleHeight(40))
		.on("click", function() { if (self.checkpoint > 0) self.checkpoint--; self.update();});

		svg.append("line")
		.classed("unselectable", true)
		.attr("x1", scaleWidth(5 + 3)).attr("y1", scaleHeight(35+buttonVertOffset))
		.attr("x2", scaleWidth(5 + 2)).attr("y2", scaleHeight(50+buttonVertOffset))
		.attr("stroke-width", 2).attr("stroke", "white");

		svg.append("line")
		.classed("unselectable", true)
		.attr("x1", scaleWidth(5+2)).attr("y1", scaleHeight(50+buttonVertOffset))
		.attr("x2", scaleWidth(5+3)).attr("y2", scaleHeight(65+buttonVertOffset))
		.attr("stroke-width", 2).attr("stroke", "white");

		svg.append("rect")
		.attr("id", "imageMenuNextButton")
		.classed("button", true)
		.attr("fill", "red")
		.attr("x", scaleWidth(95-buttonWidth)).attr("y", scaleHeight(30+buttonVertOffset))
		.attr("width", scaleWidth(buttonWidth)).attr("height", scaleHeight(40))
		.on("click", function() { if (self.checkpoint < self.maxCheckpoint) self.checkpoint++; self.update();});

		svg.append("line")
		.classed("unselectable", true)
		.attr("x1", scaleWidth(95-buttonWidth + 2)).attr("y1", scaleHeight(35+buttonVertOffset))
		.attr("x2", scaleWidth(95-buttonWidth + 3)).attr("y2", scaleHeight(50+buttonVertOffset))
		.attr("stroke-width", 2).attr("stroke", "white");

		svg.append("line")
		.classed("unselectable", true)
		.attr("x1", scaleWidth(95-buttonWidth + 3)).attr("y1", scaleHeight(50+buttonVertOffset))
		.attr("x2", scaleWidth(95-buttonWidth + 2)).attr("y2", scaleHeight(65+buttonVertOffset))
		.attr("stroke-width", 2).attr("stroke", "white");
	}

	getCheckpointHeader(checkpoint) {
		switch(checkpoint) {
		    case 0:
		        return "Ray Tracing";
		    case 1:
		        return "Renders Take a Long Time";
		    case 2:
		        return "Color Variance";
		    case 3:
		        return "Ray Bouncing and Branching";
		    case 4:
		        return "Ray Box Intersections";
		    case 5:
		        return "Ray Triangle Intersections";
		    default:
		        return "";
		}
	}

	getCheckpointLine1(checkpoint) {
		switch(checkpoint) {
		    case 0:
		        return "In computer graphics, ray tracing is a rendering technique generating";
		    case 1:
		        return "The technique is capable of producing a very ";
		    case 2:
		        return "Color variance greatly effects computational cost.";
		    case 3:
		        return "To represent reflections and refractions, rays must ";
		    case 4:
		        return "As each ray is cast, it must test for intersection with boxes ";
		    case 5:
		        return "Some boxes in that BVH contain triangles. ";
		    default:
		        return "";
		}
	}

	getCheckpointLine2(checkpoint) {
		switch(checkpoint) {
		    case 0:
		        return "an image by tracing the path of light as pixels in an image plane and";
		    case 1:
		        return "high degree of visual realism, but at a greater cost.";
		    case 2:
		        return "To get a smooth edges, when variance is high, more ";
		    case 3:
		        return "reflect and refract off surfaces. This branching ";
		    case 4:
		        return "in a bounding volume hierarchy. This greatly improves ";
		    case 5:
		        return "When a ray hits a box containing triangles, it must test for ";
		    default:
		        return "";
		}
	}
	
	getCheckpointLine3(checkpoint) {
		switch(checkpoint) {
		    case 0:
		        return "simulating the effects of its encounters with virtual objects.";
		    case 1:
		        return "(here, value encodes time to render)";
		    case 2:
		        return "rays must be cast from a particular pixel. ";
		    case 3:
		        return "can also impact render time. ";
		    case 4:
		        return "ray mesh intersection times, but isn't free.";
		    case 5:
		        return "intersections with those triangles. ";
		    default:
		        return "";
		}
	}

	getImageLocation(checkpoint) {
		switch(checkpoint) {
			case 0: return "./Data/" + this.dimensions[0] + "/render.png";
			case 1: return "./Data/" + this.dimensions[0] + "/renderTime.png";
			case 2: return "./Data/" + this.dimensions[0] + "/variance.png";
			case 3: return "./Data/" + this.dimensions[0] + "/secondaryRays.png";
			case 4: return "./Data/" + this.dimensions[0] + "/boxIntersections.png";
			case 5: return "./Data/" + this.dimensions[0] + "/objIntersections.png";
			default: return "./Data/" + this.dimensions[0] + "/render.png"
		}
	}

	updateViewport() 
	{
		if (!this.dimensions) return;
		let self = this;
		let svg = this.svg;

		/* Clear the viewport */
		svg.selectAll("g").remove();
		svg.selectAll("rect").remove();

		/* Create an overlay to prevent image drag event */
		svg.append("rect")
				.attr("class", "overlay") // class not used yet
				.attr("width", this.width + "px")
				.attr("height", this.height + "px");


		/* Create a group with an initial transform for the image being shown */
		let scale = 1;
		let translate = [0,0];
		let group = svg.append("g")
				.attr("transform", "translate(" + translate +")scale(" + scale + ")");

		svg.call(d3.zoom().scaleExtent([1.0,10]).on("zoom", () => {
				group.attr("transform", d3.event.transform)
				self.brush(d3.event.transform);
				this.lastTransform = d3.event.transform;
			}
		));

		/* Call brush once with the default transform to update other views. */
		group.attr("transform", "translate(" + this.lastTransform.x + "," + this.lastTransform.y + ") scale(" + this.lastTransform.k+")");
		this.brush(this.lastTransform);

		/* Add the image to the transformable group */
		this.image = group.append("image");
		this.image.attr("xlink:href", this.getImageLocation(this.checkpoint));
		this.image.attr("height", this.height).attr("width", this.width);
		
		/* Add a small square to show the brush region. */
		svg.append("rect")
			.attr("class", "brush")
			.attr("width", this.width / 5)
			.attr("height", this.height / 5)
			.attr("x", this.width * 2 / 5)
			.attr("y", this.height * 2/ 5);
	}

	brush(transform) {
		/* Don't ask me to explain this. ( sorry future self. ;( ) */
		let x0 = (((((2/5) * $("#image-viewport").width())  - transform.x) / transform.k) / $("#image-viewport").width())  * this.dimensions[0];
		let y0 = (((((2/5) * $("#image-viewport").height()) - transform.y) / transform.k) / $("#image-viewport").height()) * this.dimensions[1];

		let x1 = (((((3/5) * $("#image-viewport").width())  - transform.x) / transform.k) / $("#image-viewport").width())  * this.dimensions[0];
		let y1 = (((((3/5) * $("#image-viewport").height()) - transform.y) / transform.k) / $("#image-viewport").height()) * this.dimensions[1];

		x0 = Math.ceil(x0);
		y0 = Math.ceil(y0);
		x1 = Math.floor(x1);
		y1 = Math.floor(y1);

		console.log(" X0 " + x0 + " Y0 " + y0 + " X1 " + x1 + " Y1 " + y1);

		let bounds = [x0, y0, x1, y1];
		this.parallelBarView.updateBounds(bounds);
	}

	selectPixel(x, y) {
		console.log("TODO: select pixel " + x + " " + y + " in the image view. ");
	}

	resize() {
		/* Update width and height just incase canvas size has changed */
		this.width = $("#image-viewport svg").parent().width();
    	this.height = $("#image-viewport svg").parent().height();
		this.svg.attr("width", this.width).attr("height", this.height);

		this.menuWidth = $("#image-menu svg").parent().width();
    	this.menuHeight = $("#image-menu svg").parent().height();
    	this.menuSVG.attr("width", this.menuWidth).attr("height", this.menuHeight);

		this.update();
	}
}