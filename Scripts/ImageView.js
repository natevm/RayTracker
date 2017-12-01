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
		this.maxCheckpoint = 7;

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

		this.hoveredLocation = [-1, -1];
		this.selectedLocation = [-1, -1];
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

		let textData = this.getText(this.checkpoint);

		textGroup.append("text")
		.text(textData.header)
		.classed("h", true).classed("unselectable", true)
		.attr("text-anchor", "middle")
		.attr("x", scaleWidth(50)).attr("y", scaleHeight(25));


		textGroup.append("text")
		.text(textData.l1)
		.classed("p", true).classed("unselectable", true)
		.attr("text-anchor", "middle")
		.attr("x", scaleWidth(50)).attr("y", scaleHeight(50));
		textGroup.append("text")
		.text(textData.l2)
		.classed("p", true).classed("unselectable", true)
		.attr("text-anchor", "middle")
		.attr("x", scaleWidth(50)).attr("y", scaleHeight(65));
		textGroup.append("text")
		.text(textData.l3)
		.classed("p", true).classed("unselectable", true)
		.attr("text-anchor", "middle")
		.attr("x", scaleWidth(50)).attr("y", scaleHeight(80));


		/* Add next/previous buttons */
		svg.append("rect")
		.classed("button", true)
		.attr("id", "imageMenuBackButton")
		.attr("x", scaleWidth(2)).attr("y", scaleHeight(30+buttonVertOffset))
		.attr("width", scaleWidth(buttonWidth)).attr("height", scaleHeight(40))
		.on("click", function() { if (self.checkpoint > 0) self.checkpoint--; self.update();});

		svg.append("line")
		.classed("unselectable", true)
		.attr("x1", scaleWidth(2 + 3)).attr("y1", scaleHeight(35+buttonVertOffset))
		.attr("x2", scaleWidth(2 + 2)).attr("y2", scaleHeight(50+buttonVertOffset))
		.attr("stroke-width", 2).attr("stroke", "white");

		svg.append("line")
		.classed("unselectable", true)
		.attr("x1", scaleWidth(2+2)).attr("y1", scaleHeight(50+buttonVertOffset))
		.attr("x2", scaleWidth(2+3)).attr("y2", scaleHeight(65+buttonVertOffset))
		.attr("stroke-width", 2).attr("stroke", "white");

		svg.append("rect")
		.attr("id", "imageMenuNextButton")
		.classed("button", true)
		.attr("fill", "red")
		.attr("x", scaleWidth(98-buttonWidth)).attr("y", scaleHeight(30+buttonVertOffset))
		.attr("width", scaleWidth(buttonWidth)).attr("height", scaleHeight(40))
		.on("click", function() { if (self.checkpoint < self.maxCheckpoint) self.checkpoint++; self.update();});

		svg.append("line")
		.classed("unselectable", true)
		.attr("x1", scaleWidth(98-buttonWidth + 2)).attr("y1", scaleHeight(35+buttonVertOffset))
		.attr("x2", scaleWidth(98-buttonWidth + 3)).attr("y2", scaleHeight(50+buttonVertOffset))
		.attr("stroke-width", 2).attr("stroke", "white");

		svg.append("line")
		.classed("unselectable", true)
		.attr("x1", scaleWidth(98-buttonWidth + 3)).attr("y1", scaleHeight(50+buttonVertOffset))
		.attr("x2", scaleWidth(98-buttonWidth + 2)).attr("y2", scaleHeight(65+buttonVertOffset))
		.attr("stroke-width", 2).attr("stroke", "white");
	}

	getText(checkpoint) {
		switch(checkpoint) {
		    case 0:
		        return {"header": "Ray Tracing",
		        "l1": "In computer graphics, ray tracing is a rendering technique which",
		    	"l2": "generates images with a very high degree of visual realism. Unfortunately,",
		    	"l3": " an image can take a long time to render, for a wide variety of reasons."
		    };
		    case 1:
		        return {"header": "Per Pixel Render Time",
		        "l1": "Render time varies on a per pixel basis. This  ",
		    	"l2": "variation can be difficult to parallelize, since per ",
		    	"l3": "pixel parallelization can lead to load imbalancing."
		    };
		    case 2:
		        return {"header": "Color Variance",
		        "l1": "Color variance occurs when multiple objects, or detailed textures all",
		    	"l2": "lie within the boundries of a pixel. This can cause aliasing, since ",
		    	"l3": "one ray might lead to a poor estimation of a pixel's true color."
		    };
		    case 3:
		        return {"header": "Multisampling, and Antialiasing",
		        "l1": "To smooth the pixel out, multiple rays are cast for",
		    	"l2": "each pixel based on color variance. The more samples that",
		    	"l3": "are cast, the higher the computation time for that pixel"
		    };
		    case 4:
		        return {"header": "Reflections and Refractions",
		        "l1": "To represent reflections and refractions, rays must reflect ",
		    	"l2": "and refract off surfaces. This branching can also impact ",
		    	"l3": "render time, since more computation is required."
		    };
		    case 5:
		        return {"header": "Bounding Volume Hierarchies",
		        "l1": "When a ray is cast, it has to test for intersections with a",
		    	"l2": "large number of triangles. To improve performance, bounding volume",
		    	"l3": "hierarchies are used to avoid triangle intersections."
		    };
		    case 6:
		        return {"header": "Ray Triangle Intersections",
		        "l1": "Eventually, some rays must test for intersection with a couple triangles. ",
		    	"l2": "Triangle intersections can also be quite expensive.",
		    	"l3": "Rays intersect more triangles on the edge of a mesh than in the center."
		    };
		    case 7:
		        return {"header": "What's next...",
		        "l1": "Research is being conducted here at the U to improve ray tracing ",
		    	"l2": "times as well as power consumption. Hopefully, in a couple  ",
		    	"l3": "years we'll be able to achieve real time interactive ray tracing."
		    };
		    default:
		        return "";
		}
	}

	getImageLocation(checkpoint) {
		switch(checkpoint) {
			case 0: return "./Data/" + this.dimensions[0] + "/render.png";
			case 1: return "./Data/" + this.dimensions[0] + "/renderTime.png";
			case 2: return "./Data/" + this.dimensions[0] + "/variance.png";
			case 3: return "./Data/" + this.dimensions[0] + "/sampleCount.png";
			case 4: return "./Data/" + this.dimensions[0] + "/secondaryRays.png";
			case 5: return "./Data/" + this.dimensions[0] + "/boxIntersections.png";
			case 6: return "./Data/" + this.dimensions[0] + "/objIntersections.png";
			case 7: return "./Data/" + this.dimensions[0] + "/highRes.png";
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
		this.image.attr("height", this.height).attr("width", this.height);

		/* Small square to show selected pixel*/
		let xScale = d3.scaleLinear().domain([0, this.dimensions[0]]).range([0, this.height]);
		let yScale = d3.scaleLinear().domain([0, this.dimensions[0]]).range([0, this.height]);

		group.append("rect")
			.attr("class", "selectedPixel")
			.attr("width", this.width / this.dimensions[0])
			.attr("height", this.width / this.dimensions[1])
			.attr("x", xScale(this.selectedLocation[0]))
			.attr("y", yScale(this.selectedLocation[1]));

		group.append("rect")
			.attr("class", "hoveredPixel")
			.attr("width", this.width / this.dimensions[0])
			.attr("height", this.width / this.dimensions[1])
			.attr("x", xScale(this.hoveredLocation[0]))
			.attr("y", yScale(this.hoveredLocation[1]));
		
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
		
		/* Small square to show selected pixel*/
		let xScale = d3.scaleLinear().domain([0, this.dimensions[0]]).range([0, this.height]);
		let yScale = d3.scaleLinear().domain([0, this.dimensions[0]]).range([0, this.height]);
		this.selectedLocation[0] = x;
		this.selectedLocation[1] = y;

		d3.select(".selectedPixel")
			.attr("width", this.width / this.dimensions[0])
			.attr("height", this.width / this.dimensions[1])
			.attr("x", xScale(x))
			.attr("y", yScale(y));

	}

	hoverPixel(x, y) {
		/* Small square to show selected pixel*/
		let xScale = d3.scaleLinear().domain([0, this.dimensions[0]]).range([0, this.height]);
		let yScale = d3.scaleLinear().domain([0, this.dimensions[0]]).range([0, this.height]);
		this.hoveredLocation[0] = x;
		this.hoveredLocation[1] = y;

		d3.select(".hoveredPixel")
			.attr("width", this.width / this.dimensions[0])
			.attr("height", this.width / this.dimensions[1])
			.attr("x", xScale(x))
			.attr("y", yScale(y));

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