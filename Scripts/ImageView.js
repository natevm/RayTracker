/*
	Nate Morrical. 
	An "ImageView" displays an image, and allows panning, zooming, 
	and brushing pixels for analysis.
*/
class ImageView {
	constructor() {
		this.viewportDiv = d3.select("#viewport");
		this.svg = this.viewportDiv.append("svg");
		this.width = $("#viewport").width();
    	this.height = $("#viewport").height();
    	this.imgWidth = 256; 
    	this.imgHeight = 256; 
		this.svg.attr("width", this.width).attr("height", this.height);
	};

	update() {
		let self = this;
		let svg = this.svg;

		svg.selectAll("g").remove();
		svg.selectAll("rect").remove();

		let width = this.width;
		let height = this.height;

		let scale0 = 1;
		let translate0 = [0,0];

		svg.append("rect")
				.attr("class", "overlay") // class not used yet
				.attr("width", width + "px")
				.attr("height", height + "px");

		let group = svg.append("g")
				.attr("transform", "translate(" + translate0 +")scale(" + scale0 + ")");

		svg.call(d3.zoom().scaleExtent([1,10]).on("zoom", () => {
				group.attr("transform", d3.event.transform)
				self.brush(d3.event.transform);
			}
		));

		this.image = group.append("image");
		this.image.attr("xlink:href", "./Data/Images/Image.png");
		this.image.attr("height", height).attr("width", width);
		
		svg.append("rect")
			.attr("class", "brush")
			.attr("width", this.width / 5)
			.attr("height", this.height / 5)
			.attr("x", this.width * 2 / 5)
			.attr("y", this.height * 2/ 5);
	}

	brush(transform) {
		/* Don't ask me to explain this. ( sorry future self. ;( ) */
		let x0 = (((((2/5) * $("#viewport").width())  - transform.x) / transform.k) / $("#viewport").width())  * 256;
		let y0 = (((((2/5) * $("#viewport").height()) - transform.y) / transform.k) / $("#viewport").height()) * 256;

		let x1 = (((((3/5) * $("#viewport").width())  - transform.x) / transform.k) / $("#viewport").width())  * 256;
		let y1 = (((((3/5) * $("#viewport").height()) - transform.y) / transform.k) / $("#viewport").height()) * 256;

		x0 = Math.floor(x0);
		y0 = Math.floor(y0);
		x1 = Math.floor(x1);
		y1 = Math.floor(y1);

		// console.log(" x0 " + x0 + " y0 " + y0 + " x1 " + x1 + " y1 " + y1);
		d3.select("#viewport-options").text(" X0 " + x0 + " Y0 " + y0 + " X1 " + x1 + " Y1 " + y1);
	}

	resize() {
		/* Update width and height just incase canvas size has changed */
		this.width = $("#viewport svg").parent().width();
    	this.height = $("#viewport svg").parent().height();
		this.svg.attr("width", this.width).attr("height", this.height);

		this.update();
	}
}