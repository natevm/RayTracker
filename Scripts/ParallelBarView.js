/*
	Nate Morrical. 
	An "ParallelBarView" displays a set of bars visualizing the current selection through a 
	variety of characteristics.
*/
class ParallelBarView {
	constructor() {
		this.barchartDiv = d3.select("#bar-chart");
		this.svg = this.barchartDiv.append("svg");
		this.width = $("#bar-chart").width();
    	this.height = $("#bar-chart").height();
		this.svg.attr("width", this.width).attr("height", this.height);
	};

	update() {
		let self = this;
		let svg = this.svg;

		svg.selectAll("g").remove();
		svg.selectAll("rect").remove();
		
		/* concept art */
		svg.selectAll("image").remove();
		this.image = svg.append("image");
		this.image.attr("xlink:href", "./Sketches/Vis4DataSci1bars.png");
		this.image.attr("height", this.height).attr("width", this.width);

		/* TEMPORARY RECTANGLES */
		// <rect x="10" y="10" width="100" height="100"/>
		let totalRects = 7;
		let totalPixels = 16;

		for (let i = 0; i < totalRects; ++i) {
			for (let j = 0; j < totalPixels; ++j) {
				svg.append("rect")
					.attr("x", (this.width / totalRects) * i)
					.attr("y", (this.height / totalPixels) * j)
					.attr("width", this.width / (totalRects * 2))
					.attr("height", this.height / (totalPixels))
					.attr("fill", "rgba(0, 0, 0, 0.3)");

				//<line x1="20" y1="100" x2="100" y2="20" stroke-width="2" stroke="black"/>
				svg.append("line")
					.attr("x1", (this.width / totalRects) * i)
					.attr("y1", (this.height / totalPixels) * j)
					.attr("x2", ((this.width / totalRects) * i) +((this.width / totalRects) * i + 1)  )
					.attr("y2", (this.height / totalPixels) * j)
					.attr("stroke-width", 2)
					.attr("stroke", "rgba(0, 0, 0, 0.3)");
			}
		}


		// <img src="./Sketches/Vis4DataSci1bars.png"/>
		// let width = this.width;
		// let height = this.height;

		// let scale0 = 1;
		// let translate0 = [0,0];

		// svg.append("rect")
		// 		.attr("class", "overlay") // class not used yet
		// 		.attr("width", width + "px")
		// 		.attr("height", height + "px");

		// let group = svg.append("g")
		// 		.attr("transform", "translate(" + translate0 +")scale(" + scale0 + ")");

		// svg.call(d3.zoom().scaleExtent([1,10]).on("zoom", () => {
		// 		group.attr("transform", d3.event.transform)
		// 		self.brush(d3.event.transform);
		// 	}
		// ));

		// this.image = group.append("image");
		// this.image.attr("xlink:href", "./Data/Images/Image.png");
		// this.image.attr("height", height).attr("width", width);
		
		// svg.append("rect")
		// 	.attr("class", "brush")
		// 	.attr("width", this.width / 5)
		// 	.attr("height", this.height / 5)
		// 	.attr("x", this.width * 2 / 5)
		// 	.attr("y", this.height * 2/ 5);
	}

	resize() {
		/* Update width and height just incase canvas size has changed */
		this.width = $("#bar-chart svg").parent().width();
    	this.height = $("#bar-chart svg").parent().height();
		this.svg.attr("width", this.width).attr("height", this.height);

		this.update();
	}

	setData(_data){
		this.rawData = _data;
		/* TODO: read this data */
	}
}