/*
	Nate Morrical. 
	An "ParallelBarView" displays a set of bars visualizing the current selection through a 
	variety of characteristics.
*/
class ParallelBarView {
	constructor() {
		d3.select("#bar-chart").html("");
		this.barchartDiv = d3.select("#bar-chart");
		this.svg = this.barchartDiv.append("svg");
		this.width = $("#bar-chart").width();
    	this.height = $("#bar-chart").height();
		this.svg.attr("width", this.width).attr("height", this.height);
		this.totalBars = 6; // Currently only 6 characteristics

		d3.select("#pixel-analysis-options").html("");
		this.optionsDiv = d3.select("#pixel-analysis-options");
		this.optionsSvg =  this.optionsDiv.append("svg");
		this.optionsWidth = $("#pixel-analysis-options").width();
    	this.optionsHeight = $("#pixel-analysis-options").height();
		this.optionsSvg.attr("width", this.optionsWidth)
						.attr("height", this.optionsHeight);

		this.createButtons();
		this.createGroups();

		this.colorAsc = false;
		this.timeAsc = false;
		this.branchesAsc = false;
		this.samplesAsc = false;
		this.depthAsc = false;
		this.varianceAsc = false;
	};

	/* Creates all buttons */
	createButtons() {
		this.optionsSvg.selectAll("g").remove();
		let group = this.optionsSvg.append("g");

		/* Color button */
		let colorButton = this.createButton(group, 0, "Color");
		colorButton.on("click", () => {this.sortColor()});

		/* Time button */
		let timeButton = this.createButton(group, 1, "Time");
		timeButton.on("click", () => {this.sort("time")});

		/* Branches button */
		let branchesButton = this.createButton(group, 2, "Branches");
		branchesButton.on("click", () => {this.sort("branches")});
		
		/* Samples button */
		let samplesButton = this.createButton(group, 3, "Samples");
		samplesButton.on("click", () => {this.sort("samples")});
		
		/* Depth button */
		let depthButton = this.createButton(group, 4, "Depth");
		depthButton.on("click", () => {this.sort("depth")});
		
		/* Variance button */
		let varianceButton = this.createButton(group, 5, "Variance");
		varianceButton.on("click", () => {this.sort("variance")});
	}

	/* Creates a specific button */
	createButton(svgGroup, index, label) {
		let padding = 3;

		let group = svgGroup.append("g");
		group.append("rect")
			.attr("y", padding)
			.attr("height", this.optionsHeight - 2 * padding)
			.attr("x", (index * (this.optionsWidth / this.totalBars) + padding))
			.attr("width", (this.optionsWidth / this.totalBars) - 2 * padding)
			.classed("button", true);

		/**/
		group.append("text").text(label)
			.attr("x", (index + .5) * (this.optionsWidth / this.totalBars))
			.attr("y", ( 1.5 * this.optionsHeight / 4))
			.style("text-anchor", "middle")
			.classed("label", true);

		/**/
		group.append("text").text("===")
			.attr("x", (index + .5) * (this.optionsWidth / this.totalBars))
			.attr("y", (3.5 * this.optionsHeight / 4))
			.style("text-anchor", "middle")
			.classed("label", true);

		group.append("rect")
			.attr("y", padding)
			.attr("height", this.optionsHeight - 2 * padding)
			.attr("x", (index * (this.optionsWidth / this.totalBars) + padding))
			.attr("width", (this.optionsWidth / this.totalBars) - 2 * padding)
			.classed("button-overlay", true);

		return group;
	}
	
	/* Creates a group for each characteristic bar */
	createGroups() {
		this.svg.append("g").attr("id", "colorbar");
		this.svg.append("g").attr("id", "timebar");
		this.svg.append("g").attr("id", "branchesbar");
		this.svg.append("g").attr("id", "samplesbar");
		this.svg.append("g").attr("id", "depthbar");
		this.svg.append("g").attr("id", "variancebar");
		this.svg.append("g").attr("id", "selectableBar");
	}

	/* Alters bounds to be within the image */
	sanitizeBounds(bounds, dimensions) {
		bounds[0] = Math.min(bounds[0], dimensions[0]);
		bounds[1] = Math.min(bounds[1], dimensions[1]);
		bounds[2] = Math.min(bounds[2], dimensions[0]);
		bounds[3] = Math.min(bounds[3], dimensions[1]);

		bounds[0] = Math.max(bounds[0], 0);
		bounds[1] = Math.max(bounds[1], 0);
		bounds[2] = Math.max(bounds[2], 0);
		bounds[3] = Math.max(bounds[3], 0);
	}

	/* Updates the visualization to show data for pixels within the provided bounds */
	updateBounds(bounds) {
		if (this.rawData != null) {
			this.sanitizeBounds(bounds, this.rawData.dimensions);
			
			let self = this;
			let svg = this.svg;

			/* Wrangle data from bounds */
			let data=[];
			for (let y = bounds[1]; y < bounds[3]; ++y) {
				for (let x = bounds[0]; x < bounds[2]; ++x) {
					let entry = {};

					let color = {
						"red" : -1,
						"blue" : -1,
						"green" : -1
					};
					color.red = this.rawData.Color[3 * (x + y * this.rawData.dimensions[0])];
					color.blue = this.rawData.Color[3 * (x + y * this.rawData.dimensions[0]) + 1];
					color.green = this.rawData.Color[3 * (x + y * this.rawData.dimensions[0]) + 2];
					
					entry.color = color;

					entry.time = this.rawData.Render_time[x + y * this.rawData.dimensions[0]];
					entry.branches = this.rawData.Secondary_rays[x + y * this.rawData.dimensions[0]];
					entry.samples = this.rawData.Sample_count[x + y * this.rawData.dimensions[0]];
					entry.depth = this.rawData.Depth_buffer[x + y * this.rawData.dimensions[0]];
					entry.variance = this.rawData.Variance[x + y * this.rawData.dimensions[0]];
					entry.x = x;
					entry.y = y;
					data.push(entry);

				}
			}
			this.data = data;
			this.updateBars(this.data);
		}
	}

	/* Updates all bars given an array of pixel data */
	updateBars(data) {
		/* Pixel color */
		let colorBar = this.updateBar(data, "color", d3.select("#colorbar"), 6, 0);
		colorBar.style("fill", (d, i) => {return d3.rgb(d.color.red,d.color.blue,d.color.green);});

		/* Pixel render time */
		let timeBar = this.updateBar(data, "time", d3.select("#timebar"), 6, 1);
		var timeColor = d3.scaleLinear()
			.domain([this.rawData.min_render_time, this.rawData.max_render_time])
            .range(['white','red']);
		timeBar.style("fill", (d) => {return timeColor(d.time);});

		/* Total secondary rays */
		let secondaryRayBar = this.updateBar(data, "branches", d3.select("#branchesbar"), 6, 2);
		var srColor = d3.scaleLinear()
			.domain([this.rawData.min_secondary_rays, this.rawData.max_secondary_rays])
            .range(['white','red']);
		secondaryRayBar.style("fill", (d) => {return srColor(d.branches);});

		/* Total samples */
		let samplesBar = this.updateBar(data, "samples", d3.select("#samplesbar"), 6, 3);
		var scColor = d3.scaleLinear()
			.domain([this.rawData.min_sample_count, this.rawData.max_sample_count])
            .range(['white','red']);
		samplesBar.style("fill", (d) => {return scColor(d.samples);});

		/* Pixel depth */
		let depthBar = this.updateBar(data, "depth", d3.select("#depthbar"), 6, 4);
		var dColor = d3.scaleLinear()
			.domain([this.rawData.min_depth, this.rawData.max_depth])
            .range(['white','red']);
		depthBar.style("fill", (d) => {return dColor(d.depth);});

		/* Pixel variance */
		let varianceBar = this.updateBar(data, "variances", d3.select("#variancebar"), 6, 5);
		var vColor = d3.scaleLinear()
			.domain([this.rawData.min_variance, this.rawData.max_variance])
            .range(['white','red']);
		varianceBar.style("fill", (d) => {return vColor(d.variance);});

		/* Selectable bars */
		let selectableBar = this.updateSelectable(data, d3.select("#selectableBar"))
	}
	
	/* Updates a specific bar */
	updateBar(data, field, svgGroup, totalBars, index) {
		let rects = svgGroup.selectAll("rect").data(data);
		rects.exit().remove();
		let enterRects = rects.enter().append("rect");

		let barWidth = (this.width / totalBars);
		let barHeight = (1 / data.length) * this.height;
		let dataLenInv = 1 / data.length; 

		enterRects
			.attr("x", (d, i) => {return index * barWidth;})
			.attr("width", barWidth)
			.style("stroke-width", 0);

		let allRects = enterRects.merge(rects);
		allRects.attr("y", (d, i) => {return (i * dataLenInv) * this.height})
			.attr("height", barHeight);

		return allRects;
	}

	/* Similar to updateBar, this updates a selectable row overlay. When rows are selected, 
			the treeview and image view are also updated. */
	updateSelectable(data, svgGroup) {
		let self = this;
        
        /* Temporarily deselect the selected row, if selected. */
        d3.select(".row-overlay.active").classed("active", false);

		let rects = svgGroup.selectAll("rect").data(data);
		rects.exit().remove();
		let enterRects = rects.enter().append("rect");

		let barWidth = this.width;
		let barHeight = (1 / data.length) * this.height;
		let dataLenInv = 1 / data.length; 

		enterRects
			.attr("x", 0)
			.attr("width", barWidth)
			.style("stroke-width", 0);

		let allRects = enterRects.merge(rects);
		allRects.attr("y", (d, i) => {return (i * dataLenInv) * this.height})
			.attr("height", barHeight)
			.classed("row-overlay", true)
			.classed("active", (d,i) => { return (d.x == self.selectedx && d.y == self.selectedy);})
			.on("mouseover", function(d,i) {
				d3.select(this).classed("hover", true);
			})
          	.on("mouseout", function(d,i) {
          		d3.select(this).classed("hover", false);
          	})
          	.on("click", function(d,i) {
          		self.selectedx = d.x;
          		self.selectedy = d.y;
          		d3.select(".row-overlay.active").classed("active", false);
          		d3.select(this).classed("active", true);
          		if (self.imageview != null)
          			self.imageview.selectPixel(d.x, d.y);
          		/* TODO: update tree view */
          	})

		return allRects;
	}

	/* Resizes the visualization. Good for device rotation/browser scaling */
	resize() {
		/* Update width and height just incase canvas size has changed */
		this.width = $("#bar-chart svg").parent().width();
    	this.height = $("#bar-chart svg").parent().height();
    	this.optionsWidth = $("#pixel-analysis-options").width();
    	this.optionsHeight = $("#pixel-analysis-options").height();
		this.svg.attr("width", this.width).attr("height", this.height);
		this.optionsSvg.attr("width", this.optionsWidth)
						.attr("height", this.optionsHeight);
		this.createButtons();
	}

	setData(_data){
		this.rawData = _data;
	}

	rgbToHsv(r, g, b){
	    r = r/255, g = g/255, b = b/255;
	    var max = Math.max(r, g, b), min = Math.min(r, g, b);
	    var h, s, v = max;

	    var d = max - min;
	    s = max == 0 ? 0 : d / max;

	    if(max == min){
	        h = 0; // achromatic
	    }else{
	        switch(max){
	            case r: h = (g - b) / d + (g < b ? 6 : 0); break;
	            case g: h = (b - r) / d + 2; break;
	            case b: h = (r - g) / d + 4; break;
	        }
	        h /= 6;
	    }

	    return [h, s, v];
	}

	/* Sort by either hue, saturation, or value. */
	sortColor() {
		this.colorAsc = !this.colorAsc;
		this.data.sort(
            (x, y) => {
	        	return this.rgbToHsv(x.color.red, x.color.green, x.color.blue)[0] 
            		 - this.rgbToHsv(y.color.red, y.color.green, y.color.blue)[0];
          }
        );

        if (!this.colorAsc) this.data.reverse();

        this.updateBars(this.data);
	}

	/* Sort the provided numberic field */
	sort(field) {
		this[field + "Asc"] = !this[field + "Asc"];

		this.data.sort(
            (x, y) => {
            return x[field] - y[field];
          }
        );

        if (!this[field + "Asc"]) this.data.reverse();

        this.updateBars(this.data);
	}

	/* This could probably be designed a little cleaner... */
	setImageView(imageview) {
		this.imageview = imageview;
	}
}