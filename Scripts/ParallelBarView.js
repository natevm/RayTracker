/*
	Nate Morrical. 
	An "ParallelBarView" displays a set of bars visualizing the current selection through a 
	variety of characteristics.
*/
class ParallelBarView {
	constructor() {
		/* Create viewport */		
		this.barchartDiv = d3.select("#pixel-analysis-viewport");
		this.svg = this.barchartDiv.select("svg");
		this.width = $("#pixel-analysis-viewport").width();
    	this.height = $("#pixel-analysis-viewport").height();
    	this.barGroupHeight = this.height * .85;
    	this.barGroupWidth = this.width;
    	this.buttonMenuHeight = this.height * .075;
    	this.buttonMenuWidth = this.width;
		this.svg.attr("width", this.width).attr("height", this.height);
		this.totalBars = 8; // Currently only 6 characteristics
		this.lastSort = "time";
		this.colorAsc = false;
		this.timeAsc = false;
		this.branchesAsc = false;
		this.samplesAsc = false;
		this.depthAsc = false;
		this.varianceAsc = false;
		this.boxIntersectionsAsc = false;
		this.objIntersectionsAsc = false;
		this.hidden = false;

		this.createGroups();
	};


	/* Group drag events */
	dragstarted(d) { d3.select(this).raise().classed("active", true); }
	dragged(d) {
		/* Determine the column we might swap with */
		let hoveredColumn = Math.round( Math.min( Math.max( d3.event.x , 0 ) ) / d.width );
		let secondClass = d3.select(".b" + hoveredColumn).attr("class").split(" ")[0];

		/* Store that column's class for reference*/
		d.newClass = secondClass;

		/* Translate the current column. */
		let translateToZero = [-d.x, 0];
		let translate = [Math.min( Math.max( d3.event.x , 0 )), 0];
		d3.select(this).attr("transform", "translate(" + translateToZero + ") translate(" + translate +")");

		/* Add active to this */
		d3.select(this).classed("active", true);
	}
	dragended(d) {
		/* Remove the transform from the dragged bar. */
		d3.select(this).attr("transform", "");



		/* Instead, swap classes between the bars. */
		let hoveredColumn = Math.round( Math.min( Math.max( d3.event.x , 0 ) ) / d.width );
		d3.select(this).attr("class", d.newClass + " b" + d.idx  + " bargroup");
		d3.select(".b" + hoveredColumn).attr("class",  d.class + " b" + hoveredColumn + " bargroup");
		
		if (d3.event.y < d.self.buttonMenuHeight)
			d.self.sort(d.sortField);
		
		/* Update the bars */
		d.self.update();


		let otherData = d3.select(".b" + hoveredColumn).datum();

		otherData.class = d.class;
		otherData.newClass = d.class;
		d.class = d.newClass;
		d.newClass = d.class;
		
		let temp = otherData.sortField;
		otherData.sortField = d.sortField;
		d.sortField = temp;
	}


	/* Creates a group for each characteristic bar */
	createGroups() {
		let heightScale = d3.scaleLinear().domain([0, 100]).range([0, this.height - this.barGroupHeight]);
		let widthScale = d3.scaleLinear().domain([0, 100]).range([0, this.width]);
		let invWidth = 1.0 / this.totalBars;

		this.svg.selectAll("g").remove();
		let bargroupsgroup = this.svg.append("g");

		let barData = [];
		let barNames = ["colorbar", "timebar", "branchesbar", "samplesbar", "depthbar", "variancebar", "boxIntersections", "objIntersections"];
		let sortFields = ["color", "time", "branches", "samples", "depth", "variance", "boxIntersections", "objIntersections"];
		for (let i = 0; i < this.totalBars; ++i) {
			let entry = {
				"class" : barNames[i],
				"newClass" : barNames[i],
				"idx" : i,
				"width" : this.width / this.totalBars,
				"self" : this,
				"sortField" : sortFields[i],
				"x": i * widthScale(100 / this.totalBars)
			}
			barData.push(entry);
		}

		let groups = bargroupsgroup.selectAll("g").data(barData);
		let enterGroups = groups.enter().append("g");
		let allGroups = groups.merge(enterGroups)

		allGroups.each(function(d) { this.classList.add(d.class); this.classList.add("b" + d.idx); this.classList.add("bargroup"); })
				.call(d3.drag().on("start", this.dragstarted).on("drag", this.dragged).on("end", this.dragended));


		this.svg.append("g").classed("selectableBar overlay", true);
		//this.update();
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

	update() {
		if (this.hidden) return;
		if (this.data) {
			this.updateButtonMenu();
			this.updateBars(this.data);
		}
	}

	/* Updates the visualization to show data for pixels within the provided bounds */
	updateBounds(bounds) {
		if (this.rawData == null) return;

		this.sanitizeBounds(bounds, this.rawData.dimensions);
		
		let skip = 0;
		if (bounds[2] - bounds[0] > 8) skip = 2; // Helps keep the total bars under control

		let self = this;
		let svg = this.svg;

		/* Wrangle data from bounds */
		let data=[];
		for (let y = bounds[1]; y < bounds[3]; y += 1 + skip) {
			for (let x = bounds[0]; x < bounds[2]; x += 1 + skip) {
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
				entry.boxIntersections = this.rawData.Box_intersections[x + y * this.rawData.dimensions[0]];
				entry.objIntersections = this.rawData.Primitive_intersections[x + y * this.rawData.dimensions[0]];
				entry.x = x;
				entry.y = y;
				data.push(entry);

			}
		}
		this.data = data;
		//this.createGroups();
		this.updateButtonMenu();
		this.updateBars(this.data);
		this.sort(this.lastSort);
	}

	
	/* Creates a specific button */
	createButton(svgGroup, label, sortFunction, parameter) {
		if (svgGroup.empty()) return;
		let data = svgGroup.datum();
		let index = data.idx;
		d3.select(".btn" + index + ".srt").remove();
		let heightScale = d3.scaleLinear().domain([0, 100]).range([0, this.buttonMenuHeight]);
		let widthScale = d3.scaleLinear().domain([0, 100]).range([0, this.width]);
		let invWidth = 1.0 / this.totalBars;

		let verticalPad = 2;//heightScale(10);
		let horizontalPad = 2;//widthScale(1);

		let group = svgGroup.append("g");
		group.classed("btn" + index, true).classed("srt", true);
		group.append("rect")
			.attr("x", (widthScale(index * invWidth * 100) + horizontalPad))
			.attr("y", heightScale(0) + verticalPad)
			.attr("width", (this.buttonMenuWidth * invWidth) - (2.0 * horizontalPad)) 
			.attr("height", heightScale(100) - (2.0 * verticalPad))
			.attr("class", "button")
			.classed("button", true);

		group.append("text").text(label)
			.attr("x", widthScale((index + .5) * invWidth * 100))
			.attr("y", heightScale(50))
			.style("text-anchor", "middle")
			.style("alignment-baseline", "middle")
			.classed("label", true).classed("unselectable", true);

		return group;
	}

	createMoveHandle(svgGroup) {
		if (svgGroup.empty()) return;
		
		let data = svgGroup.datum();
		let index = data.idx;
		d3.select(".btn" + index + ".hndl").remove();
		let heightScale = d3.scaleLinear().domain([0, 100]).range([0, this.buttonMenuHeight]);
		let widthScale = d3.scaleLinear().domain([0, 100]).range([0, this.width]);
		let yScale = d3.scaleLinear().domain([0, 100]).range([this.buttonMenuHeight + this.barGroupHeight, this.height]);
		let invWidth = 1.0 / this.totalBars;

		let verticalPad = 2;//yScale(10);
		let horizontalPad = 2;//widthScale(1);

		let x = (widthScale(index * invWidth * 100) + horizontalPad);
		let y = (yScale(0) + verticalPad);
		let btnWidth = (this.buttonMenuWidth * invWidth) - (2.0 * horizontalPad);
		let btnHeight = heightScale(100) - (2.0 * verticalPad);

		let group = svgGroup.append("g");
		group.classed(" btn" + index, true).classed("hndl", true);
		group.append("rect").attr("x", x).attr("y", y)
			.attr("width", btnWidth).attr("height", btnHeight)
			.attr("class", "button")
			.classed("button", true);

		group.append("line")
		.classed("unselectable", true)
		.attr("x1", x + btnWidth * .25)
		.attr("y1", y + btnHeight * .5)
		.attr("x2", x + btnWidth * .75)
		.attr("y2", y + btnHeight * .5)
		.attr("stroke-width", 2).attr("stroke", "white");

		group.append("line")
		.classed("unselectable", true)
		.attr("x1", x + btnWidth * .25)
		.attr("y1", y + btnHeight * .5)
		.attr("x2", x + btnWidth * .4)
		.attr("y2", y + btnHeight * .6)
		.attr("stroke-width", 2).attr("stroke", "white");

		group.append("line")
		.classed("unselectable", true)
		.attr("x1", x + btnWidth * .25)
		.attr("y1", y + btnHeight * .5)
		.attr("x2", x + btnWidth * .4)
		.attr("y2", y + btnHeight * .4)
		.attr("stroke-width", 2).attr("stroke", "white");


		group.append("line")
		.classed("unselectable", true)
		.attr("x1", x + btnWidth * .75)
		.attr("y1", y + btnHeight * .5)
		.attr("x2", x + btnWidth * .6)
		.attr("y2", y + btnHeight * .6)
		.attr("stroke-width", 2).attr("stroke", "white");

		group.append("line")
		.classed("unselectable", true)
		.attr("x1", x + btnWidth * .75)
		.attr("y1", y + btnHeight * .5)
		.attr("x2", x + btnWidth * .6)
		.attr("y2", y + btnHeight * .4)
		.attr("stroke-width", 2).attr("stroke", "white");

		return group;
	}

	updateButtonMenu() {
		/* Create sort buttons */
		this.createButton(d3.select(".colorbar"), "Pixel", this.sort, "color");
		this.createButton(d3.select(".timebar"), "Time", this.sort, "time");
		this.createButton(d3.select(".branchesbar"), "Branches", this.sort, "branches");
		this.createButton(d3.select(".samplesbar"), "Samples", this.sort, "samples");
		this.createButton(d3.select(".depthbar"), "Depth", this.sort, "depth");
		this.createButton(d3.select(".variancebar"), "Variance", this.sort, "variance");
		this.createButton(d3.select(".boxIntersections"), "Box", this.sort, "boxIntersections");
		this.createButton(d3.select(".objIntersections"), "Obj", this.sort, "objIntersections");
		
		/* Create move handles */
		this.createMoveHandle(d3.select(".colorbar"));
		this.createMoveHandle(d3.select(".timebar"));
		this.createMoveHandle(d3.select(".branchesbar"));
		this.createMoveHandle(d3.select(".samplesbar"));
		this.createMoveHandle(d3.select(".depthbar"));
		this.createMoveHandle(d3.select(".variancebar"));
		this.createMoveHandle(d3.select(".boxIntersections"));
		this.createMoveHandle(d3.select(".objIntersections"));
	}

	/* Updates a specific bar */
	updateBar(data, field, svgGroup, totalBars) {
		if (svgGroup.empty()) return;

		let parentData = svgGroup.datum();
		let index = parentData.idx; 

		let heightScale = d3.scaleLinear().domain([0, 100]).range([this.buttonMenuHeight, this.buttonMenuHeight + this.barGroupHeight]);
		let widthScale = d3.scaleLinear().domain([0, 100])
			.range([index * (this.barGroupWidth / totalBars), (index + 1) * (this.barGroupWidth / totalBars)]);
		let dataLenInv = 1 / data.length; 

		/* Remove all exiting rectangles*/
		let rects = svgGroup.selectAll(".bar").data(data);
		rects.exit().remove();

		/* Append entering rectangles*/
		let enterRects = rects.enter().append("rect").classed("bar", true);

		/* For each enter rectangle, initialize an x offset and width */
		enterRects.style("stroke-width", 0);

		let allRects = enterRects.merge(rects);
		allRects
			.attr("height", dataLenInv * this.barGroupHeight)
			.attr("x", (d, i) => {return widthScale(0);})
			.attr("width", this.barGroupWidth / totalBars)
			.attr("y", (d, i) => {return heightScale((i * dataLenInv) * 100)});

		return allRects;
	}
	
	/* Updates all bars given an array of pixel data */
	updateBars(data) {
		/* Pixel color */
		let colorBar = this.updateBar(data, "color", d3.select(".colorbar"), this.totalBars);
		if (colorBar) colorBar.style("fill", (d, i) => {return d3.rgb(d.color.red,d.color.blue,d.color.green);});

		/* Pixel render time */
		let timeBar = this.updateBar(data, "time", d3.select(".timebar"), this.totalBars);
		var timeColor = d3.scaleLinear()
			.domain([this.rawData.min_render_time, this.rawData.max_render_time])
            .range(['black','white']);
		if (timeBar) timeBar.style("fill", (d) => {return timeColor(d.time);});

		/* Total secondary rays */
		let secondaryRayBar = this.updateBar(data, "branches", d3.select(".branchesbar"), this.totalBars);
		var srColor = d3.scaleLinear()
			.domain([this.rawData.min_secondary_rays, this.rawData.max_secondary_rays])
            .range(['black','white']);
		if (secondaryRayBar) secondaryRayBar.style("fill", (d) => {return srColor(d.branches);});

		/* Total samples */
		let samplesBar = this.updateBar(data, "samples", d3.select(".samplesbar"), this.totalBars);
		var scColor = d3.scaleLinear()
			.domain([this.rawData.min_sample_count, this.rawData.max_sample_count])
            .range(['black','white']);
		if (samplesBar) samplesBar.style("fill", (d) => {return scColor(d.samples);});

		/* Pixel depth */
		let depthBar = this.updateBar(data, "depth", d3.select(".depthbar"), this.totalBars);
		var dColor = d3.scaleLinear()
			.domain([this.rawData.min_depth, this.rawData.max_depth])
            .range(['black','white']);
		if (depthBar) depthBar.style("fill", (d) => {return dColor(d.depth);});

		/* Pixel variance */
		let varianceBar = this.updateBar(data, "variances", d3.select(".variancebar"), this.totalBars);
		var vColor = d3.scaleLinear()
			.domain([this.rawData.min_variance, this.rawData.max_variance])
            .range(['black','white']);
		if (varianceBar) varianceBar.style("fill", (d) => {return vColor(d.variance);});

		/* Ray Box interections */
		let boxBar = this.updateBar(data, "boxIntersections", d3.select(".boxIntersections"), this.totalBars);
		var vColor = d3.scaleLinear()
			.domain([this.rawData.min_box_intersections, this.rawData.max_box_intersections])
            .range(['black','white']);
		if (boxBar) boxBar.style("fill", (d) => {return vColor(d.boxIntersections);});

		/* Ray Obj interections */
		let objBar = this.updateBar(data, "objIntersections", d3.select(".objIntersections"), this.totalBars);
		var vColor = d3.scaleLinear()
			.domain([this.rawData.min_primitive_intersections, this.rawData.max_primitive_intersections])
            .range(['black','white']);
		if (objBar) objBar.style("fill", (d) => {return vColor(d.objIntersections);});

		/* Selectable bars */
		let selectableBar = this.updateSelectable(data, d3.select(".selectableBar"))
	}

	/* Similar to updateBar, this updates a selectable row overlay. When rows are selected, 
			the treeview and image view are also updated. */
	updateSelectable(data, svgGroup) {
		let self = this;
		let heightScale = d3.scaleLinear().domain([0, 100]).range([this.buttonMenuHeight, this.buttonMenuHeight + this.barGroupHeight]);
		let widthScale = d3.scaleLinear().domain([0, 100]).range([0, this.barGroupWidth]);
		let dataLenInv = 1 / data.length; 
        
        /* Temporarily deselect the selected row, if selected. */
        d3.select(".row-overlay.active").classed("active", false);

        svgGroup.append("rect").classed(".bar", true);
		let rects = svgGroup.selectAll(".bar").data(data);
		rects.exit().remove();
		let enterRects = rects.enter().append("rect");

		enterRects
			.attr("x", 0)
			.attr("width", widthScale(100))
			.style("stroke-width", 0);

		let allRects = enterRects.merge(rects);
		allRects.attr("y", (d, i) => {return heightScale((i * dataLenInv) * 100)})
			.attr("height", dataLenInv * this.barGroupHeight)
			.classed("row-overlay", true)
			.classed("active", (d,i) => { return (d.x == self.selectedx && d.y == self.selectedy);})
			.on("mouseover", function(d,i) {
				d3.select(this).classed("hover", true);
				if (self.imageview != null)
          			self.imageview.hoverPixel(d.x, d.y);
			})
          	.on("mouseout", function(d,i) {
          		d3.select(this).classed("hover", false);
          		if (self.imageview != null)
          			self.imageview.hoverPixel(-1000, -1000);
          	})
          	.on("click", function(d,i) {
          		self.selectedx = d.x;
          		self.selectedy = d.y;
          		d3.select(".row-overlay.active").classed("active", false);
          		d3.select(this).classed("active", true);
          		if (self.imageview != null)
          			self.imageview.selectPixel(d.x, d.y);
          		if (self.raytreeview != null)
					self.raytreeview.selectPixel(d.x, d.y);

          		/* TODO: update tree view */
          	});

		return allRects;
	}

	/* Resizes the visualization. Good for device rotation/browser scaling */
	resize() {
		/* Update width and height just incase canvas size has changed */
		this.width = $("#pixel-analysis-viewport svg").parent().width();
    	this.height = $("#pixel-analysis-viewport svg").parent().height();
    	this.barGroupHeight = this.height * .85;
    	this.barGroupWidth = this.width;
    	this.buttonMenuHeight = this.height * .075;
    	this.buttonMenuWidth = this.width;
		this.svg.attr("width", this.width).attr("height", this.height);
		this.createGroups();

		this.update()
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
		this.lastSort = "color";
		this.data.sort(
            (x, y) => {
            	if (x.y - y.y < 0) return 1;
            	else if (x.y - y.y > 0) return -1;
            	else if (x.x - y.x < 0) return 1;
            	else if (x.x - y.x > 0) return -1;
            	else{
            		return 0;
            	} 
	        	//return this.rgbToHsv(x.color.red, x.color.green, x.color.blue)[0] 
            	//	 - this.rgbToHsv(y.color.red, y.color.green, y.color.blue)[0];
          }
        );

        if (!this.colorAsc) this.data.reverse();

        this.updateBars(this.data);
	}

	/* Sort the provided numberic field */
	sort(field) {
		this.data.reverse();
		if (field == "color") {this.sortColor(); return;}

		this.lastSort = field;

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

	setRayTreeView(raytreeview){
		this.raytreeview = raytreeview;
	}
}