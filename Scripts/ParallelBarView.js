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

		this.totalBars = 7; // Currently only 7 characteristics
	};

	sanitize(bounds, dimensions) {
		bounds[0] = Math.min(bounds[0], dimensions[0]);
		bounds[1] = Math.min(bounds[1], dimensions[1]);
		bounds[2] = Math.min(bounds[2], dimensions[0]);
		bounds[3] = Math.min(bounds[3], dimensions[1]);

		bounds[0] = Math.max(bounds[0], 0);
		bounds[1] = Math.max(bounds[1], 0);
		bounds[2] = Math.max(bounds[2], 0);
		bounds[3] = Math.max(bounds[3], 0);
	}

	update(bounds) {
		if (this.rawData != null) {
			this.sanitize(bounds, this.rawData.dimensions);
			
			let self = this;
			let svg = this.svg;

			svg.selectAll("g").remove();
			svg.selectAll("rect").remove();
			
			this.updateColors(bounds);
			this.updateTimes(bounds);
			this.updateSecondaryRays(bounds);
			this.updateSampleCounts(bounds);
			this.updateDepths(bounds);
			this.updateVariances(bounds);
		}
	}

	updateColors(bounds) {
		let colors = [];
		for (let y = bounds[1]; y < bounds[3]; ++y) {
			for (let x = bounds[0]; x < bounds[2]; ++x) {
				let color = {
					"red" : -1,
					"blue" : -1,
					"green" : -1
				};
				color.red = this.rawData.Color[3 * (x + y * this.rawData.dimensions[0])];
				color.blue = this.rawData.Color[3 * (x + y * this.rawData.dimensions[0]) + 1];
				color.green = this.rawData.Color[3 * (x + y * this.rawData.dimensions[0]) + 2];

				colors.push(color);
			}
		}

		let group = this.svg.append("g");
		let colorGroup = group.append("g");

		let colorRects = colorGroup.selectAll("rect").data(colors);
		colorRects.exit().remove();
		let enterColorRects = colorRects.enter().append("rect");
		let allColorRects = enterColorRects.merge(colorRects);

		allColorRects
			.attr("x", (d, i) => {return 0})
			.attr("y", (d, i) => {return (i / colors.length) * this.height})
			.attr("width", (d, i) => {return 100})
			.attr("height", (d, i) => {return (1 / colors.length) * this.height})
			.style("fill", (d, i) => {return d3.rgb(d.red,d.blue,d.green);});
	}

	updateTimes(bounds) {
		var color = d3.scaleLinear()
			.domain([this.rawData.min_render_time, this.rawData.max_render_time])
            .range(['white','red']);

		let times = [];
		for (let y = bounds[1]; y < bounds[3]; ++y) {
			for (let x = bounds[0]; x < bounds[2]; ++x) {
				times.push(this.rawData.Render_time[x + y * this.rawData.dimensions[0]]);
			}
		}

		let group = this.svg.append("g");
		let timeGroup = group.append("g");

		let timeRects = timeGroup.selectAll("rect").data(times);
		timeRects.exit().remove();
		let enterTimeRects = timeRects.enter().append("rect");
		let allTimeRects = enterTimeRects.merge(timeRects);

		allTimeRects
			.attr("x", (d, i) => {return 100})
			.attr("y", (d, i) => {return (i / times.length) * this.height})
			.attr("width", (d, i) => {return 100})
			.attr("height", (d, i) => {return (1 / times.length) * this.height})
			.style("fill", (d, i) => {return color(d);});
	}

	updateSecondaryRays(bounds) {
		var color = d3.scaleLinear()
			.domain([this.rawData.min_secondary_rays, this.rawData.max_secondary_rays])
            .range(['white','red']);

		let secondaryRays = [];
		for (let y = bounds[1]; y < bounds[3]; ++y) {
			for (let x = bounds[0]; x < bounds[2]; ++x) {
				secondaryRays.push(this.rawData.Secondary_rays[x + y * this.rawData.dimensions[0]]);
			}
		}

		let group = this.svg.append("g");
		let secondaryRayGroup = group.append("g");

		let secondaryRayRects = secondaryRayGroup.selectAll("rect").data(secondaryRays);
		secondaryRayRects.exit().remove();
		let enterSecondaryRayRects = secondaryRayRects.enter().append("rect");
		let allSecondaryRayRects = enterSecondaryRayRects.merge(secondaryRayRects);

		allSecondaryRayRects
			.attr("x", (d, i) => {return 200})
			.attr("y", (d, i) => {return (i / secondaryRays.length) * this.height})
			.attr("width", (d, i) => {return 100})
			.attr("height", (d, i) => {return (1 / secondaryRays.length) * this.height})
			.style("fill", (d, i) => {return color(d);});
	}
	
	updateSampleCounts(bounds) {
		var color = d3.scaleLinear()
			.domain([this.rawData.min_sample_count, this.rawData.max_sample_count])
            .range(['white','red']);

		let sampleCounts = [];
		for (let y = bounds[1]; y < bounds[3]; ++y) {
			for (let x = bounds[0]; x < bounds[2]; ++x) {
				sampleCounts.push(this.rawData.Sample_count[x + y * this.rawData.dimensions[0]]);
			}
		}

		let group = this.svg.append("g");
		let sampleCountGroup = group.append("g");

		let sampleCountRects = sampleCountGroup.selectAll("rect").data(sampleCounts);
		sampleCountRects.exit().remove();
		let enterSampleCountRects = sampleCountRects.enter().append("rect");
		let allSampleCountRects = enterSampleCountRects.merge(sampleCountRects);

		allSampleCountRects
			.attr("x", (d, i) => {return 300})
			.attr("y", (d, i) => {return (i / sampleCounts.length) * this.height})
			.attr("width", (d, i) => {return 100})
			.attr("height", (d, i) => {return (1 / sampleCounts.length) * this.height})
			.style("fill", (d, i) => {return color(d);});
	}

	updateDepths(bounds) {
		var color = d3.scaleLinear()
			.domain([this.rawData.min_depth, this.rawData.max_depth])
            .range(['white','red']);

		let depths = [];
		for (let y = bounds[1]; y < bounds[3]; ++y) {
			for (let x = bounds[0]; x < bounds[2]; ++x) {
				depths.push(this.rawData.Depth_buffer[x + y * this.rawData.dimensions[0]]);
			}
		}

		let group = this.svg.append("g");
		let depthGroup = group.append("g");

		let depthRects = depthGroup.selectAll("rect").data(depths);
		depthRects.exit().remove();
		let enterDepthRects = depthRects.enter().append("rect");
		let allDepthRects = enterDepthRects.merge(depthRects);

		allDepthRects
			.attr("x", (d, i) => {return 400})
			.attr("y", (d, i) => {return (i / depths.length) * this.height})
			.attr("width", (d, i) => {return 100})
			.attr("height", (d, i) => {return (1 / depths.length) * this.height})
			.style("fill", (d, i) => {return color(d);});
	}


	updateVariances(bounds){
		var color = d3.scaleLinear()
			.domain([this.rawData.min_variance, this.rawData.max_variance])
            .range(['white','red']);

		let variances = [];
		for (let y = bounds[1]; y < bounds[3]; ++y) {
			for (let x = bounds[0]; x < bounds[2]; ++x) {
				variances.push(this.rawData.Variance[x + y * this.rawData.dimensions[0]]);
			}
		}

		let group = this.svg.append("g");
		let varianceGroup = group.append("g");

		let varianceRects = varianceGroup.selectAll("rect").data(variances);
		varianceRects.exit().remove();
		let enterVarianceRects = varianceRects.enter().append("rect");
		let allVarianceRects = enterVarianceRects.merge(varianceRects);

		allVarianceRects
			.attr("x", (d, i) => {return 500})
			.attr("y", (d, i) => {return (i / variances.length) * this.height})
			.attr("width", (d, i) => {return 100})
			.attr("height", (d, i) => {return (1 / variances.length) * this.height})
			.style("fill", (d, i) => {return color(d);});
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