/*
*	Justin Jensen
*	A RayTreeView displays the ray tree diagram for one or more samples of a single pixel.
*/

class RayTreeView {
	constructor(){
		this.treeDiv = d3.select("#pixel-analysis-viewport");
		this.svg = this.treeDiv.select("svg");
		this.width = $("#pixel-analysis-viewport").width();
    	this.height = $("#pixel-analysis-viewport").height();
    	this.svg.attr("width", this.width).attr("height", this.height);

    	this.hidden = true;
	};
	
	update(_sampleIdx){
		if (this.hidden) return;
		if (this.data == null) {
			this.svg.selectAll("text").remove();

			this.svg.append("text")
				.text("Select a pixel to see it's samples here.")
				.classed("h", true)
				.attr("text-anchor", "middle")
				.attr("x", this.width / 2.0)
				.attr("y", this.height / 2.0);
			return;
		} 

		var vWidth = this.width;
    	var vHeight = this.height;
    	
    	// Prepare our physical space
    	var g = this.svg.append('g').attr('transform', 'translate(' + vWidth/2 + ',' + vHeight/2 + ')');


		// Declare d3 layout
        var vLayout = d3.tree().size([2 * Math.PI, Math.min(vWidth, vHeight)/2 - 30]); // margin!

		// Layout + Data
        var vRoot = d3.hierarchy(this.data, function(d) {
        		if (!d.ch && d.length != 0) return d; 
        		if (d.ch.length == 0) return null;
        		return d.ch;
        	});
        var vNodes = vRoot.descendants();
        var vLinks = vLayout(vRoot).links();

        // Draw on screen

		var link = g.selectAll(".link")
		    .data(vLinks)
		    .enter().append("line")
		      .attr("class", "link")
		      .attr("stroke","#ccc")
		      .attr("x1", function(d) { return d3.pointRadial(d.source.x,d.source.y)[0]; })
		      .attr("y1", function(d) { return d3.pointRadial(d.source.x,d.source.y)[1]; })
		      .attr("x2", function(d) { return d3.pointRadial(d.target.x,d.target.y)[0]; })
		      .attr("y2", function(d) { return d3.pointRadial(d.target.x,d.target.y)[1]; }) ;

		let rectSize = 15;

        g.selectAll('rect').data(vNodes).enter().append('rect')
        	.classed("treePixel", true)
            .attr('x', (d) => {return d3.pointRadial(d.x, d.y)[0] - rectSize/2;} ).attr("y", (d) => {return d3.pointRadial(d.x, d.y)[1] - rectSize/2;})
            .attr('width', rectSize).attr("height", rectSize)
            .attr("fill", (d)=>{ 
            	console.log(d); 
            	if (d.data.c) 
            		return d3.rgb(d.data.c[0]*255, d.data.c[1]*255, d.data.c[2]*255); 
            	else return "";
            });

            //.attr("transform", function (d) { return "translate(" + d3.pointRadial(d.x, d.y) + ")"; });

		// /* concept art TEMPORARY */
		// this.svg.selectAll("image").remove();
		// this.image = this.svg.append("image");
		// this.image.attr("xlink:href", "./Sketches/Vis4DataSci1tree.png");
		// this.image.attr("height", this.height).attr("width", this.width);
	}
	
	resize(){
		//get the width and height of the div and copy those to the svg
		this.width = $("#pixel-analysis-viewport svg").parent().width();
    	this.height = $("#pixel-analysis-viewport svg").parent().height();

    	this.svg.attr("width", this.width)
				.attr("height", this.height);

		this.svg.selectAll("g").remove();
		
		this.update();
	}
	
	selectPixel(x, y){
		let self = this;
		/* For now, assume the image is 128 by 128, and 512 pixels per file in row major */
		let pixelIdx = x + y * 128;
		let fileIdx = Math.floor(pixelIdx / 512);
		let dataIdx = pixelIdx % 512;

		d3.json("./Data/128/raydata/raydata" + fileIdx + ".json", function(error, data) {
			if (!error) {
				self.data = data.RayData[dataIdx];
				self.update();
			} else {
				console.log(error);
			}
		});
	}
}
















