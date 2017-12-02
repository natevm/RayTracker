/*
	Nate Morrical. 
	An "PixelAnalysisView" displays a set of bars visualizing the current selection through a 
	variety of characteristics.
*/
class PixelAnalysisView {
	constructor() {
		/* Create menu */

		d3.select("#pixel-analysis-menu").html("");
		this.menuDiv = d3.select("#pixel-analysis-menu");
		this.menuSVG =  this.menuDiv.append("svg");
		this.menuWidth = $("#pixel-analysis-menu").width();
    	this.menuHeight = $("#pixel-analysis-menu").height();
		this.menuSVG.attr("width", this.menuWidth)
						.attr("height", this.menuHeight);

		this.checkpoint = 0;
		this.maxCheckpoint = 1;

		/* Create viewport */
		d3.select("#pixel-analysis-viewport").html("");
		this.viewportSVG = d3.select("#pixel-analysis-viewport").append("svg");
		this.viewportWidth = $("#pixel-analysis-viewport").width();
    	this.viewportHeight = $("#pixel-analysis-viewport").height();
		this.viewportSVG.attr("width", this.width).attr("height", this.height);

		this.parallelBarView = new ParallelBarView();
		this.rayTreeView = new RayTreeView();

		this.mode = 0;
	};

	update() {
		let buttonVertOffset = 10;
		let buttonWidth = 10;


		let self = this;
		let svg = this.menuSVG;
		let scaleHeight = d3.scaleLinear().domain([0,100]).range([0, this.menuHeight]);
		let scaleWidth = d3.scaleLinear().domain([0,100]).range([0, this.menuWidth]);

		/* Clear the menu */
		svg.selectAll("*").remove();

		if (this.mode == 0) {
			/* Add text for the current checkpoint */
			let textGroup = svg.append("g").classed("textGroup", true);
			let text = this.getText(this.checkpoint);

			textGroup.append("text")
			.text(text.header)
			.classed("h", true).classed("unselectable", true)
			.attr("text-anchor", "middle")
			.attr("x", scaleWidth(50)).attr("y", scaleHeight(25));


			textGroup.append("text")
			.text(text.l1)
			.classed("p", true).classed("unselectable", true)
			.attr("text-anchor", "middle")
			.attr("x", scaleWidth(50)).attr("y", scaleHeight(50));
			textGroup.append("text")
			.text(text.l2)
			.classed("p", true).classed("unselectable", true)
			.attr("text-anchor", "middle")
			.attr("x", scaleWidth(50)).attr("y", scaleHeight(65));
			textGroup.append("text")
			.text(text.l3)
			.classed("p", true).classed("unselectable", true)
			.attr("text-anchor", "middle")
			.attr("x", scaleWidth(50)).attr("y", scaleHeight(80));


			/* Add next/previous buttons */
			svg.append("rect")
			.classed("button", true).classed("toggled", (this.checkpoint==0))
			.attr("id", "imageMenuBackButton")
			.attr("x", scaleWidth(1)).attr("y", scaleHeight(buttonVertOffset))
			.attr("width", scaleWidth(buttonWidth)).attr("height", scaleHeight(80))
			.on("click", function() { if (self.checkpoint > 0) self.setCheckpoint(self.checkpoint-1);});

			svg.append("text").text("Bars").classed("p", true).classed("unselectable", true)
				.attr("x", scaleWidth(1) + .5 * scaleWidth(buttonWidth)).attr("y", scaleHeight(buttonVertOffset) + .5 * scaleHeight(80))
				.attr("text-anchor", "middle").attr("alignment-baseline", "middle");

			svg.append("rect")
			.attr("id", "imageMenuNextButton")
			.classed("button", true).classed("toggled", (this.checkpoint==1))
			.attr("fill", "red")
			.attr("x", scaleWidth(99-buttonWidth)).attr("y", scaleHeight(buttonVertOffset))
			.attr("width", scaleWidth(buttonWidth)).attr("height", scaleHeight(80))
			.on("click", function() { if (self.checkpoint < self.maxCheckpoint) self.setCheckpoint(self.checkpoint+1);});

			svg.append("text").text("Tree").classed("p", true).classed("unselectable", true)
				.attr("x", scaleWidth(99-buttonWidth) + .5 * scaleWidth(buttonWidth)).attr("y", scaleHeight(buttonVertOffset) + .5 * scaleHeight(80))
				.attr("text-anchor", "middle").attr("alignment-baseline", "middle");

		} else {
			this.createButton("Parallel Bar View", 0, 0, svg, 0);
			this.createButton("Ray Tree View", 0, 1, svg, 1);
		}

		if (this.checkpoint == 0) {
			this.parallelBarView.createGroups();
			this.parallelBarView.update();
		} else if (this.checkpoint == 1) {
			this.rayTreeView.update();
		}
	}

	createButton(name, row, column, parent, checkpoint) {
		let self = this;

		let totalRows = 1;
		let totalColumns = 2;

		let yScale = d3.scaleLinear().domain([0, 100]).range([0, this.menuWidth]);
		let xScale = d3.scaleLinear().domain([0, 100]).range([0, this.menuHeight]);
		let hScale = d3.scaleLinear().domain([0, 100]).range([0, this.menuWidth]);
		let wScale = d3.scaleLinear().domain([0, 100]).range([0, this.menuHeight]);

		let width = hScale(100 / totalColumns);
		let height = wScale(100 / totalRows);
		let x = yScale(100 / totalColumns) * column;
		let y = xScale(100 / totalRows) * row;

		parent.append("rect").classed("button", true).classed("toggled", (this.checkpoint == checkpoint))
			.attr("x", x).attr("y", y).attr("width", width).attr("height",height)
			.on("click", function() {
				self.setCheckpoint(checkpoint);
			});

		parent.append("text").text(name).classed("unselectable", true)
			.attr("x", x + width / 2.0).attr("y", y + height / 2.0).classed("p", true)
			.attr("text-anchor", "middle")
			.attr("alignment-baseline", "middle");
	}
	setCheckpoint(checkpoint) {
		this.viewportSVG.selectAll("*").remove();
		this.checkpoint = checkpoint;

		this.parallelBarView.hidden = (this.checkpoint != 0);
		this.rayTreeView.hidden = (this.checkpoint != 1);

		this.update();
	}

	getText(checkpoint) {
		switch(checkpoint) {
		    case 0:
		        return {"header": "Pixel Data",
		        "l1": "Hover over the rows and click to select a pixel.",
		    	"l2": "Click the button above a column to sort on that attribute.",
		    	"l3": "Clcik and drag the button below the column to move the column."
		    };
		    case 1:
		        return {"header": "Pixel Samples",
		        "l1": "The tree below shows how rays bounce and split for ",
		    	"l2": "the selected pixel. Colors from the leaves are combined",
		    	"l3": "to create the final pixel color in the center of the tree."
		    };
		    default:
		        return "";
		}
	}

	/* Resizes the visualization. Good for device rotation/browser scaling */
	resize() {
		/* Update width and height just incase canvas size has changed */
    	this.menuWidth = $("#pixel-analysis-menu").width();
    	this.menuHeight = $("#pixel-analysis-menu").height();
		this.menuSVG.attr("width", this.menuWidth)
						.attr("height", this.menuHeight);

		if (this.checkpoint == 0) {
			this.rayTreeView.resize();
			this.parallelBarView.resize();
		} else if (this.checkpoint = 1) {
			this.parallelBarView.resize();
			
			this.rayTreeView.resize();
		}
						
		this.update()
	}

	setData(_data){
		this.rawData = _data;
		this.parallelBarView.setData(this.rawData);
	}

	toggleExploreMode() {
		this.mode = 1;
		this.update();
	}

	toggleStoryMode() {
		this.mode = 0;
		this.update();
	}
}