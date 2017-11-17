/*
*	Justin Jensen
*	A RayTreeView displays the ray tree diagram for one or more samples of a single pixel.
*/

class RayTreeView {
	constructor(){
		this.treeDiv = d3.select("#tree-view");
		this.svg = this.treeDiv.append("svg");
		this.width = $("#tree-view").width();
		this.height = $("#tree-view").height();
		this.svg.attr("width", this.width)
				.attr("height", this.height);
		//create a group in the svg for easy transformation
		this.svg.append("g").attr("id", "tree");
	};
	
	update(_sampleIdx){
		//clear the svg
		this.svg.selectAll("g#tree > g.node").remove();
		this.svg.selectAll("g#tree > path.link").remove();
		
		//if _sampleIdx is undefined
		if(_sampleIdx === undefined){
			//set it to -1
			_sampleIdx = -1;
		}
		//if _sampleIdx is -1
		if(_sampleIdx == -1){
			//we'll just leave it cleared
			return;
		}
		
		//extract the data for _sampleIdx
		let treeMap = d3.tree().size([this.width,this.height]);
		
		/*let root = d3.stratify()
			.id(function(d, i){
				return i;
			})
			.parentId(function(d){
				return 
			})*/
		//root = d3.hierarchy(treeData, function(d) { return d.children; });
		let root = d3.hierarchy(this.data[_sampleIdx][0], function(d){
			return d.children;
		});
		
		let theTreeData = treeMap(root);
		let nodes = theTreeData.descendants(),
			links = theTreeData.descendants().slice(1);
			
		
		//draw the tree
		this.node = d3.select("g#tree").selectAll("g.node")
			.data(nodes);
		let nodeEnter = this.node.enter().append('g')
			.attr('class', 'node')
			.attr('transform', function(d){
				return `translate(${d.y}, ${d.x})`;
			});
		nodeEnter.append('circle')
			.attr('class', 'node')
			.attr('r', 6);
		this.node = this.node.merge(nodeEnter);
		
		this.link = d3.select('g#tree').selectAll('path.link')
			.data(links);
		let linkEnter = this.link.enter().insert('path', 'g')
			.attr('class', 'link')
			.attr('d', function(d){
				return diagonal(d, d.parent);
			})
			.attr("fill", "transparent")
			.attr("stroke", "white");
		this.link = this.link.merge(linkEnter);
		
		function diagonal(s, d){
			return `M ${s.y} ${s.x}
					C ${(s.y + d.y) / 2} ${s.x},
					  ${(s.y + d.y) / 2} ${d.x},
					  ${d.y} ${d.x}`;
		}
	}
	
	resize(){
		//get the width and height of the div and copy those to the svg
		let bounds = this.treeDiv.node().getBoundingClientRect();
		this.width = bounds.width;
		this.height = bounds.height;
		this.svg.attr("width", this.width)
				.attr("height", this.height);
		
		this.update();
	}
	
	setData(_data){
		this.rawData = _data;
		//wrangle the data into a format that D3 trees understand
		/*
		The format should look something like this:
		{
			"name":"theName",
			"children":[
				{
					"name":"the first child",
					"children":[]
				},
				{
					"name":"the second child",
					"children":[]
				}
			]
		}
		*/
		
		this.data = [];
		//for each pixel
		for(let i = 0; i < this.rawData.length; i++){
			let pixelData = this.rawData[i];
			//for each sample
			let pixel = [];
			for(let s = 0; s < pixelData.length; s++){
				pixel.push(r_processData(pixelData[s]));
			}
			this.data.push(pixel);
		}
		
		var rayTypes = [
			'camera',
			'reflection',
			'refraction',
			'shadow'
		];
		function r_processData(rawNode){
			let newNode = {};
			//add "name" to newNode
			newNode.name = "TODO";
			//add node data
			//TODO
			//add link data
			newNode.rayType = +rawNode['t'];
			//add "children"[] to newNode
			newNode.children = [];
			//for each object in rawNode["ch"]
			for(let i = 0; i < rawNode['ch'].length; i++){
				//call r_processData on the rawNode's child
				//add the returned object to "children"
				newNode.children.push(r_processData(rawNode['ch'][i]));
			}
			//return newNode
			return newNode;
		}
	}
}
















