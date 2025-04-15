import * as d3 from "https://cdn.jsdelivr.net/npm/d3@7/+esm";

// Declare the chart dimensions and margins.
const width = 928;
const height = 800;
// const marginTop = 20;
// const marginRight = 20;
// const marginBottom = 30;
// const marginLeft = 40;

const typeToSymbol = {
  'peak': d3.symbol().type(d3.symbolTriangle),
  'gene': d3.symbol().type(d3.symbolCircle)
}

const build = async (grnPath) => {
  const data = await d3.json(grnPath);

  // Specify the color scale.
  const color = d3.scaleOrdinal(d3.schemeCategory10);
  const organSpecific = color(true)
  const universal = color(false)

  // The force simulation mutates links and nodes, so create a copy
  // so that re-evaluating this cell produces the same result.
  const allNodes = data.nodes.map(node => ({...node, id: String(node.id)})).filter(node => node.condition_sensitivity > 0.05);
  const allLinks = Object.entries(data.edges).flatMap(([source, targets]) => targets.map(t => ({source: String(source), target: String(t)})));

  // Filter out links without both nodes and nodes with no links iteratively
  let nodes = allNodes;
  let links = allLinks;
  let nodeLength = nodes.length + 1
  while (nodeLength !== nodes.length) {
    nodeLength = nodes.length
    const nodeIdSet = new Set(nodes.map(n => n.id));
    links = links.filter(edge => nodeIdSet.has(edge.source) && nodeIdSet.has(edge.target));
    links.forEach(edge => nodeIdSet.delete(edge.source) && nodeIdSet.delete(edge.target));
    nodes = nodes.filter(node => !nodeIdSet.has(node.id));
  }

  // Create a simulation with several forces.
  const simulation = d3.forceSimulation(nodes)
      .force("link", d3.forceLink(links).id(d => d.id))
      .force("charge", d3.forceManyBody())
      .force("x", d3.forceX())
      .force("y", d3.forceY());

  const container = d3.create("div").style("display", "flex")

  // Create the SVG container.
  const outerSvg = container.append("svg")
      .attr("width", width)
      .attr("height", height)
      .attr("style", "max-width: 100%; height: auto; flex: 1;");


  const textContainer = container.append("div").attr("style", 
    "position: fixed; top: 0; right: 0; padding: 0.5rem; width: 15rem; background-color: #fff; border: 1px solid #000; border-top: 0; border-right: 0; border-radius: 0 0 0 3px;");

  // Add search functionality
  const searchForm = textContainer.append("form");
  const searchText = searchForm.append("input").attr("type", "text");
  searchForm.append("button").attr("type", "submit").text("Go");
  searchForm.on("submit", (evt) => {
    evt.preventDefault();
    selectNode(searchText.property("value"));
    return false;
  });
  
  // Add display text for showing selected node information
  let selected = '';
  const displayText = textContainer.append("xhtml:p").style("margin", "0.5rem 0");
  displayText.html(" <br> <br> <br> <br>");

  const svg = outerSvg.append("g");

  // Add an arrow for each link, and a circle/triangle for each node.
  svg.append("defs").append("marker")
    .attr("id", "arrow")
    .attr("viewBox", "0 -5 10 10")
    .attr("refX", 20)
    .attr("refY", 0)
    .attr("markerWidth", 5)
    .attr("markerHeight", 5)
    .attr("orient", "auto")
  .append("svg:path")
    .attr("d", "M0,-5L10,0L0,5");

  const link = svg.append("g")
      .attr("stroke", "#999")
      .attr("stroke-opacity", 0.6)
    .selectAll("line")
    .data(links)
    .join("line")
      .attr("stroke-width", 1)
      // .attr("marker-height", 3)
      .attr("marker-end", "url(#arrow)");

  const node = svg.append("g")
      .attr("stroke", "#fff")
      .attr("stroke-width", 1.5)
    .selectAll("path")
    .data(nodes)
    .join("path")
      .attr("cursor", "pointer")
      .attr("id", d => `id${d.id}`)
      .attr("fill", d => color(d.organ_specificity > d.universality))
      .attr("d", d => typeToSymbol[d.node_type].size(100)());

  // Add a legend for colors
  textContainer.append("svg")
    .attr("width", 12)
    .attr("height", 12)
    .append("rect")
      .attr("width", 12)
      .attr("height", 12)
      .attr("rx", 2)
      .attr("fill", organSpecific)
  textContainer.append("span").style("margin", "0 15px 0 5px").text("Organ-specific")
  textContainer.append("svg")
    .attr("width", 12)
    .attr("height", 12)
    .append("rect")
      .attr("width", 12)
      .attr("height", 12)
      .attr("rx", 2)
      .attr("fill", universal)
  textContainer.append("span").style("margin", "0 5px").text("Universal")

  // Add a legend for shapes
  textContainer.append("br")
  textContainer.append("svg")
    .attr("width", 12)
    .attr("height", 12)
    .attr("viewBox", "-8 -10 16 16")
    .append("path")
      .attr("d", typeToSymbol['peak']())
      .attr("fill", "black")
  textContainer.append("span").style("margin", "0 15px 0 5px").text("Peak")
  textContainer.append("svg")
    .attr("width", 12)
    .attr("height", 12)
    .attr("viewBox", "-6 -6 12 12")
    .append("path")
      .attr("d", typeToSymbol['gene']())
      .attr("fill", "black")
  textContainer.append("span").style("margin", "0 15px 0 5px").text("Gene")

  // Add a drag behavior.
  node.call(d3.drag()
        .on("start", dragstarted)
        .on("drag", dragged)
        .on("end", dragended));

  // Add tooltips for nodes
  const tooltip = container.append("p")
    .style("opacity", 0)
    .style("position", "absolute")
    .style("margin", 0)
    .style("background-color", "#fff")
    .style("border", "1px solid black")
    .style("border-radius", "3px")
    .style("padding", "0.25rem 0.5rem");
  
  node.on("mouseover", (evt, d) => {
    tooltip.style("opacity", 1).style("top", `${evt.pageY - 10}px`).style("left", `${evt.pageX + 10}px`).text(d.id);
  })
  node.on("mouseleave", (_evt) => {
    tooltip.style("opacity", 0)
  })

  // On click node, select that node
  node.on("click", (_evt, d) => selectNode(d.id))

  // On click edge, if one node is selected, select the other node
  link.on("click", (_evt, d) => {
    if (selected === d.source.id) { 
      selectNode(d.target.id) 
    } else if (selected === d.target.id) {
      selectNode(d.source.id)
    }
})

  // Set the position attributes of links and nodes each time the simulation ticks.
  simulation.on("tick", () => {
    link
        .attr("x1", d => d.source.x)
        .attr("y1", d => d.source.y)
        .attr("x2", d => d.target.x)
        .attr("y2", d => d.target.y);

    node
      .attr("transform", d => `translate(${d.x},${d.y})`);
  });

  // Reheat the simulation when drag starts, and fix the subject position.
  function dragstarted(event) {
    if (!event.active) simulation.alphaTarget(0.3).restart();
    event.subject.fx = event.subject.x;
    event.subject.fy = event.subject.y;
  }

  // Update the subject (dragged node) position during drag.
  function dragged(event) {
    event.subject.fx = event.x;
    event.subject.fy = event.y;
  }

  // Restore the target alpha so the simulation cools after dragging ends.
  // Unfix the subject position now that it’s no longer being dragged.
  function dragended(event) {
    if (!event.active) simulation.alphaTarget(0);
    event.subject.fx = null;
    event.subject.fy = null;
  }

  // Allow for panning and zooming
  let currTransform = [0, 0, width]

  svg.attr("transform", transform(currTransform))

  outerSvg.call(d3.drag()
    .on("drag", (evt) => {
      const scale = currTransform[2] / height
      currTransform = [currTransform[0] - scale * evt.dx, currTransform[1] - scale * evt.dy, currTransform[2]]
      svg.attr("transform", transform(currTransform))
    })
  );

  let zoom = d3.zoom()
  outerSvg.call(zoom
    .on("zoom", ({ sourceEvent }) => {
      if (sourceEvent.type === "wheel") {
        currTransform = [currTransform[0], currTransform[1], currTransform[2] - currTransform[2] * sourceEvent.wheelDelta / height]
        svg.attr("transform", transform(currTransform));
      }
    })
  );

  // Zoom to a node with the given data, show data in displayText
  function selectNode(nodeId) {
    node.attr("stroke", "#fff")
    const d = nodes.find(n => n.id === nodeId)
    selected = nodeId
    if (!d) {
      displayText.html(`Gene/Peak '${nodeId}' not found.`)
      return
    }
    svg.select(`#id${d.id}`).attr("stroke", "red")
    console.log(d)
    displayText.html(`${d.node_type === 'peak' ? 'Peak' : 'Gene'}: ${d.id}<br/>
      Condition sensitivity score: ${d.condition_sensitivity.toFixed(5)}<br/>
      Organ specificity score: ${d.organ_specificity.toFixed(5)}<br/>
      Universality score: ${d.universality.toFixed(5)}`);

    const endTransform = [d.x, d.y, width / 5];
    const i = d3.interpolateZoom(currTransform, endTransform)

    svg.transition()
      .duration(i.duration)
      .attrTween("transform", () => t => transform(currTransform = i(t)));
    d3.zoomTransform(outerSvg, endTransform)
  }

  // Used for transforming the view to a specific xy coordinate and with view radius r
  function transform([x, y, r]) {
      return `
        translate(${width / 2}, ${height / 2})
        scale(${height / r})
        translate(${-x}, ${-y})
      `;
    }

  // When this cell is re-run, stop the previous simulation. (This doesn’t
  // really matter since the target alpha is zero and the simulation will
  // stop naturally, but it’s a good practice.)
  // invalidation.then(() => simulation.stop());

  return container.node();
}

export { build };