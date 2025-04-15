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

  // The force simulation mutates links and nodes, so create a copy
  // so that re-evaluating this cell produces the same result.
  const allNodes = data.nodes.map(node => ({...node, id: String(node.id)})).filter(node => node.condition_sensitivity > 0.05);
  const allLinks = Object.entries(data.edges).flatMap(([source, targets]) => targets.map(t => ({source: String(source), target: String(t)})));
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

  const container = d3.create("div")

  // Create the SVG container.
  const outerSvg = container.append("svg")
      .attr("width", width)
      .attr("height", height)
      .attr("viewBox", [-width, -height, width * 2, height * 2])
      .attr("style", "max-width: 100%; height: auto;");


  const searchText = container.append("input").attr("type", "text");
  const searchButton = container.append("button").text("Go");

  const displayText = container.append("xhtml:body");


  // Add a zoom behavior to container
  const zoom = d3.zoom().scaleExtent([1, 5]).translateExtent([[-width, -height], [width * 2, height * 2]])
  // outerSvg.call(zoom
  //   .on("zoom", (evt) => {
  //     const prevViewBox = outerSvg.attr("viewBox").split(",").map(s => Number.parseInt(s))
  //     const scale = evt.transform.k
  //     console.log(prevViewBox)
  //     console.log(evt.transform)
  //     outerSvg.attr("viewBox", [[prevViewBox[0] - 2 * evt.dx, prevViewBox[1] - 2 * evt.dy], prevViewBox.slice(2, 4)])
  //     // console.log(evt.transform)
  //     // outerSvg.attr("transform", `translate(${evt.transform.x},${evt.transform.y}) scale(${evt.transform.k})`);
  //   })
  // );

  // Add a drag behavior to container
  // outerSvg.call(d3.drag()
  //   .on("drag", (evt) => {
  //     console.log(evt)
  //     if (evt.sourceEvent.type !== "mousemove") {
  //       // outerSvg.attr("transform", `translate(${evt.transform.x},${evt.transform.y}) scale(${evt.transform.k})`);
  //     }
  //   })
  // );

  // outerSvg.call(d3.drag()
  //   .on("drag", (evt) => {
  //     const prevViewBox = outerSvg.attr("viewBox").split(",").map(s => Number.parseInt(s))
  //     outerSvg.attr("viewBox", [[prevViewBox[0] - 2 * evt.dx, prevViewBox[1] - 2 * evt.dy], prevViewBox.slice(2, 4)])
  //   })
  // );

  const svg = outerSvg.append("g");

  // Add an arrow for each link, and a circle/triangle for each node.
  svg.append("defs").append("marker")
    .attr("id", "arrow")
    .attr("viewBox", "0 -5 10 10")
    .attr("refX", 20)
    .attr("refY", 0)
    .attr("markerWidth", 2)
    .attr("markerHeight", 2)
    .attr("orient", "auto")
  .append("svg:path")
    .attr("d", "M0,-5L10,0L0,5");

  const link = svg.append("g")
      .attr("stroke", "#999")
      .attr("stroke-opacity", 0.6)
    .selectAll("line")
    .data(links)
    .join("line")
      .attr("stroke-width", 3)
      .attr("marker-end", "url(#arrow)");

  const node = svg.append("g")
      .attr("stroke", "#fff")
      .attr("stroke-width", 1.5)
    .selectAll("path")
    .data(nodes)
    .join("path")
      .attr("id", d => `id${d.id}`)
      .attr("fill", d => color(d.organ_specificity > d.universality))
      .attr("d", d => typeToSymbol[d.node_type].size(100)());

  // Add a drag behavior.
  node.call(d3.drag()
        .on("start", dragstarted)
        .on("drag", dragged)
        .on("end", dragended));

  node.on("click", (_evt, d) => selectNode(d.id))

  searchButton.on("click", (_evt) => selectNode(searchText.property("value")))
  
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

  // Zoom to a node with the given data, show data in displayText
  function selectNode(nodeId) {
    const d = nodes.find(n => n.id === nodeId)
    if (!d) {
      displayText.html(`Gene/Peak ${nodeId} not found.`)
      return
    }
    node.attr("stroke", "#fff")
    svg.select(`#id${d.id}`).attr("stroke", "red")
    displayText.html(d.node_type === 'peak' ? 
        `Peak: ${d.id}<br/>Condition sensitivity score: ${d.condition_sensitivity.toFixed(5)}<br/>Organ specificity score: ${d.organ_specificity.toFixed(5)}<br/>Universality score: ${d.universality.toFixed(5)}` : 
        `Gene: ${d.id}<br/>Condition sensitivity score: ${d.condition_sensitivity.toFixed(5)}<br/>Organ specificity score: ${d.organ_specificity.toFixed(5)}<br/>Universality score: ${d.universality.toFixed(5)}`);

    svg.transition().duration(700).call(
      zoom.transform,
      d3.zoomIdentity.translate((width / 2), (height / 2)).scale(4).translate(-d.x, -d.y),
    )
  }

  // When this cell is re-run, stop the previous simulation. (This doesn’t
  // really matter since the target alpha is zero and the simulation will
  // stop naturally, but it’s a good practice.)
  // invalidation.then(() => simulation.stop());

  return container.node();
}

export { build };