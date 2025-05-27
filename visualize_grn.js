import * as d3 from "https://cdn.jsdelivr.net/npm/d3@7/+esm";

// Declare the chart dimensions and margins.
let width = window.innerWidth;
let height = window.innerHeight;

// Used for transforming the view to a specific xy coordinate and with view radius r
function transform([x, y, r]) {
  return `
    translate(${width / 2}, ${height / 2})
    scale(${height / r})
    translate(${-x}, ${-y})
  `;
}

function dataToSymbol(data) {
  if (data.node_type === 'peak') {
    return d3.symbol().type(d3.symbolTriangle)
  } else if (data.node_type === 'gene') {
    if (data.tf_coding) {
      return d3.symbol().type(d3.symbolSquare)
    } else {
      return d3.symbol().type(d3.symbolCircle)
    }
  }
}

function dataToColor(data) {
  const red = data.Heart.organ_specificity
  const green = data.Lung.organ_specificity
  const blue = data.Liver.organ_specificity
  const scalePositive = [red, green, blue].map(col => col + 1 / 3).map(col => Math.max(col, 0.0))
  const max = Math.max(...scalePositive)
  const hex = scalePositive.map(col => col / max).map(col => Math.round(col * 255).toString(16).padStart(2, "0")).join("")
  return `#${hex}`;
}

async function fetchWithRetry(url, options = undefined, nRetries = 5) {
  const statusCodes = []
  for (let i = 0; i < nRetries; i += 1) {
    const res = await fetch(url, options);
    if (res.ok) {
      return await res.json();
    } else {
      statusCodes.push(res.status)
      await new Promise(resolve => setTimeout(resolve, 1000))
    }
  }
  return await Promise.reject(`Failed to fetch, even with retry. Status codes: ${statusCodes.join(", ")}`);
}


const build = (grnPath) => {
  // Start fetching grn data, but build the static parts and display in the meantime, then load the grn data into the svg
  const res = fetchWithRetry(grnPath);

  const container = d3.create("div").style("display", "flex")

  // Create the SVG container.
  const outerSvg = container.append("svg")
      .attr("width", width)
      .attr("height", height)
      .attr("style", "max-width: 100%; flex: 1;");


  const textContainer = container.append("div").attr("style", 
    "position: fixed; top: 0; right: 0; padding: 0.5rem; width: 15rem; background-color: #fff; border: 1px solid #000; border-top: 0; border-right: 0; border-radius: 0 0 0 3px;");

  // Add search functionality
  const searchForm = textContainer.append("form");
  const searchText = searchForm.append("input").attr("type", "text");
  searchForm.append("button").attr("type", "submit").text("Go");
  
  // Add display text for showing selected node information
  const displayText = textContainer.append("xhtml:p").style("margin", "0.5rem 0");
  const scoreText = textContainer.append("xhtml:p").style("margin", "0.5rem 0");
  const displayTextBreaks = " <br> <br> <br> <br> <br> <br>"
  const emptyDisplayText = `Click around or search to see more.${displayTextBreaks}`
  displayText.html(emptyDisplayText);

  const svg = outerSvg.append("g");

  const loadingText = svg.append("text").text("Loading...")

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

  // Add a legend for colors
  textContainer.append("p").style("margin", "0").style("text-align", "center").text("Heart")
  
  const colorKey = textContainer.append("div")
    .style("margin", "0")
    .style("justify-self", "center")

  const sideStep = 0.02
  const vertColorStep = sideStep / Math.sqrt(3)  // Length = sqrt(3)x
  const horColorStep = sideStep / 2  // Length = 2x
  for (let heart = 1 + vertColorStep; heart > -vertColorStep; heart -= vertColorStep) {
    const row = colorKey.append("div")
      .style("display", "block")
      .style("height", "2px")
      .style("margin", 0)
      .style("padding", 0)
    const addColor = (h, l) => {
      const liver = 1 - (h + l)
      let fillColor = '#fff'
      if (liver >= -horColorStep) {
        fillColor = dataToColor({
          Heart: { organ_specificity: h - 1/3 },
          Lung: { organ_specificity: l - 1/3 },
          Liver: { organ_specificity: liver - 1/3 }
        })
      }
      const rgb = [fillColor.slice(1, 3), fillColor.slice(3, 5), fillColor.slice(5, 7)].map(n => Number.parseInt(n, 16))
      const sumScaled = Math.pow(Math.round((rgb[0] + rgb[1] + rgb[2]) / 255) + 1, 2)
      const colorsWithinBuffer = []
      for (let offsetR = -sumScaled; offsetR <= sumScaled; offsetR += 1) {
        for (let offsetG = -sumScaled; offsetG <= sumScaled; offsetG += 1) {
          for (let offsetB = -sumScaled; offsetB <= sumScaled; offsetB += 1) {
            colorsWithinBuffer.push((rgb[0] + offsetR).toString(16).padStart(2, "0") + 
                                    (rgb[1] + offsetG).toString(16).padStart(2, "0") + 
                                    (rgb[2] + offsetB).toString(16).padStart(2, "0"))
          }
        }
      }
      row.append("svg")
      .attr("width", 2)
      .attr("height", 2)
      .style("margin", 0)
      .style("vertical-align", "top")
      .append("rect")
        .attr("width", 2)
        .attr("height", 2)
        .attr("fill", fillColor)
        .style("stroke-width", 2)
        .attr("class", `${colorsWithinBuffer.map(c => `color${c}`).join(" ")}`)
    }
    for (let i = 0; i < heart; i += 2 * horColorStep) {
      addColor(1, 1)
    }
    for (let lung = 1 - heart; lung >= -horColorStep; lung -= horColorStep) {
      addColor(heart, lung)
    }
  }

  textContainer.append("span").style("margin", "0 170px 0 0").text("Lung")
  textContainer.append("span").text("Liver")

  // Add a legend for shapes
  textContainer.append("br")
  textContainer.append("svg")
    .attr("width", 12)
    .attr("height", 12)
    .attr("viewBox", "-8 -10 16 16")
    .append("path")
      .attr("d", dataToSymbol({node_type: 'peak'})())
      .attr("fill", "black")
  textContainer.append("span").style("margin", "0 15px 0 5px").text("Peak")
  textContainer.append("svg")
    .attr("width", 12)
    .attr("height", 12)
    .attr("viewBox", "-6 -6 12 12")
    .append("path")
      .attr("d", dataToSymbol({node_type: 'gene', tf_coding: false})())
      .attr("fill", "black")
  textContainer.append("span").style("margin", "0 15px 0 5px").text("Gene")
  textContainer.append("svg")
    .attr("width", 12)
    .attr("height", 12)
    .attr("viewBox", "-6 -6 12 12")
    .append("path")
      .attr("d", dataToSymbol({node_type: 'gene', tf_coding: true})())
      .attr("fill", "black")
  textContainer.append("span").style("margin", "0 15px 0 5px").text("TF")

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

  // Allow for window resizing
  window.addEventListener("resize", () => {
    width = window.innerWidth;
    height = window.innerHeight;
    outerSvg
      .attr("width", width)
      .attr("height", height);
    svg.attr("transform", transform(currTransform))
  });

  res.then(data => {
    // The force simulation mutates links and nodes, so create a copy
    // so that re-evaluating this cell produces the same result.
    const allNodes = data.nodes.map(node => ({...node, id: String(node.id)})).filter(node => 
      (node.Heart.condition_sensitivity + node.Liver.condition_sensitivity + node.Lung.condition_sensitivity) > 0.05 || 
      (node.node_type === 'gene' && node.tf_coding && (node.Heart.condition_sensitivity + node.Liver.condition_sensitivity + node.Lung.condition_sensitivity) > 0.03));
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
        .force("link", d3.forceLink(links).id(d => d.id).distance(_ => 20))
        .force("charge", d3.forceManyBody().strength(node => node.node_type === 'gene' && node.tf_coding ? -80 : -30))
        .force("x", d3.forceX())
        .force("y", d3.forceY());

    // Add links
    const link = svg.append("g")
        .attr("stroke", "#999")
        .attr("stroke-opacity", 0.6)
      .selectAll("line")
      .data(links)
      .join("line")
        .attr("id", d => `index${d.index}`)
        .attr("stroke-width", 1)
        .attr("marker-end", "url(#arrow)");

    // Add nodes
    const node = svg.append("g")
        .attr("stroke", "#000")
        .attr("stroke-width", 0.5)
      .selectAll("path")
      .data(nodes)
      .join("path")
        .attr("cursor", "pointer")
        .attr("id", d => `id${d.id}`)
        .attr("fill", d => dataToColor(d))
        .attr("d", d => dataToSymbol(d).size(100)());

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
    node.on("click", (evt, d) => {
      evt.stopPropagation()
      selectNode(d.id)
    })

    // On click edge, if one node is selected, select the other node
    let selected = '';
    link.on("click", (_evt, d) => {
      if (selected === d.source.id) { 
        selectNode(d.target.id) 
      } else if (selected === d.target.id) {
        selectNode(d.source.id)
      }
    })

    outerSvg.on("click", () => deselectNode())

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
      colorKey.selectAll("rect").style("stroke", "none")
      node.attr("stroke-width", 0.5)
      node.attr("opacity", 0.1)
      link.attr("opacity", 0.1)
      const d = nodes.find(n => n.id === nodeId)
      const outEdges = links.filter(e => e.source.id === nodeId)
      const inEdges = links.filter(e => e.target.id === nodeId)
      const selectedEdges = [...outEdges, ...inEdges]
      if (d.node_type === 'gene') {
        // Map back to closest TF
        selectedEdges.push(...inEdges.flatMap(inE => links.filter(e => e.target.id === inE.source.id)))
      }
      selected = nodeId
      if (!d) {
        deselectNode()
        displayText.html(`Gene/Peak '${nodeId}' not found. ${displayTextBreaks}`)
        scoreText.html("")
        return
      }
      fetchWithRetry('/plot', 
        { method: 'POST', body: JSON.stringify(d), headers: { 'Content-Type': 'application/json'}}
      ).then((res) => {
        scoreText.html(`<image src="${res.path}" style="width: 100%; min-height: 180px">`)
      }).catch((_err) => {
        scoreText.html(`Condition sensitivity: ${['Heart', 'Lung', 'Liver'].map(organ => d[organ].condition_sensitivity.toFixed(2)).join(" ")}<br/>
        Organ specificity: ${['Heart', 'Lung', 'Liver'].map(organ => d[organ].organ_specificity.toFixed(2)).join(" ")}<br/>
        Universality: ${d.universality.toFixed(5)}`);
      })
      svg.select(`#id${d.id}`).attr("stroke-width", 1.5)
      const selectedColor = svg.select(`#id${d.id}`).attr("fill").substring(1);
      colorKey.select(`.color${selectedColor}`).style("stroke", "black")
      selectedEdges.forEach(e => {
        svg.select(`#index${e.index}`).attr("opacity", 1)
        svg.select(`#id${e.source.id}`).attr("opacity", 1)
        svg.select(`#id${e.target.id}`).attr("opacity", 1)
      })

      displayText.html(`${d.node_type === 'peak' ? 'Peak' : 'Gene'}: ${d.id}<br/>
        ${d.chr}:${d.start}-${d.end}<br/>
        ${d.node_type === 'gene' ? `<a target="_blank" href="https://www.genecards.org/cgi-bin/carddisp.pl?gene=${d.id}">View on GeneCards</a>` : 
          `<a target="_blank" href="https://screen.wenglab.org/search?assembly=GRCh38&chromosome=${d.chr}&start=${d.start}&end=${d.end}">View in SCREEN</a>`}`)

      const endTransform = [d.x, d.y, width / 5];
      const i = d3.interpolateZoom(currTransform, endTransform)

      svg.transition()
        .duration(i.duration)
        .attrTween("transform", () => t => transform(currTransform = i(t)));
      d3.zoomTransform(outerSvg, endTransform)
    }

    function deselectNode() {
      node.attr("stroke-width", 0.5)
      node.attr("opacity", 1)
      link.attr("opacity", 1)
      displayText.html(emptyDisplayText)
      colorKey.selectAll("rect").style("stroke", "none")
    }

    // Add search listener
    searchForm.on("submit", (evt) => {
      evt.preventDefault();
      selectNode(searchText.property("value"));
      return false;
    });

    // remove loading text
    loadingText.remove();
  }).catch(err => window.alert(err));

  return container.node();
}

export { build };