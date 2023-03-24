let main
let tris
let heightMap
let seaDrop
let jagged
let numPoints = 4096
let points = []
let delaunay
let voronoi
let vorVertices // array of vertices of voronoi cells,[[x1, y1], [x2, y2], [x3, y3]...]
let vorVertexToIndex // dictionary that maps indices of voronoi verticies to its coordinates
let vorAdjacent // array of indices that correspond with other voronoi vertices adjacent to a given vertex (can be connected by a edge on a voronoi cell)
let vorEdges // array of edges in the map, where each edge is represented as an array [v0, v1, leftC, rightC] where v0 and v1 are the indices of the verticies on either end of the edge, and leftC and rightC are the indices of input points associated with voronoi cells on either side of the edge
let vorTriangles // array for indexes of vertices of the delaunay triangle that each vertex in vorVertices can be found in [[x1, y1, z1], [x2, y2, z2], ...]
let heights
let heightMax = 0
let numLloydRelaxation = 2
let bufferSize = 15
let downhill
let waterFlow
let slope
let erosionRates
let riverSegments

function setup(){
  createCanvas(2200, 400)
  main = createGraphics(400, 400)
  tris = createGraphics(main.width, main.height)
  heightMap = createGraphics(main.width, main.height)
  seaDrop = createGraphics(main.width, main.height)
  jagged = createGraphics(main.width, main.height)
  colorMode(RGB)
  noStroke()
  noLoop()
}

function draw(){
  makeRandomPoints(numPoints)
  delaunay = d3.Delaunay.from(points)
  voronoi = delaunay.voronoi([-bufferSize, -bufferSize, main.width+bufferSize, main.height+bufferSize])
  drawRandomPoints(points)
  for (let i=0; i<numLloydRelaxation; i++){
    main.clear()
    points = computeCentroids(voronoi)
    delaunay = d3.Delaunay.from(points)
    voronoi = delaunay.voronoi([-bufferSize, -bufferSize, main.width+bufferSize, main.height+bufferSize])
    drawVoronoiCells(voronoi)
    drawDelaunayTriangles(delaunay)
    drawRandomPoints(points)
  }
  makeMapRepresentation()
  tris.copy(main, 0, 0, main.width, main.height, 0, 0, main.width, main.height)
  image(tris, 1800, 0)
  makeHills(200, 5) // num of hills, radius of hill
  addRandomSlope()
  colorHeights()
  heightMap.copy(main, 0, 0, main.width, main.height, 0, 0, main.width, main.height)
  image(heightMap, 1350, 0)
  setSeaLevelToMedian()
  colorHeights()
  seaDrop.copy(main, 0, 0, main.width, main.height, 0, 0, main.width, main.height)
  image(seaDrop, 900, 0)
  // erode for the first time
  fillDepressions()
  getDownhill()
  getWaterFlow()
  getSlope()
  getErosion()
  erode(100, 5) // intensity of erosion, num times erosion happens
  setSeaLevelToMedian()
  fillDepressions()

  colorHeights()
  jagged.copy(main, 0, 0, main.width, main.height, 0, 0, main.width, main.height)
  image(jagged, 450, 0)
  cleanCoastlines(5) // num times cleaning happens

  colorHeights()
  drawRivers(d3.max(waterFlow)) // water flow required at a certain vor vertex before it is rendered
  image(main, 0, 0)
  drawSeaLevel()
}

function mergeSegments(segments){
  // segments is a bunch of pairs of x,y coords [[a, b], [c, d]] is one element
  let done = []
  let paths = []
  let workingPath = null
  while (true){
    if (workingPath == null){
      for (let i=0; i<segments.length; i++){
        if (done[i] == true){
          continue
        }
        workingPath = [segments[i][0], segments[i][1]]
        done[i] = true
        break
      }
      if (workingPath == null){
        break
      }
    }
    let changed = false
    for (let j=0; j<segments.length; j++){
      if (done[j] == true){
        continue
      }
      if (segments[j][0] == workingPath[0]){
        workingPath.unshift(segments[j][1])
      }else if (segments[j][1] == workingPath[0]){
        workingPath.unshift(segments[j][0])
      }else if (segments[j][0] == workingPath[workingPath.length-1]){
        workingPath.push(segments[j][1])
      }else if (segments[j][1] == workingPath[workingPath.length-1]){
        workingPath.push(segments[j][0])
      }else{
        continue
      }
      done[j] = true
      changed = true
      break
    }
    if (changed == false){
      paths.push(workingPath)
      workingPath = null
    }
  }
  console.log(paths)
  return paths
}

function drawRivers(waterFlowRequired){
  console.log(`the water flow required is ${waterFlowRequired}`)
  console.log("this runs")
  riverSegments = []
  let above = 0
  for (let i=0; i<heights.length; i++){
    if (heights[i] > 0){
      above++
    }
  }
  waterFlowRequired = (above/heights.length) * waterFlowRequired
  console.log(`the water flow required is ${waterFlowRequired}`)
  for (let i=0; i<vorVertices.length; i++){
    let vx =  vorVertices[i][0]
    let vy =  vorVertices[i][1]
    if (vx < main.width*0.05 || vx > main.width*0.95 || vy < main.height*0.05 || vy > main.height*0.95){
      continue
    }
    if ((waterFlow[i] > waterFlowRequired) && (heights[i] > 0) && (heights[downhill[i]] >= 0)){
      let up = vorVertices[i]
      let down = vorVertices[downhill[i]]
      riverSegments.push([up, down]) // pushes a pair of x,y coords [[a, b], [c, d]]
    }
  }
  console.log(riverSegments)

  /*for (let i=0; i<riverSegments.length; i++){
    stroke("blue")
    strokeWeight(2)
    line(riverSegments[i][0][0], riverSegments[i][0][1], riverSegments[i][1][0], riverSegments[i][1][1])
  }*/

  riverSegments = mergeSegments(riverSegments)
  console.log("river segments")
  console.log(riverSegments)

  for (let i=0; i<riverSegments.length; i++){
    if (riverSegments[i].length >= 15){
      for (let j=0; j<riverSegments[i].length-1; j++){
        main.stroke("blue")
        main.strokeWeight(2)
        main.line(riverSegments[i][j][0], riverSegments[i][j][1], riverSegments[i][j+1][0], riverSegments[i][j+1][1])
      }
    }
  }
}

function cleanCoastlines(times){
  for (let i=0; i<times; i++){
    let newHeights = [...heights]
    for (let j=0; j<heights.length; j++){
      let nWater = 0
      let nLand = 0
      let maxWaterHeight = -99999
      let minLandHeight = 99999
      let vheight = heights[j]
      let neighbors = vorAdjacent[j]
      if (neighbors.length < 3){
        continue
      }
      for (let k=0; k<neighbors.length; k++){
        if (heights[neighbors[k]] < 0){
          nWater += 1
          if (maxWaterHeight < heights[neighbors[k]]){
            maxWaterHeight = heights[neighbors[k]]
          }
        }else if (heights[neighbors[k]] > 0){
          nLand += 1
          if (minLandHeight > heights[neighbors[k]]){
            minLandHeight = heights[neighbors[k]]
          }
        }
      }
      if ((vheight > 0) && (nLand <= 1)){
        newHeights[j] = (heights[j]+maxWaterHeight)/2
      }
      if ((vheight < 0) && (nWater <= 1)){
        newHeights[j] = (heights[j]+minLandHeight)/2
      }
    }
    heights = newHeights
  }
}

function erode(intensity, times){
  fillDepressions()
  for (let i=0; i<times; i++){
    getErosion()
    let maxErosion = d3.max(erosionRates)
    // console.log(maxErosion)
    for (let j=0; j<heights.length; j++){
      heights[j] -= intensity * (erosionRates[j]/maxErosion)
    }
    fillDepressions()
  }
}

function getErosion(){
  //getDownhill()
  //getWaterFlow()
  getSlope()
  erosionRates = []
  let riverErosion, soilErosion, totalErosion
  for (let i=0; i<heights.length; i++){
    riverErosion = Math.sqrt(waterFlow[i]) * slope[i]
    soilErosion = slope[i] * slope[i]
    totalErosion = 2000 * riverErosion + soilErosion
    totalErosion = totalErosion > 200 ? 200 : totalErosion
    erosionRates.push(totalErosion)
  }
}

function getSlope(){
  slope = []
  for (let i=0; i<heights.length; i++){
    let neighbors = vorAdjacent[i]
    if (neighbors.length < 3){
      // if it's near the edge just assume it's slope is 0 (not fully correct, but makes it easier) -- and you also don't need to erode the edge vor indexes as much
      slope.push(0)
      continue
    }
    let p0x = vorVertices[neighbors[0]][0]
    let p0y = vorVertices[neighbors[0]][1]
    let p1x = vorVertices[neighbors[1]][0]
    let p1y = vorVertices[neighbors[1]][1]
    let p2x = vorVertices[neighbors[2]][0]
    let p2y = vorVertices[neighbors[2]][1]

    let x1 = p1x - p0x
    let x2 = p2x - p0x
    let y1 = p1y - p0y
    let y2 = p2y - p0y

    // get the determinant to figure out how much the slope should scale (the determinant tells you how much a linear transformation stretches / squishes shapes)
    let det = (x1*y2) - (x2*y1)

    let h1 = heights[neighbors[1]] - heights[neighbors[0]]
    let h2 = heights[neighbors[2]] - heights[neighbors[0]]

    // calculate the final slope vector at a point. you do this by taking the determinant (used as a simplified 2d cross product) to get a vector that's perpendicular to the height vector and direction vector
    if (det == 0){
      slope.push(0)
      continue
    }
    let xComponent = (y2 * h1 - y1 * h2) / det
    let yComponent = (-x2 * h1 + x1 * h2) / det

    slope.push(Math.sqrt(xComponent*xComponent + yComponent*yComponent))
  }
}

function getWaterFlow(){
  waterFlow = []
  let indices = []
  for (let i=0; i<heights.length; i++){
    waterFlow.push(1/heights.length)
    indices.push(i)
  }
  indices.sort((a, b) => heights[a] - heights[b]) // sort the indices by height. it takes two elemets from indices, a and b, and calculates diff in height. if height[a] - height[b] is positive, a comes in front of b
  for (let i=0; i<indices.length; i++){
    let waterToAdd = waterFlow[indices[i]]
    let downhillIndex = downhill[indices[i]]
    if (downhillIndex >= 0){
      waterFlow[downhillIndex] += waterToAdd
    }
  }
  //console.log(waterFlow)
}

function getDownhill(){
  downhill = []
  for (let i=0; i<heights.length; i++){
    let neighbors = vorAdjacent[i]
    if (neighbors.length < 3){
      downhill.push(-2)
      continue
    }
    let minHeight = heights[i]
    let minIndex
    for (let j=0; j<neighbors.length; j++){
      if (heights[neighbors[j]] < minHeight){
        minHeight = heights[neighbors[j]]
        minIndex = neighbors[j]
      }
    }
    if (minHeight == heights[i]){
      downhill.push(-1)
    }else{
      downhill.push(minIndex)
    }
  }
}

function fillDepressions(){
  let newHeights = []
  for (let i=0; i<heights.length; i++){
    if (vorAdjacent[i].length < 3){
      newHeights.push(heights[i])
    }else{
      newHeights.push(Infinity)
    }
  }
  let changed = true
  while (changed == true){
    changed = false
    for (let i=0; i<newHeights.length; i++){
      if (newHeights[i] == heights[i]){
        continue
      }
      let neighbors = vorAdjacent[i] // array of indices of adjacent voronoi vertices
      for (let j=0; j<neighbors.length; j++){
        let neighborIndex = neighbors[j]
        if (newHeights[neighborIndex]+Number.EPSILON< newHeights[i]){
          newHeights[i] = Math.max(heights[i], newHeights[neighborIndex]+Number.EPSILON)
          changed = true
        }
      }
    }
  }
  heights = newHeights.slice()
}

function setSeaLevelToMedian(){
  let median = d3.median(heights)
  for (let i=0; i<heights.length; i++){
    heights[i] = heights[i] - median
  }
}

function drawSeaLevel(){
  for (let i=0; i<vorEdges.length; i++){
    let startI = vorEdges[i][0]
    let endI = vorEdges[i][1]
    if ((heights[startI]>=0 && heights[endI]<0) || (heights[startI]<0 && heights[endI>=0])){
      //console.log("this runs!")
      let leftPoint = vorEdges[i][2]
      let rightPoint = vorEdges[i][3]
      let leftPx, leftPy, rightPx, rightPy
      leftPx = delaunay.points[leftPoint*2]
      leftPy = delaunay.points[leftPoint*2+1]
      rightPx = delaunay.points[rightPoint*2]
      rightPy = delaunay.points[rightPoint*2+1]
      if (leftPx >= 0 && leftPx <= main.width && leftPy >= 0 && leftPy <= main.height && rightPx >= 0 && rightPx <= main.width && rightPy >= 0 && rightPy <= main.height){
        stroke("black")
        strokeWeight(2)
        line(leftPx, leftPy, rightPx, rightPy)
      }
    }
  }
}

function addRandomSlope(){
  let dir = map(Math.random(), 0, 1, 0, 360)
  let slopeAngle = map(Math.random(), 0, 1, 30, 60)
  if (Math.random() < .5){
    slopeAngle = slopeAngle*-1
  }
  let slopeMultiplier = 2
  let centerX = main.width / 2
  let centerY = main.height / 2
  let dirRadians = radians(dir)
  let dirVector = {x: cos(dirRadians), y: sin(dirRadians)}

  for (let i=0; i<vorVertices.length; i++){
    let heightInc
    let vx = vorVertices[i][0]
    let vy = vorVertices[i][1]
    let dx = vx - centerX
    let dy = vy - centerY
    let distance = dist(centerX, centerY, vx, vy)
    let centerToVertexVector = {x: dx / distance, y: dy / distance}
    let dotProduct = dirVector.x * centerToVertexVector.x + dirVector.y * centerToVertexVector.y
    heightInc = slopeMultiplier * tan(radians(slopeAngle)) * distance * dotProduct
    heights[i] += heightInc
  }

  /*for (let i=0; i<vorVertices.length; i++){
    let heightInc
    let vx = vorVertices[i][0]
    let vy = vorVertices[i][1]
    let dx = vx - centerX
    let dy = vy - centerY
    let angle = atan(dy/dx) // angle from center to vertex
    if (degrees(angle)-90 < dir && degrees(angle)+90 > dir){
      heightInc = tan(slopeAngle) * dx
    }else{
      heightInc = -tan(slopeAngle) * dx
    }
    heights[i] += heightInc
  }*/
}

function colorHeights(){
  heightMax = Math.max(...heights)
  for (let i=0; i<vorTriangles.length; i++){
    //console.log(heightMax)
    let triIndex1 = vorTriangles[i][0]
    let triIndex2 = vorTriangles[i][1]
    let triIndex3 = vorTriangles[i][2]
    let t1x = delaunay.points[triIndex1*2]
    let t1y = delaunay.points[triIndex1*2+1]
    let t2x = delaunay.points[triIndex2*2]
    let t2y = delaunay.points[triIndex2*2+1]
    let t3x = delaunay.points[triIndex3*2]
    let t3y = delaunay.points[triIndex3*2+1]
    if (t1x >= 0 && t1x <= main.width && t1y >= 0 && t1y <= main.height && t2x >= 0 && t2x <= main.width && t2y >= 0 && t2y <= main.height && t3x >= 0 && t3x <= main.width && t3y >= 0 && t3y <= main.height){
      let normalizedHeight = heights[i]/300
      if (heights[i] == 0){
        main.fill(0, 128, 128)
        main.triangle(t1x, t1y, t2x, t2y, t3x, t3y)
      }
      if (heights[i]>0){
        let r = Math.min(255, 255*normalizedHeight)
        let g = Math.min(255, 128+(128*normalizedHeight))
        let b = Math.max(0, 128-(128*normalizedHeight))
        //console.log(r, g, b)
        main.fill(r, g, b)
        main.triangle(t1x, t1y, t2x, t2y, t3x, t3y)
      }
      if (heights[i]<0){
        let h = -1*normalizedHeight
        let r = Math.min(100, 100*h)
        let g = Math.max(50, 50-(50*h))
        let b = Math.min(100, 100+(100*h))
        main.fill(r, g, b)
        main.triangle(t1x, t1y, t2x, t2y, t3x, t3y)
      }
    }
  }
}

function makeHills(n, r){ // vorVertices, num hills, radius of each hill
  let hillLoc = []
  let x, y
  for (let i=0; i<n; i++){
    x = map(Math.random(), 0, 1, 0, main.width)
    y = map(Math.random(), 0, 1, 0, main.height)
    hillLoc.push([x, y])
  }
  for (let i=0; i<vorVertices.length; i++){
    for (let j=0; j<hillLoc.length; j++){
      let vx = vorVertices[i][0]
      let vy = vorVertices[i][1]
      let hx = hillLoc[j][0]
      let hy = hillLoc[j][1]
      let d = dist(vx, vy, hx, hy)
      //console.log(d)
      let heightContribution = Math.pow(Math.exp(-d/(2*r*r)), 4) * 100
      heights[i] += heightContribution
    }
  }
}

function makeMapRepresentation(){
  vorVertexToIndex = {}
  vorVertices = []
  vorAdjacent = []
  vorEdges = []
  vorTriangles = []
  heights = []

  for (let i=0; i<voronoi.circumcenters.length; i+=2){
    let x = voronoi.circumcenters[i]
    let y = voronoi.circumcenters[i+1]
    let vertex = [x,y]
    vorVertices.push(vertex)
    vorVertexToIndex[`${x},${y}`] = i/2
    heights.push(0)
  }

  for (let i=0; i<vorVertices.length; i++){
    vorAdjacent.push([])
    vorTriangles.push([])
  }

  for (let i = 0; i < delaunay.halfedges.length; i++) {
    let j = delaunay.halfedges[i]
    if (j < i) continue
    let triangle1index = Math.floor(i/3)
    let triangle2index = Math.floor(j/3)
    let vx1 = voronoi.circumcenters[triangle1index*2]
    let vy1 = voronoi.circumcenters[triangle1index*2+1]
    let vx2 = voronoi.circumcenters[triangle2index*2]
    let vy2 = voronoi.circumcenters[triangle2index*2+1]
    let index1 = vorVertexToIndex[`${vx1},${vy1}`]
    let index2 = vorVertexToIndex[`${vx2},${vy2}`]
    let leftC = delaunay.triangles[i]
    let rightC = delaunay.triangles[j]
    if (index1 != undefined && index2 != undefined){
      vorEdges.push([index1, index2, leftC, rightC])
      vorEdges.push([index2, index1, rightC, leftC])
      vorAdjacent[index1].push(index2)
      vorAdjacent[index2].push(index1)
    }
    if (index1 != undefined && vorTriangles[index1].length == 0) {
      let t1 = delaunay.triangles[triangle1index*3]
      let t2 = delaunay.triangles[triangle1index*3+1]
      let t3 = delaunay.triangles[triangle1index*3+2]
      vorTriangles[index1] = [t1, t2, t3]
    }
    if (index2 != undefined && vorTriangles[index2].length == 0) {
      let t1 = delaunay.triangles[triangle2index*3]
      let t2 = delaunay.triangles[triangle2index*3+1]
      let t3 = delaunay.triangles[triangle2index*3+2]
      vorTriangles[index2] = [t1, t2, t3]
    }
  }
}

function computeCentroids(v){
  let centroids = []
  for (let polygon of v.cellPolygons()){
    let A = 0
    let Cx = 0
    let Cy = 0
    let xi, xn, yi, yn
    for (let i=0; i<polygon.length; i++){
      xi = polygon[i][0]
      yi = polygon[i][1]
      xn = polygon[(i+1)%polygon.length][0]
      yn = polygon[(i+1)%polygon.length][1]
      Cx += (xi+xn)*((xi*yn)-(xn*yi))
      Cy += (yi+yn)*((xi*yn)-(xn*yi))
      A += ((xi*yn)-(xn*yi))
    }
    A = A/2
    Cx = (1/(6*A))*Cx
    Cy = (1/(6*A))*Cy
    centroids.push([Cx, Cy])
  }
  return centroids
}

function drawVoronoiCells(v){
  for (let polygon of v.cellPolygons()){
    main.noFill()
    main.stroke("blue")
    main.strokeWeight(.5)
    main.beginShape()
    for (let [x, y] of polygon){
      main.vertex(x, y)
    }
    main.endShape(CLOSE)
  }
}

function drawDelaunayTriangles(d){
  let v1, v2, v3
  let x1, y1, x2, y2, x3, y3
  for (let i=0; i<d.triangles.length; i+=3){
    v1 = d.triangles[i]
    v2 = d.triangles[i+1]
    v3 = d.triangles[i+2]
    x1 = d.points[v1*2]
    y1 = d.points[v1*2+1]
    x2 = d.points[v2*2]
    y2 = d.points[v2*2+1]
    x3 = d.points[v3*2]
    y3 = d.points[v3*2+1]
    if (x1 >= 0 && x1 <= main.width && y1 >= 0 && y1 <= main.height && x2 >= 0 && x2 <= main.width && y2 >= 0 && y2 <= main.height && x3 >= 0 && x3 <= main.width && y3 >= 0 && y3 <= main.height){
      main.noFill()
      main.stroke("black")
      main.strokeWeight(.5)
      main.triangle(x1, y1, x2, y2, x3, y3)
    }
  }
}

function drawRandomPoints(p){
  let x, y
  for (let i=0; i<p.length; i++){
    x = p[i][0]
    y = p[i][1]
    main.stroke("white")
    main.strokeWeight(.5)
    main.point(x, y)
  }
}

function makeRandomPoints(num){
  let x, y
  for (let i=0; i<num; i++){
    x = map(Math.random(), 0, 1, -bufferSize, main.width+bufferSize)
    y = map(Math.random(), 0, 1, -bufferSize, main.height+bufferSize)
    points.push([x, y])
  }
}