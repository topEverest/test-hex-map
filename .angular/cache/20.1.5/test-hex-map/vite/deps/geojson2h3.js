import {
  __commonJS,
  __esm,
  __export,
  __toCommonJS
} from "./chunk-DZYXDVEG.js";

// node_modules/geojson2h3/node_modules/h3-js/dist/browser/h3-js.es.js
var h3_js_es_exports = {};
__export(h3_js_es_exports, {
  UNITS: () => UNITS,
  cellArea: () => cellArea,
  compact: () => compact,
  degsToRads: () => degsToRads,
  edgeLength: () => edgeLength,
  exactEdgeLength: () => exactEdgeLength,
  experimentalH3ToLocalIj: () => experimentalH3ToLocalIj,
  experimentalLocalIjToH3: () => experimentalLocalIjToH3,
  geoToH3: () => geoToH3,
  getDestinationH3IndexFromUnidirectionalEdge: () => getDestinationH3IndexFromUnidirectionalEdge,
  getH3IndexesFromUnidirectionalEdge: () => getH3IndexesFromUnidirectionalEdge,
  getH3UnidirectionalEdge: () => getH3UnidirectionalEdge,
  getH3UnidirectionalEdgeBoundary: () => getH3UnidirectionalEdgeBoundary,
  getH3UnidirectionalEdgesFromHexagon: () => getH3UnidirectionalEdgesFromHexagon,
  getOriginH3IndexFromUnidirectionalEdge: () => getOriginH3IndexFromUnidirectionalEdge,
  getPentagonIndexes: () => getPentagonIndexes,
  getRes0Indexes: () => getRes0Indexes,
  h3Distance: () => h3Distance,
  h3GetBaseCell: () => h3GetBaseCell,
  h3GetFaces: () => h3GetFaces,
  h3GetResolution: () => h3GetResolution,
  h3IndexToSplitLong: () => h3IndexToSplitLong,
  h3IndexesAreNeighbors: () => h3IndexesAreNeighbors,
  h3IsPentagon: () => h3IsPentagon,
  h3IsResClassIII: () => h3IsResClassIII,
  h3IsValid: () => h3IsValid,
  h3Line: () => h3Line,
  h3SetToMultiPolygon: () => h3SetToMultiPolygon,
  h3ToCenterChild: () => h3ToCenterChild,
  h3ToChildren: () => h3ToChildren,
  h3ToGeo: () => h3ToGeo,
  h3ToGeoBoundary: () => h3ToGeoBoundary,
  h3ToParent: () => h3ToParent,
  h3UnidirectionalEdgeIsValid: () => h3UnidirectionalEdgeIsValid,
  hexArea: () => hexArea,
  hexRing: () => hexRing,
  kRing: () => kRing,
  kRingDistances: () => kRingDistances,
  numHexagons: () => numHexagons,
  pointDist: () => pointDist,
  polyfill: () => polyfill,
  radsToDegs: () => radsToDegs,
  splitLongToh3Index: () => splitLongToh3Index,
  uncompact: () => uncompact
});
function validateRes(res) {
  if (typeof res !== "number" || res < 0 || res > 15 || Math.floor(res) !== res) {
    throw new Error("Invalid resolution: " + res);
  }
}
function h3IndexToSplitLong(h3Index) {
  if (Array.isArray(h3Index) && h3Index.length === 2 && Number.isInteger(h3Index[0]) && Number.isInteger(h3Index[1])) {
    return h3Index;
  }
  if (typeof h3Index !== "string" || INVALID_HEXIDECIMAL_CHAR.test(h3Index)) {
    return [0, 0];
  }
  var upper = parseInt(h3Index.substring(0, h3Index.length - 8), BASE_16);
  var lower = parseInt(h3Index.substring(h3Index.length - 8), BASE_16);
  return [lower, upper];
}
function hexFrom32Bit(num) {
  if (num >= 0) {
    return num.toString(BASE_16);
  }
  num = num & 2147483647;
  var tempStr = zeroPad(8, num.toString(BASE_16));
  var topNum = (parseInt(tempStr[0], BASE_16) + 8).toString(BASE_16);
  tempStr = topNum + tempStr.substring(1);
  return tempStr;
}
function splitLongToh3Index(lower, upper) {
  return hexFrom32Bit(upper) + zeroPad(8, hexFrom32Bit(lower));
}
function zeroPad(fullLen, numStr) {
  var numZeroes = fullLen - numStr.length;
  var outStr = "";
  for (var i = 0; i < numZeroes; i++) {
    outStr += "0";
  }
  outStr = outStr + numStr;
  return outStr;
}
function polygonArrayToGeofence(polygonArray, geofence, isGeoJson) {
  var numVerts = polygonArray.length;
  var geoCoordArray = libh3._calloc(numVerts, SZ_GEOCOORD);
  var latIndex = isGeoJson ? 1 : 0;
  var lngIndex = isGeoJson ? 0 : 1;
  for (var i = 0; i < numVerts * 2; i += 2) {
    libh3.HEAPF64.set([polygonArray[i / 2][latIndex], polygonArray[i / 2][lngIndex]].map(degsToRads), geoCoordArray / SZ_DBL + i);
  }
  libh3.HEAPU32.set([numVerts, geoCoordArray], geofence / SZ_INT);
  return geofence;
}
function coordinatesToGeoPolygon(coordinates, isGeoJson) {
  var numHoles = coordinates.length - 1;
  var geoPolygon = libh3._calloc(SZ_GEOPOLYGON);
  var geofenceOffset = 0;
  var numHolesOffset = geofenceOffset + SZ_GEOFENCE;
  var holesOffset = numHolesOffset + SZ_INT;
  polygonArrayToGeofence(coordinates[0], geoPolygon + geofenceOffset, isGeoJson);
  var holes;
  if (numHoles > 0) {
    holes = libh3._calloc(numHoles, SZ_GEOFENCE);
    for (var i = 0; i < numHoles; i++) {
      polygonArrayToGeofence(coordinates[i + 1], holes + SZ_GEOFENCE * i, isGeoJson);
    }
  }
  libh3.setValue(geoPolygon + numHolesOffset, numHoles, "i32");
  libh3.setValue(geoPolygon + holesOffset, holes, "i32");
  return geoPolygon;
}
function destroyGeoPolygon(geoPolygon) {
  var geofenceOffset = 0;
  var numHolesOffset = geofenceOffset + SZ_GEOFENCE;
  var holesOffset = numHolesOffset + SZ_INT;
  var geofenceArrayOffset = SZ_INT;
  libh3._free(libh3.getValue(geoPolygon + geofenceOffset + geofenceArrayOffset, "i8*"));
  var numHoles = libh3.getValue(geoPolygon + numHolesOffset, "i32");
  if (numHoles > 0) {
    var holes = libh3.getValue(geoPolygon + holesOffset, "i32");
    for (var i = 0; i < numHoles; i++) {
      libh3._free(libh3.getValue(holes + SZ_GEOFENCE * i + geofenceArrayOffset, "i8*"));
    }
    libh3._free(holes);
  }
  libh3._free(geoPolygon);
}
function readLong(invocation) {
  var upper = libh3.getTempRet0();
  return [invocation, upper];
}
function readH3Index(invocation) {
  var ref = readLong(invocation);
  var lower = ref[0];
  var upper = ref[1];
  return upper ? splitLongToh3Index(lower, upper) : null;
}
function readH3IndexFromPointer(cAddress, offset) {
  if (offset === void 0) offset = 0;
  var lower = libh3.getValue(cAddress + SZ_INT * offset * 2, "i32");
  var upper = libh3.getValue(cAddress + SZ_INT * (offset * 2 + 1), "i32");
  return upper ? splitLongToh3Index(lower, upper) : null;
}
function storeH3Index(h3Index, cAddress, offset) {
  libh3.HEAPU32.set(h3IndexToSplitLong(h3Index), cAddress / SZ_INT + 2 * offset);
}
function readArrayOfHexagons(cAddress, maxCount) {
  var out = [];
  for (var i = 0; i < maxCount; i++) {
    var h3Index = readH3IndexFromPointer(cAddress, i);
    if (h3Index !== null) {
      out.push(h3Index);
    }
  }
  return out;
}
function storeArrayOfHexagons(cAddress, hexagons) {
  var count = hexagons.length;
  for (var i = 0; i < count; i++) {
    storeH3Index(hexagons[i], cAddress, i);
  }
}
function storeGeoCoord(lat, lng) {
  var geoCoord = libh3._calloc(1, SZ_GEOCOORD);
  libh3.HEAPF64.set([lat, lng].map(degsToRads), geoCoord / SZ_DBL);
  return geoCoord;
}
function readSingleCoord(cAddress) {
  return radsToDegs(libh3.getValue(cAddress, "double"));
}
function readGeoCoord(cAddress) {
  return [readSingleCoord(cAddress), readSingleCoord(cAddress + SZ_DBL)];
}
function readGeoCoordGeoJson(cAddress) {
  return [readSingleCoord(cAddress + SZ_DBL), readSingleCoord(cAddress)];
}
function readGeoBoundary(geoBoundary, geoJsonCoords, closedLoop) {
  var numVerts = libh3.getValue(geoBoundary, "i32");
  var vertsPos = geoBoundary + SZ_DBL;
  var out = [];
  var readCoord = geoJsonCoords ? readGeoCoordGeoJson : readGeoCoord;
  for (var i = 0; i < numVerts * 2; i += 2) {
    out.push(readCoord(vertsPos + SZ_DBL * i));
  }
  if (closedLoop) {
    out.push(out[0]);
  }
  return out;
}
function readMultiPolygon(polygon, formatAsGeoJson) {
  var output = [];
  var readCoord = formatAsGeoJson ? readGeoCoordGeoJson : readGeoCoord;
  var loops;
  var loop;
  var coords;
  var coord;
  while (polygon) {
    output.push(loops = []);
    loop = libh3.getValue(polygon, "i8*");
    while (loop) {
      loops.push(coords = []);
      coord = libh3.getValue(loop, "i8*");
      while (coord) {
        coords.push(readCoord(coord));
        coord = libh3.getValue(coord + SZ_DBL * 2, "i8*");
      }
      if (formatAsGeoJson) {
        coords.push(coords[0]);
      }
      loop = libh3.getValue(loop + SZ_PTR * 2, "i8*");
    }
    polygon = libh3.getValue(polygon + SZ_PTR * 2, "i8*");
  }
  return output;
}
function readCoordIJ(cAddress) {
  return {
    i: libh3.getValue(cAddress, "i32"),
    j: libh3.getValue(cAddress + SZ_INT, "i32")
  };
}
function storeCoordIJ(cAddress, ref) {
  var i = ref.i;
  var j = ref.j;
  libh3.setValue(cAddress, i, "i32");
  libh3.setValue(cAddress + SZ_INT, j, "i32");
}
function readArrayOfPositiveIntegers(cAddress, count) {
  var out = [];
  for (var i = 0; i < count; i++) {
    var int = libh3.getValue(cAddress + SZ_INT * i, "i32");
    if (int >= 0) {
      out.push(int);
    }
  }
  return out;
}
function h3IsValid(h3Index) {
  var ref = h3IndexToSplitLong(h3Index);
  var lower = ref[0];
  var upper = ref[1];
  return Boolean(H3.h3IsValid(lower, upper));
}
function h3IsPentagon(h3Index) {
  var ref = h3IndexToSplitLong(h3Index);
  var lower = ref[0];
  var upper = ref[1];
  return Boolean(H3.h3IsPentagon(lower, upper));
}
function h3IsResClassIII(h3Index) {
  var ref = h3IndexToSplitLong(h3Index);
  var lower = ref[0];
  var upper = ref[1];
  return Boolean(H3.h3IsResClassIII(lower, upper));
}
function h3GetBaseCell(h3Index) {
  var ref = h3IndexToSplitLong(h3Index);
  var lower = ref[0];
  var upper = ref[1];
  return H3.h3GetBaseCell(lower, upper);
}
function h3GetFaces(h3Index) {
  var ref = h3IndexToSplitLong(h3Index);
  var lower = ref[0];
  var upper = ref[1];
  var count = H3.maxFaceCount(lower, upper);
  var faces = libh3._malloc(SZ_INT * count);
  H3.h3GetFaces(lower, upper, faces);
  var out = readArrayOfPositiveIntegers(faces, count);
  libh3._free(faces);
  return out;
}
function h3GetResolution(h3Index) {
  var ref = h3IndexToSplitLong(h3Index);
  var lower = ref[0];
  var upper = ref[1];
  if (!H3.h3IsValid(lower, upper)) {
    return -1;
  }
  return H3.h3GetResolution(lower, upper);
}
function geoToH3(lat, lng, res) {
  var latlng = libh3._malloc(SZ_GEOCOORD);
  libh3.HEAPF64.set([lat, lng].map(degsToRads), latlng / SZ_DBL);
  var h3Index = readH3Index(H3.geoToH3(latlng, res));
  libh3._free(latlng);
  return h3Index;
}
function h3ToGeo(h3Index) {
  var latlng = libh3._malloc(SZ_GEOCOORD);
  var ref = h3IndexToSplitLong(h3Index);
  var lower = ref[0];
  var upper = ref[1];
  H3.h3ToGeo(lower, upper, latlng);
  var out = readGeoCoord(latlng);
  libh3._free(latlng);
  return out;
}
function h3ToGeoBoundary(h3Index, formatAsGeoJson) {
  var geoBoundary = libh3._malloc(SZ_GEOBOUNDARY);
  var ref = h3IndexToSplitLong(h3Index);
  var lower = ref[0];
  var upper = ref[1];
  H3.h3ToGeoBoundary(lower, upper, geoBoundary);
  var out = readGeoBoundary(geoBoundary, formatAsGeoJson, formatAsGeoJson);
  libh3._free(geoBoundary);
  return out;
}
function h3ToParent(h3Index, res) {
  var ref = h3IndexToSplitLong(h3Index);
  var lower = ref[0];
  var upper = ref[1];
  return readH3Index(H3.h3ToParent(lower, upper, res));
}
function h3ToChildren(h3Index, res) {
  if (!h3IsValid(h3Index)) {
    return [];
  }
  var ref = h3IndexToSplitLong(h3Index);
  var lower = ref[0];
  var upper = ref[1];
  var maxCount = H3.maxH3ToChildrenSize(lower, upper, res);
  var hexagons = libh3._calloc(maxCount, SZ_H3INDEX);
  H3.h3ToChildren(lower, upper, res, hexagons);
  var out = readArrayOfHexagons(hexagons, maxCount);
  libh3._free(hexagons);
  return out;
}
function h3ToCenterChild(h3Index, res) {
  var ref = h3IndexToSplitLong(h3Index);
  var lower = ref[0];
  var upper = ref[1];
  return readH3Index(H3.h3ToCenterChild(lower, upper, res));
}
function kRing(h3Index, ringSize) {
  var ref = h3IndexToSplitLong(h3Index);
  var lower = ref[0];
  var upper = ref[1];
  var maxCount = H3.maxKringSize(ringSize);
  var hexagons = libh3._calloc(maxCount, SZ_H3INDEX);
  H3.kRing(lower, upper, ringSize, hexagons);
  var out = readArrayOfHexagons(hexagons, maxCount);
  libh3._free(hexagons);
  return out;
}
function kRingDistances(h3Index, ringSize) {
  var ref = h3IndexToSplitLong(h3Index);
  var lower = ref[0];
  var upper = ref[1];
  var maxCount = H3.maxKringSize(ringSize);
  var kRings = libh3._calloc(maxCount, SZ_H3INDEX);
  var distances = libh3._calloc(maxCount, SZ_INT);
  H3.kRingDistances(lower, upper, ringSize, kRings, distances);
  var out = [];
  for (var i = 0; i < ringSize + 1; i++) {
    out.push([]);
  }
  for (var i$1 = 0; i$1 < maxCount * 2; i$1 += 2) {
    var hexLower = libh3.getValue(kRings + SZ_INT * i$1, "i32");
    var hexUpper = libh3.getValue(kRings + SZ_INT * (i$1 + 1), "i32");
    var index = libh3.getValue(distances + SZ_INT * (i$1 / 2), "i32");
    if (hexLower !== 0 || hexUpper !== 0) {
      out[index].push(splitLongToh3Index(hexLower, hexUpper));
    }
  }
  libh3._free(kRings);
  libh3._free(distances);
  return out;
}
function hexRing(h3Index, ringSize) {
  var maxCount = ringSize === 0 ? 1 : 6 * ringSize;
  var hexagons = libh3._calloc(maxCount, SZ_H3INDEX);
  var retVal = H3.hexRing.apply(H3, h3IndexToSplitLong(h3Index).concat([ringSize], [hexagons]));
  if (retVal !== 0) {
    libh3._free(hexagons);
    throw new Error("Failed to get hexRing (encountered a pentagon?)");
  }
  var out = readArrayOfHexagons(hexagons, maxCount);
  libh3._free(hexagons);
  return out;
}
function polyfill(coordinates, res, isGeoJson) {
  validateRes(res);
  isGeoJson = Boolean(isGeoJson);
  if (coordinates.length === 0 || coordinates[0].length === 0) {
    return [];
  }
  if (typeof coordinates[0][0] === "number") {
    coordinates = [coordinates];
  }
  var geoPolygon = coordinatesToGeoPolygon(coordinates, isGeoJson);
  var arrayLen = H3.maxPolyfillSize(geoPolygon, res);
  var hexagons = libh3._calloc(arrayLen, SZ_H3INDEX);
  H3.polyfill(geoPolygon, res, hexagons);
  var out = readArrayOfHexagons(hexagons, arrayLen);
  libh3._free(hexagons);
  destroyGeoPolygon(geoPolygon);
  return out;
}
function h3SetToMultiPolygon(h3Indexes, formatAsGeoJson) {
  if (!h3Indexes || !h3Indexes.length) {
    return [];
  }
  var indexCount = h3Indexes.length;
  var set = libh3._calloc(indexCount, SZ_H3INDEX);
  storeArrayOfHexagons(set, h3Indexes);
  var polygon = libh3._calloc(SZ_LINKED_GEOPOLYGON);
  var originalPolygon = polygon;
  H3.h3SetToLinkedGeo(set, indexCount, polygon);
  var multiPolygon = readMultiPolygon(polygon, formatAsGeoJson);
  H3.destroyLinkedPolygon(originalPolygon);
  libh3._free(originalPolygon);
  libh3._free(set);
  return multiPolygon;
}
function compact(h3Set) {
  if (!h3Set || !h3Set.length) {
    return [];
  }
  var count = h3Set.length;
  var set = libh3._calloc(count, SZ_H3INDEX);
  storeArrayOfHexagons(set, h3Set);
  var compactedSet = libh3._calloc(count, SZ_H3INDEX);
  var retVal = H3.compact(set, compactedSet, count);
  if (retVal !== 0) {
    libh3._free(set);
    libh3._free(compactedSet);
    throw new Error("Failed to compact, malformed input data (duplicate hexagons?)");
  }
  var out = readArrayOfHexagons(compactedSet, count);
  libh3._free(set);
  libh3._free(compactedSet);
  return out;
}
function uncompact(compactedSet, res) {
  validateRes(res);
  if (!compactedSet || !compactedSet.length) {
    return [];
  }
  var count = compactedSet.length;
  var set = libh3._calloc(count, SZ_H3INDEX);
  storeArrayOfHexagons(set, compactedSet);
  var maxUncompactedNum = H3.maxUncompactSize(set, count, res);
  var uncompactedSet = libh3._calloc(maxUncompactedNum, SZ_H3INDEX);
  var retVal = H3.uncompact(set, count, uncompactedSet, maxUncompactedNum, res);
  if (retVal !== 0) {
    libh3._free(set);
    libh3._free(uncompactedSet);
    throw new Error("Failed to uncompact (bad resolution?)");
  }
  var out = readArrayOfHexagons(uncompactedSet, maxUncompactedNum);
  libh3._free(set);
  libh3._free(uncompactedSet);
  return out;
}
function h3IndexesAreNeighbors(origin, destination) {
  var ref = h3IndexToSplitLong(origin);
  var oLower = ref[0];
  var oUpper = ref[1];
  var ref$1 = h3IndexToSplitLong(destination);
  var dLower = ref$1[0];
  var dUpper = ref$1[1];
  return Boolean(H3.h3IndexesAreNeighbors(oLower, oUpper, dLower, dUpper));
}
function getH3UnidirectionalEdge(origin, destination) {
  var ref = h3IndexToSplitLong(origin);
  var oLower = ref[0];
  var oUpper = ref[1];
  var ref$1 = h3IndexToSplitLong(destination);
  var dLower = ref$1[0];
  var dUpper = ref$1[1];
  return readH3Index(H3.getH3UnidirectionalEdge(oLower, oUpper, dLower, dUpper));
}
function getOriginH3IndexFromUnidirectionalEdge(edgeIndex) {
  var ref = h3IndexToSplitLong(edgeIndex);
  var lower = ref[0];
  var upper = ref[1];
  return readH3Index(H3.getOriginH3IndexFromUnidirectionalEdge(lower, upper));
}
function getDestinationH3IndexFromUnidirectionalEdge(edgeIndex) {
  var ref = h3IndexToSplitLong(edgeIndex);
  var lower = ref[0];
  var upper = ref[1];
  return readH3Index(H3.getDestinationH3IndexFromUnidirectionalEdge(lower, upper));
}
function h3UnidirectionalEdgeIsValid(edgeIndex) {
  var ref = h3IndexToSplitLong(edgeIndex);
  var lower = ref[0];
  var upper = ref[1];
  return Boolean(H3.h3UnidirectionalEdgeIsValid(lower, upper));
}
function getH3IndexesFromUnidirectionalEdge(edgeIndex) {
  var ref = h3IndexToSplitLong(edgeIndex);
  var lower = ref[0];
  var upper = ref[1];
  var count = 2;
  var hexagons = libh3._calloc(count, SZ_H3INDEX);
  H3.getH3IndexesFromUnidirectionalEdge(lower, upper, hexagons);
  var out = readArrayOfHexagons(hexagons, count);
  libh3._free(hexagons);
  return out;
}
function getH3UnidirectionalEdgesFromHexagon(h3Index) {
  var ref = h3IndexToSplitLong(h3Index);
  var lower = ref[0];
  var upper = ref[1];
  var count = 6;
  var edges = libh3._calloc(count, SZ_H3INDEX);
  H3.getH3UnidirectionalEdgesFromHexagon(lower, upper, edges);
  var out = readArrayOfHexagons(edges, count);
  libh3._free(edges);
  return out;
}
function getH3UnidirectionalEdgeBoundary(edgeIndex, formatAsGeoJson) {
  var geoBoundary = libh3._malloc(SZ_GEOBOUNDARY);
  var ref = h3IndexToSplitLong(edgeIndex);
  var lower = ref[0];
  var upper = ref[1];
  H3.getH3UnidirectionalEdgeBoundary(lower, upper, geoBoundary);
  var out = readGeoBoundary(geoBoundary, formatAsGeoJson);
  libh3._free(geoBoundary);
  return out;
}
function h3Distance(origin, destination) {
  var ref = h3IndexToSplitLong(origin);
  var oLower = ref[0];
  var oUpper = ref[1];
  var ref$1 = h3IndexToSplitLong(destination);
  var dLower = ref$1[0];
  var dUpper = ref$1[1];
  return H3.h3Distance(oLower, oUpper, dLower, dUpper);
}
function h3Line(origin, destination) {
  var ref = h3IndexToSplitLong(origin);
  var oLower = ref[0];
  var oUpper = ref[1];
  var ref$1 = h3IndexToSplitLong(destination);
  var dLower = ref$1[0];
  var dUpper = ref$1[1];
  var count = H3.h3LineSize(oLower, oUpper, dLower, dUpper);
  if (count < 0) {
    throw new Error("Line cannot be calculated");
  }
  var hexagons = libh3._calloc(count, SZ_H3INDEX);
  H3.h3Line(oLower, oUpper, dLower, dUpper, hexagons);
  var out = readArrayOfHexagons(hexagons, count);
  libh3._free(hexagons);
  return out;
}
function experimentalH3ToLocalIj(origin, destination) {
  var ij = libh3._malloc(SZ_COORDIJ);
  var retVal = H3.experimentalH3ToLocalIj.apply(H3, h3IndexToSplitLong(origin).concat(h3IndexToSplitLong(destination), [ij]));
  var coords = readCoordIJ(ij);
  libh3._free(ij);
  switch (retVal) {
    case 0:
      return coords;
    case 1:
      throw new Error("Incompatible origin and index.");
    case 2:
    default:
      throw new Error("Local IJ coordinates undefined for this origin and index pair. The index may be too far from the origin.");
    case 3:
    case 4:
    case 5:
      throw new Error("Encountered possible pentagon distortion");
  }
}
function experimentalLocalIjToH3(origin, coords) {
  if (!coords || typeof coords.i !== "number" || typeof coords.j !== "number") {
    throw new Error("Coordinates must be provided as an {i, j} object");
  }
  var ij = libh3._malloc(SZ_COORDIJ);
  var out = libh3._malloc(SZ_H3INDEX);
  storeCoordIJ(ij, coords);
  var retVal = H3.experimentalLocalIjToH3.apply(H3, h3IndexToSplitLong(origin).concat([ij], [out]));
  var h3Index = readH3IndexFromPointer(out);
  libh3._free(ij);
  libh3._free(out);
  if (retVal !== 0) {
    throw new Error("Index not defined for this origin and IJ coordinates pair. IJ coordinates may be too far from origin, or a pentagon distortion was encountered.");
  }
  return h3Index;
}
function pointDist(latlng1, latlng2, unit) {
  var coord1 = storeGeoCoord(latlng1[0], latlng1[1]);
  var coord2 = storeGeoCoord(latlng2[0], latlng2[1]);
  var result;
  switch (unit) {
    case UNITS.m:
      result = H3.pointDistM(coord1, coord2);
      break;
    case UNITS.km:
      result = H3.pointDistKm(coord1, coord2);
      break;
    case UNITS.rads:
      result = H3.pointDistRads(coord1, coord2);
      break;
    default:
      result = null;
  }
  libh3._free(coord1);
  libh3._free(coord2);
  if (result === null) {
    throw new Error("Unknown unit: " + unit);
  }
  return result;
}
function cellArea(h3Index, unit) {
  var ref = h3IndexToSplitLong(h3Index);
  var lower = ref[0];
  var upper = ref[1];
  switch (unit) {
    case UNITS.m2:
      return H3.cellAreaM2(lower, upper);
    case UNITS.km2:
      return H3.cellAreaKm2(lower, upper);
    case UNITS.rads2:
      return H3.cellAreaRads2(lower, upper);
    default:
      throw new Error("Unknown unit: " + unit);
  }
}
function exactEdgeLength(edge, unit) {
  var ref = h3IndexToSplitLong(edge);
  var lower = ref[0];
  var upper = ref[1];
  switch (unit) {
    case UNITS.m:
      return H3.exactEdgeLengthM(lower, upper);
    case UNITS.km:
      return H3.exactEdgeLengthKm(lower, upper);
    case UNITS.rads:
      return H3.exactEdgeLengthRads(lower, upper);
    default:
      throw new Error("Unknown unit: " + unit);
  }
}
function hexArea(res, unit) {
  validateRes(res);
  switch (unit) {
    case UNITS.m2:
      return H3.hexAreaM2(res);
    case UNITS.km2:
      return H3.hexAreaKm2(res);
    default:
      throw new Error("Unknown unit: " + unit);
  }
}
function edgeLength(res, unit) {
  validateRes(res);
  switch (unit) {
    case UNITS.m:
      return H3.edgeLengthM(res);
    case UNITS.km:
      return H3.edgeLengthKm(res);
    default:
      throw new Error("Unknown unit: " + unit);
  }
}
function numHexagons(res) {
  validateRes(res);
  var ref = readLong(H3.numHexagons(res));
  var lower = ref[0];
  var upper = ref[1];
  if (!upper) {
    return lower;
  }
  return upper * Math.pow(2, 32) + lower;
}
function getRes0Indexes() {
  var count = H3.res0IndexCount();
  var hexagons = libh3._malloc(SZ_H3INDEX * count);
  H3.getRes0Indexes(hexagons);
  var out = readArrayOfHexagons(hexagons, count);
  libh3._free(hexagons);
  return out;
}
function getPentagonIndexes(res) {
  validateRes(res);
  var count = H3.pentagonIndexCount();
  var hexagons = libh3._malloc(SZ_H3INDEX * count);
  H3.getPentagonIndexes(res, hexagons);
  var out = readArrayOfHexagons(hexagons, count);
  libh3._free(hexagons);
  return out;
}
function degsToRads(deg) {
  return deg * Math.PI / 180;
}
function radsToDegs(rad) {
  return rad * 180 / Math.PI;
}
var libh3, NUMBER, BOOLEAN, H3_LOWER, H3_UPPER, RESOLUTION, POINTER, BINDINGS, H3, BASE_16, SZ_INT, SZ_PTR, SZ_DBL, SZ_H3INDEX, SZ_GEOCOORD, SZ_GEOBOUNDARY, SZ_GEOPOLYGON, SZ_GEOFENCE, SZ_LINKED_GEOPOLYGON, SZ_COORDIJ, UNITS, INVALID_HEXIDECIMAL_CHAR;
var init_h3_js_es = __esm({
  "node_modules/geojson2h3/node_modules/h3-js/dist/browser/h3-js.es.js"() {
    libh3 = function(libh32) {
      libh32 = libh32 || {};
      var Module = typeof libh32 !== "undefined" ? libh32 : {};
      var moduleOverrides = {};
      var key;
      for (key in Module) {
        if (Module.hasOwnProperty(key)) {
          moduleOverrides[key] = Module[key];
        }
      }
      var arguments_ = [];
      var scriptDirectory = "";
      function locateFile(path) {
        if (Module["locateFile"]) {
          return Module["locateFile"](path, scriptDirectory);
        }
        return scriptDirectory + path;
      }
      var readAsync;
      {
        if (document.currentScript) {
          scriptDirectory = document.currentScript.src;
        }
        if (scriptDirectory.indexOf("blob:") !== 0) {
          scriptDirectory = scriptDirectory.substr(0, scriptDirectory.lastIndexOf("/") + 1);
        } else {
          scriptDirectory = "";
        }
        readAsync = function readAsync2(url, onload, onerror) {
          var xhr = new XMLHttpRequest();
          xhr.open("GET", url, true);
          xhr.responseType = "arraybuffer";
          xhr.onload = function xhr_onload() {
            if (xhr.status == 200 || xhr.status == 0 && xhr.response) {
              onload(xhr.response);
              return;
            }
            var data = tryParseAsDataURI(url);
            if (data) {
              onload(data.buffer);
              return;
            }
            onerror();
          };
          xhr.onerror = onerror;
          xhr.send(null);
        };
      }
      var out = Module["print"] || console.log.bind(console);
      var err = Module["printErr"] || console.warn.bind(console);
      for (key in moduleOverrides) {
        if (moduleOverrides.hasOwnProperty(key)) {
          Module[key] = moduleOverrides[key];
        }
      }
      moduleOverrides = null;
      if (Module["arguments"]) {
        arguments_ = Module["arguments"];
      }
      var tempRet0 = 0;
      var setTempRet0 = function(value) {
        tempRet0 = value;
      };
      var getTempRet0 = function() {
        return tempRet0;
      };
      var GLOBAL_BASE = 8;
      function setValue(ptr, value, type, noSafe) {
        type = type || "i8";
        if (type.charAt(type.length - 1) === "*") {
          type = "i32";
        }
        switch (type) {
          case "i1":
            HEAP8[ptr >> 0] = value;
            break;
          case "i8":
            HEAP8[ptr >> 0] = value;
            break;
          case "i16":
            HEAP16[ptr >> 1] = value;
            break;
          case "i32":
            HEAP32[ptr >> 2] = value;
            break;
          case "i64":
            tempI64 = [value >>> 0, (tempDouble = value, +Math_abs(tempDouble) >= 1 ? tempDouble > 0 ? (Math_min(+Math_floor(tempDouble / 4294967296), 4294967295) | 0) >>> 0 : ~~+Math_ceil((tempDouble - +(~~tempDouble >>> 0)) / 4294967296) >>> 0 : 0)], HEAP32[ptr >> 2] = tempI64[0], HEAP32[ptr + 4 >> 2] = tempI64[1];
            break;
          case "float":
            HEAPF32[ptr >> 2] = value;
            break;
          case "double":
            HEAPF64[ptr >> 3] = value;
            break;
          default:
            abort("invalid type for setValue: " + type);
        }
      }
      function getValue(ptr, type, noSafe) {
        type = type || "i8";
        if (type.charAt(type.length - 1) === "*") {
          type = "i32";
        }
        switch (type) {
          case "i1":
            return HEAP8[ptr >> 0];
          case "i8":
            return HEAP8[ptr >> 0];
          case "i16":
            return HEAP16[ptr >> 1];
          case "i32":
            return HEAP32[ptr >> 2];
          case "i64":
            return HEAP32[ptr >> 2];
          case "float":
            return HEAPF32[ptr >> 2];
          case "double":
            return HEAPF64[ptr >> 3];
          default:
            abort("invalid type for getValue: " + type);
        }
        return null;
      }
      var ABORT = false;
      function assert(condition, text) {
        if (!condition) {
          abort("Assertion failed: " + text);
        }
      }
      function getCFunc(ident) {
        var func = Module["_" + ident];
        assert(func, "Cannot call unknown function " + ident + ", make sure it is exported");
        return func;
      }
      function ccall(ident, returnType, argTypes, args, opts) {
        var toC = {
          "string": function(str) {
            var ret2 = 0;
            if (str !== null && str !== void 0 && str !== 0) {
              var len = (str.length << 2) + 1;
              ret2 = stackAlloc(len);
              stringToUTF8(str, ret2, len);
            }
            return ret2;
          },
          "array": function(arr) {
            var ret2 = stackAlloc(arr.length);
            writeArrayToMemory(arr, ret2);
            return ret2;
          }
        };
        function convertReturnValue(ret2) {
          if (returnType === "string") {
            return UTF8ToString(ret2);
          }
          if (returnType === "boolean") {
            return Boolean(ret2);
          }
          return ret2;
        }
        var func = getCFunc(ident);
        var cArgs = [];
        var stack = 0;
        if (args) {
          for (var i = 0; i < args.length; i++) {
            var converter = toC[argTypes[i]];
            if (converter) {
              if (stack === 0) {
                stack = stackSave();
              }
              cArgs[i] = converter(args[i]);
            } else {
              cArgs[i] = args[i];
            }
          }
        }
        var ret = func.apply(null, cArgs);
        ret = convertReturnValue(ret);
        if (stack !== 0) {
          stackRestore(stack);
        }
        return ret;
      }
      function cwrap(ident, returnType, argTypes, opts) {
        argTypes = argTypes || [];
        var numericArgs = argTypes.every(function(type) {
          return type === "number";
        });
        var numericRet = returnType !== "string";
        if (numericRet && numericArgs && !opts) {
          return getCFunc(ident);
        }
        return function() {
          return ccall(ident, returnType, argTypes, arguments, opts);
        };
      }
      var UTF8Decoder = typeof TextDecoder !== "undefined" ? new TextDecoder("utf8") : void 0;
      function UTF8ArrayToString(u8Array, idx, maxBytesToRead) {
        var endIdx = idx + maxBytesToRead;
        var endPtr = idx;
        while (u8Array[endPtr] && !(endPtr >= endIdx)) {
          ++endPtr;
        }
        if (endPtr - idx > 16 && u8Array.subarray && UTF8Decoder) {
          return UTF8Decoder.decode(u8Array.subarray(idx, endPtr));
        } else {
          var str = "";
          while (idx < endPtr) {
            var u0 = u8Array[idx++];
            if (!(u0 & 128)) {
              str += String.fromCharCode(u0);
              continue;
            }
            var u1 = u8Array[idx++] & 63;
            if ((u0 & 224) == 192) {
              str += String.fromCharCode((u0 & 31) << 6 | u1);
              continue;
            }
            var u2 = u8Array[idx++] & 63;
            if ((u0 & 240) == 224) {
              u0 = (u0 & 15) << 12 | u1 << 6 | u2;
            } else {
              u0 = (u0 & 7) << 18 | u1 << 12 | u2 << 6 | u8Array[idx++] & 63;
            }
            if (u0 < 65536) {
              str += String.fromCharCode(u0);
            } else {
              var ch = u0 - 65536;
              str += String.fromCharCode(55296 | ch >> 10, 56320 | ch & 1023);
            }
          }
        }
        return str;
      }
      function UTF8ToString(ptr, maxBytesToRead) {
        return ptr ? UTF8ArrayToString(HEAPU8, ptr, maxBytesToRead) : "";
      }
      function stringToUTF8Array(str, outU8Array, outIdx, maxBytesToWrite) {
        if (!(maxBytesToWrite > 0)) {
          return 0;
        }
        var startIdx = outIdx;
        var endIdx = outIdx + maxBytesToWrite - 1;
        for (var i = 0; i < str.length; ++i) {
          var u = str.charCodeAt(i);
          if (u >= 55296 && u <= 57343) {
            var u1 = str.charCodeAt(++i);
            u = 65536 + ((u & 1023) << 10) | u1 & 1023;
          }
          if (u <= 127) {
            if (outIdx >= endIdx) {
              break;
            }
            outU8Array[outIdx++] = u;
          } else if (u <= 2047) {
            if (outIdx + 1 >= endIdx) {
              break;
            }
            outU8Array[outIdx++] = 192 | u >> 6;
            outU8Array[outIdx++] = 128 | u & 63;
          } else if (u <= 65535) {
            if (outIdx + 2 >= endIdx) {
              break;
            }
            outU8Array[outIdx++] = 224 | u >> 12;
            outU8Array[outIdx++] = 128 | u >> 6 & 63;
            outU8Array[outIdx++] = 128 | u & 63;
          } else {
            if (outIdx + 3 >= endIdx) {
              break;
            }
            outU8Array[outIdx++] = 240 | u >> 18;
            outU8Array[outIdx++] = 128 | u >> 12 & 63;
            outU8Array[outIdx++] = 128 | u >> 6 & 63;
            outU8Array[outIdx++] = 128 | u & 63;
          }
        }
        outU8Array[outIdx] = 0;
        return outIdx - startIdx;
      }
      function stringToUTF8(str, outPtr, maxBytesToWrite) {
        return stringToUTF8Array(str, HEAPU8, outPtr, maxBytesToWrite);
      }
      var UTF16Decoder = typeof TextDecoder !== "undefined" ? new TextDecoder("utf-16le") : void 0;
      function writeArrayToMemory(array, buffer2) {
        HEAP8.set(array, buffer2);
      }
      function alignUp(x, multiple) {
        if (x % multiple > 0) {
          x += multiple - x % multiple;
        }
        return x;
      }
      var buffer, HEAP8, HEAPU8, HEAP16, HEAPU16, HEAP32, HEAPU32, HEAPF32, HEAPF64;
      function updateGlobalBufferAndViews(buf) {
        buffer = buf;
        Module["HEAP8"] = HEAP8 = new Int8Array(buf);
        Module["HEAP16"] = HEAP16 = new Int16Array(buf);
        Module["HEAP32"] = HEAP32 = new Int32Array(buf);
        Module["HEAPU8"] = HEAPU8 = new Uint8Array(buf);
        Module["HEAPU16"] = HEAPU16 = new Uint16Array(buf);
        Module["HEAPU32"] = HEAPU32 = new Uint32Array(buf);
        Module["HEAPF32"] = HEAPF32 = new Float32Array(buf);
        Module["HEAPF64"] = HEAPF64 = new Float64Array(buf);
      }
      var DYNAMIC_BASE = 5266928, DYNAMICTOP_PTR = 24016;
      var INITIAL_TOTAL_MEMORY = Module["TOTAL_MEMORY"] || 33554432;
      if (Module["buffer"]) {
        buffer = Module["buffer"];
      } else {
        buffer = new ArrayBuffer(INITIAL_TOTAL_MEMORY);
      }
      INITIAL_TOTAL_MEMORY = buffer.byteLength;
      updateGlobalBufferAndViews(buffer);
      HEAP32[DYNAMICTOP_PTR >> 2] = DYNAMIC_BASE;
      function callRuntimeCallbacks(callbacks) {
        while (callbacks.length > 0) {
          var callback = callbacks.shift();
          if (typeof callback == "function") {
            callback();
            continue;
          }
          var func = callback.func;
          if (typeof func === "number") {
            if (callback.arg === void 0) {
              Module["dynCall_v"](func);
            } else {
              Module["dynCall_vi"](func, callback.arg);
            }
          } else {
            func(callback.arg === void 0 ? null : callback.arg);
          }
        }
      }
      var __ATPRERUN__ = [];
      var __ATINIT__ = [];
      var __ATMAIN__ = [];
      var __ATPOSTRUN__ = [];
      function preRun() {
        if (Module["preRun"]) {
          if (typeof Module["preRun"] == "function") {
            Module["preRun"] = [Module["preRun"]];
          }
          while (Module["preRun"].length) {
            addOnPreRun(Module["preRun"].shift());
          }
        }
        callRuntimeCallbacks(__ATPRERUN__);
      }
      function initRuntime() {
        callRuntimeCallbacks(__ATINIT__);
      }
      function preMain() {
        callRuntimeCallbacks(__ATMAIN__);
      }
      function postRun() {
        if (Module["postRun"]) {
          if (typeof Module["postRun"] == "function") {
            Module["postRun"] = [Module["postRun"]];
          }
          while (Module["postRun"].length) {
            addOnPostRun(Module["postRun"].shift());
          }
        }
        callRuntimeCallbacks(__ATPOSTRUN__);
      }
      function addOnPreRun(cb) {
        __ATPRERUN__.unshift(cb);
      }
      function addOnPostRun(cb) {
        __ATPOSTRUN__.unshift(cb);
      }
      var Math_abs = Math.abs;
      var Math_ceil = Math.ceil;
      var Math_floor = Math.floor;
      var Math_min = Math.min;
      var runDependencies = 0;
      var runDependencyWatcher = null;
      var dependenciesFulfilled = null;
      function addRunDependency(id) {
        runDependencies++;
        if (Module["monitorRunDependencies"]) {
          Module["monitorRunDependencies"](runDependencies);
        }
      }
      function removeRunDependency(id) {
        runDependencies--;
        if (Module["monitorRunDependencies"]) {
          Module["monitorRunDependencies"](runDependencies);
        }
        if (runDependencies == 0) {
          if (runDependencyWatcher !== null) {
            clearInterval(runDependencyWatcher);
            runDependencyWatcher = null;
          }
          if (dependenciesFulfilled) {
            var callback = dependenciesFulfilled;
            dependenciesFulfilled = null;
            callback();
          }
        }
      }
      Module["preloadedImages"] = {};
      Module["preloadedAudios"] = {};
      var memoryInitializer = null;
      var dataURIPrefix = "data:application/octet-stream;base64,";
      function isDataURI(filename) {
        return String.prototype.startsWith ? filename.startsWith(dataURIPrefix) : filename.indexOf(dataURIPrefix) === 0;
      }
      var tempDouble;
      var tempI64;
      memoryInitializer = "data:application/octet-stream;base64,AAAAAAAAAAACAAAAAwAAAAEAAAAFAAAABAAAAAYAAAAAAAAAAAAAAAAAAAABAAAAAgAAAAMAAAAEAAAABQAAAAYAAAABAAAABAAAAAMAAAAGAAAABQAAAAIAAAAAAAAAAgAAAAMAAAABAAAABAAAAAYAAAAAAAAABQAAAAMAAAAGAAAABAAAAAUAAAAAAAAAAQAAAAIAAAAEAAAABQAAAAYAAAAAAAAAAgAAAAMAAAABAAAABQAAAAIAAAAAAAAAAQAAAAMAAAAGAAAABAAAAAYAAAAAAAAABQAAAAIAAAABAAAABAAAAAMAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAEAAAAAAAAAAQAAAAAAAAAFAAAAAAAAAAAAAAAAAAAAAgAAAAMAAAAAAAAAAAAAAAIAAAAAAAAAAQAAAAMAAAADAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAQAAAAEAAAABgAAAAAAAAAFAAAAAAAAAAAAAAAEAAAABQAAAAAAAAAAAAAAAAAAAAIAAAAAAAAABgAAAAAAAAAGAAAAAAAAAAAAAAAAAAAAAAAAAAEAAAACAAAAAwAAAAQAAAAFAAAABgAAAAEAAAACAAAAAwAAAAQAAAAFAAAABgAAAAAAAAACAAAAAwAAAAQAAAAFAAAABgAAAAAAAAABAAAAAwAAAAQAAAAFAAAABgAAAAAAAAABAAAAAgAAAAQAAAAFAAAABgAAAAAAAAABAAAAAgAAAAMAAAAFAAAABgAAAAAAAAABAAAAAgAAAAMAAAAEAAAABgAAAAAAAAABAAAAAgAAAAMAAAAEAAAABQAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAQAAAAAAAAADAAAAAAAAAAEAAAAAAAAAAAAAAAAAAAACAAAAAgAAAAAAAAAAAAAABgAAAAAAAAADAAAAAgAAAAMAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAABAAAAAUAAAAEAAAAAAAAAAEAAAAAAAAAAAAAAAUAAAAFAAAAAAAAAAAAAAAAAAAABgAAAAAAAAAEAAAAAAAAAAYAAAAAAAAAAAAAAAAAAAAAAAAAAQAAAAUAAAACAAAABAAAAAMAAAAIAAAAAQAAAAcAAAAGAAAACQAAAAAAAAADAAAAAgAAAAIAAAAGAAAACgAAAAsAAAAAAAAAAQAAAAUAAAADAAAADQAAAAEAAAAHAAAABAAAAAwAAAAAAAAABAAAAH8AAAAPAAAACAAAAAMAAAAAAAAADAAAAAUAAAACAAAAEgAAAAoAAAAIAAAAAAAAABAAAAAGAAAADgAAAAsAAAARAAAAAQAAAAkAAAACAAAABwAAABUAAAAJAAAAEwAAAAMAAAANAAAAAQAAAAgAAAAFAAAAFgAAABAAAAAEAAAAAAAAAA8AAAAJAAAAEwAAAA4AAAAUAAAAAQAAAAcAAAAGAAAACgAAAAsAAAAYAAAAFwAAAAUAAAACAAAAEgAAAAsAAAARAAAAFwAAABkAAAACAAAABgAAAAoAAAAMAAAAHAAAAA0AAAAaAAAABAAAAA8AAAADAAAADQAAABoAAAAVAAAAHQAAAAMAAAAMAAAABwAAAA4AAAB/AAAAEQAAABsAAAAJAAAAFAAAAAYAAAAPAAAAFgAAABwAAAAfAAAABAAAAAgAAAAMAAAAEAAAABIAAAAhAAAAHgAAAAgAAAAFAAAAFgAAABEAAAALAAAADgAAAAYAAAAjAAAAGQAAABsAAAASAAAAGAAAAB4AAAAgAAAABQAAAAoAAAAQAAAAEwAAACIAAAAUAAAAJAAAAAcAAAAVAAAACQAAABQAAAAOAAAAEwAAAAkAAAAoAAAAGwAAACQAAAAVAAAAJgAAABMAAAAiAAAADQAAAB0AAAAHAAAAFgAAABAAAAApAAAAIQAAAA8AAAAIAAAAHwAAABcAAAAYAAAACwAAAAoAAAAnAAAAJQAAABkAAAAYAAAAfwAAACAAAAAlAAAACgAAABcAAAASAAAAGQAAABcAAAARAAAACwAAAC0AAAAnAAAAIwAAABoAAAAqAAAAHQAAACsAAAAMAAAAHAAAAA0AAAAbAAAAKAAAACMAAAAuAAAADgAAABQAAAARAAAAHAAAAB8AAAAqAAAALAAAAAwAAAAPAAAAGgAAAB0AAAArAAAAJgAAAC8AAAANAAAAGgAAABUAAAAeAAAAIAAAADAAAAAyAAAAEAAAABIAAAAhAAAAHwAAACkAAAAsAAAANQAAAA8AAAAWAAAAHAAAACAAAAAeAAAAGAAAABIAAAA0AAAAMgAAACUAAAAhAAAAHgAAADEAAAAwAAAAFgAAABAAAAApAAAAIgAAABMAAAAmAAAAFQAAADYAAAAkAAAAMwAAACMAAAAuAAAALQAAADgAAAARAAAAGwAAABkAAAAkAAAAFAAAACIAAAATAAAANwAAACgAAAA2AAAAJQAAACcAAAA0AAAAOQAAABgAAAAXAAAAIAAAACYAAAB/AAAAIgAAADMAAAAdAAAALwAAABUAAAAnAAAAJQAAABkAAAAXAAAAOwAAADkAAAAtAAAAKAAAABsAAAAkAAAAFAAAADwAAAAuAAAANwAAACkAAAAxAAAANQAAAD0AAAAWAAAAIQAAAB8AAAAqAAAAOgAAACsAAAA+AAAAHAAAACwAAAAaAAAAKwAAAD4AAAAvAAAAQAAAABoAAAAqAAAAHQAAACwAAAA1AAAAOgAAAEEAAAAcAAAAHwAAACoAAAAtAAAAJwAAACMAAAAZAAAAPwAAADsAAAA4AAAALgAAADwAAAA4AAAARAAAABsAAAAoAAAAIwAAAC8AAAAmAAAAKwAAAB0AAABFAAAAMwAAAEAAAAAwAAAAMQAAAB4AAAAhAAAAQwAAAEIAAAAyAAAAMQAAAH8AAAA9AAAAQgAAACEAAAAwAAAAKQAAADIAAAAwAAAAIAAAAB4AAABGAAAAQwAAADQAAAAzAAAARQAAADYAAABHAAAAJgAAAC8AAAAiAAAANAAAADkAAABGAAAASgAAACAAAAAlAAAAMgAAADUAAAA9AAAAQQAAAEsAAAAfAAAAKQAAACwAAAA2AAAARwAAADcAAABJAAAAIgAAADMAAAAkAAAANwAAACgAAAA2AAAAJAAAAEgAAAA8AAAASQAAADgAAABEAAAAPwAAAE0AAAAjAAAALgAAAC0AAAA5AAAAOwAAAEoAAABOAAAAJQAAACcAAAA0AAAAOgAAAH8AAAA+AAAATAAAACwAAABBAAAAKgAAADsAAAA/AAAATgAAAE8AAAAnAAAALQAAADkAAAA8AAAASAAAAEQAAABQAAAAKAAAADcAAAAuAAAAPQAAADUAAAAxAAAAKQAAAFEAAABLAAAAQgAAAD4AAAArAAAAOgAAACoAAABSAAAAQAAAAEwAAAA/AAAAfwAAADgAAAAtAAAATwAAADsAAABNAAAAQAAAAC8AAAA+AAAAKwAAAFQAAABFAAAAUgAAAEEAAAA6AAAANQAAACwAAABWAAAATAAAAEsAAABCAAAAQwAAAFEAAABVAAAAMQAAADAAAAA9AAAAQwAAAEIAAAAyAAAAMAAAAFcAAABVAAAARgAAAEQAAAA4AAAAPAAAAC4AAABaAAAATQAAAFAAAABFAAAAMwAAAEAAAAAvAAAAWQAAAEcAAABUAAAARgAAAEMAAAA0AAAAMgAAAFMAAABXAAAASgAAAEcAAABZAAAASQAAAFsAAAAzAAAARQAAADYAAABIAAAAfwAAAEkAAAA3AAAAUAAAADwAAABYAAAASQAAAFsAAABIAAAAWAAAADYAAABHAAAANwAAAEoAAABOAAAAUwAAAFwAAAA0AAAAOQAAAEYAAABLAAAAQQAAAD0AAAA1AAAAXgAAAFYAAABRAAAATAAAAFYAAABSAAAAYAAAADoAAABBAAAAPgAAAE0AAAA/AAAARAAAADgAAABdAAAATwAAAFoAAABOAAAASgAAADsAAAA5AAAAXwAAAFwAAABPAAAATwAAAE4AAAA/AAAAOwAAAF0AAABfAAAATQAAAFAAAABEAAAASAAAADwAAABjAAAAWgAAAFgAAABRAAAAVQAAAF4AAABlAAAAPQAAAEIAAABLAAAAUgAAAGAAAABUAAAAYgAAAD4AAABMAAAAQAAAAFMAAAB/AAAASgAAAEYAAABkAAAAVwAAAFwAAABUAAAARQAAAFIAAABAAAAAYQAAAFkAAABiAAAAVQAAAFcAAABlAAAAZgAAAEIAAABDAAAAUQAAAFYAAABMAAAASwAAAEEAAABoAAAAYAAAAF4AAABXAAAAUwAAAGYAAABkAAAAQwAAAEYAAABVAAAAWAAAAEgAAABbAAAASQAAAGMAAABQAAAAaQAAAFkAAABhAAAAWwAAAGcAAABFAAAAVAAAAEcAAABaAAAATQAAAFAAAABEAAAAagAAAF0AAABjAAAAWwAAAEkAAABZAAAARwAAAGkAAABYAAAAZwAAAFwAAABTAAAATgAAAEoAAABsAAAAZAAAAF8AAABdAAAATwAAAFoAAABNAAAAbQAAAF8AAABqAAAAXgAAAFYAAABRAAAASwAAAGsAAABoAAAAZQAAAF8AAABcAAAATwAAAE4AAABtAAAAbAAAAF0AAABgAAAAaAAAAGIAAABuAAAATAAAAFYAAABSAAAAYQAAAH8AAABiAAAAVAAAAGcAAABZAAAAbwAAAGIAAABuAAAAYQAAAG8AAABSAAAAYAAAAFQAAABjAAAAUAAAAGkAAABYAAAAagAAAFoAAABxAAAAZAAAAGYAAABTAAAAVwAAAGwAAAByAAAAXAAAAGUAAABmAAAAawAAAHAAAABRAAAAVQAAAF4AAABmAAAAZQAAAFcAAABVAAAAcgAAAHAAAABkAAAAZwAAAFsAAABhAAAAWQAAAHQAAABpAAAAbwAAAGgAAABrAAAAbgAAAHMAAABWAAAAXgAAAGAAAABpAAAAWAAAAGcAAABbAAAAcQAAAGMAAAB0AAAAagAAAF0AAABjAAAAWgAAAHUAAABtAAAAcQAAAGsAAAB/AAAAZQAAAF4AAABzAAAAaAAAAHAAAABsAAAAZAAAAF8AAABcAAAAdgAAAHIAAABtAAAAbQAAAGwAAABdAAAAXwAAAHUAAAB2AAAAagAAAG4AAABiAAAAaAAAAGAAAAB3AAAAbwAAAHMAAABvAAAAYQAAAG4AAABiAAAAdAAAAGcAAAB3AAAAcAAAAGsAAABmAAAAZQAAAHgAAABzAAAAcgAAAHEAAABjAAAAdAAAAGkAAAB1AAAAagAAAHkAAAByAAAAcAAAAGQAAABmAAAAdgAAAHgAAABsAAAAcwAAAG4AAABrAAAAaAAAAHgAAAB3AAAAcAAAAHQAAABnAAAAdwAAAG8AAABxAAAAaQAAAHkAAAB1AAAAfwAAAG0AAAB2AAAAcQAAAHkAAABqAAAAdgAAAHgAAABsAAAAcgAAAHUAAAB5AAAAbQAAAHcAAABvAAAAcwAAAG4AAAB5AAAAdAAAAHgAAAB4AAAAcwAAAHIAAABwAAAAeQAAAHcAAAB2AAAAeQAAAHQAAAB4AAAAdwAAAHUAAABxAAAAdgAAAAAAAAAAAAAAAAAAAAUAAAAAAAAAAAAAAAEAAAAFAAAAAQAAAAAAAAAAAAAAAQAAAAAAAAABAAAAAAAAAAEAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAUAAAAAAAAAAAAAAAUAAAAAAAAAAAAAAAIAAAAFAAAAAQAAAAAAAAD/////AQAAAAAAAAADAAAABAAAAAIAAAAAAAAAAAAAAAEAAAAAAAAAAQAAAAAAAAABAAAAAAAAAAAAAAAAAAAAAwAAAAUAAAAFAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAABQAAAAAAAAAAAAAABQAAAAAAAAAAAAAAAAAAAAUAAAABAAAAAAAAAAAAAAABAAAAAwAAAAAAAAAAAAAAAQAAAAAAAAAAAAAAAQAAAAMAAAAAAAAAAAAAAAEAAAAAAAAAAwAAAAMAAAADAAAAAAAAAAAAAAAAAAAAAAAAAAUAAAAAAAAAAAAAAAMAAAAFAAAAAQAAAAAAAAAAAAAAAQAAAAAAAAABAAAAAAAAAAEAAAAAAAAA/////wMAAAAAAAAABQAAAAIAAAAAAAAAAAAAAAUAAAAAAAAAAAAAAAQAAAAFAAAAAQAAAAAAAAAAAAAAAAAAAAAAAAAAAAAABQAAAAAAAAAAAAAAAwAAAAMAAAADAAAAAwAAAAAAAAADAAAAAAAAAAAAAAAAAAAAAwAAAAUAAAAFAAAAAAAAAAAAAAADAAAAAwAAAAMAAAAAAAAAAAAAAAAAAAAAAAAAAwAAAAMAAAADAAAAAAAAAAMAAAAAAAAAAAAAAAAAAAAAAAAAAwAAAAUAAAAFAAAAAAAAAAAAAAAAAAAAAQAAAAAAAAABAAAAAAAAAAEAAAAAAAAAAwAAAAMAAAADAAAAAAAAAAMAAAAAAAAAAAAAAP////8DAAAAAAAAAAUAAAACAAAAAAAAAAAAAAAAAAAAAAAAAAMAAAAAAAAAAAAAAAMAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAUAAAAAAAAAAAAAAAMAAAAAAAAAAAAAAAAAAAADAAAAAwAAAAAAAAAAAAAAAQAAAAAAAAABAAAAAAAAAAEAAAAAAAAAAAAAAAEAAAADAAAAAAAAAAAAAAABAAAAAAAAAAMAAAADAAAAAwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAABQAAAAAAAAAAAAAAAwAAAAMAAAADAAAAAwAAAAAAAAADAAAAAAAAAAAAAAABAAAAAwAAAAAAAAAAAAAAAQAAAAAAAAADAAAAAwAAAAMAAAADAAAAAAAAAAMAAAAAAAAAAAAAAAMAAAAAAAAAAwAAAAAAAAADAAAAAAAAAAAAAAAAAAAAAwAAAAAAAAAAAAAAAwAAAAAAAAADAAAAAAAAAAAAAAAAAAAAAwAAAAMAAAAAAAAA/////wMAAAAAAAAABQAAAAIAAAAAAAAAAAAAAAMAAAAAAAAAAAAAAAMAAAADAAAAAAAAAAAAAAADAAAAAAAAAAAAAAADAAAAAwAAAAAAAAAAAAAAAAAAAAAAAAADAAAABQAAAAUAAAAAAAAAAAAAAAAAAAAAAAAAAwAAAAUAAAAFAAAAAAAAAAAAAAADAAAAAwAAAAMAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAEAAAADAAAAAAAAAAAAAAABAAAAAAAAAAAAAAADAAAAAAAAAAAAAAADAAAAAwAAAAAAAAAAAAAAAAAAAAMAAAAAAAAAAwAAAAAAAAAAAAAAAwAAAAMAAAADAAAAAAAAAAMAAAAAAAAAAAAAAAMAAAADAAAAAwAAAAAAAAADAAAAAAAAAAAAAAD/////AwAAAAAAAAAFAAAAAgAAAAAAAAAAAAAAAAAAAAAAAAADAAAAAAAAAAAAAAADAAAAAAAAAAMAAAAAAAAAAAAAAAAAAAADAAAAAwAAAAAAAAAAAAAAAwAAAAAAAAADAAAAAAAAAAMAAAAAAAAAAwAAAAMAAAADAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAADAAAAAAAAAAMAAAAAAAAAAwAAAAAAAAAAAAAAAwAAAAAAAAAAAAAAAwAAAAMAAAAAAAAAAwAAAAMAAAADAAAAAAAAAAAAAAADAAAAAAAAAAAAAAAAAAAAAwAAAAAAAAADAAAAAAAAAAAAAAD/////AwAAAAAAAAAFAAAAAgAAAAAAAAAAAAAAAwAAAAMAAAADAAAAAwAAAAMAAAAAAAAAAAAAAAMAAAADAAAAAwAAAAMAAAADAAAAAAAAAAAAAAADAAAAAwAAAAMAAAADAAAAAAAAAAMAAAAAAAAAAwAAAAMAAAADAAAAAwAAAAAAAAADAAAAAAAAAP////8DAAAAAAAAAAUAAAACAAAAAAAAAAAAAAAAAAAAAAAAAAMAAAAAAAAAAAAAAAMAAAAAAAAAAwAAAAMAAAADAAAAAAAAAAMAAAAAAAAAAAAAAAMAAAAAAAAAAAAAAAAAAAADAAAAAwAAAAAAAAADAAAAAAAAAAAAAAADAAAAAwAAAAAAAAAAAAAAAwAAAAMAAAADAAAAAAAAAAAAAAAAAAAAAAAAAAMAAAAAAAAAAAAAAAMAAAADAAAAAAAAAAAAAAAAAAAAAwAAAAAAAAAAAAAAAwAAAAMAAAAAAAAAAAAAAAAAAAADAAAAAAAAAAMAAAAAAAAAAAAAAP////8DAAAAAAAAAAUAAAACAAAAAAAAAAAAAAADAAAAAwAAAAMAAAAAAAAAAAAAAAMAAAAAAAAAAwAAAAMAAAADAAAAAAAAAAAAAAADAAAAAAAAAAAAAAAAAAAAAwAAAAAAAAAAAAAAAwAAAAAAAAADAAAAAAAAAAAAAAAAAAAAAwAAAAMAAAAAAAAAAAAAAAAAAAADAAAAAAAAAAUAAAAAAAAAAAAAAAMAAAADAAAAAwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAQAAAAMAAAABAAAAAAAAAAEAAAAAAAAAAAAAAAEAAAADAAAAAQAAAAAAAAABAAAAAAAAAAAAAAADAAAAAAAAAAMAAAAAAAAAAwAAAAAAAAAAAAAAAwAAAAAAAAADAAAAAAAAAAMAAAAAAAAA/////wMAAAAAAAAABQAAAAIAAAAAAAAAAAAAAAAAAAADAAAAAAAAAAAAAAADAAAAAwAAAAAAAAAAAAAAAAAAAAMAAAAAAAAAAwAAAAAAAAAAAAAAAwAAAAAAAAAAAAAAAwAAAAMAAAAAAAAAAAAAAAMAAAADAAAAAwAAAAMAAAADAAAAAAAAAAAAAAAAAAAAAAAAAAMAAAAAAAAABQAAAAAAAAAAAAAAAwAAAAMAAAADAAAAAwAAAAMAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAQAAAAAAAAADAAAAAwAAAAMAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAADAAAAAAAAAAUAAAAAAAAAAAAAAAUAAAAAAAAAAAAAAAUAAAAFAAAAAAAAAAAAAAAAAAAAAwAAAAAAAAAAAAAAAwAAAAMAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAABAAAAAAAAAAAAAAAAAAAAAwAAAAAAAAADAAAAAAAAAAAAAAD/////AwAAAAAAAAAFAAAAAgAAAAAAAAAAAAAAAwAAAAMAAAADAAAAAAAAAAAAAAADAAAAAAAAAAUAAAAAAAAAAAAAAAUAAAAFAAAAAAAAAAAAAAAAAAAAAQAAAAMAAAABAAAAAAAAAAEAAAAAAAAAAwAAAAMAAAADAAAAAAAAAAAAAAADAAAAAAAAAAMAAAADAAAAAwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAQAAAAMAAAABAAAAAAAAAAEAAAAAAAAAAwAAAAMAAAADAAAAAwAAAAMAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAQAAAAAAAAAAAAAAAQAAAAAAAAADAAAABQAAAAEAAAAAAAAA/////wMAAAAAAAAABQAAAAIAAAAAAAAAAAAAAAUAAAAAAAAAAAAAAAUAAAAFAAAAAAAAAAAAAAAAAAAAAQAAAAAAAAAEAAAABQAAAAEAAAAAAAAAAwAAAAMAAAADAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAwAAAAAAAAAFAAAAAAAAAAAAAAAAAAAAAAAAAAMAAAAAAAAABQAAAAAAAAAAAAAAAAAAAAEAAAAAAAAAAgAAAAUAAAABAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAQAAAAAAAAAAAAAAAQAAAAMAAAABAAAAAAAAAAEAAAAAAAAABQAAAAAAAAAAAAAABQAAAAUAAAAAAAAAAAAAAP////8BAAAAAAAAAAMAAAAEAAAAAgAAAAAAAAAAAAAAAQAAAAAAAAAAAAAABQAAAAEAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAABAAAAAAAAAAUAAAAAAAAAAAAAAAUAAAAFAAAAAAAAAAAAAAAAAAAAAQAAAAAAAAABAAAABQAAAAEAAAAAAAAAAAAAAAEAAAABAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACAAAAAQAAAAEAAAAAAAAAAAAAAAAAAAAAAAAAAQAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAIAAAABAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAgAAAAAAAAAAAAAAAQAAAP//////////AQAAAAEAAAABAAAAAAAAAAAAAAAAAAAAAAAAAAEAAAAAAAAAAAAAAAEAAAAAAAAAAAAAAAAAAAACAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAEAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAIAAAAAAAAAAQAAAAAAAAAAAAAAAAAAAAAAAAABAAAAAAAAAAEAAAAAAAAAAAAAAAAAAAAAAAAAAQAAAAAAAAABAAAAAQAAAAAAAAAAAAAAAAAAAAMAAAABAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAADAAAAAQAAAAEAAAAAAAAAAAAAAAAAAAAAAAAACwAAAAIAAAAAAAAAAAAAAAEAAAACAAAABgAAAAQAAAABAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAABgAAAAAAAAABAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAEAAAAAAAAAAAAAAAAAAAACAAAAAAAAAAEAAAABAAAAAAAAAAAAAAAAAAAABwAAAAAAAAAAAAAAAQAAAAAAAAAAAAAAAAAAAAIAAAAAAAAAAAAAAAEAAAAAAAAAAAAAAAAAAAAAAAAAAQAAAAEAAAAAAAAAAAAAAAAAAAAAAAAABgAAAAAAAAAAAAAAAQAAAAAAAAAAAAAAAAAAAAoAAAACAAAAAAAAAAAAAAABAAAAAQAAAAUAAAAGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAsAAAABAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAEAAAAAQAAAAEAAAAAAAAAAAAAAAAAAAAAAAAAAwAAAAAAAAABAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAQAAAAEAAAAAAAAAAAAAAAAAAAAEAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAABQAAAAAAAAABAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAQAAAAAAAAAAAAAAAAAAAAAAAAAHAAAAAAAAAAEAAAAAAAAAAAAAAAAAAAAAAAAACwAAAAEAAAABAAAAAAAAAAAAAAAAAAAAAAAAAAcAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAKAAAAAQAAAAAAAAAAAAAAAAAAAAAAAAAAAAAADAAAAAIAAAAAAAAAAAAAAAEAAAADAAAABwAAAAYAAAABAAAAAAAAAAEAAAAAAAAAAAAAAAAAAAAHAAAAAQAAAAAAAAABAAAAAAAAAAAAAAAAAAAABAAAAAAAAAAAAAAAAQAAAAAAAAAAAAAAAAAAAAMAAAAAAAAAAAAAAAEAAAAAAAAAAAAAAAAAAAADAAAAAAAAAAEAAAABAAAAAAAAAAAAAAAAAAAABAAAAAAAAAABAAAAAAAAAAAAAAAAAAAAAAAAAAYAAAABAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAALAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACAAAAAAAAAAAAAAAAQAAAAAAAAAAAAAAAAAAAAUAAAAAAAAAAAAAAAEAAAAAAAAAAAAAAAAAAAAOAAAAAgAAAAAAAAAAAAAAAQAAAAAAAAAJAAAABQAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAwAAAABAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAKAAAAAQAAAAEAAAAAAAAAAAAAAAAAAAAAAAAABAAAAAAAAAABAAAAAQAAAAAAAAAAAAAAAAAAAAwAAAABAAAAAQAAAAAAAAAAAAAAAAAAAAAAAAAHAAAAAQAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACwAAAAAAAAABAAAAAAAAAAAAAAAAAAAAAAAAAAoAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAANAAAAAgAAAAAAAAAAAAAAAQAAAAQAAAAIAAAACgAAAAAAAAAAAAAAAQAAAAAAAAAAAAAAAAAAAAsAAAAAAAAAAAAAAAEAAAAAAAAAAAAAAAAAAAAJAAAAAAAAAAEAAAAAAAAAAAAAAAAAAAAAAAAACAAAAAAAAAABAAAAAAAAAAAAAAAAAAAAAAAAAAYAAAACAAAAAAAAAAAAAAABAAAACwAAAA8AAAAIAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACQAAAAAAAAAAAAAAAQAAAAAAAAAAAAAAAAAAAA4AAAABAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAFAAAAAQAAAAAAAAABAAAAAAAAAAAAAAAAAAAAEAAAAAAAAAABAAAAAQAAAAAAAAAAAAAAAAAAAAgAAAABAAAAAAAAAAEAAAAAAAAAAAAAAAAAAAAFAAAAAQAAAAAAAAAAAAAAAAAAAAAAAAAAAAAADAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAcAAAACAAAAAAAAAAAAAAABAAAADAAAABAAAAAMAAAAAAAAAAEAAAAAAAAAAAAAAAAAAAAAAAAACgAAAAAAAAABAAAAAAAAAAAAAAAAAAAAAAAAAAkAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAANAAAAAQAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAEAAAAAAAAAAAAAAAAQAAAAAAAAAAAAAAAAAAAA8AAAAAAAAAAQAAAAEAAAAAAAAAAAAAAAAAAAAPAAAAAAAAAAEAAAAAAAAAAAAAAAAAAAAAAAAAEAAAAAAAAAABAAAAAAAAAAAAAAAAAAAAAAAAAA4AAAABAAAAAQAAAAAAAAAAAAAAAAAAAAAAAAANAAAAAQAAAAEAAAAAAAAAAAAAAAAAAAAAAAAABQAAAAIAAAAAAAAAAAAAAAEAAAAKAAAAEwAAAAgAAAABAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAOAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACQAAAAEAAAAAAAAAAQAAAAAAAAAAAAAAAAAAAA4AAAAAAAAAAAAAAAEAAAAAAAAAAAAAAAAAAAARAAAAAAAAAAAAAAABAAAAAAAAAAAAAAAAAAAADAAAAAAAAAAAAAAAAQAAAAAAAAAAAAAAAAAAABAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAARAAAAAAAAAAEAAAABAAAAAAAAAAAAAAAAAAAADwAAAAAAAAAAAAAAAQAAAAAAAAAAAAAAAAAAABAAAAABAAAAAAAAAAEAAAAAAAAAAAAAAAAAAAAJAAAAAQAAAAAAAAAAAAAAAAAAAAAAAAAAAAAADwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA0AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAIAAAAAgAAAAAAAAAAAAAAAQAAAA0AAAARAAAADQAAAAAAAAABAAAAAAAAAAAAAAAAAAAAAAAAABEAAAABAAAAAAAAAAEAAAAAAAAAAAAAAAAAAAATAAAAAAAAAAEAAAAAAAAAAAAAAAAAAAAAAAAADgAAAAAAAAABAAAAAAAAAAAAAAAAAAAAAAAAABMAAAAAAAAAAQAAAAEAAAAAAAAAAAAAAAAAAAARAAAAAAAAAAEAAAAAAAAAAAAAAAAAAAAAAAAADQAAAAAAAAAAAAAAAQAAAAAAAAAAAAAAAAAAABEAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAQAAAAAQAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACQAAAAIAAAAAAAAAAAAAAAEAAAAOAAAAEgAAAA8AAAABAAAAAAAAAAEAAAAAAAAAAAAAAAAAAAAPAAAAAQAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAEgAAAAAAAAABAAAAAQAAAAAAAAAAAAAAAAAAABIAAAAAAAAAAAAAAAEAAAAAAAAAAAAAAAAAAAATAAAAAAAAAAAAAAABAAAAAAAAAAAAAAAAAAAAEQAAAAEAAAAAAAAAAAAAAAAAAAAAAAAAAAAAABMAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAASAAAAAAAAAAEAAAAAAAAAAAAAAAAAAAAAAAAAEgAAAAEAAAAAAAAAAQAAAAAAAAAAAAAAAAAAABMAAAACAAAAAAAAAAAAAAABAAAA//////////8TAAAAAQAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAEgAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAABMAAAABAAAAAAAAAAEAAAAAAAAAAAAAAAAAAAASAAAAAQAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAQAAAAAAAAABIAAAAAAAAAGAAAAAAAAAAhAAAAAAAAAB4AAAAAAAAAIAAAAAMAAAAxAAAAAQAAADAAAAADAAAAMgAAAAMAAAAIAAAAAAAAAAUAAAAFAAAACgAAAAUAAAAWAAAAAAAAABAAAAAAAAAAEgAAAAAAAAApAAAAAQAAACEAAAAAAAAAHgAAAAAAAAAEAAAAAAAAAAAAAAAFAAAAAgAAAAUAAAAPAAAAAQAAAAgAAAAAAAAABQAAAAUAAAAfAAAAAQAAABYAAAAAAAAAEAAAAAAAAAACAAAAAAAAAAYAAAAAAAAADgAAAAAAAAAKAAAAAAAAAAsAAAAAAAAAEQAAAAMAAAAYAAAAAQAAABcAAAADAAAAGQAAAAMAAAAAAAAAAAAAAAEAAAAFAAAACQAAAAUAAAAFAAAAAAAAAAIAAAAAAAAABgAAAAAAAAASAAAAAQAAAAoAAAAAAAAACwAAAAAAAAAEAAAAAQAAAAMAAAAFAAAABwAAAAUAAAAIAAAAAQAAAAAAAAAAAAAAAQAAAAUAAAAQAAAAAQAAAAUAAAAAAAAAAgAAAAAAAAAHAAAAAAAAABUAAAAAAAAAJgAAAAAAAAAJAAAAAAAAABMAAAAAAAAAIgAAAAMAAAAOAAAAAQAAABQAAAADAAAAJAAAAAMAAAADAAAAAAAAAA0AAAAFAAAAHQAAAAUAAAABAAAAAAAAAAcAAAAAAAAAFQAAAAAAAAAGAAAAAQAAAAkAAAAAAAAAEwAAAAAAAAAEAAAAAgAAAAwAAAAFAAAAGgAAAAUAAAAAAAAAAQAAAAMAAAAAAAAADQAAAAUAAAACAAAAAQAAAAEAAAAAAAAABwAAAAAAAAAaAAAAAAAAACoAAAAAAAAAOgAAAAAAAAAdAAAAAAAAACsAAAAAAAAAPgAAAAMAAAAmAAAAAQAAAC8AAAADAAAAQAAAAAMAAAAMAAAAAAAAABwAAAAFAAAALAAAAAUAAAANAAAAAAAAABoAAAAAAAAAKgAAAAAAAAAVAAAAAQAAAB0AAAAAAAAAKwAAAAAAAAAEAAAAAwAAAA8AAAAFAAAAHwAAAAUAAAADAAAAAQAAAAwAAAAAAAAAHAAAAAUAAAAHAAAAAQAAAA0AAAAAAAAAGgAAAAAAAAAfAAAAAAAAACkAAAAAAAAAMQAAAAAAAAAsAAAAAAAAADUAAAAAAAAAPQAAAAMAAAA6AAAAAQAAAEEAAAADAAAASwAAAAMAAAAPAAAAAAAAABYAAAAFAAAAIQAAAAUAAAAcAAAAAAAAAB8AAAAAAAAAKQAAAAAAAAAqAAAAAQAAACwAAAAAAAAANQAAAAAAAAAEAAAABAAAAAgAAAAFAAAAEAAAAAUAAAAMAAAAAQAAAA8AAAAAAAAAFgAAAAUAAAAaAAAAAQAAABwAAAAAAAAAHwAAAAAAAAAyAAAAAAAAADAAAAAAAAAAMQAAAAMAAAAgAAAAAAAAAB4AAAADAAAAIQAAAAMAAAAYAAAAAwAAABIAAAADAAAAEAAAAAMAAABGAAAAAAAAAEMAAAAAAAAAQgAAAAMAAAA0AAAAAwAAADIAAAAAAAAAMAAAAAAAAAAlAAAAAwAAACAAAAAAAAAAHgAAAAMAAABTAAAAAAAAAFcAAAADAAAAVQAAAAMAAABKAAAAAwAAAEYAAAAAAAAAQwAAAAAAAAA5AAAAAQAAADQAAAADAAAAMgAAAAAAAAAZAAAAAAAAABcAAAAAAAAAGAAAAAMAAAARAAAAAAAAAAsAAAADAAAACgAAAAMAAAAOAAAAAwAAAAYAAAADAAAAAgAAAAMAAAAtAAAAAAAAACcAAAAAAAAAJQAAAAMAAAAjAAAAAwAAABkAAAAAAAAAFwAAAAAAAAAbAAAAAwAAABEAAAAAAAAACwAAAAMAAAA/AAAAAAAAADsAAAADAAAAOQAAAAMAAAA4AAAAAwAAAC0AAAAAAAAAJwAAAAAAAAAuAAAAAwAAACMAAAADAAAAGQAAAAAAAAAkAAAAAAAAABQAAAAAAAAADgAAAAMAAAAiAAAAAAAAABMAAAADAAAACQAAAAMAAAAmAAAAAwAAABUAAAADAAAABwAAAAMAAAA3AAAAAAAAACgAAAAAAAAAGwAAAAMAAAA2AAAAAwAAACQAAAAAAAAAFAAAAAAAAAAzAAAAAwAAACIAAAAAAAAAEwAAAAMAAABIAAAAAAAAADwAAAADAAAALgAAAAMAAABJAAAAAwAAADcAAAAAAAAAKAAAAAAAAABHAAAAAwAAADYAAAADAAAAJAAAAAAAAABAAAAAAAAAAC8AAAAAAAAAJgAAAAMAAAA+AAAAAAAAACsAAAADAAAAHQAAAAMAAAA6AAAAAwAAACoAAAADAAAAGgAAAAMAAABUAAAAAAAAAEUAAAAAAAAAMwAAAAMAAABSAAAAAwAAAEAAAAAAAAAALwAAAAAAAABMAAAAAwAAAD4AAAAAAAAAKwAAAAMAAABhAAAAAAAAAFkAAAADAAAARwAAAAMAAABiAAAAAwAAAFQAAAAAAAAARQAAAAAAAABgAAAAAwAAAFIAAAADAAAAQAAAAAAAAABLAAAAAAAAAEEAAAAAAAAAOgAAAAMAAAA9AAAAAAAAADUAAAADAAAALAAAAAMAAAAxAAAAAwAAACkAAAADAAAAHwAAAAMAAABeAAAAAAAAAFYAAAAAAAAATAAAAAMAAABRAAAAAwAAAEsAAAAAAAAAQQAAAAAAAABCAAAAAwAAAD0AAAAAAAAANQAAAAMAAABrAAAAAAAAAGgAAAADAAAAYAAAAAMAAABlAAAAAwAAAF4AAAAAAAAAVgAAAAAAAABVAAAAAwAAAFEAAAADAAAASwAAAAAAAAA5AAAAAAAAADsAAAAAAAAAPwAAAAMAAABKAAAAAAAAAE4AAAADAAAATwAAAAMAAABTAAAAAwAAAFwAAAADAAAAXwAAAAMAAAAlAAAAAAAAACcAAAADAAAALQAAAAMAAAA0AAAAAAAAADkAAAAAAAAAOwAAAAAAAABGAAAAAwAAAEoAAAAAAAAATgAAAAMAAAAYAAAAAAAAABcAAAADAAAAGQAAAAMAAAAgAAAAAwAAACUAAAAAAAAAJwAAAAMAAAAyAAAAAwAAADQAAAAAAAAAOQAAAAAAAAAuAAAAAAAAADwAAAAAAAAASAAAAAMAAAA4AAAAAAAAAEQAAAADAAAAUAAAAAMAAAA/AAAAAwAAAE0AAAADAAAAWgAAAAMAAAAbAAAAAAAAACgAAAADAAAANwAAAAMAAAAjAAAAAAAAAC4AAAAAAAAAPAAAAAAAAAAtAAAAAwAAADgAAAAAAAAARAAAAAMAAAAOAAAAAAAAABQAAAADAAAAJAAAAAMAAAARAAAAAwAAABsAAAAAAAAAKAAAAAMAAAAZAAAAAwAAACMAAAAAAAAALgAAAAAAAABHAAAAAAAAAFkAAAAAAAAAYQAAAAMAAABJAAAAAAAAAFsAAAADAAAAZwAAAAMAAABIAAAAAwAAAFgAAAADAAAAaQAAAAMAAAAzAAAAAAAAAEUAAAADAAAAVAAAAAMAAAA2AAAAAAAAAEcAAAAAAAAAWQAAAAAAAAA3AAAAAwAAAEkAAAAAAAAAWwAAAAMAAAAmAAAAAAAAAC8AAAADAAAAQAAAAAMAAAAiAAAAAwAAADMAAAAAAAAARQAAAAMAAAAkAAAAAwAAADYAAAAAAAAARwAAAAAAAABgAAAAAAAAAGgAAAAAAAAAawAAAAMAAABiAAAAAAAAAG4AAAADAAAAcwAAAAMAAABhAAAAAwAAAG8AAAADAAAAdwAAAAMAAABMAAAAAAAAAFYAAAADAAAAXgAAAAMAAABSAAAAAAAAAGAAAAAAAAAAaAAAAAAAAABUAAAAAwAAAGIAAAAAAAAAbgAAAAMAAAA6AAAAAAAAAEEAAAADAAAASwAAAAMAAAA+AAAAAwAAAEwAAAAAAAAAVgAAAAMAAABAAAAAAwAAAFIAAAAAAAAAYAAAAAAAAABVAAAAAAAAAFcAAAAAAAAAUwAAAAMAAABlAAAAAAAAAGYAAAADAAAAZAAAAAMAAABrAAAAAwAAAHAAAAADAAAAcgAAAAMAAABCAAAAAAAAAEMAAAADAAAARgAAAAMAAABRAAAAAAAAAFUAAAAAAAAAVwAAAAAAAABeAAAAAwAAAGUAAAAAAAAAZgAAAAMAAAAxAAAAAAAAADAAAAADAAAAMgAAAAMAAAA9AAAAAwAAAEIAAAAAAAAAQwAAAAMAAABLAAAAAwAAAFEAAAAAAAAAVQAAAAAAAABfAAAAAAAAAFwAAAAAAAAAUwAAAAAAAABPAAAAAAAAAE4AAAAAAAAASgAAAAMAAAA/AAAAAQAAADsAAAADAAAAOQAAAAMAAABtAAAAAAAAAGwAAAAAAAAAZAAAAAUAAABdAAAAAQAAAF8AAAAAAAAAXAAAAAAAAABNAAAAAQAAAE8AAAAAAAAATgAAAAAAAAB1AAAABAAAAHYAAAAFAAAAcgAAAAUAAABqAAAAAQAAAG0AAAAAAAAAbAAAAAAAAABaAAAAAQAAAF0AAAABAAAAXwAAAAAAAABaAAAAAAAAAE0AAAAAAAAAPwAAAAAAAABQAAAAAAAAAEQAAAAAAAAAOAAAAAMAAABIAAAAAQAAADwAAAADAAAALgAAAAMAAABqAAAAAAAAAF0AAAAAAAAATwAAAAUAAABjAAAAAQAAAFoAAAAAAAAATQAAAAAAAABYAAAAAQAAAFAAAAAAAAAARAAAAAAAAAB1AAAAAwAAAG0AAAAFAAAAXwAAAAUAAABxAAAAAQAAAGoAAAAAAAAAXQAAAAAAAABpAAAAAQAAAGMAAAABAAAAWgAAAAAAAABpAAAAAAAAAFgAAAAAAAAASAAAAAAAAABnAAAAAAAAAFsAAAAAAAAASQAAAAMAAABhAAAAAQAAAFkAAAADAAAARwAAAAMAAABxAAAAAAAAAGMAAAAAAAAAUAAAAAUAAAB0AAAAAQAAAGkAAAAAAAAAWAAAAAAAAABvAAAAAQAAAGcAAAAAAAAAWwAAAAAAAAB1AAAAAgAAAGoAAAAFAAAAWgAAAAUAAAB5AAAAAQAAAHEAAAAAAAAAYwAAAAAAAAB3AAAAAQAAAHQAAAABAAAAaQAAAAAAAAB3AAAAAAAAAG8AAAAAAAAAYQAAAAAAAABzAAAAAAAAAG4AAAAAAAAAYgAAAAMAAABrAAAAAQAAAGgAAAADAAAAYAAAAAMAAAB5AAAAAAAAAHQAAAAAAAAAZwAAAAUAAAB4AAAAAQAAAHcAAAAAAAAAbwAAAAAAAABwAAAAAQAAAHMAAAAAAAAAbgAAAAAAAAB1AAAAAQAAAHEAAAAFAAAAaQAAAAUAAAB2AAAAAQAAAHkAAAAAAAAAdAAAAAAAAAByAAAAAQAAAHgAAAABAAAAdwAAAAAAAAByAAAAAAAAAHAAAAAAAAAAawAAAAAAAABkAAAAAAAAAGYAAAAAAAAAZQAAAAMAAABTAAAAAQAAAFcAAAADAAAAVQAAAAMAAAB2AAAAAAAAAHgAAAAAAAAAcwAAAAUAAABsAAAAAQAAAHIAAAAAAAAAcAAAAAAAAABcAAAAAQAAAGQAAAAAAAAAZgAAAAAAAAB1AAAAAAAAAHkAAAAFAAAAdwAAAAUAAABtAAAAAQAAAHYAAAAAAAAAeAAAAAAAAABfAAAAAQAAAGwAAAABAAAAcgAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAEAAAAAAAAAAQAAAAAAAAAAAAAAAQAAAAEAAAABAAAAAAAAAAAAAAABAAAAAAAAAAEAAAABAAAAAQAAAAAAAAAAAAAAAAAAAAAAAAB+ogX28rbpPxqumpJv+fM/165tC4ns9D+XaEnTqUsEQFrOtNlC4PA/3U+0XG6P9b9TdUUBxTTjP4PUp8ex1ty/B1rD/EN43z+lcDi6LLrZP/a45NWEHMY/oJ5ijLDZ+j/xw3rjxWPjP2B8A46ioQdAotff3wla2z+FMSpA1jj+v6b5Y1mtPbS/cIu8K0F457/2esiyJpDNv98k5Ts2NeA/pvljWa09tD88ClUJ60MDQPZ6yLImkM0/4ONKxa0UBcD2uOTVhBzGv5G7JRxGave/8cN648Vj47+HCwtkjAXIv6LX398JWtu/qyheaCAL9D9TdUUBxTTjv4gyTxslhwVAB1rD/EN4378EH/28teoFwH6iBfbytum/F6ztFYdK/r/Xrm0Liez0vwcS6wNGWeO/Ws602ULg8L9TCtRLiLT8P8pi5RexJsw/BlIKPVwR5T95Wyu0/QjnP5PjoT7YYcu/mBhKZ6zrwj8wRYS7NebuP3qW6geh+Ls/SLrixebL3r+pcyymN9XrPwmkNHp7xec/GWNMZVAA17+82s+x2BLiPwn2ytbJ9ek/LgEH1sMS1j8yp/2LhTfeP+SnWwtQBbu/d38gkp5X7z8ytsuHaADGPzUYObdf1+m/7IauECWhwz+cjSACjzniP76Z+wUhN9K/1+GEKzup67+/GYr/04baPw6idWOvsuc/ZedTWsRa5b/EJQOuRzi0v/OncYhHPes/h49PixY53j+i8wWfC03Nvw2idWOvsue/ZedTWsRa5T/EJQOuRzi0P/KncYhHPeu/iY9PixY53r+i8wWfC03NP9anWwtQBbs/d38gkp5X778ytsuHaADGvzUYObdf1+k/74auECWhw7+cjSACjzniv8CZ+wUhN9I/1uGEKzup6z+/GYr/04bavwmkNHp7xee/F2NMZVAA1z+82s+x2BLivwr2ytbJ9em/KwEH1sMS1r8yp/2LhTfev81i5RexJsy/BlIKPVwR5b95Wyu0/Qjnv5DjoT7YYcs/nBhKZ6zrwr8wRYS7Nebuv3OW6geh+Lu/SLrixebL3j+pcyymN9Xrv8rHIFfWehZAMBwUdlo0DECTUc17EOb2PxpVB1SWChdAzjbhb9pTDUDQhmdvECX5P9FlMKCC9+g/IIAzjELgE0DajDngMv8GQFhWDmDPjNs/y1guLh96EkAxPi8k7DIEQJCc4URlhRhA3eLKKLwkEECqpNAyTBD/P6xpjXcDiwVAFtl//cQm4z+Ibt3XKiYTQM7mCLUb3QdAoM1t8yVv7D8aLZv2Nk8UQEAJPV5nQwxAtSsfTCoE9z9TPjXLXIIWQBVanC5W9AtAYM3d7Adm9j++5mQz1FoWQBUThyaVBghAwH5muQsV7T89Q1qv82MUQJoWGOfNuBdAzrkClkmwDkDQjKq77t37Py+g0dtitsE/ZwAMTwVPEUBojepluNwBQGYbtuW+t9w/HNWIJs6MEkDTNuQUSlgEQKxktPP5TcQ/ixbLB8JjEUCwuWjXMQYCQAS/R09FkRdAowpiZjhhDkB7LmlczD/7P01iQmhhsAVAnrtTwDy84z/Z6jfQ2TgTQChOCXMnWwpAhrW3daoz8z/HYJvVPI4VQLT3ik5FcA5Angi7LOZd+z+NNVzDy5gXQBXdvVTFUA1AYNMgOeYe+T8+qHXGCwkXQKQTOKwa5AJA8gFVoEMW0T+FwzJyttIRQAEAAAD/////BwAAAP////8xAAAA/////1cBAAD/////YQkAAP////+nQQAA/////5HLAQD/////95AMAP/////B9lcAAAAAAAAAAAAAAAAAAgAAAP////8OAAAA/////2IAAAD/////rgIAAP/////CEgAA/////06DAAD/////IpcDAP/////uIRkA/////4LtrwAAAAAAAAAAAAAAAAAAAAAAAgAAAP//////////AQAAAAMAAAD//////////////////////////////////////////////////////////////////////////wEAAAAAAAAAAgAAAP///////////////wMAAAD//////////////////////////////////////////////////////////////////////////wEAAAAAAAAAAgAAAP///////////////wMAAAD//////////////////////////////////////////////////////////////////////////wEAAAAAAAAAAgAAAP///////////////wMAAAD//////////////////////////////////////////////////////////wIAAAD//////////wEAAAAAAAAA/////////////////////wMAAAD/////////////////////////////////////////////////////AwAAAP////////////////////8AAAAA/////////////////////wEAAAD///////////////8CAAAA////////////////////////////////AwAAAP////////////////////8AAAAA////////////////AgAAAAEAAAD/////////////////////////////////////////////////////AwAAAP////////////////////8AAAAA////////////////AgAAAAEAAAD/////////////////////////////////////////////////////AwAAAP////////////////////8AAAAA////////////////AgAAAAEAAAD/////////////////////////////////////////////////////AwAAAP////////////////////8AAAAA////////////////AgAAAAEAAAD/////////////////////////////////////////////////////AQAAAAIAAAD///////////////8AAAAA/////////////////////wMAAAD/////////////////////////////////////////////////////AQAAAAIAAAD///////////////8AAAAA/////////////////////wMAAAD/////////////////////////////////////////////////////AQAAAAIAAAD///////////////8AAAAA/////////////////////wMAAAD/////////////////////////////////////////////////////AQAAAAIAAAD///////////////8AAAAA/////////////////////wMAAAD///////////////////////////////8CAAAA////////////////AQAAAP////////////////////8AAAAA/////////////////////wMAAAD/////////////////////////////////////////////////////AwAAAP////////////////////8AAAAAAQAAAP//////////AgAAAP//////////////////////////////////////////////////////////AwAAAP///////////////wIAAAAAAAAAAQAAAP//////////////////////////////////////////////////////////////////////////AwAAAP///////////////wIAAAAAAAAAAQAAAP//////////////////////////////////////////////////////////////////////////AwAAAP///////////////wIAAAAAAAAAAQAAAP//////////////////////////////////////////////////////////////////////////AwAAAAEAAAD//////////wIAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAEAAAAAgAAAAAAAAACAAAAAQAAAAEAAAACAAAAAgAAAAAAAAAFAAAABQAAAAAAAAACAAAAAgAAAAMAAAABAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACAAAAAAAAAAIAAAABAAAAAgAAAAIAAAACAAAAAAAAAAUAAAAGAAAAAAAAAAIAAAACAAAAAwAAAAIAAAAAAAAAAAAAAAAAAAAAAAAAAQAAAAIAAAAAAAAAAgAAAAEAAAADAAAAAgAAAAIAAAAAAAAABQAAAAcAAAAAAAAAAgAAAAIAAAADAAAAAwAAAAAAAAAAAAAAAAAAAAAAAAACAAAAAgAAAAAAAAACAAAAAQAAAAQAAAACAAAAAgAAAAAAAAAFAAAACAAAAAAAAAACAAAAAgAAAAMAAAAEAAAAAAAAAAAAAAAAAAAAAAAAAAMAAAACAAAAAAAAAAIAAAABAAAAAAAAAAIAAAACAAAAAAAAAAUAAAAJAAAAAAAAAAIAAAACAAAAAwAAAAUAAAAAAAAAAAAAAAAAAAAAAAAACgAAAAIAAAACAAAAAAAAAAMAAAAOAAAAAgAAAAAAAAACAAAAAwAAAAAAAAAAAAAAAgAAAAIAAAADAAAABgAAAAAAAAAAAAAAAAAAAAAAAAALAAAAAgAAAAIAAAAAAAAAAwAAAAoAAAACAAAAAAAAAAIAAAADAAAAAQAAAAAAAAACAAAAAgAAAAMAAAAHAAAAAAAAAAAAAAAAAAAAAAAAAAwAAAACAAAAAgAAAAAAAAADAAAACwAAAAIAAAAAAAAAAgAAAAMAAAACAAAAAAAAAAIAAAACAAAAAwAAAAgAAAAAAAAAAAAAAAAAAAAAAAAADQAAAAIAAAACAAAAAAAAAAMAAAAMAAAAAgAAAAAAAAACAAAAAwAAAAMAAAAAAAAAAgAAAAIAAAADAAAACQAAAAAAAAAAAAAAAAAAAAAAAAAOAAAAAgAAAAIAAAAAAAAAAwAAAA0AAAACAAAAAAAAAAIAAAADAAAABAAAAAAAAAACAAAAAgAAAAMAAAAKAAAAAAAAAAAAAAAAAAAAAAAAAAUAAAACAAAAAgAAAAAAAAADAAAABgAAAAIAAAAAAAAAAgAAAAMAAAAPAAAAAAAAAAIAAAACAAAAAwAAAAsAAAAAAAAAAAAAAAAAAAAAAAAABgAAAAIAAAACAAAAAAAAAAMAAAAHAAAAAgAAAAAAAAACAAAAAwAAABAAAAAAAAAAAgAAAAIAAAADAAAADAAAAAAAAAAAAAAAAAAAAAAAAAAHAAAAAgAAAAIAAAAAAAAAAwAAAAgAAAACAAAAAAAAAAIAAAADAAAAEQAAAAAAAAACAAAAAgAAAAMAAAANAAAAAAAAAAAAAAAAAAAAAAAAAAgAAAACAAAAAgAAAAAAAAADAAAACQAAAAIAAAAAAAAAAgAAAAMAAAASAAAAAAAAAAIAAAACAAAAAwAAAA4AAAAAAAAAAAAAAAAAAAAAAAAACQAAAAIAAAACAAAAAAAAAAMAAAAFAAAAAgAAAAAAAAACAAAAAwAAABMAAAAAAAAAAgAAAAIAAAADAAAADwAAAAAAAAAAAAAAAAAAAAAAAAAQAAAAAgAAAAAAAAACAAAAAQAAABMAAAACAAAAAgAAAAAAAAAFAAAACgAAAAAAAAACAAAAAgAAAAMAAAAQAAAAAAAAAAAAAAAAAAAAAAAAABEAAAACAAAAAAAAAAIAAAABAAAADwAAAAIAAAACAAAAAAAAAAUAAAALAAAAAAAAAAIAAAACAAAAAwAAABEAAAAAAAAAAAAAAAAAAAAAAAAAEgAAAAIAAAAAAAAAAgAAAAEAAAAQAAAAAgAAAAIAAAAAAAAABQAAAAwAAAAAAAAAAgAAAAIAAAADAAAAEgAAAAAAAAAAAAAAAAAAAAAAAAATAAAAAgAAAAAAAAACAAAAAQAAABEAAAACAAAAAgAAAAAAAAAFAAAADQAAAAAAAAACAAAAAgAAAAMAAAATAAAAAAAAAAAAAAAAAAAAAAAAAA8AAAACAAAAAAAAAAIAAAABAAAAEgAAAAIAAAACAAAAAAAAAAUAAAAOAAAAAAAAAAIAAAACAAAAAwAAAAIAAAABAAAAAAAAAAEAAAACAAAAAAAAAAAAAAACAAAAAQAAAAAAAAABAAAAAgAAAAEAAAAAAAAAAgAAAAAAAAAFAAAABAAAAAAAAAABAAAABQAAAAAAAAAAAAAABQAAAAQAAAAAAAAAAQAAAAUAAAAEAAAAAAAAAAUAAAAAAAAAAgAAAAEAAAAAAAAAAQAAAAIAAAAAAAAAAAAAAAIAAAABAAAAAAAAAAEAAAACAAAAAQAAAAAAAAACAAAAAgAAAAAAAAABAAAAAAAAAAAAAAAFAAAABAAAAAAAAAABAAAABQAAAAAAAAAAAAAABQAAAAQAAAAAAAAAAQAAAAUAAAAEAAAAAAAAAAUAAAAFAAAAAAAAAAEAAAAAAAAAAAAAAMuhRbbsNlBBYqHW9OmHIkF9XBuqnS31QAK37uYhNMhAOSo3UUupm0DC+6pc6JxvQHV9eseEEEJAzURsCyqlFEB8BQ4NMJjnPyy3tBoS97o/xawXQznRjj89J2K2CZxhP6vX43RIIDQ/S8isgygEBz+LvFHQkmzaPjFFFO7wMq4+AADMLkTtjkIAAOgkJqxhQgAAU7B0MjRCAADwpBcVB0IAAACYP2HaQQAAAIn/Ja5BzczM4Eg6gUHNzMxMU7BTQTMzMzNfgCZBAAAAAEi3+UAAAAAAwGPNQDMzMzMzy6BAmpmZmZkxc0AzMzMzM/NFQDMzMzMzMxlAzczMzMzM7D+ygXSx2U6RQKimJOvQKnpA23hmONTHY0A/AGcxyudNQNb3K647mzZA+S56rrwWIUAm4kUQ+9UJQKre9hGzh/M/BLvoy9WG3T+LmqMf8VHGP2m3nYNV37A/gbFHcyeCmT+cBPWBckiDP61tZACjKW0/q2RbYVUYVj8uDypVyLNAP6jGS5cA5zBBwcqhBdCNGUEGEhQ/JVEDQT6WPnRbNO1AB/AWSJgT1kDfUWNCNLDAQNk+5C33OqlAchWL34QSk0DKvtDIrNV8QNF0G3kFzGVASSeWhBl6UED+/0mNGuk4QGjA/dm/1CJALPLPMql6DEDSHoDrwpP1P2jouzWST+A/egAAAAAAAABKAwAAAAAAAPoWAAAAAAAAyqAAAAAAAAB6ZQQAAAAAAErGHgAAAAAA+mvXAAAAAADK8+MFAAAAAHqqOykAAAAASqmhIAEAAAD6oGvkBwAAAMpm8T43AAAAes+ZuIIBAABKrDQMkwoAAPq1cFUFSgAAyvkUViUGAgAAAAAAAwAAAAYAAAACAAAABQAAAAEAAAAEAAAAAAAAAAAAAAAFAAAAAwAAAAEAAAAGAAAABAAAAAIAAAAAAAAAAAAAAP////8AAAAAAAAAAAAAAAAAAAAAAAAAAP////////////////////////////////////8AAAAA/////wAAAAAAAAAAAAAAAAEAAAAAAAAAAAAAAP////8AAAAAAAAAAAEAAAABAAAAAAAAAAAAAAD/////AAAAAAUAAAAAAAAAAAAAAAAAAAAAAAAA/////wUAAAAFAAAAAAAAAAAAAAAAAAAAAAAAAP////8AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAD/////////////////////////////////////AAAAAAEAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAABAAAAAAAAAAAAAAAAAAAAAQAAAAAAAAAAAAAABQAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAUAAAAAAAAABQAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA/////////////////////////////////////wAAAAABAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAQAAAAAAAAAAAAAAAAAAAAEAAAAAAAAAAAAAAAUAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAABAAAAAAAAAAUAAAABAAAAAQAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAP////////////////////////////////////8AAAAAAQAAAAEAAAABAAAAAQAAAAEAAAABAAAAAAAAAAEAAAAAAAAAAAAAAAAAAAABAAAAAAAAAAAAAAABAAAAAAAAAAAAAAABAAAAAQAAAAEAAAAAAAAAAQAAAAAAAAAFAAAAAQAAAAEAAAAAAAAAAAAAAAEAAAABAAAAAAAAAAEAAAABAAAAAQAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAEBAAAAAAABAAEAAAEBAAAAAAABAAAAAQAAAAEAAQAAAAAAAAAAAAAAAAAAAAAEAAAABAAAAAAAAAACAAAAAQAAAAMAAAAOAAAABgAAAAsAAAACAAAABwAAAAEAAAAYAAAABQAAAAoAAAABAAAABgAAAAAAAAAmAAAABwAAAAwAAAADAAAACAAAAAIAAAAxAAAACQAAAA4AAAAAAAAABQAAAAQAAAA6AAAACAAAAA0AAAAEAAAACQAAAAMAAAA/AAAACwAAAAYAAAAPAAAACgAAABAAAABIAAAADAAAAAcAAAAQAAAACwAAABEAAABTAAAACgAAAAUAAAATAAAADgAAAA8AAABhAAAADQAAAAgAAAARAAAADAAAABIAAABrAAAADgAAAAkAAAASAAAADQAAABMAAAB1AAAADwAAABMAAAARAAAAEgAAABAAAAAHAAAABwAAAAEAAAACAAAABAAAAAMAAAAAAAAAAAAAAAcAAAADAAAAAQAAAAIAAAAFAAAABAAAAAAAAAAAAAAAYWxnb3MuYwBfcG9seWZpbGxJbnRlcm5hbABhZGphY2VudEZhY2VEaXJbdG1wRmlqay5mYWNlXVtmaWprLmZhY2VdID09IEtJAGZhY2VpamsuYwBfZmFjZUlqa1BlbnRUb0dlb0JvdW5kYXJ5AGFkamFjZW50RmFjZURpcltjZW50ZXJJSksuZmFjZV1bZmFjZTJdID09IEtJAF9mYWNlSWprVG9HZW9Cb3VuZGFyeQBwb2x5Z29uLT5uZXh0ID09IE5VTEwAbGlua2VkR2VvLmMAYWRkTmV3TGlua2VkUG9seWdvbgBuZXh0ICE9IE5VTEwAbG9vcCAhPSBOVUxMAGFkZE5ld0xpbmtlZExvb3AAcG9seWdvbi0+Zmlyc3QgPT0gTlVMTABhZGRMaW5rZWRMb29wAGNvb3JkICE9IE5VTEwAYWRkTGlua2VkQ29vcmQAbG9vcC0+Zmlyc3QgPT0gTlVMTABpbm5lckxvb3BzICE9IE5VTEwAbm9ybWFsaXplTXVsdGlQb2x5Z29uAGJib3hlcyAhPSBOVUxMAGNhbmRpZGF0ZXMgIT0gTlVMTABmaW5kUG9seWdvbkZvckhvbGUAY2FuZGlkYXRlQkJveGVzICE9IE5VTEwAcmV2RGlyICE9IElOVkFMSURfRElHSVQAbG9jYWxpai5jAGgzVG9Mb2NhbElqawBiYXNlQ2VsbCAhPSBvcmlnaW5CYXNlQ2VsbAAhKG9yaWdpbk9uUGVudCAmJiBpbmRleE9uUGVudCkAcGVudGFnb25Sb3RhdGlvbnMgPj0gMABkaXJlY3Rpb25Sb3RhdGlvbnMgPj0gMABiYXNlQ2VsbCA9PSBvcmlnaW5CYXNlQ2VsbABiYXNlQ2VsbCAhPSBJTlZBTElEX0JBU0VfQ0VMTABsb2NhbElqa1RvSDMAIV9pc0Jhc2VDZWxsUGVudGFnb24oYmFzZUNlbGwpAGJhc2VDZWxsUm90YXRpb25zID49IDAAd2l0aGluUGVudGFnb25Sb3RhdGlvbnMgPj0gMABncmFwaC0+YnVja2V0cyAhPSBOVUxMAHZlcnRleEdyYXBoLmMAaW5pdFZlcnRleEdyYXBoAG5vZGUgIT0gTlVMTABhZGRWZXJ0ZXhOb2Rl";
      var tempDoublePtr = 24032;
      function demangle(func) {
        return func;
      }
      function demangleAll(text) {
        var regex = /\b__Z[\w\d_]+/g;
        return text.replace(regex, function(x) {
          var y = demangle(x);
          return x === y ? x : y + " [" + x + "]";
        });
      }
      function jsStackTrace() {
        var err2 = new Error();
        if (!err2.stack) {
          try {
            throw new Error(0);
          } catch (e) {
            err2 = e;
          }
          if (!err2.stack) {
            return "(no stack trace available)";
          }
        }
        return err2.stack.toString();
      }
      function stackTrace() {
        var js = jsStackTrace();
        if (Module["extraStackTrace"]) {
          js += "\n" + Module["extraStackTrace"]();
        }
        return demangleAll(js);
      }
      function ___assert_fail(condition, filename, line, func) {
        abort("Assertion failed: " + UTF8ToString(condition) + ", at: " + [filename ? UTF8ToString(filename) : "unknown filename", line, func ? UTF8ToString(func) : "unknown function"]);
      }
      function _emscripten_get_heap_size() {
        return HEAP8.length;
      }
      function _emscripten_memcpy_big(dest, src, num) {
        HEAPU8.set(HEAPU8.subarray(src, src + num), dest);
      }
      function ___setErrNo(value) {
        if (Module["___errno_location"]) {
          HEAP32[Module["___errno_location"]() >> 2] = value;
        }
        return value;
      }
      function abortOnCannotGrowMemory(requestedSize) {
        abort("OOM");
      }
      function emscripten_realloc_buffer(size) {
        try {
          var newBuffer = new ArrayBuffer(size);
          if (newBuffer.byteLength != size) {
            return;
          }
          new Int8Array(newBuffer).set(HEAP8);
          _emscripten_replace_memory(newBuffer);
          updateGlobalBufferAndViews(newBuffer);
          return 1;
        } catch (e) {
        }
      }
      function _emscripten_resize_heap(requestedSize) {
        var oldSize = _emscripten_get_heap_size();
        var PAGE_MULTIPLE = 16777216;
        var LIMIT = 2147483648 - PAGE_MULTIPLE;
        if (requestedSize > LIMIT) {
          return false;
        }
        var MIN_TOTAL_MEMORY = 16777216;
        var newSize = Math.max(oldSize, MIN_TOTAL_MEMORY);
        while (newSize < requestedSize) {
          if (newSize <= 536870912) {
            newSize = alignUp(2 * newSize, PAGE_MULTIPLE);
          } else {
            newSize = Math.min(alignUp((3 * newSize + 2147483648) / 4, PAGE_MULTIPLE), LIMIT);
          }
        }
        var replacement = emscripten_realloc_buffer(newSize);
        if (!replacement) {
          return false;
        }
        return true;
      }
      var decodeBase64 = typeof atob === "function" ? atob : function(input) {
        var keyStr = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/=";
        var output = "";
        var chr1, chr2, chr3;
        var enc1, enc2, enc3, enc4;
        var i = 0;
        input = input.replace(/[^A-Za-z0-9\+\/\=]/g, "");
        do {
          enc1 = keyStr.indexOf(input.charAt(i++));
          enc2 = keyStr.indexOf(input.charAt(i++));
          enc3 = keyStr.indexOf(input.charAt(i++));
          enc4 = keyStr.indexOf(input.charAt(i++));
          chr1 = enc1 << 2 | enc2 >> 4;
          chr2 = (enc2 & 15) << 4 | enc3 >> 2;
          chr3 = (enc3 & 3) << 6 | enc4;
          output = output + String.fromCharCode(chr1);
          if (enc3 !== 64) {
            output = output + String.fromCharCode(chr2);
          }
          if (enc4 !== 64) {
            output = output + String.fromCharCode(chr3);
          }
        } while (i < input.length);
        return output;
      };
      function intArrayFromBase64(s) {
        try {
          var decoded = decodeBase64(s);
          var bytes = new Uint8Array(decoded.length);
          for (var i = 0; i < decoded.length; ++i) {
            bytes[i] = decoded.charCodeAt(i);
          }
          return bytes;
        } catch (_) {
          throw new Error("Converting base64 string to bytes failed.");
        }
      }
      function tryParseAsDataURI(filename) {
        if (!isDataURI(filename)) {
          return;
        }
        return intArrayFromBase64(filename.slice(dataURIPrefix.length));
      }
      var asmGlobalArg = {
        "Math": Math,
        "Int8Array": Int8Array,
        "Int32Array": Int32Array,
        "Uint8Array": Uint8Array,
        "Float32Array": Float32Array,
        "Float64Array": Float64Array
      };
      var asmLibraryArg = {
        "a": abort,
        "b": setTempRet0,
        "c": getTempRet0,
        "d": ___assert_fail,
        "e": ___setErrNo,
        "f": _emscripten_get_heap_size,
        "g": _emscripten_memcpy_big,
        "h": _emscripten_resize_heap,
        "i": abortOnCannotGrowMemory,
        "j": demangle,
        "k": demangleAll,
        "l": emscripten_realloc_buffer,
        "m": jsStackTrace,
        "n": stackTrace,
        "o": tempDoublePtr,
        "p": DYNAMICTOP_PTR
      };
      var asm = (
        /** @suppress {uselessCode} */
        function(global, env, buffer2) {
          "almost asm";
          var a = new global.Int8Array(buffer2), b = new global.Int32Array(buffer2), c = new global.Uint8Array(buffer2), d = new global.Float32Array(buffer2), e = new global.Float64Array(buffer2), g = env.p | 0, p = global.Math.floor, q = global.Math.abs, r = global.Math.sqrt, s = global.Math.pow, t = global.Math.cos, u = global.Math.sin, v = global.Math.tan, w = global.Math.acos, x = global.Math.asin, y = global.Math.atan, z = global.Math.atan2, A = global.Math.ceil, B = global.Math.imul, C = global.Math.min, D = global.Math.clz32, F = env.b, G = env.c, H = env.d, I = env.e, J = env.f, K = env.g, L = env.h, M = env.i, S = 24048;
          function V(newBuffer) {
            a = new Int8Array(newBuffer);
            c = new Uint8Array(newBuffer);
            b = new Int32Array(newBuffer);
            d = new Float32Array(newBuffer);
            e = new Float64Array(newBuffer);
            buffer2 = newBuffer;
            return true;
          }
          function W(a2) {
            a2 = a2 | 0;
            var b2 = 0;
            b2 = S;
            S = S + a2 | 0;
            S = S + 15 & -16;
            return b2 | 0;
          }
          function X() {
            return S | 0;
          }
          function Y(a2) {
            a2 = a2 | 0;
            S = a2;
          }
          function Z(a2, b2) {
            a2 = a2 | 0;
            b2 = b2 | 0;
            S = a2;
          }
          function _(a2) {
            a2 = a2 | 0;
            return (B(a2 * 3 | 0, a2 + 1 | 0) | 0) + 1 | 0;
          }
          function $(a2, b2, c2, d2) {
            a2 = a2 | 0;
            b2 = b2 | 0;
            c2 = c2 | 0;
            d2 = d2 | 0;
            var e2 = 0, f = 0;
            if (!(ba(a2, b2, c2, d2, 0) | 0)) {
              return;
            }
            f = (B(c2 * 3 | 0, c2 + 1 | 0) | 0) + 1 | 0;
            hd(d2 | 0, 0, f << 3 | 0) | 0;
            e2 = Yc(f, 4) | 0;
            if (!e2) {
              return;
            }
            ca(a2, b2, c2, d2, e2, f, 0);
            Xc(e2);
            return;
          }
          function aa(a2, b2, c2, d2, e2) {
            a2 = a2 | 0;
            b2 = b2 | 0;
            c2 = c2 | 0;
            d2 = d2 | 0;
            e2 = e2 | 0;
            var f = 0;
            if (!(ba(a2, b2, c2, d2, e2) | 0)) {
              return;
            }
            f = (B(c2 * 3 | 0, c2 + 1 | 0) | 0) + 1 | 0;
            hd(d2 | 0, 0, f << 3 | 0) | 0;
            if (e2 | 0) {
              hd(e2 | 0, 0, f << 2 | 0) | 0;
              ca(a2, b2, c2, d2, e2, f, 0);
              return;
            }
            e2 = Yc(f, 4) | 0;
            if (!e2) {
              return;
            }
            ca(a2, b2, c2, d2, e2, f, 0);
            Xc(e2);
            return;
          }
          function ba(a2, c2, d2, e2, f) {
            a2 = a2 | 0;
            c2 = c2 | 0;
            d2 = d2 | 0;
            e2 = e2 | 0;
            f = f | 0;
            var g2 = 0, h = 0, i = 0, j = 0, k = 0, l = 0, m = 0, n = 0, o = 0;
            o = S;
            S = S + 16 | 0;
            n = o;
            g2 = e2;
            b[g2 >> 2] = a2;
            b[g2 + 4 >> 2] = c2;
            g2 = (f | 0) != 0;
            if (g2) {
              b[f >> 2] = 0;
            }
            if (Fb(a2, c2) | 0) {
              n = 1;
              S = o;
              return n | 0;
            }
            b[n >> 2] = 0;
            a: do {
              if ((d2 | 0) >= 1) {
                if (g2) {
                  k = 0;
                  l = 1;
                  m = 1;
                  h = 0;
                  g2 = a2;
                  while (1) {
                    if (!(h | k)) {
                      g2 = da(g2, c2, 4, n) | 0;
                      c2 = G() | 0;
                      if ((g2 | 0) == 0 & (c2 | 0) == 0) {
                        g2 = 2;
                        break a;
                      }
                      if (Fb(g2, c2) | 0) {
                        g2 = 1;
                        break a;
                      }
                    }
                    g2 = da(g2, c2, b[16 + (k << 2) >> 2] | 0, n) | 0;
                    c2 = G() | 0;
                    if ((g2 | 0) == 0 & (c2 | 0) == 0) {
                      g2 = 2;
                      break a;
                    }
                    a2 = e2 + (m << 3) | 0;
                    b[a2 >> 2] = g2;
                    b[a2 + 4 >> 2] = c2;
                    b[f + (m << 2) >> 2] = l;
                    h = h + 1 | 0;
                    a2 = (h | 0) == (l | 0);
                    i = k + 1 | 0;
                    j = (i | 0) == 6;
                    if (Fb(g2, c2) | 0) {
                      g2 = 1;
                      break a;
                    }
                    l = l + (j & a2 & 1) | 0;
                    if ((l | 0) > (d2 | 0)) {
                      g2 = 0;
                      break;
                    } else {
                      k = a2 ? j ? 0 : i : k;
                      m = m + 1 | 0;
                      h = a2 ? 0 : h;
                    }
                  }
                } else {
                  k = 0;
                  l = 1;
                  m = 1;
                  h = 0;
                  g2 = a2;
                  while (1) {
                    if (!(h | k)) {
                      g2 = da(g2, c2, 4, n) | 0;
                      c2 = G() | 0;
                      if ((g2 | 0) == 0 & (c2 | 0) == 0) {
                        g2 = 2;
                        break a;
                      }
                      if (Fb(g2, c2) | 0) {
                        g2 = 1;
                        break a;
                      }
                    }
                    g2 = da(g2, c2, b[16 + (k << 2) >> 2] | 0, n) | 0;
                    c2 = G() | 0;
                    if ((g2 | 0) == 0 & (c2 | 0) == 0) {
                      g2 = 2;
                      break a;
                    }
                    a2 = e2 + (m << 3) | 0;
                    b[a2 >> 2] = g2;
                    b[a2 + 4 >> 2] = c2;
                    h = h + 1 | 0;
                    a2 = (h | 0) == (l | 0);
                    i = k + 1 | 0;
                    j = (i | 0) == 6;
                    if (Fb(g2, c2) | 0) {
                      g2 = 1;
                      break a;
                    }
                    l = l + (j & a2 & 1) | 0;
                    if ((l | 0) > (d2 | 0)) {
                      g2 = 0;
                      break;
                    } else {
                      k = a2 ? j ? 0 : i : k;
                      m = m + 1 | 0;
                      h = a2 ? 0 : h;
                    }
                  }
                }
              } else {
                g2 = 0;
              }
            } while (0);
            n = g2;
            S = o;
            return n | 0;
          }
          function ca(a2, c2, d2, e2, f, g2, h) {
            a2 = a2 | 0;
            c2 = c2 | 0;
            d2 = d2 | 0;
            e2 = e2 | 0;
            f = f | 0;
            g2 = g2 | 0;
            h = h | 0;
            var i = 0, j = 0, k = 0, l = 0, m = 0, n = 0, o = 0;
            m = S;
            S = S + 16 | 0;
            l = m;
            if ((a2 | 0) == 0 & (c2 | 0) == 0) {
              S = m;
              return;
            }
            i = bd(a2 | 0, c2 | 0, g2 | 0, ((g2 | 0) < 0) << 31 >> 31 | 0) | 0;
            G() | 0;
            j = e2 + (i << 3) | 0;
            n = j;
            o = b[n >> 2] | 0;
            n = b[n + 4 >> 2] | 0;
            k = (o | 0) == (a2 | 0) & (n | 0) == (c2 | 0);
            if (!((o | 0) == 0 & (n | 0) == 0 | k)) {
              do {
                i = (i + 1 | 0) % (g2 | 0) | 0;
                j = e2 + (i << 3) | 0;
                o = j;
                n = b[o >> 2] | 0;
                o = b[o + 4 >> 2] | 0;
                k = (n | 0) == (a2 | 0) & (o | 0) == (c2 | 0);
              } while (!((n | 0) == 0 & (o | 0) == 0 | k));
            }
            i = f + (i << 2) | 0;
            if (k ? (b[i >> 2] | 0) <= (h | 0) : 0) {
              S = m;
              return;
            }
            o = j;
            b[o >> 2] = a2;
            b[o + 4 >> 2] = c2;
            b[i >> 2] = h;
            if ((h | 0) >= (d2 | 0)) {
              S = m;
              return;
            }
            o = h + 1 | 0;
            b[l >> 2] = 0;
            n = da(a2, c2, 2, l) | 0;
            ca(n, G() | 0, d2, e2, f, g2, o);
            b[l >> 2] = 0;
            n = da(a2, c2, 3, l) | 0;
            ca(n, G() | 0, d2, e2, f, g2, o);
            b[l >> 2] = 0;
            n = da(a2, c2, 1, l) | 0;
            ca(n, G() | 0, d2, e2, f, g2, o);
            b[l >> 2] = 0;
            n = da(a2, c2, 5, l) | 0;
            ca(n, G() | 0, d2, e2, f, g2, o);
            b[l >> 2] = 0;
            n = da(a2, c2, 4, l) | 0;
            ca(n, G() | 0, d2, e2, f, g2, o);
            b[l >> 2] = 0;
            n = da(a2, c2, 6, l) | 0;
            ca(n, G() | 0, d2, e2, f, g2, o);
            S = m;
            return;
          }
          function da(a2, c2, d2, e2) {
            a2 = a2 | 0;
            c2 = c2 | 0;
            d2 = d2 | 0;
            e2 = e2 | 0;
            var f = 0, g2 = 0, h = 0, i = 0, j = 0, k = 0, l = 0, m = 0, n = 0, o = 0;
            if ((b[e2 >> 2] | 0) > 0) {
              f = 0;
              do {
                d2 = Pa(d2) | 0;
                f = f + 1 | 0;
              } while ((f | 0) < (b[e2 >> 2] | 0));
            }
            i = cd(a2 | 0, c2 | 0, 45) | 0;
            G() | 0;
            j = i & 127;
            g2 = Lb(a2, c2) | 0;
            f = cd(a2 | 0, c2 | 0, 52) | 0;
            G() | 0;
            f = f & 15;
            a: do {
              if (!f) {
                h = 6;
              } else {
                while (1) {
                  m = (15 - f | 0) * 3 | 0;
                  n = cd(a2 | 0, c2 | 0, m | 0) | 0;
                  G() | 0;
                  n = n & 7;
                  o = (Rb(f) | 0) == 0;
                  f = f + -1 | 0;
                  l = dd(7, 0, m | 0) | 0;
                  c2 = c2 & ~(G() | 0);
                  m = dd(b[(o ? 464 : 48) + (n * 28 | 0) + (d2 << 2) >> 2] | 0, 0, m | 0) | 0;
                  k = G() | 0;
                  d2 = b[(o ? 672 : 256) + (n * 28 | 0) + (d2 << 2) >> 2] | 0;
                  a2 = m | a2 & ~l;
                  c2 = k | c2;
                  if (!d2) {
                    d2 = 0;
                    break a;
                  }
                  if (!f) {
                    h = 6;
                    break;
                  }
                }
              }
            } while (0);
            if ((h | 0) == 6) {
              o = b[880 + (j * 28 | 0) + (d2 << 2) >> 2] | 0;
              n = dd(o | 0, 0, 45) | 0;
              a2 = n | a2;
              c2 = G() | 0 | c2 & -1040385;
              d2 = b[4304 + (j * 28 | 0) + (d2 << 2) >> 2] | 0;
              if ((o & 127 | 0) == 127) {
                o = dd(b[880 + (j * 28 | 0) + 20 >> 2] | 0, 0, 45) | 0;
                c2 = G() | 0 | c2 & -1040385;
                d2 = b[4304 + (j * 28 | 0) + 20 >> 2] | 0;
                a2 = Nb(o | a2, c2) | 0;
                c2 = G() | 0;
                b[e2 >> 2] = (b[e2 >> 2] | 0) + 1;
              }
            }
            h = cd(a2 | 0, c2 | 0, 45) | 0;
            G() | 0;
            h = h & 127;
            b: do {
              if (!(la(h) | 0)) {
                if ((d2 | 0) > 0) {
                  f = 0;
                  do {
                    a2 = Nb(a2, c2) | 0;
                    c2 = G() | 0;
                    f = f + 1 | 0;
                  } while ((f | 0) != (d2 | 0));
                }
              } else {
                c: do {
                  if ((Lb(a2, c2) | 0) == 1) {
                    if ((j | 0) != (h | 0)) {
                      if (ra(h, b[7728 + (j * 28 | 0) >> 2] | 0) | 0) {
                        a2 = Pb(a2, c2) | 0;
                        g2 = 1;
                        c2 = G() | 0;
                        break;
                      } else {
                        a2 = Nb(a2, c2) | 0;
                        g2 = 1;
                        c2 = G() | 0;
                        break;
                      }
                    }
                    switch (g2 | 0) {
                      case 5: {
                        a2 = Pb(a2, c2) | 0;
                        c2 = G() | 0;
                        b[e2 >> 2] = (b[e2 >> 2] | 0) + 5;
                        g2 = 0;
                        break c;
                      }
                      case 3: {
                        a2 = Nb(a2, c2) | 0;
                        c2 = G() | 0;
                        b[e2 >> 2] = (b[e2 >> 2] | 0) + 1;
                        g2 = 0;
                        break c;
                      }
                      default: {
                        n = 0;
                        o = 0;
                        F(n | 0);
                        return o | 0;
                      }
                    }
                  } else {
                    g2 = 0;
                  }
                } while (0);
                if ((d2 | 0) > 0) {
                  f = 0;
                  do {
                    a2 = Mb(a2, c2) | 0;
                    c2 = G() | 0;
                    f = f + 1 | 0;
                  } while ((f | 0) != (d2 | 0));
                }
                if ((j | 0) != (h | 0)) {
                  if (!(ma(h) | 0)) {
                    if ((g2 | 0) != 0 | (Lb(a2, c2) | 0) != 5) {
                      break;
                    }
                    b[e2 >> 2] = (b[e2 >> 2] | 0) + 1;
                    break;
                  }
                  switch (i & 127) {
                    case 8:
                    case 118:
                      break b;
                    default:
                  }
                  if ((Lb(a2, c2) | 0) != 3) {
                    b[e2 >> 2] = (b[e2 >> 2] | 0) + 1;
                  }
                }
              }
            } while (0);
            b[e2 >> 2] = ((b[e2 >> 2] | 0) + d2 | 0) % 6 | 0;
            n = c2;
            o = a2;
            F(n | 0);
            return o | 0;
          }
          function ea(a2, c2, d2, e2) {
            a2 = a2 | 0;
            c2 = c2 | 0;
            d2 = d2 | 0;
            e2 = e2 | 0;
            var f = 0, g2 = 0, h = 0, i = 0, j = 0, k = 0, l = 0, m = 0;
            m = S;
            S = S + 16 | 0;
            l = m;
            if (!d2) {
              l = e2;
              b[l >> 2] = a2;
              b[l + 4 >> 2] = c2;
              l = 0;
              S = m;
              return l | 0;
            }
            b[l >> 2] = 0;
            a: do {
              if (!(Fb(a2, c2) | 0)) {
                g2 = (d2 | 0) > 0;
                if (g2) {
                  f = 0;
                  k = a2;
                  do {
                    k = da(k, c2, 4, l) | 0;
                    c2 = G() | 0;
                    if ((k | 0) == 0 & (c2 | 0) == 0) {
                      a2 = 2;
                      break a;
                    }
                    f = f + 1 | 0;
                    if (Fb(k, c2) | 0) {
                      a2 = 1;
                      break a;
                    }
                  } while ((f | 0) < (d2 | 0));
                  j = e2;
                  b[j >> 2] = k;
                  b[j + 4 >> 2] = c2;
                  j = d2 + -1 | 0;
                  if (g2) {
                    g2 = 0;
                    h = 1;
                    f = k;
                    a2 = c2;
                    do {
                      f = da(f, a2, 2, l) | 0;
                      a2 = G() | 0;
                      if ((f | 0) == 0 & (a2 | 0) == 0) {
                        a2 = 2;
                        break a;
                      }
                      i = e2 + (h << 3) | 0;
                      b[i >> 2] = f;
                      b[i + 4 >> 2] = a2;
                      h = h + 1 | 0;
                      if (Fb(f, a2) | 0) {
                        a2 = 1;
                        break a;
                      }
                      g2 = g2 + 1 | 0;
                    } while ((g2 | 0) < (d2 | 0));
                    i = 0;
                    g2 = h;
                    do {
                      f = da(f, a2, 3, l) | 0;
                      a2 = G() | 0;
                      if ((f | 0) == 0 & (a2 | 0) == 0) {
                        a2 = 2;
                        break a;
                      }
                      h = e2 + (g2 << 3) | 0;
                      b[h >> 2] = f;
                      b[h + 4 >> 2] = a2;
                      g2 = g2 + 1 | 0;
                      if (Fb(f, a2) | 0) {
                        a2 = 1;
                        break a;
                      }
                      i = i + 1 | 0;
                    } while ((i | 0) < (d2 | 0));
                    h = 0;
                    do {
                      f = da(f, a2, 1, l) | 0;
                      a2 = G() | 0;
                      if ((f | 0) == 0 & (a2 | 0) == 0) {
                        a2 = 2;
                        break a;
                      }
                      i = e2 + (g2 << 3) | 0;
                      b[i >> 2] = f;
                      b[i + 4 >> 2] = a2;
                      g2 = g2 + 1 | 0;
                      if (Fb(f, a2) | 0) {
                        a2 = 1;
                        break a;
                      }
                      h = h + 1 | 0;
                    } while ((h | 0) < (d2 | 0));
                    h = 0;
                    do {
                      f = da(f, a2, 5, l) | 0;
                      a2 = G() | 0;
                      if ((f | 0) == 0 & (a2 | 0) == 0) {
                        a2 = 2;
                        break a;
                      }
                      i = e2 + (g2 << 3) | 0;
                      b[i >> 2] = f;
                      b[i + 4 >> 2] = a2;
                      g2 = g2 + 1 | 0;
                      if (Fb(f, a2) | 0) {
                        a2 = 1;
                        break a;
                      }
                      h = h + 1 | 0;
                    } while ((h | 0) < (d2 | 0));
                    h = 0;
                    do {
                      f = da(f, a2, 4, l) | 0;
                      a2 = G() | 0;
                      if ((f | 0) == 0 & (a2 | 0) == 0) {
                        a2 = 2;
                        break a;
                      }
                      i = e2 + (g2 << 3) | 0;
                      b[i >> 2] = f;
                      b[i + 4 >> 2] = a2;
                      g2 = g2 + 1 | 0;
                      if (Fb(f, a2) | 0) {
                        a2 = 1;
                        break a;
                      }
                      h = h + 1 | 0;
                    } while ((h | 0) < (d2 | 0));
                    h = 0;
                    while (1) {
                      f = da(f, a2, 6, l) | 0;
                      a2 = G() | 0;
                      if ((f | 0) == 0 & (a2 | 0) == 0) {
                        a2 = 2;
                        break a;
                      }
                      if ((h | 0) != (j | 0)) {
                        i = e2 + (g2 << 3) | 0;
                        b[i >> 2] = f;
                        b[i + 4 >> 2] = a2;
                        if (!(Fb(f, a2) | 0)) {
                          g2 = g2 + 1 | 0;
                        } else {
                          a2 = 1;
                          break a;
                        }
                      }
                      h = h + 1 | 0;
                      if ((h | 0) >= (d2 | 0)) {
                        h = k;
                        g2 = c2;
                        break;
                      }
                    }
                  } else {
                    h = k;
                    f = k;
                    g2 = c2;
                    a2 = c2;
                  }
                } else {
                  h = e2;
                  b[h >> 2] = a2;
                  b[h + 4 >> 2] = c2;
                  h = a2;
                  f = a2;
                  g2 = c2;
                  a2 = c2;
                }
                a2 = ((h | 0) != (f | 0) | (g2 | 0) != (a2 | 0)) & 1;
              } else {
                a2 = 1;
              }
            } while (0);
            l = a2;
            S = m;
            return l | 0;
          }
          function fa(a2, c2) {
            a2 = a2 | 0;
            c2 = c2 | 0;
            var d2 = 0, e2 = 0, f = 0, g2 = 0, h = 0, i = 0;
            g2 = S;
            S = S + 48 | 0;
            f = g2 + 8 | 0;
            e2 = g2;
            i = a2;
            h = b[i + 4 >> 2] | 0;
            d2 = e2;
            b[d2 >> 2] = b[i >> 2];
            b[d2 + 4 >> 2] = h;
            vc(e2, f);
            f = ya(f, c2) | 0;
            c2 = b[e2 >> 2] | 0;
            e2 = b[a2 + 8 >> 2] | 0;
            if ((e2 | 0) <= 0) {
              i = c2;
              h = (f | 0) < (i | 0);
              i = h ? i : f;
              i = i + 12 | 0;
              S = g2;
              return i | 0;
            }
            d2 = b[a2 + 12 >> 2] | 0;
            a2 = 0;
            do {
              c2 = (b[d2 + (a2 << 3) >> 2] | 0) + c2 | 0;
              a2 = a2 + 1 | 0;
            } while ((a2 | 0) < (e2 | 0));
            i = (f | 0) < (c2 | 0);
            i = i ? c2 : f;
            i = i + 12 | 0;
            S = g2;
            return i | 0;
          }
          function ga(a2, c2, d2) {
            a2 = a2 | 0;
            c2 = c2 | 0;
            d2 = d2 | 0;
            var e2 = 0, f = 0, g2 = 0, h = 0, i = 0, j = 0;
            i = S;
            S = S + 48 | 0;
            e2 = i + 8 | 0;
            f = i;
            if (!(ha(a2, c2, d2) | 0)) {
              S = i;
              return;
            }
            j = a2;
            g2 = b[j + 4 >> 2] | 0;
            h = f;
            b[h >> 2] = b[j >> 2];
            b[h + 4 >> 2] = g2;
            vc(f, e2);
            h = ya(e2, c2) | 0;
            c2 = b[f >> 2] | 0;
            g2 = b[a2 + 8 >> 2] | 0;
            if ((g2 | 0) > 0) {
              f = b[a2 + 12 >> 2] | 0;
              e2 = 0;
              do {
                c2 = (b[f + (e2 << 3) >> 2] | 0) + c2 | 0;
                e2 = e2 + 1 | 0;
              } while ((e2 | 0) != (g2 | 0));
            }
            c2 = (h | 0) < (c2 | 0) ? c2 : h;
            if ((c2 | 0) <= -12) {
              S = i;
              return;
            }
            j = c2 + 11 | 0;
            hd(d2 | 0, 0, (((j | 0) > 0 ? j : 0) << 3) + 8 | 0) | 0;
            S = i;
            return;
          }
          function ha(a2, c2, d2) {
            a2 = a2 | 0;
            c2 = c2 | 0;
            d2 = d2 | 0;
            var e2 = 0, f = 0, g2 = 0, h = 0, i = 0, j = 0, k = 0, l = 0, m = 0, n = 0, o = 0, p2 = 0, q2 = 0, r2 = 0, s2 = 0, t2 = 0, u2 = 0, v2 = 0, w2 = 0, x2 = 0, y2 = 0, z2 = 0, A2 = 0, B2 = 0, C2 = 0, D2 = 0, E = 0, F2 = 0, I2 = 0, J2 = 0;
            J2 = S;
            S = S + 112 | 0;
            D2 = J2 + 80 | 0;
            j = J2 + 72 | 0;
            E = J2;
            F2 = J2 + 56 | 0;
            k = a2 + 8 | 0;
            I2 = Wc((b[k >> 2] << 5) + 32 | 0) | 0;
            if (!I2) {
              H(22848, 22448, 800, 22456);
            }
            wc(a2, I2);
            g2 = a2;
            e2 = b[g2 + 4 >> 2] | 0;
            i = j;
            b[i >> 2] = b[g2 >> 2];
            b[i + 4 >> 2] = e2;
            vc(j, D2);
            i = ya(D2, c2) | 0;
            e2 = b[j >> 2] | 0;
            g2 = b[k >> 2] | 0;
            if ((g2 | 0) > 0) {
              h = b[a2 + 12 >> 2] | 0;
              f = 0;
              do {
                e2 = (b[h + (f << 3) >> 2] | 0) + e2 | 0;
                f = f + 1 | 0;
              } while ((f | 0) != (g2 | 0));
            }
            i = (i | 0) < (e2 | 0) ? e2 : i;
            C2 = i + 12 | 0;
            f = Yc(C2, 8) | 0;
            l = Yc(C2, 8) | 0;
            b[D2 >> 2] = 0;
            A2 = a2;
            B2 = b[A2 + 4 >> 2] | 0;
            e2 = j;
            b[e2 >> 2] = b[A2 >> 2];
            b[e2 + 4 >> 2] = B2;
            e2 = ia(j, C2, c2, D2, f, l) | 0;
            if (e2 | 0) {
              Xc(f);
              Xc(l);
              Xc(I2);
              I2 = e2;
              S = J2;
              return I2 | 0;
            }
            a: do {
              if ((b[k >> 2] | 0) > 0) {
                g2 = a2 + 12 | 0;
                e2 = 0;
                while (1) {
                  h = ia((b[g2 >> 2] | 0) + (e2 << 3) | 0, C2, c2, D2, f, l) | 0;
                  e2 = e2 + 1 | 0;
                  if (h | 0) {
                    break;
                  }
                  if ((e2 | 0) >= (b[k >> 2] | 0)) {
                    break a;
                  }
                }
                Xc(f);
                Xc(l);
                Xc(I2);
                I2 = h;
                S = J2;
                return I2 | 0;
              }
            } while (0);
            if ((i | 0) > -12) {
              hd(l | 0, 0, ((C2 | 0) > 1 ? C2 : 1) << 3 | 0) | 0;
            }
            b: do {
              if ((b[D2 >> 2] | 0) > 0) {
                B2 = ((C2 | 0) < 0) << 31 >> 31;
                v2 = f;
                w2 = l;
                x2 = f;
                y2 = f;
                z2 = l;
                A2 = f;
                e2 = f;
                r2 = f;
                s2 = l;
                t2 = l;
                u2 = l;
                f = l;
                c: while (1) {
                  q2 = b[D2 >> 2] | 0;
                  o = 0;
                  p2 = 0;
                  g2 = 0;
                  while (1) {
                    h = E;
                    i = h + 56 | 0;
                    do {
                      b[h >> 2] = 0;
                      h = h + 4 | 0;
                    } while ((h | 0) < (i | 0));
                    c2 = v2 + (o << 3) | 0;
                    j = b[c2 >> 2] | 0;
                    c2 = b[c2 + 4 >> 2] | 0;
                    if (ba(j, c2, 1, E, 0) | 0) {
                      h = E;
                      i = h + 56 | 0;
                      do {
                        b[h >> 2] = 0;
                        h = h + 4 | 0;
                      } while ((h | 0) < (i | 0));
                      h = Yc(7, 4) | 0;
                      if (h | 0) {
                        ca(j, c2, 1, E, h, 7, 0);
                        Xc(h);
                      }
                    }
                    n = 0;
                    do {
                      m = E + (n << 3) | 0;
                      l = b[m >> 2] | 0;
                      m = b[m + 4 >> 2] | 0;
                      d: do {
                        if (!((l | 0) == 0 & (m | 0) == 0)) {
                          j = bd(l | 0, m | 0, C2 | 0, B2 | 0) | 0;
                          G() | 0;
                          h = d2 + (j << 3) | 0;
                          i = h;
                          c2 = b[i >> 2] | 0;
                          i = b[i + 4 >> 2] | 0;
                          if (!((c2 | 0) == 0 & (i | 0) == 0)) {
                            k = 0;
                            while (1) {
                              if ((k | 0) > (C2 | 0)) {
                                break c;
                              }
                              if ((c2 | 0) == (l | 0) & (i | 0) == (m | 0)) {
                                break d;
                              }
                              j = (j + 1 | 0) % (C2 | 0) | 0;
                              h = d2 + (j << 3) | 0;
                              i = h;
                              c2 = b[i >> 2] | 0;
                              i = b[i + 4 >> 2] | 0;
                              if ((c2 | 0) == 0 & (i | 0) == 0) {
                                break;
                              } else {
                                k = k + 1 | 0;
                              }
                            }
                          }
                          if (!((l | 0) == 0 & (m | 0) == 0)) {
                            Vb(l, m, F2);
                            if (xc(a2, I2, F2) | 0) {
                              k = h;
                              b[k >> 2] = l;
                              b[k + 4 >> 2] = m;
                              k = w2 + (g2 << 3) | 0;
                              b[k >> 2] = l;
                              b[k + 4 >> 2] = m;
                              g2 = g2 + 1 | 0;
                            }
                          }
                        }
                      } while (0);
                      n = n + 1 | 0;
                    } while (n >>> 0 < 7);
                    p2 = p2 + 1 | 0;
                    if ((p2 | 0) >= (q2 | 0)) {
                      break;
                    } else {
                      o = o + 1 | 0;
                    }
                  }
                  if ((q2 | 0) > 0) {
                    hd(x2 | 0, 0, q2 << 3 | 0) | 0;
                  }
                  b[D2 >> 2] = g2;
                  if ((g2 | 0) > 0) {
                    l = f;
                    m = u2;
                    n = A2;
                    o = t2;
                    p2 = s2;
                    q2 = w2;
                    f = r2;
                    u2 = e2;
                    t2 = y2;
                    s2 = x2;
                    r2 = l;
                    e2 = m;
                    A2 = z2;
                    z2 = n;
                    y2 = o;
                    x2 = p2;
                    w2 = v2;
                    v2 = q2;
                  } else {
                    break b;
                  }
                }
                Xc(y2);
                Xc(z2);
                Xc(I2);
                I2 = -1;
                S = J2;
                return I2 | 0;
              } else {
                e2 = l;
              }
            } while (0);
            Xc(I2);
            Xc(f);
            Xc(e2);
            I2 = 0;
            S = J2;
            return I2 | 0;
          }
          function ia(a2, c2, d2, f, g2, h) {
            a2 = a2 | 0;
            c2 = c2 | 0;
            d2 = d2 | 0;
            f = f | 0;
            g2 = g2 | 0;
            h = h | 0;
            var i = 0, j = 0, k = 0, l = 0, m = 0, n = 0, o = 0, p2 = 0, q2 = 0, r2 = 0, s2 = 0, t2 = 0, u2 = 0, v2 = 0, w2 = 0, x2 = 0, y2 = 0, z2 = 0, A2 = 0, B2 = 0, C2 = 0, D2 = 0, E = 0;
            C2 = S;
            S = S + 48 | 0;
            y2 = C2 + 32 | 0;
            z2 = C2 + 16 | 0;
            A2 = C2;
            i = b[a2 >> 2] | 0;
            if ((i | 0) <= 0) {
              B2 = 0;
              S = C2;
              return B2 | 0;
            }
            t2 = a2 + 4 | 0;
            u2 = y2 + 8 | 0;
            v2 = z2 + 8 | 0;
            w2 = A2 + 8 | 0;
            x2 = ((c2 | 0) < 0) << 31 >> 31;
            s2 = 0;
            a: while (1) {
              j = b[t2 >> 2] | 0;
              q2 = j + (s2 << 4) | 0;
              b[y2 >> 2] = b[q2 >> 2];
              b[y2 + 4 >> 2] = b[q2 + 4 >> 2];
              b[y2 + 8 >> 2] = b[q2 + 8 >> 2];
              b[y2 + 12 >> 2] = b[q2 + 12 >> 2];
              if ((s2 | 0) == (i + -1 | 0)) {
                b[z2 >> 2] = b[j >> 2];
                b[z2 + 4 >> 2] = b[j + 4 >> 2];
                b[z2 + 8 >> 2] = b[j + 8 >> 2];
                b[z2 + 12 >> 2] = b[j + 12 >> 2];
              } else {
                q2 = j + (s2 + 1 << 4) | 0;
                b[z2 >> 2] = b[q2 >> 2];
                b[z2 + 4 >> 2] = b[q2 + 4 >> 2];
                b[z2 + 8 >> 2] = b[q2 + 8 >> 2];
                b[z2 + 12 >> 2] = b[q2 + 12 >> 2];
              }
              q2 = za(y2, z2, d2) | 0;
              b: do {
                if ((q2 | 0) > 0) {
                  r2 = +(q2 | 0);
                  p2 = 0;
                  c: while (1) {
                    E = +(q2 - p2 | 0);
                    D2 = +(p2 | 0);
                    e[A2 >> 3] = +e[y2 >> 3] * E / r2 + +e[z2 >> 3] * D2 / r2;
                    e[w2 >> 3] = +e[u2 >> 3] * E / r2 + +e[v2 >> 3] * D2 / r2;
                    n = Sb(A2, d2) | 0;
                    o = G() | 0;
                    j = bd(n | 0, o | 0, c2 | 0, x2 | 0) | 0;
                    G() | 0;
                    i = h + (j << 3) | 0;
                    k = i;
                    l = b[k >> 2] | 0;
                    k = b[k + 4 >> 2] | 0;
                    d: do {
                      if ((l | 0) == 0 & (k | 0) == 0) {
                        B2 = 14;
                      } else {
                        m = 0;
                        while (1) {
                          if ((m | 0) > (c2 | 0)) {
                            i = 1;
                            break d;
                          }
                          if ((l | 0) == (n | 0) & (k | 0) == (o | 0)) {
                            i = 7;
                            break d;
                          }
                          j = (j + 1 | 0) % (c2 | 0) | 0;
                          i = h + (j << 3) | 0;
                          k = i;
                          l = b[k >> 2] | 0;
                          k = b[k + 4 >> 2] | 0;
                          if ((l | 0) == 0 & (k | 0) == 0) {
                            B2 = 14;
                            break;
                          } else {
                            m = m + 1 | 0;
                          }
                        }
                      }
                    } while (0);
                    if ((B2 | 0) == 14) {
                      B2 = 0;
                      if ((n | 0) == 0 & (o | 0) == 0) {
                        i = 7;
                      } else {
                        b[i >> 2] = n;
                        b[i + 4 >> 2] = o;
                        i = b[f >> 2] | 0;
                        m = g2 + (i << 3) | 0;
                        b[m >> 2] = n;
                        b[m + 4 >> 2] = o;
                        b[f >> 2] = i + 1;
                        i = 0;
                      }
                    }
                    switch (i & 7) {
                      case 7:
                      case 0:
                        break;
                      default:
                        break c;
                    }
                    p2 = p2 + 1 | 0;
                    if ((q2 | 0) <= (p2 | 0)) {
                      B2 = 8;
                      break b;
                    }
                  }
                  if (i | 0) {
                    i = -1;
                    B2 = 20;
                    break a;
                  }
                } else {
                  B2 = 8;
                }
              } while (0);
              if ((B2 | 0) == 8) {
                B2 = 0;
              }
              s2 = s2 + 1 | 0;
              i = b[a2 >> 2] | 0;
              if ((s2 | 0) >= (i | 0)) {
                i = 0;
                B2 = 20;
                break;
              }
            }
            if ((B2 | 0) == 20) {
              S = C2;
              return i | 0;
            }
            return 0;
          }
          function ja(a2, c2, d2) {
            a2 = a2 | 0;
            c2 = c2 | 0;
            d2 = d2 | 0;
            var e2 = 0, f = 0, g2 = 0, h = 0, i = 0, j = 0, k = 0;
            k = S;
            S = S + 176 | 0;
            j = k;
            if ((c2 | 0) < 1) {
              Mc(d2, 0, 0);
              S = k;
              return;
            }
            h = a2;
            h = cd(b[h >> 2] | 0, b[h + 4 >> 2] | 0, 52) | 0;
            G() | 0;
            Mc(d2, (c2 | 0) > 6 ? c2 : 6, h & 15);
            h = 0;
            do {
              e2 = a2 + (h << 3) | 0;
              Wb(b[e2 >> 2] | 0, b[e2 + 4 >> 2] | 0, j);
              e2 = b[j >> 2] | 0;
              if ((e2 | 0) > 0) {
                i = 0;
                do {
                  g2 = j + 8 + (i << 4) | 0;
                  i = i + 1 | 0;
                  e2 = j + 8 + (((i | 0) % (e2 | 0) | 0) << 4) | 0;
                  f = Rc(d2, e2, g2) | 0;
                  if (!f) {
                    Qc(d2, g2, e2) | 0;
                  } else {
                    Pc(d2, f) | 0;
                  }
                  e2 = b[j >> 2] | 0;
                } while ((i | 0) < (e2 | 0));
              }
              h = h + 1 | 0;
            } while ((h | 0) != (c2 | 0));
            S = k;
            return;
          }
          function ka(a2, c2, d2) {
            a2 = a2 | 0;
            c2 = c2 | 0;
            d2 = d2 | 0;
            var e2 = 0, f = 0, g2 = 0, h = 0;
            g2 = S;
            S = S + 32 | 0;
            e2 = g2;
            f = g2 + 16 | 0;
            ja(a2, c2, f);
            b[d2 >> 2] = 0;
            b[d2 + 4 >> 2] = 0;
            b[d2 + 8 >> 2] = 0;
            a2 = Oc(f) | 0;
            if (!a2) {
              kc(d2) | 0;
              Nc(f);
              S = g2;
              return;
            }
            do {
              c2 = hc(d2) | 0;
              do {
                ic(c2, a2) | 0;
                h = a2 + 16 | 0;
                b[e2 >> 2] = b[h >> 2];
                b[e2 + 4 >> 2] = b[h + 4 >> 2];
                b[e2 + 8 >> 2] = b[h + 8 >> 2];
                b[e2 + 12 >> 2] = b[h + 12 >> 2];
                Pc(f, a2) | 0;
                a2 = Sc(f, e2) | 0;
              } while ((a2 | 0) != 0);
              a2 = Oc(f) | 0;
            } while ((a2 | 0) != 0);
            kc(d2) | 0;
            Nc(f);
            S = g2;
            return;
          }
          function la(a2) {
            a2 = a2 | 0;
            return b[7728 + (a2 * 28 | 0) + 16 >> 2] | 0;
          }
          function ma(a2) {
            a2 = a2 | 0;
            return (a2 | 0) == 4 | (a2 | 0) == 117 | 0;
          }
          function na(a2) {
            a2 = a2 | 0;
            return b[11152 + ((b[a2 >> 2] | 0) * 216 | 0) + ((b[a2 + 4 >> 2] | 0) * 72 | 0) + ((b[a2 + 8 >> 2] | 0) * 24 | 0) + (b[a2 + 12 >> 2] << 3) >> 2] | 0;
          }
          function oa(a2) {
            a2 = a2 | 0;
            return b[11152 + ((b[a2 >> 2] | 0) * 216 | 0) + ((b[a2 + 4 >> 2] | 0) * 72 | 0) + ((b[a2 + 8 >> 2] | 0) * 24 | 0) + (b[a2 + 12 >> 2] << 3) + 4 >> 2] | 0;
          }
          function pa(a2, c2) {
            a2 = a2 | 0;
            c2 = c2 | 0;
            a2 = 7728 + (a2 * 28 | 0) | 0;
            b[c2 >> 2] = b[a2 >> 2];
            b[c2 + 4 >> 2] = b[a2 + 4 >> 2];
            b[c2 + 8 >> 2] = b[a2 + 8 >> 2];
            b[c2 + 12 >> 2] = b[a2 + 12 >> 2];
            return;
          }
          function qa(a2, c2) {
            a2 = a2 | 0;
            c2 = c2 | 0;
            var d2 = 0, e2 = 0;
            if (c2 >>> 0 > 20) {
              c2 = -1;
              return c2 | 0;
            }
            do {
              if ((b[11152 + (c2 * 216 | 0) >> 2] | 0) != (a2 | 0)) {
                if ((b[11152 + (c2 * 216 | 0) + 8 >> 2] | 0) != (a2 | 0)) {
                  if ((b[11152 + (c2 * 216 | 0) + 16 >> 2] | 0) != (a2 | 0)) {
                    if ((b[11152 + (c2 * 216 | 0) + 24 >> 2] | 0) != (a2 | 0)) {
                      if ((b[11152 + (c2 * 216 | 0) + 32 >> 2] | 0) != (a2 | 0)) {
                        if ((b[11152 + (c2 * 216 | 0) + 40 >> 2] | 0) != (a2 | 0)) {
                          if ((b[11152 + (c2 * 216 | 0) + 48 >> 2] | 0) != (a2 | 0)) {
                            if ((b[11152 + (c2 * 216 | 0) + 56 >> 2] | 0) != (a2 | 0)) {
                              if ((b[11152 + (c2 * 216 | 0) + 64 >> 2] | 0) != (a2 | 0)) {
                                if ((b[11152 + (c2 * 216 | 0) + 72 >> 2] | 0) != (a2 | 0)) {
                                  if ((b[11152 + (c2 * 216 | 0) + 80 >> 2] | 0) != (a2 | 0)) {
                                    if ((b[11152 + (c2 * 216 | 0) + 88 >> 2] | 0) != (a2 | 0)) {
                                      if ((b[11152 + (c2 * 216 | 0) + 96 >> 2] | 0) != (a2 | 0)) {
                                        if ((b[11152 + (c2 * 216 | 0) + 104 >> 2] | 0) != (a2 | 0)) {
                                          if ((b[11152 + (c2 * 216 | 0) + 112 >> 2] | 0) != (a2 | 0)) {
                                            if ((b[11152 + (c2 * 216 | 0) + 120 >> 2] | 0) != (a2 | 0)) {
                                              if ((b[11152 + (c2 * 216 | 0) + 128 >> 2] | 0) != (a2 | 0)) {
                                                if ((b[11152 + (c2 * 216 | 0) + 136 >> 2] | 0) == (a2 | 0)) {
                                                  a2 = 2;
                                                  d2 = 1;
                                                  e2 = 2;
                                                } else {
                                                  if ((b[11152 + (c2 * 216 | 0) + 144 >> 2] | 0) == (a2 | 0)) {
                                                    a2 = 0;
                                                    d2 = 2;
                                                    e2 = 0;
                                                    break;
                                                  }
                                                  if ((b[11152 + (c2 * 216 | 0) + 152 >> 2] | 0) == (a2 | 0)) {
                                                    a2 = 0;
                                                    d2 = 2;
                                                    e2 = 1;
                                                    break;
                                                  }
                                                  if ((b[11152 + (c2 * 216 | 0) + 160 >> 2] | 0) == (a2 | 0)) {
                                                    a2 = 0;
                                                    d2 = 2;
                                                    e2 = 2;
                                                    break;
                                                  }
                                                  if ((b[11152 + (c2 * 216 | 0) + 168 >> 2] | 0) == (a2 | 0)) {
                                                    a2 = 1;
                                                    d2 = 2;
                                                    e2 = 0;
                                                    break;
                                                  }
                                                  if ((b[11152 + (c2 * 216 | 0) + 176 >> 2] | 0) == (a2 | 0)) {
                                                    a2 = 1;
                                                    d2 = 2;
                                                    e2 = 1;
                                                    break;
                                                  }
                                                  if ((b[11152 + (c2 * 216 | 0) + 184 >> 2] | 0) == (a2 | 0)) {
                                                    a2 = 1;
                                                    d2 = 2;
                                                    e2 = 2;
                                                    break;
                                                  }
                                                  if ((b[11152 + (c2 * 216 | 0) + 192 >> 2] | 0) == (a2 | 0)) {
                                                    a2 = 2;
                                                    d2 = 2;
                                                    e2 = 0;
                                                    break;
                                                  }
                                                  if ((b[11152 + (c2 * 216 | 0) + 200 >> 2] | 0) == (a2 | 0)) {
                                                    a2 = 2;
                                                    d2 = 2;
                                                    e2 = 1;
                                                    break;
                                                  }
                                                  if ((b[11152 + (c2 * 216 | 0) + 208 >> 2] | 0) == (a2 | 0)) {
                                                    a2 = 2;
                                                    d2 = 2;
                                                    e2 = 2;
                                                    break;
                                                  } else {
                                                    a2 = -1;
                                                  }
                                                  return a2 | 0;
                                                }
                                              } else {
                                                a2 = 2;
                                                d2 = 1;
                                                e2 = 1;
                                              }
                                            } else {
                                              a2 = 2;
                                              d2 = 1;
                                              e2 = 0;
                                            }
                                          } else {
                                            a2 = 1;
                                            d2 = 1;
                                            e2 = 2;
                                          }
                                        } else {
                                          a2 = 1;
                                          d2 = 1;
                                          e2 = 1;
                                        }
                                      } else {
                                        a2 = 1;
                                        d2 = 1;
                                        e2 = 0;
                                      }
                                    } else {
                                      a2 = 0;
                                      d2 = 1;
                                      e2 = 2;
                                    }
                                  } else {
                                    a2 = 0;
                                    d2 = 1;
                                    e2 = 1;
                                  }
                                } else {
                                  a2 = 0;
                                  d2 = 1;
                                  e2 = 0;
                                }
                              } else {
                                a2 = 2;
                                d2 = 0;
                                e2 = 2;
                              }
                            } else {
                              a2 = 2;
                              d2 = 0;
                              e2 = 1;
                            }
                          } else {
                            a2 = 2;
                            d2 = 0;
                            e2 = 0;
                          }
                        } else {
                          a2 = 1;
                          d2 = 0;
                          e2 = 2;
                        }
                      } else {
                        a2 = 1;
                        d2 = 0;
                        e2 = 1;
                      }
                    } else {
                      a2 = 1;
                      d2 = 0;
                      e2 = 0;
                    }
                  } else {
                    a2 = 0;
                    d2 = 0;
                    e2 = 2;
                  }
                } else {
                  a2 = 0;
                  d2 = 0;
                  e2 = 1;
                }
              } else {
                a2 = 0;
                d2 = 0;
                e2 = 0;
              }
            } while (0);
            c2 = b[11152 + (c2 * 216 | 0) + (d2 * 72 | 0) + (a2 * 24 | 0) + (e2 << 3) + 4 >> 2] | 0;
            return c2 | 0;
          }
          function ra(a2, c2) {
            a2 = a2 | 0;
            c2 = c2 | 0;
            if ((b[7728 + (a2 * 28 | 0) + 20 >> 2] | 0) == (c2 | 0)) {
              c2 = 1;
              return c2 | 0;
            }
            c2 = (b[7728 + (a2 * 28 | 0) + 24 >> 2] | 0) == (c2 | 0);
            return c2 | 0;
          }
          function sa(a2, c2) {
            a2 = a2 | 0;
            c2 = c2 | 0;
            return b[880 + (a2 * 28 | 0) + (c2 << 2) >> 2] | 0;
          }
          function ta(a2, c2) {
            a2 = a2 | 0;
            c2 = c2 | 0;
            if ((b[880 + (a2 * 28 | 0) >> 2] | 0) == (c2 | 0)) {
              c2 = 0;
              return c2 | 0;
            }
            if ((b[880 + (a2 * 28 | 0) + 4 >> 2] | 0) == (c2 | 0)) {
              c2 = 1;
              return c2 | 0;
            }
            if ((b[880 + (a2 * 28 | 0) + 8 >> 2] | 0) == (c2 | 0)) {
              c2 = 2;
              return c2 | 0;
            }
            if ((b[880 + (a2 * 28 | 0) + 12 >> 2] | 0) == (c2 | 0)) {
              c2 = 3;
              return c2 | 0;
            }
            if ((b[880 + (a2 * 28 | 0) + 16 >> 2] | 0) == (c2 | 0)) {
              c2 = 4;
              return c2 | 0;
            }
            if ((b[880 + (a2 * 28 | 0) + 20 >> 2] | 0) == (c2 | 0)) {
              c2 = 5;
              return c2 | 0;
            } else {
              return ((b[880 + (a2 * 28 | 0) + 24 >> 2] | 0) == (c2 | 0) ? 6 : 7) | 0;
            }
            return 0;
          }
          function ua() {
            return 122;
          }
          function va(a2) {
            a2 = a2 | 0;
            var c2 = 0, d2 = 0, e2 = 0;
            c2 = 0;
            do {
              dd(c2 | 0, 0, 45) | 0;
              e2 = G() | 0 | 134225919;
              d2 = a2 + (c2 << 3) | 0;
              b[d2 >> 2] = -1;
              b[d2 + 4 >> 2] = e2;
              c2 = c2 + 1 | 0;
            } while ((c2 | 0) != 122);
            return;
          }
          function wa(a2) {
            a2 = a2 | 0;
            return +e[a2 + 16 >> 3] < +e[a2 + 24 >> 3] | 0;
          }
          function xa(a2, b2) {
            a2 = a2 | 0;
            b2 = b2 | 0;
            var c2 = 0, d2 = 0, f = 0;
            c2 = +e[b2 >> 3];
            if (!(c2 >= +e[a2 + 8 >> 3])) {
              b2 = 0;
              return b2 | 0;
            }
            if (!(c2 <= +e[a2 >> 3])) {
              b2 = 0;
              return b2 | 0;
            }
            d2 = +e[a2 + 16 >> 3];
            c2 = +e[a2 + 24 >> 3];
            f = +e[b2 + 8 >> 3];
            b2 = f >= c2;
            a2 = f <= d2 & 1;
            if (d2 < c2) {
              if (b2) {
                a2 = 1;
              }
            } else if (!b2) {
              a2 = 0;
            }
            b2 = (a2 | 0) != 0;
            return b2 | 0;
          }
          function ya(a2, c2) {
            a2 = a2 | 0;
            c2 = c2 | 0;
            var d2 = 0, f = 0, g2 = 0, h = 0, i = 0, j = 0, k = 0, l = 0;
            i = S;
            S = S + 288 | 0;
            d2 = i + 264 | 0;
            f = i + 96 | 0;
            g2 = i;
            h = g2;
            j = h + 96 | 0;
            do {
              b[h >> 2] = 0;
              h = h + 4 | 0;
            } while ((h | 0) < (j | 0));
            _b(c2, g2);
            h = g2;
            j = b[h >> 2] | 0;
            h = b[h + 4 >> 2] | 0;
            Vb(j, h, d2);
            Wb(j, h, f);
            k = +jb(d2, f + 8 | 0);
            e[d2 >> 3] = +e[a2 >> 3];
            h = d2 + 8 | 0;
            e[h >> 3] = +e[a2 + 16 >> 3];
            e[f >> 3] = +e[a2 + 8 >> 3];
            j = f + 8 | 0;
            e[j >> 3] = +e[a2 + 24 >> 3];
            l = +jb(d2, f);
            j = ~~+A(+(l * l / +ed(+ +q(+((+e[h >> 3] - +e[j >> 3]) / (+e[d2 >> 3] - +e[f >> 3]))), 3) / (k * (k * 2.59807621135) * 0.8)));
            S = i;
            return ((j | 0) == 0 ? 1 : j) | 0;
          }
          function za(a2, c2, d2) {
            a2 = a2 | 0;
            c2 = c2 | 0;
            d2 = d2 | 0;
            var e2 = 0, f = 0, g2 = 0, h = 0, i = 0, j = 0, k = 0;
            i = S;
            S = S + 288 | 0;
            e2 = i + 264 | 0;
            f = i + 96 | 0;
            g2 = i;
            h = g2;
            j = h + 96 | 0;
            do {
              b[h >> 2] = 0;
              h = h + 4 | 0;
            } while ((h | 0) < (j | 0));
            _b(d2, g2);
            j = g2;
            h = b[j >> 2] | 0;
            j = b[j + 4 >> 2] | 0;
            Vb(h, j, e2);
            Wb(h, j, f);
            k = +jb(e2, f + 8 | 0);
            j = ~~+A(+(+jb(a2, c2) / (k * 2)));
            S = i;
            return ((j | 0) == 0 ? 1 : j) | 0;
          }
          function Aa(a2, c2, d2, e2) {
            a2 = a2 | 0;
            c2 = c2 | 0;
            d2 = d2 | 0;
            e2 = e2 | 0;
            b[a2 >> 2] = c2;
            b[a2 + 4 >> 2] = d2;
            b[a2 + 8 >> 2] = e2;
            return;
          }
          function Ba(a2, c2) {
            a2 = a2 | 0;
            c2 = c2 | 0;
            var d2 = 0, f = 0, g2 = 0, h = 0, i = 0, j = 0, k = 0, l = 0, m = 0, n = 0, o = 0;
            n = c2 + 8 | 0;
            b[n >> 2] = 0;
            k = +e[a2 >> 3];
            i = +q(+k);
            l = +e[a2 + 8 >> 3];
            j = +q(+l) / 0.8660254037844386;
            i = i + j * 0.5;
            d2 = ~~i;
            a2 = ~~j;
            i = i - +(d2 | 0);
            j = j - +(a2 | 0);
            do {
              if (i < 0.5) {
                if (i < 0.3333333333333333) {
                  b[c2 >> 2] = d2;
                  if (j < (i + 1) * 0.5) {
                    b[c2 + 4 >> 2] = a2;
                    break;
                  } else {
                    a2 = a2 + 1 | 0;
                    b[c2 + 4 >> 2] = a2;
                    break;
                  }
                } else {
                  o = 1 - i;
                  a2 = (!(j < o) & 1) + a2 | 0;
                  b[c2 + 4 >> 2] = a2;
                  if (o <= j & j < i * 2) {
                    d2 = d2 + 1 | 0;
                    b[c2 >> 2] = d2;
                    break;
                  } else {
                    b[c2 >> 2] = d2;
                    break;
                  }
                }
              } else {
                if (!(i < 0.6666666666666666)) {
                  d2 = d2 + 1 | 0;
                  b[c2 >> 2] = d2;
                  if (j < i * 0.5) {
                    b[c2 + 4 >> 2] = a2;
                    break;
                  } else {
                    a2 = a2 + 1 | 0;
                    b[c2 + 4 >> 2] = a2;
                    break;
                  }
                }
                if (j < 1 - i) {
                  b[c2 + 4 >> 2] = a2;
                  if (i * 2 + -1 < j) {
                    b[c2 >> 2] = d2;
                    break;
                  }
                } else {
                  a2 = a2 + 1 | 0;
                  b[c2 + 4 >> 2] = a2;
                }
                d2 = d2 + 1 | 0;
                b[c2 >> 2] = d2;
              }
            } while (0);
            do {
              if (k < 0) {
                if (!(a2 & 1)) {
                  m = (a2 | 0) / 2 | 0;
                  m = _c(d2 | 0, ((d2 | 0) < 0) << 31 >> 31 | 0, m | 0, ((m | 0) < 0) << 31 >> 31 | 0) | 0;
                  d2 = ~~(+(d2 | 0) - (+(m >>> 0) + 4294967296 * +(G() | 0)) * 2);
                  b[c2 >> 2] = d2;
                  break;
                } else {
                  m = (a2 + 1 | 0) / 2 | 0;
                  m = _c(d2 | 0, ((d2 | 0) < 0) << 31 >> 31 | 0, m | 0, ((m | 0) < 0) << 31 >> 31 | 0) | 0;
                  d2 = ~~(+(d2 | 0) - ((+(m >>> 0) + 4294967296 * +(G() | 0)) * 2 + 1));
                  b[c2 >> 2] = d2;
                  break;
                }
              }
            } while (0);
            m = c2 + 4 | 0;
            if (l < 0) {
              d2 = d2 - ((a2 << 1 | 1 | 0) / 2 | 0) | 0;
              b[c2 >> 2] = d2;
              a2 = 0 - a2 | 0;
              b[m >> 2] = a2;
            }
            f = a2 - d2 | 0;
            if ((d2 | 0) < 0) {
              g2 = 0 - d2 | 0;
              b[m >> 2] = f;
              b[n >> 2] = g2;
              b[c2 >> 2] = 0;
              a2 = f;
              d2 = 0;
            } else {
              g2 = 0;
            }
            if ((a2 | 0) < 0) {
              d2 = d2 - a2 | 0;
              b[c2 >> 2] = d2;
              g2 = g2 - a2 | 0;
              b[n >> 2] = g2;
              b[m >> 2] = 0;
              a2 = 0;
            }
            h = d2 - g2 | 0;
            f = a2 - g2 | 0;
            if ((g2 | 0) < 0) {
              b[c2 >> 2] = h;
              b[m >> 2] = f;
              b[n >> 2] = 0;
              a2 = f;
              d2 = h;
              g2 = 0;
            }
            f = (a2 | 0) < (d2 | 0) ? a2 : d2;
            f = (g2 | 0) < (f | 0) ? g2 : f;
            if ((f | 0) <= 0) {
              return;
            }
            b[c2 >> 2] = d2 - f;
            b[m >> 2] = a2 - f;
            b[n >> 2] = g2 - f;
            return;
          }
          function Ca(a2) {
            a2 = a2 | 0;
            var c2 = 0, d2 = 0, e2 = 0, f = 0, g2 = 0, h = 0;
            c2 = b[a2 >> 2] | 0;
            h = a2 + 4 | 0;
            d2 = b[h >> 2] | 0;
            if ((c2 | 0) < 0) {
              d2 = d2 - c2 | 0;
              b[h >> 2] = d2;
              g2 = a2 + 8 | 0;
              b[g2 >> 2] = (b[g2 >> 2] | 0) - c2;
              b[a2 >> 2] = 0;
              c2 = 0;
            }
            if ((d2 | 0) < 0) {
              c2 = c2 - d2 | 0;
              b[a2 >> 2] = c2;
              g2 = a2 + 8 | 0;
              f = (b[g2 >> 2] | 0) - d2 | 0;
              b[g2 >> 2] = f;
              b[h >> 2] = 0;
              d2 = 0;
            } else {
              f = a2 + 8 | 0;
              g2 = f;
              f = b[f >> 2] | 0;
            }
            if ((f | 0) < 0) {
              c2 = c2 - f | 0;
              b[a2 >> 2] = c2;
              d2 = d2 - f | 0;
              b[h >> 2] = d2;
              b[g2 >> 2] = 0;
              f = 0;
            }
            e2 = (d2 | 0) < (c2 | 0) ? d2 : c2;
            e2 = (f | 0) < (e2 | 0) ? f : e2;
            if ((e2 | 0) <= 0) {
              return;
            }
            b[a2 >> 2] = c2 - e2;
            b[h >> 2] = d2 - e2;
            b[g2 >> 2] = f - e2;
            return;
          }
          function Da(a2, c2) {
            a2 = a2 | 0;
            c2 = c2 | 0;
            var d2 = 0, f = 0;
            f = b[a2 + 8 >> 2] | 0;
            d2 = +((b[a2 + 4 >> 2] | 0) - f | 0);
            e[c2 >> 3] = +((b[a2 >> 2] | 0) - f | 0) - d2 * 0.5;
            e[c2 + 8 >> 3] = d2 * 0.8660254037844386;
            return;
          }
          function Ea(a2, c2, d2) {
            a2 = a2 | 0;
            c2 = c2 | 0;
            d2 = d2 | 0;
            b[d2 >> 2] = (b[c2 >> 2] | 0) + (b[a2 >> 2] | 0);
            b[d2 + 4 >> 2] = (b[c2 + 4 >> 2] | 0) + (b[a2 + 4 >> 2] | 0);
            b[d2 + 8 >> 2] = (b[c2 + 8 >> 2] | 0) + (b[a2 + 8 >> 2] | 0);
            return;
          }
          function Fa(a2, c2, d2) {
            a2 = a2 | 0;
            c2 = c2 | 0;
            d2 = d2 | 0;
            b[d2 >> 2] = (b[a2 >> 2] | 0) - (b[c2 >> 2] | 0);
            b[d2 + 4 >> 2] = (b[a2 + 4 >> 2] | 0) - (b[c2 + 4 >> 2] | 0);
            b[d2 + 8 >> 2] = (b[a2 + 8 >> 2] | 0) - (b[c2 + 8 >> 2] | 0);
            return;
          }
          function Ga(a2, c2) {
            a2 = a2 | 0;
            c2 = c2 | 0;
            var d2 = 0, e2 = 0;
            d2 = B(b[a2 >> 2] | 0, c2) | 0;
            b[a2 >> 2] = d2;
            d2 = a2 + 4 | 0;
            e2 = B(b[d2 >> 2] | 0, c2) | 0;
            b[d2 >> 2] = e2;
            a2 = a2 + 8 | 0;
            c2 = B(b[a2 >> 2] | 0, c2) | 0;
            b[a2 >> 2] = c2;
            return;
          }
          function Ha(a2) {
            a2 = a2 | 0;
            var c2 = 0, d2 = 0, e2 = 0, f = 0, g2 = 0, h = 0, i = 0;
            h = b[a2 >> 2] | 0;
            i = (h | 0) < 0;
            e2 = (b[a2 + 4 >> 2] | 0) - (i ? h : 0) | 0;
            g2 = (e2 | 0) < 0;
            f = (g2 ? 0 - e2 | 0 : 0) + ((b[a2 + 8 >> 2] | 0) - (i ? h : 0)) | 0;
            d2 = (f | 0) < 0;
            a2 = d2 ? 0 : f;
            c2 = (g2 ? 0 : e2) - (d2 ? f : 0) | 0;
            f = (i ? 0 : h) - (g2 ? e2 : 0) - (d2 ? f : 0) | 0;
            d2 = (c2 | 0) < (f | 0) ? c2 : f;
            d2 = (a2 | 0) < (d2 | 0) ? a2 : d2;
            e2 = (d2 | 0) > 0;
            a2 = a2 - (e2 ? d2 : 0) | 0;
            c2 = c2 - (e2 ? d2 : 0) | 0;
            a: do {
              switch (f - (e2 ? d2 : 0) | 0) {
                case 0:
                  switch (c2 | 0) {
                    case 0: {
                      i = (a2 | 0) == 0 ? 0 : (a2 | 0) == 1 ? 1 : 7;
                      return i | 0;
                    }
                    case 1: {
                      i = (a2 | 0) == 0 ? 2 : (a2 | 0) == 1 ? 3 : 7;
                      return i | 0;
                    }
                    default:
                      break a;
                  }
                case 1:
                  switch (c2 | 0) {
                    case 0: {
                      i = (a2 | 0) == 0 ? 4 : (a2 | 0) == 1 ? 5 : 7;
                      return i | 0;
                    }
                    case 1: {
                      if (!a2) {
                        a2 = 6;
                      } else {
                        break a;
                      }
                      return a2 | 0;
                    }
                    default:
                      break a;
                  }
                default:
              }
            } while (0);
            i = 7;
            return i | 0;
          }
          function Ia(a2) {
            a2 = a2 | 0;
            var c2 = 0, d2 = 0, e2 = 0, f = 0, g2 = 0, h = 0, i = 0;
            h = a2 + 8 | 0;
            d2 = b[h >> 2] | 0;
            c2 = (b[a2 >> 2] | 0) - d2 | 0;
            i = a2 + 4 | 0;
            d2 = (b[i >> 2] | 0) - d2 | 0;
            e2 = Vc(+((c2 * 3 | 0) - d2 | 0) / 7) | 0;
            b[a2 >> 2] = e2;
            c2 = Vc(+((d2 << 1) + c2 | 0) / 7) | 0;
            b[i >> 2] = c2;
            b[h >> 2] = 0;
            d2 = c2 - e2 | 0;
            if ((e2 | 0) < 0) {
              g2 = 0 - e2 | 0;
              b[i >> 2] = d2;
              b[h >> 2] = g2;
              b[a2 >> 2] = 0;
              c2 = d2;
              e2 = 0;
              d2 = g2;
            } else {
              d2 = 0;
            }
            if ((c2 | 0) < 0) {
              e2 = e2 - c2 | 0;
              b[a2 >> 2] = e2;
              d2 = d2 - c2 | 0;
              b[h >> 2] = d2;
              b[i >> 2] = 0;
              c2 = 0;
            }
            g2 = e2 - d2 | 0;
            f = c2 - d2 | 0;
            if ((d2 | 0) < 0) {
              b[a2 >> 2] = g2;
              b[i >> 2] = f;
              b[h >> 2] = 0;
              c2 = f;
              f = g2;
              d2 = 0;
            } else {
              f = e2;
            }
            e2 = (c2 | 0) < (f | 0) ? c2 : f;
            e2 = (d2 | 0) < (e2 | 0) ? d2 : e2;
            if ((e2 | 0) <= 0) {
              return;
            }
            b[a2 >> 2] = f - e2;
            b[i >> 2] = c2 - e2;
            b[h >> 2] = d2 - e2;
            return;
          }
          function Ja(a2) {
            a2 = a2 | 0;
            var c2 = 0, d2 = 0, e2 = 0, f = 0, g2 = 0, h = 0, i = 0;
            h = a2 + 8 | 0;
            d2 = b[h >> 2] | 0;
            c2 = (b[a2 >> 2] | 0) - d2 | 0;
            i = a2 + 4 | 0;
            d2 = (b[i >> 2] | 0) - d2 | 0;
            e2 = Vc(+((c2 << 1) + d2 | 0) / 7) | 0;
            b[a2 >> 2] = e2;
            c2 = Vc(+((d2 * 3 | 0) - c2 | 0) / 7) | 0;
            b[i >> 2] = c2;
            b[h >> 2] = 0;
            d2 = c2 - e2 | 0;
            if ((e2 | 0) < 0) {
              g2 = 0 - e2 | 0;
              b[i >> 2] = d2;
              b[h >> 2] = g2;
              b[a2 >> 2] = 0;
              c2 = d2;
              e2 = 0;
              d2 = g2;
            } else {
              d2 = 0;
            }
            if ((c2 | 0) < 0) {
              e2 = e2 - c2 | 0;
              b[a2 >> 2] = e2;
              d2 = d2 - c2 | 0;
              b[h >> 2] = d2;
              b[i >> 2] = 0;
              c2 = 0;
            }
            g2 = e2 - d2 | 0;
            f = c2 - d2 | 0;
            if ((d2 | 0) < 0) {
              b[a2 >> 2] = g2;
              b[i >> 2] = f;
              b[h >> 2] = 0;
              c2 = f;
              f = g2;
              d2 = 0;
            } else {
              f = e2;
            }
            e2 = (c2 | 0) < (f | 0) ? c2 : f;
            e2 = (d2 | 0) < (e2 | 0) ? d2 : e2;
            if ((e2 | 0) <= 0) {
              return;
            }
            b[a2 >> 2] = f - e2;
            b[i >> 2] = c2 - e2;
            b[h >> 2] = d2 - e2;
            return;
          }
          function Ka(a2) {
            a2 = a2 | 0;
            var c2 = 0, d2 = 0, e2 = 0, f = 0, g2 = 0, h = 0, i = 0;
            c2 = b[a2 >> 2] | 0;
            h = a2 + 4 | 0;
            d2 = b[h >> 2] | 0;
            i = a2 + 8 | 0;
            e2 = b[i >> 2] | 0;
            f = d2 + (c2 * 3 | 0) | 0;
            b[a2 >> 2] = f;
            d2 = e2 + (d2 * 3 | 0) | 0;
            b[h >> 2] = d2;
            c2 = (e2 * 3 | 0) + c2 | 0;
            b[i >> 2] = c2;
            e2 = d2 - f | 0;
            if ((f | 0) < 0) {
              c2 = c2 - f | 0;
              b[h >> 2] = e2;
              b[i >> 2] = c2;
              b[a2 >> 2] = 0;
              d2 = e2;
              e2 = 0;
            } else {
              e2 = f;
            }
            if ((d2 | 0) < 0) {
              e2 = e2 - d2 | 0;
              b[a2 >> 2] = e2;
              c2 = c2 - d2 | 0;
              b[i >> 2] = c2;
              b[h >> 2] = 0;
              d2 = 0;
            }
            g2 = e2 - c2 | 0;
            f = d2 - c2 | 0;
            if ((c2 | 0) < 0) {
              b[a2 >> 2] = g2;
              b[h >> 2] = f;
              b[i >> 2] = 0;
              e2 = g2;
              c2 = 0;
            } else {
              f = d2;
            }
            d2 = (f | 0) < (e2 | 0) ? f : e2;
            d2 = (c2 | 0) < (d2 | 0) ? c2 : d2;
            if ((d2 | 0) <= 0) {
              return;
            }
            b[a2 >> 2] = e2 - d2;
            b[h >> 2] = f - d2;
            b[i >> 2] = c2 - d2;
            return;
          }
          function La(a2) {
            a2 = a2 | 0;
            var c2 = 0, d2 = 0, e2 = 0, f = 0, g2 = 0, h = 0, i = 0;
            f = b[a2 >> 2] | 0;
            h = a2 + 4 | 0;
            c2 = b[h >> 2] | 0;
            i = a2 + 8 | 0;
            d2 = b[i >> 2] | 0;
            e2 = (c2 * 3 | 0) + f | 0;
            f = d2 + (f * 3 | 0) | 0;
            b[a2 >> 2] = f;
            b[h >> 2] = e2;
            c2 = (d2 * 3 | 0) + c2 | 0;
            b[i >> 2] = c2;
            d2 = e2 - f | 0;
            if ((f | 0) < 0) {
              c2 = c2 - f | 0;
              b[h >> 2] = d2;
              b[i >> 2] = c2;
              b[a2 >> 2] = 0;
              f = 0;
            } else {
              d2 = e2;
            }
            if ((d2 | 0) < 0) {
              f = f - d2 | 0;
              b[a2 >> 2] = f;
              c2 = c2 - d2 | 0;
              b[i >> 2] = c2;
              b[h >> 2] = 0;
              d2 = 0;
            }
            g2 = f - c2 | 0;
            e2 = d2 - c2 | 0;
            if ((c2 | 0) < 0) {
              b[a2 >> 2] = g2;
              b[h >> 2] = e2;
              b[i >> 2] = 0;
              f = g2;
              c2 = 0;
            } else {
              e2 = d2;
            }
            d2 = (e2 | 0) < (f | 0) ? e2 : f;
            d2 = (c2 | 0) < (d2 | 0) ? c2 : d2;
            if ((d2 | 0) <= 0) {
              return;
            }
            b[a2 >> 2] = f - d2;
            b[h >> 2] = e2 - d2;
            b[i >> 2] = c2 - d2;
            return;
          }
          function Ma(a2, c2) {
            a2 = a2 | 0;
            c2 = c2 | 0;
            var d2 = 0, e2 = 0, f = 0, g2 = 0, h = 0, i = 0;
            if ((c2 + -1 | 0) >>> 0 >= 6) {
              return;
            }
            f = (b[15472 + (c2 * 12 | 0) >> 2] | 0) + (b[a2 >> 2] | 0) | 0;
            b[a2 >> 2] = f;
            i = a2 + 4 | 0;
            e2 = (b[15472 + (c2 * 12 | 0) + 4 >> 2] | 0) + (b[i >> 2] | 0) | 0;
            b[i >> 2] = e2;
            h = a2 + 8 | 0;
            c2 = (b[15472 + (c2 * 12 | 0) + 8 >> 2] | 0) + (b[h >> 2] | 0) | 0;
            b[h >> 2] = c2;
            d2 = e2 - f | 0;
            if ((f | 0) < 0) {
              c2 = c2 - f | 0;
              b[i >> 2] = d2;
              b[h >> 2] = c2;
              b[a2 >> 2] = 0;
              e2 = 0;
            } else {
              d2 = e2;
              e2 = f;
            }
            if ((d2 | 0) < 0) {
              e2 = e2 - d2 | 0;
              b[a2 >> 2] = e2;
              c2 = c2 - d2 | 0;
              b[h >> 2] = c2;
              b[i >> 2] = 0;
              d2 = 0;
            }
            g2 = e2 - c2 | 0;
            f = d2 - c2 | 0;
            if ((c2 | 0) < 0) {
              b[a2 >> 2] = g2;
              b[i >> 2] = f;
              b[h >> 2] = 0;
              e2 = g2;
              c2 = 0;
            } else {
              f = d2;
            }
            d2 = (f | 0) < (e2 | 0) ? f : e2;
            d2 = (c2 | 0) < (d2 | 0) ? c2 : d2;
            if ((d2 | 0) <= 0) {
              return;
            }
            b[a2 >> 2] = e2 - d2;
            b[i >> 2] = f - d2;
            b[h >> 2] = c2 - d2;
            return;
          }
          function Na(a2) {
            a2 = a2 | 0;
            var c2 = 0, d2 = 0, e2 = 0, f = 0, g2 = 0, h = 0, i = 0;
            f = b[a2 >> 2] | 0;
            h = a2 + 4 | 0;
            c2 = b[h >> 2] | 0;
            i = a2 + 8 | 0;
            d2 = b[i >> 2] | 0;
            e2 = c2 + f | 0;
            f = d2 + f | 0;
            b[a2 >> 2] = f;
            b[h >> 2] = e2;
            c2 = d2 + c2 | 0;
            b[i >> 2] = c2;
            d2 = e2 - f | 0;
            if ((f | 0) < 0) {
              c2 = c2 - f | 0;
              b[h >> 2] = d2;
              b[i >> 2] = c2;
              b[a2 >> 2] = 0;
              e2 = 0;
            } else {
              d2 = e2;
              e2 = f;
            }
            if ((d2 | 0) < 0) {
              e2 = e2 - d2 | 0;
              b[a2 >> 2] = e2;
              c2 = c2 - d2 | 0;
              b[i >> 2] = c2;
              b[h >> 2] = 0;
              d2 = 0;
            }
            g2 = e2 - c2 | 0;
            f = d2 - c2 | 0;
            if ((c2 | 0) < 0) {
              b[a2 >> 2] = g2;
              b[h >> 2] = f;
              b[i >> 2] = 0;
              e2 = g2;
              c2 = 0;
            } else {
              f = d2;
            }
            d2 = (f | 0) < (e2 | 0) ? f : e2;
            d2 = (c2 | 0) < (d2 | 0) ? c2 : d2;
            if ((d2 | 0) <= 0) {
              return;
            }
            b[a2 >> 2] = e2 - d2;
            b[h >> 2] = f - d2;
            b[i >> 2] = c2 - d2;
            return;
          }
          function Oa(a2) {
            a2 = a2 | 0;
            var c2 = 0, d2 = 0, e2 = 0, f = 0, g2 = 0, h = 0, i = 0;
            c2 = b[a2 >> 2] | 0;
            h = a2 + 4 | 0;
            e2 = b[h >> 2] | 0;
            i = a2 + 8 | 0;
            d2 = b[i >> 2] | 0;
            f = e2 + c2 | 0;
            b[a2 >> 2] = f;
            e2 = d2 + e2 | 0;
            b[h >> 2] = e2;
            c2 = d2 + c2 | 0;
            b[i >> 2] = c2;
            d2 = e2 - f | 0;
            if ((f | 0) < 0) {
              c2 = c2 - f | 0;
              b[h >> 2] = d2;
              b[i >> 2] = c2;
              b[a2 >> 2] = 0;
              e2 = 0;
            } else {
              d2 = e2;
              e2 = f;
            }
            if ((d2 | 0) < 0) {
              e2 = e2 - d2 | 0;
              b[a2 >> 2] = e2;
              c2 = c2 - d2 | 0;
              b[i >> 2] = c2;
              b[h >> 2] = 0;
              d2 = 0;
            }
            g2 = e2 - c2 | 0;
            f = d2 - c2 | 0;
            if ((c2 | 0) < 0) {
              b[a2 >> 2] = g2;
              b[h >> 2] = f;
              b[i >> 2] = 0;
              e2 = g2;
              c2 = 0;
            } else {
              f = d2;
            }
            d2 = (f | 0) < (e2 | 0) ? f : e2;
            d2 = (c2 | 0) < (d2 | 0) ? c2 : d2;
            if ((d2 | 0) <= 0) {
              return;
            }
            b[a2 >> 2] = e2 - d2;
            b[h >> 2] = f - d2;
            b[i >> 2] = c2 - d2;
            return;
          }
          function Pa(a2) {
            a2 = a2 | 0;
            switch (a2 | 0) {
              case 1: {
                a2 = 5;
                break;
              }
              case 5: {
                a2 = 4;
                break;
              }
              case 4: {
                a2 = 6;
                break;
              }
              case 6: {
                a2 = 2;
                break;
              }
              case 2: {
                a2 = 3;
                break;
              }
              case 3: {
                a2 = 1;
                break;
              }
              default:
            }
            return a2 | 0;
          }
          function Qa(a2) {
            a2 = a2 | 0;
            switch (a2 | 0) {
              case 1: {
                a2 = 3;
                break;
              }
              case 3: {
                a2 = 2;
                break;
              }
              case 2: {
                a2 = 6;
                break;
              }
              case 6: {
                a2 = 4;
                break;
              }
              case 4: {
                a2 = 5;
                break;
              }
              case 5: {
                a2 = 1;
                break;
              }
              default:
            }
            return a2 | 0;
          }
          function Ra(a2) {
            a2 = a2 | 0;
            var c2 = 0, d2 = 0, e2 = 0, f = 0, g2 = 0, h = 0, i = 0;
            c2 = b[a2 >> 2] | 0;
            h = a2 + 4 | 0;
            d2 = b[h >> 2] | 0;
            i = a2 + 8 | 0;
            e2 = b[i >> 2] | 0;
            f = d2 + (c2 << 1) | 0;
            b[a2 >> 2] = f;
            d2 = e2 + (d2 << 1) | 0;
            b[h >> 2] = d2;
            c2 = (e2 << 1) + c2 | 0;
            b[i >> 2] = c2;
            e2 = d2 - f | 0;
            if ((f | 0) < 0) {
              c2 = c2 - f | 0;
              b[h >> 2] = e2;
              b[i >> 2] = c2;
              b[a2 >> 2] = 0;
              d2 = e2;
              e2 = 0;
            } else {
              e2 = f;
            }
            if ((d2 | 0) < 0) {
              e2 = e2 - d2 | 0;
              b[a2 >> 2] = e2;
              c2 = c2 - d2 | 0;
              b[i >> 2] = c2;
              b[h >> 2] = 0;
              d2 = 0;
            }
            g2 = e2 - c2 | 0;
            f = d2 - c2 | 0;
            if ((c2 | 0) < 0) {
              b[a2 >> 2] = g2;
              b[h >> 2] = f;
              b[i >> 2] = 0;
              e2 = g2;
              c2 = 0;
            } else {
              f = d2;
            }
            d2 = (f | 0) < (e2 | 0) ? f : e2;
            d2 = (c2 | 0) < (d2 | 0) ? c2 : d2;
            if ((d2 | 0) <= 0) {
              return;
            }
            b[a2 >> 2] = e2 - d2;
            b[h >> 2] = f - d2;
            b[i >> 2] = c2 - d2;
            return;
          }
          function Sa(a2) {
            a2 = a2 | 0;
            var c2 = 0, d2 = 0, e2 = 0, f = 0, g2 = 0, h = 0, i = 0;
            f = b[a2 >> 2] | 0;
            h = a2 + 4 | 0;
            c2 = b[h >> 2] | 0;
            i = a2 + 8 | 0;
            d2 = b[i >> 2] | 0;
            e2 = (c2 << 1) + f | 0;
            f = d2 + (f << 1) | 0;
            b[a2 >> 2] = f;
            b[h >> 2] = e2;
            c2 = (d2 << 1) + c2 | 0;
            b[i >> 2] = c2;
            d2 = e2 - f | 0;
            if ((f | 0) < 0) {
              c2 = c2 - f | 0;
              b[h >> 2] = d2;
              b[i >> 2] = c2;
              b[a2 >> 2] = 0;
              f = 0;
            } else {
              d2 = e2;
            }
            if ((d2 | 0) < 0) {
              f = f - d2 | 0;
              b[a2 >> 2] = f;
              c2 = c2 - d2 | 0;
              b[i >> 2] = c2;
              b[h >> 2] = 0;
              d2 = 0;
            }
            g2 = f - c2 | 0;
            e2 = d2 - c2 | 0;
            if ((c2 | 0) < 0) {
              b[a2 >> 2] = g2;
              b[h >> 2] = e2;
              b[i >> 2] = 0;
              f = g2;
              c2 = 0;
            } else {
              e2 = d2;
            }
            d2 = (e2 | 0) < (f | 0) ? e2 : f;
            d2 = (c2 | 0) < (d2 | 0) ? c2 : d2;
            if ((d2 | 0) <= 0) {
              return;
            }
            b[a2 >> 2] = f - d2;
            b[h >> 2] = e2 - d2;
            b[i >> 2] = c2 - d2;
            return;
          }
          function Ta(a2, c2) {
            a2 = a2 | 0;
            c2 = c2 | 0;
            var d2 = 0, e2 = 0, f = 0, g2 = 0, h = 0, i = 0;
            h = (b[a2 >> 2] | 0) - (b[c2 >> 2] | 0) | 0;
            i = (h | 0) < 0;
            e2 = (b[a2 + 4 >> 2] | 0) - (b[c2 + 4 >> 2] | 0) - (i ? h : 0) | 0;
            g2 = (e2 | 0) < 0;
            f = (i ? 0 - h | 0 : 0) + (b[a2 + 8 >> 2] | 0) - (b[c2 + 8 >> 2] | 0) + (g2 ? 0 - e2 | 0 : 0) | 0;
            a2 = (f | 0) < 0;
            c2 = a2 ? 0 : f;
            d2 = (g2 ? 0 : e2) - (a2 ? f : 0) | 0;
            f = (i ? 0 : h) - (g2 ? e2 : 0) - (a2 ? f : 0) | 0;
            a2 = (d2 | 0) < (f | 0) ? d2 : f;
            a2 = (c2 | 0) < (a2 | 0) ? c2 : a2;
            e2 = (a2 | 0) > 0;
            c2 = c2 - (e2 ? a2 : 0) | 0;
            d2 = d2 - (e2 ? a2 : 0) | 0;
            a2 = f - (e2 ? a2 : 0) | 0;
            a2 = (a2 | 0) > -1 ? a2 : 0 - a2 | 0;
            d2 = (d2 | 0) > -1 ? d2 : 0 - d2 | 0;
            c2 = (c2 | 0) > -1 ? c2 : 0 - c2 | 0;
            c2 = (d2 | 0) > (c2 | 0) ? d2 : c2;
            return ((a2 | 0) > (c2 | 0) ? a2 : c2) | 0;
          }
          function Ua(a2, c2) {
            a2 = a2 | 0;
            c2 = c2 | 0;
            var d2 = 0;
            d2 = b[a2 + 8 >> 2] | 0;
            b[c2 >> 2] = (b[a2 >> 2] | 0) - d2;
            b[c2 + 4 >> 2] = (b[a2 + 4 >> 2] | 0) - d2;
            return;
          }
          function Va(a2, c2) {
            a2 = a2 | 0;
            c2 = c2 | 0;
            var d2 = 0, e2 = 0, f = 0, g2 = 0, h = 0, i = 0;
            e2 = b[a2 >> 2] | 0;
            b[c2 >> 2] = e2;
            a2 = b[a2 + 4 >> 2] | 0;
            h = c2 + 4 | 0;
            b[h >> 2] = a2;
            i = c2 + 8 | 0;
            b[i >> 2] = 0;
            d2 = a2 - e2 | 0;
            if ((e2 | 0) < 0) {
              a2 = 0 - e2 | 0;
              b[h >> 2] = d2;
              b[i >> 2] = a2;
              b[c2 >> 2] = 0;
              e2 = 0;
            } else {
              d2 = a2;
              a2 = 0;
            }
            if ((d2 | 0) < 0) {
              e2 = e2 - d2 | 0;
              b[c2 >> 2] = e2;
              a2 = a2 - d2 | 0;
              b[i >> 2] = a2;
              b[h >> 2] = 0;
              d2 = 0;
            }
            g2 = e2 - a2 | 0;
            f = d2 - a2 | 0;
            if ((a2 | 0) < 0) {
              b[c2 >> 2] = g2;
              b[h >> 2] = f;
              b[i >> 2] = 0;
              d2 = f;
              f = g2;
              a2 = 0;
            } else {
              f = e2;
            }
            e2 = (d2 | 0) < (f | 0) ? d2 : f;
            e2 = (a2 | 0) < (e2 | 0) ? a2 : e2;
            if ((e2 | 0) <= 0) {
              return;
            }
            b[c2 >> 2] = f - e2;
            b[h >> 2] = d2 - e2;
            b[i >> 2] = a2 - e2;
            return;
          }
          function Wa(a2) {
            a2 = a2 | 0;
            var c2 = 0, d2 = 0, e2 = 0, f = 0;
            c2 = a2 + 8 | 0;
            f = b[c2 >> 2] | 0;
            d2 = f - (b[a2 >> 2] | 0) | 0;
            b[a2 >> 2] = d2;
            e2 = a2 + 4 | 0;
            a2 = (b[e2 >> 2] | 0) - f | 0;
            b[e2 >> 2] = a2;
            b[c2 >> 2] = 0 - (a2 + d2);
            return;
          }
          function Xa(a2) {
            a2 = a2 | 0;
            var c2 = 0, d2 = 0, e2 = 0, f = 0, g2 = 0, h = 0, i = 0;
            d2 = b[a2 >> 2] | 0;
            c2 = 0 - d2 | 0;
            b[a2 >> 2] = c2;
            h = a2 + 8 | 0;
            b[h >> 2] = 0;
            i = a2 + 4 | 0;
            e2 = b[i >> 2] | 0;
            f = e2 + d2 | 0;
            if ((d2 | 0) > 0) {
              b[i >> 2] = f;
              b[h >> 2] = d2;
              b[a2 >> 2] = 0;
              c2 = 0;
              e2 = f;
            } else {
              d2 = 0;
            }
            if ((e2 | 0) < 0) {
              g2 = c2 - e2 | 0;
              b[a2 >> 2] = g2;
              d2 = d2 - e2 | 0;
              b[h >> 2] = d2;
              b[i >> 2] = 0;
              f = g2 - d2 | 0;
              c2 = 0 - d2 | 0;
              if ((d2 | 0) < 0) {
                b[a2 >> 2] = f;
                b[i >> 2] = c2;
                b[h >> 2] = 0;
                e2 = c2;
                d2 = 0;
              } else {
                e2 = 0;
                f = g2;
              }
            } else {
              f = c2;
            }
            c2 = (e2 | 0) < (f | 0) ? e2 : f;
            c2 = (d2 | 0) < (c2 | 0) ? d2 : c2;
            if ((c2 | 0) <= 0) {
              return;
            }
            b[a2 >> 2] = f - c2;
            b[i >> 2] = e2 - c2;
            b[h >> 2] = d2 - c2;
            return;
          }
          function Ya(a2, b2, c2) {
            a2 = a2 | 0;
            b2 = b2 | 0;
            c2 = c2 | 0;
            var d2 = 0, e2 = 0;
            d2 = S;
            S = S + 16 | 0;
            e2 = d2;
            Za(a2, b2, c2, e2);
            Ba(e2, c2 + 4 | 0);
            S = d2;
            return;
          }
          function Za(a2, c2, d2, f) {
            a2 = a2 | 0;
            c2 = c2 | 0;
            d2 = d2 | 0;
            f = f | 0;
            var g2 = 0, h = 0, i = 0, j = 0, k = 0;
            k = S;
            S = S + 32 | 0;
            h = k;
            Jc(a2, h);
            b[d2 >> 2] = 0;
            g2 = +Ic(15888, h);
            i = +Ic(15912, h);
            if (i < g2) {
              b[d2 >> 2] = 1;
              g2 = i;
            }
            i = +Ic(15936, h);
            if (i < g2) {
              b[d2 >> 2] = 2;
              g2 = i;
            }
            i = +Ic(15960, h);
            if (i < g2) {
              b[d2 >> 2] = 3;
              g2 = i;
            }
            i = +Ic(15984, h);
            if (i < g2) {
              b[d2 >> 2] = 4;
              g2 = i;
            }
            i = +Ic(16008, h);
            if (i < g2) {
              b[d2 >> 2] = 5;
              g2 = i;
            }
            i = +Ic(16032, h);
            if (i < g2) {
              b[d2 >> 2] = 6;
              g2 = i;
            }
            i = +Ic(16056, h);
            if (i < g2) {
              b[d2 >> 2] = 7;
              g2 = i;
            }
            i = +Ic(16080, h);
            if (i < g2) {
              b[d2 >> 2] = 8;
              g2 = i;
            }
            i = +Ic(16104, h);
            if (i < g2) {
              b[d2 >> 2] = 9;
              g2 = i;
            }
            i = +Ic(16128, h);
            if (i < g2) {
              b[d2 >> 2] = 10;
              g2 = i;
            }
            i = +Ic(16152, h);
            if (i < g2) {
              b[d2 >> 2] = 11;
              g2 = i;
            }
            i = +Ic(16176, h);
            if (i < g2) {
              b[d2 >> 2] = 12;
              g2 = i;
            }
            i = +Ic(16200, h);
            if (i < g2) {
              b[d2 >> 2] = 13;
              g2 = i;
            }
            i = +Ic(16224, h);
            if (i < g2) {
              b[d2 >> 2] = 14;
              g2 = i;
            }
            i = +Ic(16248, h);
            if (i < g2) {
              b[d2 >> 2] = 15;
              g2 = i;
            }
            i = +Ic(16272, h);
            if (i < g2) {
              b[d2 >> 2] = 16;
              g2 = i;
            }
            i = +Ic(16296, h);
            if (i < g2) {
              b[d2 >> 2] = 17;
              g2 = i;
            }
            i = +Ic(16320, h);
            if (i < g2) {
              b[d2 >> 2] = 18;
              g2 = i;
            }
            i = +Ic(16344, h);
            if (i < g2) {
              b[d2 >> 2] = 19;
              g2 = i;
            }
            i = +w(+(1 - g2 * 0.5));
            if (i < 1e-16) {
              b[f >> 2] = 0;
              b[f + 4 >> 2] = 0;
              b[f + 8 >> 2] = 0;
              b[f + 12 >> 2] = 0;
              S = k;
              return;
            }
            d2 = b[d2 >> 2] | 0;
            g2 = +e[16368 + (d2 * 24 | 0) >> 3];
            g2 = +gb(g2 - +gb(+lb(15568 + (d2 << 4) | 0, a2)));
            if (!(Rb(c2) | 0)) {
              j = g2;
            } else {
              j = +gb(g2 + -0.3334731722518321);
            }
            g2 = +v(+i) / 0.381966011250105;
            if ((c2 | 0) > 0) {
              h = 0;
              do {
                g2 = g2 * 2.6457513110645907;
                h = h + 1 | 0;
              } while ((h | 0) != (c2 | 0));
            }
            i = +t(+j) * g2;
            e[f >> 3] = i;
            j = +u(+j) * g2;
            e[f + 8 >> 3] = j;
            S = k;
            return;
          }
          function _a(a2, c2, d2, f, g2) {
            a2 = a2 | 0;
            c2 = c2 | 0;
            d2 = d2 | 0;
            f = f | 0;
            g2 = g2 | 0;
            var h = 0, i = 0;
            h = +Fc(a2);
            if (h < 1e-16) {
              c2 = 15568 + (c2 << 4) | 0;
              b[g2 >> 2] = b[c2 >> 2];
              b[g2 + 4 >> 2] = b[c2 + 4 >> 2];
              b[g2 + 8 >> 2] = b[c2 + 8 >> 2];
              b[g2 + 12 >> 2] = b[c2 + 12 >> 2];
              return;
            }
            i = +z(+ +e[a2 + 8 >> 3], + +e[a2 >> 3]);
            if ((d2 | 0) > 0) {
              a2 = 0;
              do {
                h = h / 2.6457513110645907;
                a2 = a2 + 1 | 0;
              } while ((a2 | 0) != (d2 | 0));
            }
            if (!f) {
              h = +y(+(h * 0.381966011250105));
              if (Rb(d2) | 0) {
                i = +gb(i + 0.3334731722518321);
              }
            } else {
              h = h / 3;
              d2 = (Rb(d2) | 0) == 0;
              h = +y(+((d2 ? h : h / 2.6457513110645907) * 0.381966011250105));
            }
            mb(15568 + (c2 << 4) | 0, +gb(+e[16368 + (c2 * 24 | 0) >> 3] - i), h, g2);
            return;
          }
          function $a(a2, c2, d2) {
            a2 = a2 | 0;
            c2 = c2 | 0;
            d2 = d2 | 0;
            var e2 = 0, f = 0;
            e2 = S;
            S = S + 16 | 0;
            f = e2;
            Da(a2 + 4 | 0, f);
            _a(f, b[a2 >> 2] | 0, c2, 0, d2);
            S = e2;
            return;
          }
          function ab(a2, c2, d2, f, g2) {
            a2 = a2 | 0;
            c2 = c2 | 0;
            d2 = d2 | 0;
            f = f | 0;
            g2 = g2 | 0;
            var h = 0, i = 0, j = 0, k = 0, l = 0, m = 0, n = 0, o = 0, p2 = 0, q2 = 0, r2 = 0, s2 = 0, t2 = 0, u2 = 0, v2 = 0, w2 = 0, x2 = 0, y2 = 0, z2 = 0, A2 = 0, B2 = 0, C2 = 0, D2 = 0, E = 0, F2 = 0, G2 = 0, I2 = 0, J2 = 0;
            G2 = S;
            S = S + 272 | 0;
            h = G2 + 256 | 0;
            u2 = G2 + 240 | 0;
            D2 = G2;
            E = G2 + 224 | 0;
            F2 = G2 + 208 | 0;
            v2 = G2 + 176 | 0;
            w2 = G2 + 160 | 0;
            x2 = G2 + 192 | 0;
            y2 = G2 + 144 | 0;
            z2 = G2 + 128 | 0;
            A2 = G2 + 112 | 0;
            B2 = G2 + 96 | 0;
            C2 = G2 + 80 | 0;
            b[h >> 2] = c2;
            b[u2 >> 2] = b[a2 >> 2];
            b[u2 + 4 >> 2] = b[a2 + 4 >> 2];
            b[u2 + 8 >> 2] = b[a2 + 8 >> 2];
            b[u2 + 12 >> 2] = b[a2 + 12 >> 2];
            bb(u2, h, D2);
            b[g2 >> 2] = 0;
            u2 = f + d2 + ((f | 0) == 5 & 1) | 0;
            if ((u2 | 0) <= (d2 | 0)) {
              S = G2;
              return;
            }
            k = b[h >> 2] | 0;
            l = E + 4 | 0;
            m = v2 + 4 | 0;
            n = d2 + 5 | 0;
            o = 16848 + (k << 2) | 0;
            p2 = 16928 + (k << 2) | 0;
            q2 = z2 + 8 | 0;
            r2 = A2 + 8 | 0;
            s2 = B2 + 8 | 0;
            t2 = F2 + 4 | 0;
            j = d2;
            a: while (1) {
              i = D2 + (((j | 0) % 5 | 0) << 4) | 0;
              b[F2 >> 2] = b[i >> 2];
              b[F2 + 4 >> 2] = b[i + 4 >> 2];
              b[F2 + 8 >> 2] = b[i + 8 >> 2];
              b[F2 + 12 >> 2] = b[i + 12 >> 2];
              do {
              } while ((cb(F2, k, 0, 1) | 0) == 2);
              if ((j | 0) > (d2 | 0) & (Rb(c2) | 0) != 0) {
                b[v2 >> 2] = b[F2 >> 2];
                b[v2 + 4 >> 2] = b[F2 + 4 >> 2];
                b[v2 + 8 >> 2] = b[F2 + 8 >> 2];
                b[v2 + 12 >> 2] = b[F2 + 12 >> 2];
                Da(l, w2);
                f = b[v2 >> 2] | 0;
                h = b[17008 + (f * 80 | 0) + (b[E >> 2] << 2) >> 2] | 0;
                b[v2 >> 2] = b[18608 + (f * 80 | 0) + (h * 20 | 0) >> 2];
                i = b[18608 + (f * 80 | 0) + (h * 20 | 0) + 16 >> 2] | 0;
                if ((i | 0) > 0) {
                  a2 = 0;
                  do {
                    Na(m);
                    a2 = a2 + 1 | 0;
                  } while ((a2 | 0) < (i | 0));
                }
                i = 18608 + (f * 80 | 0) + (h * 20 | 0) + 4 | 0;
                b[x2 >> 2] = b[i >> 2];
                b[x2 + 4 >> 2] = b[i + 4 >> 2];
                b[x2 + 8 >> 2] = b[i + 8 >> 2];
                Ga(x2, (b[o >> 2] | 0) * 3 | 0);
                Ea(m, x2, m);
                Ca(m);
                Da(m, y2);
                I2 = +(b[p2 >> 2] | 0);
                e[z2 >> 3] = I2 * 3;
                e[q2 >> 3] = 0;
                J2 = I2 * -1.5;
                e[A2 >> 3] = J2;
                e[r2 >> 3] = I2 * 2.598076211353316;
                e[B2 >> 3] = J2;
                e[s2 >> 3] = I2 * -2.598076211353316;
                switch (b[17008 + ((b[v2 >> 2] | 0) * 80 | 0) + (b[F2 >> 2] << 2) >> 2] | 0) {
                  case 1: {
                    a2 = A2;
                    f = z2;
                    break;
                  }
                  case 3: {
                    a2 = B2;
                    f = A2;
                    break;
                  }
                  case 2: {
                    a2 = z2;
                    f = B2;
                    break;
                  }
                  default: {
                    a2 = 12;
                    break a;
                  }
                }
                Gc(w2, y2, f, a2, C2);
                _a(C2, b[v2 >> 2] | 0, k, 1, g2 + 8 + (b[g2 >> 2] << 4) | 0);
                b[g2 >> 2] = (b[g2 >> 2] | 0) + 1;
              }
              if ((j | 0) < (n | 0)) {
                Da(t2, v2);
                _a(v2, b[F2 >> 2] | 0, k, 1, g2 + 8 + (b[g2 >> 2] << 4) | 0);
                b[g2 >> 2] = (b[g2 >> 2] | 0) + 1;
              }
              b[E >> 2] = b[F2 >> 2];
              b[E + 4 >> 2] = b[F2 + 4 >> 2];
              b[E + 8 >> 2] = b[F2 + 8 >> 2];
              b[E + 12 >> 2] = b[F2 + 12 >> 2];
              j = j + 1 | 0;
              if ((j | 0) >= (u2 | 0)) {
                a2 = 3;
                break;
              }
            }
            if ((a2 | 0) == 3) {
              S = G2;
              return;
            } else if ((a2 | 0) == 12) {
              H(22474, 22521, 581, 22531);
            }
          }
          function bb(a2, c2, d2) {
            a2 = a2 | 0;
            c2 = c2 | 0;
            d2 = d2 | 0;
            var e2 = 0, f = 0, g2 = 0, h = 0, i = 0, j = 0;
            j = S;
            S = S + 128 | 0;
            e2 = j + 64 | 0;
            f = j;
            g2 = e2;
            h = 20208;
            i = g2 + 60 | 0;
            do {
              b[g2 >> 2] = b[h >> 2];
              g2 = g2 + 4 | 0;
              h = h + 4 | 0;
            } while ((g2 | 0) < (i | 0));
            g2 = f;
            h = 20272;
            i = g2 + 60 | 0;
            do {
              b[g2 >> 2] = b[h >> 2];
              g2 = g2 + 4 | 0;
              h = h + 4 | 0;
            } while ((g2 | 0) < (i | 0));
            i = (Rb(b[c2 >> 2] | 0) | 0) == 0;
            e2 = i ? e2 : f;
            f = a2 + 4 | 0;
            Ra(f);
            Sa(f);
            if (Rb(b[c2 >> 2] | 0) | 0) {
              La(f);
              b[c2 >> 2] = (b[c2 >> 2] | 0) + 1;
            }
            b[d2 >> 2] = b[a2 >> 2];
            c2 = d2 + 4 | 0;
            Ea(f, e2, c2);
            Ca(c2);
            b[d2 + 16 >> 2] = b[a2 >> 2];
            c2 = d2 + 20 | 0;
            Ea(f, e2 + 12 | 0, c2);
            Ca(c2);
            b[d2 + 32 >> 2] = b[a2 >> 2];
            c2 = d2 + 36 | 0;
            Ea(f, e2 + 24 | 0, c2);
            Ca(c2);
            b[d2 + 48 >> 2] = b[a2 >> 2];
            c2 = d2 + 52 | 0;
            Ea(f, e2 + 36 | 0, c2);
            Ca(c2);
            b[d2 + 64 >> 2] = b[a2 >> 2];
            d2 = d2 + 68 | 0;
            Ea(f, e2 + 48 | 0, d2);
            Ca(d2);
            S = j;
            return;
          }
          function cb(a2, c2, d2, e2) {
            a2 = a2 | 0;
            c2 = c2 | 0;
            d2 = d2 | 0;
            e2 = e2 | 0;
            var f = 0, g2 = 0, h = 0, i = 0, j = 0, k = 0, l = 0, m = 0, n = 0, o = 0, p2 = 0;
            p2 = S;
            S = S + 32 | 0;
            n = p2 + 12 | 0;
            i = p2;
            o = a2 + 4 | 0;
            m = b[16928 + (c2 << 2) >> 2] | 0;
            l = (e2 | 0) != 0;
            m = l ? m * 3 | 0 : m;
            f = b[o >> 2] | 0;
            k = a2 + 8 | 0;
            h = b[k >> 2] | 0;
            if (l) {
              g2 = a2 + 12 | 0;
              e2 = b[g2 >> 2] | 0;
              f = h + f + e2 | 0;
              if ((f | 0) == (m | 0)) {
                o = 1;
                S = p2;
                return o | 0;
              } else {
                j = g2;
              }
            } else {
              j = a2 + 12 | 0;
              e2 = b[j >> 2] | 0;
              f = h + f + e2 | 0;
            }
            if ((f | 0) <= (m | 0)) {
              o = 0;
              S = p2;
              return o | 0;
            }
            do {
              if ((e2 | 0) > 0) {
                e2 = b[a2 >> 2] | 0;
                if ((h | 0) > 0) {
                  g2 = 18608 + (e2 * 80 | 0) + 60 | 0;
                  e2 = a2;
                  break;
                }
                e2 = 18608 + (e2 * 80 | 0) + 40 | 0;
                if (!d2) {
                  g2 = e2;
                  e2 = a2;
                } else {
                  Aa(n, m, 0, 0);
                  Fa(o, n, i);
                  Oa(i);
                  Ea(i, n, o);
                  g2 = e2;
                  e2 = a2;
                }
              } else {
                g2 = 18608 + ((b[a2 >> 2] | 0) * 80 | 0) + 20 | 0;
                e2 = a2;
              }
            } while (0);
            b[e2 >> 2] = b[g2 >> 2];
            f = g2 + 16 | 0;
            if ((b[f >> 2] | 0) > 0) {
              e2 = 0;
              do {
                Na(o);
                e2 = e2 + 1 | 0;
              } while ((e2 | 0) < (b[f >> 2] | 0));
            }
            a2 = g2 + 4 | 0;
            b[n >> 2] = b[a2 >> 2];
            b[n + 4 >> 2] = b[a2 + 4 >> 2];
            b[n + 8 >> 2] = b[a2 + 8 >> 2];
            c2 = b[16848 + (c2 << 2) >> 2] | 0;
            Ga(n, l ? c2 * 3 | 0 : c2);
            Ea(o, n, o);
            Ca(o);
            if (l) {
              e2 = ((b[k >> 2] | 0) + (b[o >> 2] | 0) + (b[j >> 2] | 0) | 0) == (m | 0) ? 1 : 2;
            } else {
              e2 = 2;
            }
            o = e2;
            S = p2;
            return o | 0;
          }
          function db(a2, b2) {
            a2 = a2 | 0;
            b2 = b2 | 0;
            var c2 = 0;
            do {
              c2 = cb(a2, b2, 0, 1) | 0;
            } while ((c2 | 0) == 2);
            return c2 | 0;
          }
          function eb(a2, c2, d2, f, g2) {
            a2 = a2 | 0;
            c2 = c2 | 0;
            d2 = d2 | 0;
            f = f | 0;
            g2 = g2 | 0;
            var h = 0, i = 0, j = 0, k = 0, l = 0, m = 0, n = 0, o = 0, p2 = 0, q2 = 0, r2 = 0, s2 = 0, t2 = 0, u2 = 0, v2 = 0, w2 = 0, x2 = 0, y2 = 0, z2 = 0, A2 = 0, B2 = 0, C2 = 0, D2 = 0;
            B2 = S;
            S = S + 240 | 0;
            h = B2 + 224 | 0;
            x2 = B2 + 208 | 0;
            y2 = B2;
            z2 = B2 + 192 | 0;
            A2 = B2 + 176 | 0;
            s2 = B2 + 160 | 0;
            t2 = B2 + 144 | 0;
            u2 = B2 + 128 | 0;
            v2 = B2 + 112 | 0;
            w2 = B2 + 96 | 0;
            b[h >> 2] = c2;
            b[x2 >> 2] = b[a2 >> 2];
            b[x2 + 4 >> 2] = b[a2 + 4 >> 2];
            b[x2 + 8 >> 2] = b[a2 + 8 >> 2];
            b[x2 + 12 >> 2] = b[a2 + 12 >> 2];
            fb(x2, h, y2);
            b[g2 >> 2] = 0;
            r2 = f + d2 + ((f | 0) == 6 & 1) | 0;
            if ((r2 | 0) <= (d2 | 0)) {
              S = B2;
              return;
            }
            k = b[h >> 2] | 0;
            l = d2 + 6 | 0;
            m = 16928 + (k << 2) | 0;
            n = t2 + 8 | 0;
            o = u2 + 8 | 0;
            p2 = v2 + 8 | 0;
            q2 = z2 + 4 | 0;
            i = 0;
            j = d2;
            f = -1;
            a: while (1) {
              h = (j | 0) % 6 | 0;
              a2 = y2 + (h << 4) | 0;
              b[z2 >> 2] = b[a2 >> 2];
              b[z2 + 4 >> 2] = b[a2 + 4 >> 2];
              b[z2 + 8 >> 2] = b[a2 + 8 >> 2];
              b[z2 + 12 >> 2] = b[a2 + 12 >> 2];
              a2 = i;
              i = cb(z2, k, 0, 1) | 0;
              if ((j | 0) > (d2 | 0) & (Rb(c2) | 0) != 0 ? (a2 | 0) != 1 ? (b[z2 >> 2] | 0) != (f | 0) : 0 : 0) {
                Da(y2 + (((h + 5 | 0) % 6 | 0) << 4) + 4 | 0, A2);
                Da(y2 + (h << 4) + 4 | 0, s2);
                C2 = +(b[m >> 2] | 0);
                e[t2 >> 3] = C2 * 3;
                e[n >> 3] = 0;
                D2 = C2 * -1.5;
                e[u2 >> 3] = D2;
                e[o >> 3] = C2 * 2.598076211353316;
                e[v2 >> 3] = D2;
                e[p2 >> 3] = C2 * -2.598076211353316;
                h = b[x2 >> 2] | 0;
                switch (b[17008 + (h * 80 | 0) + (((f | 0) == (h | 0) ? b[z2 >> 2] | 0 : f) << 2) >> 2] | 0) {
                  case 1: {
                    a2 = u2;
                    f = t2;
                    break;
                  }
                  case 3: {
                    a2 = v2;
                    f = u2;
                    break;
                  }
                  case 2: {
                    a2 = t2;
                    f = v2;
                    break;
                  }
                  default: {
                    a2 = 8;
                    break a;
                  }
                }
                Gc(A2, s2, f, a2, w2);
                if (!(Hc(A2, w2) | 0) ? !(Hc(s2, w2) | 0) : 0) {
                  _a(w2, b[x2 >> 2] | 0, k, 1, g2 + 8 + (b[g2 >> 2] << 4) | 0);
                  b[g2 >> 2] = (b[g2 >> 2] | 0) + 1;
                }
              }
              if ((j | 0) < (l | 0)) {
                Da(q2, A2);
                _a(A2, b[z2 >> 2] | 0, k, 1, g2 + 8 + (b[g2 >> 2] << 4) | 0);
                b[g2 >> 2] = (b[g2 >> 2] | 0) + 1;
              }
              j = j + 1 | 0;
              if ((j | 0) >= (r2 | 0)) {
                a2 = 3;
                break;
              } else {
                f = b[z2 >> 2] | 0;
              }
            }
            if ((a2 | 0) == 3) {
              S = B2;
              return;
            } else if ((a2 | 0) == 8) {
              H(22557, 22521, 746, 22602);
            }
          }
          function fb(a2, c2, d2) {
            a2 = a2 | 0;
            c2 = c2 | 0;
            d2 = d2 | 0;
            var e2 = 0, f = 0, g2 = 0, h = 0, i = 0, j = 0;
            j = S;
            S = S + 160 | 0;
            e2 = j + 80 | 0;
            f = j;
            g2 = e2;
            h = 20336;
            i = g2 + 72 | 0;
            do {
              b[g2 >> 2] = b[h >> 2];
              g2 = g2 + 4 | 0;
              h = h + 4 | 0;
            } while ((g2 | 0) < (i | 0));
            g2 = f;
            h = 20416;
            i = g2 + 72 | 0;
            do {
              b[g2 >> 2] = b[h >> 2];
              g2 = g2 + 4 | 0;
              h = h + 4 | 0;
            } while ((g2 | 0) < (i | 0));
            i = (Rb(b[c2 >> 2] | 0) | 0) == 0;
            e2 = i ? e2 : f;
            f = a2 + 4 | 0;
            Ra(f);
            Sa(f);
            if (Rb(b[c2 >> 2] | 0) | 0) {
              La(f);
              b[c2 >> 2] = (b[c2 >> 2] | 0) + 1;
            }
            b[d2 >> 2] = b[a2 >> 2];
            c2 = d2 + 4 | 0;
            Ea(f, e2, c2);
            Ca(c2);
            b[d2 + 16 >> 2] = b[a2 >> 2];
            c2 = d2 + 20 | 0;
            Ea(f, e2 + 12 | 0, c2);
            Ca(c2);
            b[d2 + 32 >> 2] = b[a2 >> 2];
            c2 = d2 + 36 | 0;
            Ea(f, e2 + 24 | 0, c2);
            Ca(c2);
            b[d2 + 48 >> 2] = b[a2 >> 2];
            c2 = d2 + 52 | 0;
            Ea(f, e2 + 36 | 0, c2);
            Ca(c2);
            b[d2 + 64 >> 2] = b[a2 >> 2];
            c2 = d2 + 68 | 0;
            Ea(f, e2 + 48 | 0, c2);
            Ca(c2);
            b[d2 + 80 >> 2] = b[a2 >> 2];
            d2 = d2 + 84 | 0;
            Ea(f, e2 + 60 | 0, d2);
            Ca(d2);
            S = j;
            return;
          }
          function gb(a2) {
            a2 = +a2;
            var b2 = 0;
            b2 = a2 < 0 ? a2 + 6.283185307179586 : a2;
            return +(!(a2 >= 6.283185307179586) ? b2 : b2 + -6.283185307179586);
          }
          function hb(a2, b2) {
            a2 = a2 | 0;
            b2 = b2 | 0;
            if (!(+q(+(+e[a2 >> 3] - +e[b2 >> 3])) < 17453292519943298e-27)) {
              b2 = 0;
              return b2 | 0;
            }
            b2 = +q(+(+e[a2 + 8 >> 3] - +e[b2 + 8 >> 3])) < 17453292519943298e-27;
            return b2 | 0;
          }
          function ib(a2, b2) {
            a2 = a2 | 0;
            b2 = b2 | 0;
            var c2 = 0, d2 = 0, f = 0, g2 = 0;
            f = +e[b2 >> 3];
            d2 = +e[a2 >> 3];
            g2 = +u(+((f - d2) * 0.5));
            c2 = +u(+((+e[b2 + 8 >> 3] - +e[a2 + 8 >> 3]) * 0.5));
            c2 = g2 * g2 + c2 * (+t(+f) * +t(+d2) * c2);
            return +(+z(+ +r(+c2), + +r(+(1 - c2))) * 2);
          }
          function jb(a2, b2) {
            a2 = a2 | 0;
            b2 = b2 | 0;
            var c2 = 0, d2 = 0, f = 0, g2 = 0;
            f = +e[b2 >> 3];
            d2 = +e[a2 >> 3];
            g2 = +u(+((f - d2) * 0.5));
            c2 = +u(+((+e[b2 + 8 >> 3] - +e[a2 + 8 >> 3]) * 0.5));
            c2 = g2 * g2 + c2 * (+t(+f) * +t(+d2) * c2);
            return +(+z(+ +r(+c2), + +r(+(1 - c2))) * 2 * 6371.007180918475);
          }
          function kb(a2, b2) {
            a2 = a2 | 0;
            b2 = b2 | 0;
            var c2 = 0, d2 = 0, f = 0, g2 = 0;
            f = +e[b2 >> 3];
            d2 = +e[a2 >> 3];
            g2 = +u(+((f - d2) * 0.5));
            c2 = +u(+((+e[b2 + 8 >> 3] - +e[a2 + 8 >> 3]) * 0.5));
            c2 = g2 * g2 + c2 * (+t(+f) * +t(+d2) * c2);
            return +(+z(+ +r(+c2), + +r(+(1 - c2))) * 2 * 6371.007180918475 * 1e3);
          }
          function lb(a2, b2) {
            a2 = a2 | 0;
            b2 = b2 | 0;
            var c2 = 0, d2 = 0, f = 0, g2 = 0, h = 0;
            g2 = +e[b2 >> 3];
            d2 = +t(+g2);
            f = +e[b2 + 8 >> 3] - +e[a2 + 8 >> 3];
            h = d2 * +u(+f);
            c2 = +e[a2 >> 3];
            return + +z(+h, +(+u(+g2) * +t(+c2) - +t(+f) * (d2 * +u(+c2))));
          }
          function mb(a2, c2, d2, f) {
            a2 = a2 | 0;
            c2 = +c2;
            d2 = +d2;
            f = f | 0;
            var g2 = 0, h = 0, i = 0, j = 0;
            if (d2 < 1e-16) {
              b[f >> 2] = b[a2 >> 2];
              b[f + 4 >> 2] = b[a2 + 4 >> 2];
              b[f + 8 >> 2] = b[a2 + 8 >> 2];
              b[f + 12 >> 2] = b[a2 + 12 >> 2];
              return;
            }
            h = c2 < 0 ? c2 + 6.283185307179586 : c2;
            h = !(c2 >= 6.283185307179586) ? h : h + -6.283185307179586;
            do {
              if (h < 1e-16) {
                c2 = +e[a2 >> 3] + d2;
                e[f >> 3] = c2;
                g2 = f;
              } else {
                g2 = +q(+(h + -3.141592653589793)) < 1e-16;
                c2 = +e[a2 >> 3];
                if (g2) {
                  c2 = c2 - d2;
                  e[f >> 3] = c2;
                  g2 = f;
                  break;
                }
                i = +t(+d2);
                d2 = +u(+d2);
                c2 = i * +u(+c2) + +t(+h) * (d2 * +t(+c2));
                c2 = c2 > 1 ? 1 : c2;
                c2 = +x(+(c2 < -1 ? -1 : c2));
                e[f >> 3] = c2;
                if (+q(+(c2 + -1.5707963267948966)) < 1e-16) {
                  e[f >> 3] = 1.5707963267948966;
                  e[f + 8 >> 3] = 0;
                  return;
                }
                if (+q(+(c2 + 1.5707963267948966)) < 1e-16) {
                  e[f >> 3] = -1.5707963267948966;
                  e[f + 8 >> 3] = 0;
                  return;
                }
                j = +t(+c2);
                h = d2 * +u(+h) / j;
                d2 = +e[a2 >> 3];
                c2 = (i - +u(+c2) * +u(+d2)) / +t(+d2) / j;
                i = h > 1 ? 1 : h;
                c2 = c2 > 1 ? 1 : c2;
                c2 = +e[a2 + 8 >> 3] + +z(+(i < -1 ? -1 : i), +(c2 < -1 ? -1 : c2));
                if (c2 > 3.141592653589793) {
                  do {
                    c2 = c2 + -6.283185307179586;
                  } while (c2 > 3.141592653589793);
                }
                if (c2 < -3.141592653589793) {
                  do {
                    c2 = c2 + 6.283185307179586;
                  } while (c2 < -3.141592653589793);
                }
                e[f + 8 >> 3] = c2;
                return;
              }
            } while (0);
            if (+q(+(c2 + -1.5707963267948966)) < 1e-16) {
              e[g2 >> 3] = 1.5707963267948966;
              e[f + 8 >> 3] = 0;
              return;
            }
            if (+q(+(c2 + 1.5707963267948966)) < 1e-16) {
              e[g2 >> 3] = -1.5707963267948966;
              e[f + 8 >> 3] = 0;
              return;
            }
            c2 = +e[a2 + 8 >> 3];
            if (c2 > 3.141592653589793) {
              do {
                c2 = c2 + -6.283185307179586;
              } while (c2 > 3.141592653589793);
            }
            if (c2 < -3.141592653589793) {
              do {
                c2 = c2 + 6.283185307179586;
              } while (c2 < -3.141592653589793);
            }
            e[f + 8 >> 3] = c2;
            return;
          }
          function nb(a2) {
            a2 = a2 | 0;
            return + +e[20496 + (a2 << 3) >> 3];
          }
          function ob(a2) {
            a2 = a2 | 0;
            return + +e[20624 + (a2 << 3) >> 3];
          }
          function pb(a2) {
            a2 = a2 | 0;
            return + +e[20752 + (a2 << 3) >> 3];
          }
          function qb(a2) {
            a2 = a2 | 0;
            return + +e[20880 + (a2 << 3) >> 3];
          }
          function rb(a2) {
            a2 = a2 | 0;
            var c2 = 0;
            c2 = 21008 + (a2 << 3) | 0;
            a2 = b[c2 >> 2] | 0;
            F(b[c2 + 4 >> 2] | 0);
            return a2 | 0;
          }
          function sb(a2, b2, c2) {
            a2 = a2 | 0;
            b2 = b2 | 0;
            c2 = c2 | 0;
            var d2 = 0, f = 0, g2 = 0, h = 0, i = 0, j = 0, k = 0, l = 0, m = 0, n = 0;
            n = +e[b2 >> 3];
            l = +e[a2 >> 3];
            j = +u(+((n - l) * 0.5));
            g2 = +e[b2 + 8 >> 3];
            k = +e[a2 + 8 >> 3];
            h = +u(+((g2 - k) * 0.5));
            i = +t(+l);
            m = +t(+n);
            h = j * j + h * (m * i * h);
            h = +z(+ +r(+h), + +r(+(1 - h))) * 2;
            j = +e[c2 >> 3];
            n = +u(+((j - n) * 0.5));
            d2 = +e[c2 + 8 >> 3];
            g2 = +u(+((d2 - g2) * 0.5));
            f = +t(+j);
            g2 = n * n + g2 * (m * f * g2);
            g2 = +z(+ +r(+g2), + +r(+(1 - g2))) * 2;
            j = +u(+((l - j) * 0.5));
            d2 = +u(+((k - d2) * 0.5));
            d2 = j * j + d2 * (i * f * d2);
            d2 = +z(+ +r(+d2), + +r(+(1 - d2))) * 2;
            f = (h + g2 + d2) * 0.5;
            return +(+y(+ +r(+(+v(+(f * 0.5)) * +v(+((f - h) * 0.5)) * +v(+((f - g2) * 0.5)) * +v(+((f - d2) * 0.5))))) * 4);
          }
          function tb(a2, c2) {
            a2 = a2 | 0;
            c2 = c2 | 0;
            var d2 = 0, e2 = 0, f = 0, g2 = 0, h = 0;
            g2 = S;
            S = S + 192 | 0;
            e2 = g2 + 168 | 0;
            f = g2;
            Vb(a2, c2, e2);
            Wb(a2, c2, f);
            c2 = b[f >> 2] | 0;
            if ((c2 | 0) <= 0) {
              d2 = 0;
              S = g2;
              return +d2;
            }
            d2 = +sb(f + 8 | 0, f + 8 + (((c2 | 0) != 1 & 1) << 4) | 0, e2) + 0;
            if ((c2 | 0) == 1) {
              S = g2;
              return +d2;
            }
            a2 = 1;
            do {
              h = a2;
              a2 = a2 + 1 | 0;
              d2 = d2 + +sb(f + 8 + (h << 4) | 0, f + 8 + (((a2 | 0) % (c2 | 0) | 0) << 4) | 0, e2);
            } while ((a2 | 0) < (c2 | 0));
            S = g2;
            return +d2;
          }
          function ub(a2, c2) {
            a2 = a2 | 0;
            c2 = c2 | 0;
            var d2 = 0, e2 = 0, f = 0, g2 = 0, h = 0;
            g2 = S;
            S = S + 192 | 0;
            e2 = g2 + 168 | 0;
            f = g2;
            Vb(a2, c2, e2);
            Wb(a2, c2, f);
            c2 = b[f >> 2] | 0;
            if ((c2 | 0) > 0) {
              d2 = +sb(f + 8 | 0, f + 8 + (((c2 | 0) != 1 & 1) << 4) | 0, e2) + 0;
              if ((c2 | 0) != 1) {
                a2 = 1;
                do {
                  h = a2;
                  a2 = a2 + 1 | 0;
                  d2 = d2 + +sb(f + 8 + (h << 4) | 0, f + 8 + (((a2 | 0) % (c2 | 0) | 0) << 4) | 0, e2);
                } while ((a2 | 0) < (c2 | 0));
              }
            } else {
              d2 = 0;
            }
            S = g2;
            return +(d2 * 6371.007180918475 * 6371.007180918475);
          }
          function vb(a2, c2) {
            a2 = a2 | 0;
            c2 = c2 | 0;
            var d2 = 0, e2 = 0, f = 0, g2 = 0, h = 0;
            g2 = S;
            S = S + 192 | 0;
            e2 = g2 + 168 | 0;
            f = g2;
            Vb(a2, c2, e2);
            Wb(a2, c2, f);
            c2 = b[f >> 2] | 0;
            if ((c2 | 0) > 0) {
              d2 = +sb(f + 8 | 0, f + 8 + (((c2 | 0) != 1 & 1) << 4) | 0, e2) + 0;
              if ((c2 | 0) != 1) {
                a2 = 1;
                do {
                  h = a2;
                  a2 = a2 + 1 | 0;
                  d2 = d2 + +sb(f + 8 + (h << 4) | 0, f + 8 + (((a2 | 0) % (c2 | 0) | 0) << 4) | 0, e2);
                } while ((a2 | 0) < (c2 | 0));
              }
            } else {
              d2 = 0;
            }
            S = g2;
            return +(d2 * 6371.007180918475 * 6371.007180918475 * 1e3 * 1e3);
          }
          function wb(a2, c2) {
            a2 = a2 | 0;
            c2 = c2 | 0;
            var d2 = 0, f = 0, g2 = 0, h = 0, i = 0, j = 0, k = 0, l = 0;
            i = S;
            S = S + 176 | 0;
            h = i;
            gc(a2, c2, h);
            a2 = b[h >> 2] | 0;
            if ((a2 | 0) <= 1) {
              g2 = 0;
              S = i;
              return +g2;
            }
            c2 = a2 + -1 | 0;
            a2 = 0;
            d2 = 0;
            f = +e[h + 8 >> 3];
            g2 = +e[h + 16 >> 3];
            do {
              a2 = a2 + 1 | 0;
              k = f;
              f = +e[h + 8 + (a2 << 4) >> 3];
              l = +u(+((f - k) * 0.5));
              j = g2;
              g2 = +e[h + 8 + (a2 << 4) + 8 >> 3];
              j = +u(+((g2 - j) * 0.5));
              j = l * l + j * (+t(+f) * +t(+k) * j);
              d2 = d2 + +z(+ +r(+j), + +r(+(1 - j))) * 2;
            } while ((a2 | 0) < (c2 | 0));
            S = i;
            return +d2;
          }
          function xb(a2, c2) {
            a2 = a2 | 0;
            c2 = c2 | 0;
            var d2 = 0, f = 0, g2 = 0, h = 0, i = 0, j = 0, k = 0, l = 0;
            i = S;
            S = S + 176 | 0;
            h = i;
            gc(a2, c2, h);
            a2 = b[h >> 2] | 0;
            if ((a2 | 0) <= 1) {
              g2 = 0;
              S = i;
              return +g2;
            }
            c2 = a2 + -1 | 0;
            a2 = 0;
            d2 = 0;
            f = +e[h + 8 >> 3];
            g2 = +e[h + 16 >> 3];
            do {
              a2 = a2 + 1 | 0;
              k = f;
              f = +e[h + 8 + (a2 << 4) >> 3];
              l = +u(+((f - k) * 0.5));
              j = g2;
              g2 = +e[h + 8 + (a2 << 4) + 8 >> 3];
              j = +u(+((g2 - j) * 0.5));
              j = l * l + j * (+t(+k) * +t(+f) * j);
              d2 = d2 + +z(+ +r(+j), + +r(+(1 - j))) * 2;
            } while ((a2 | 0) != (c2 | 0));
            l = d2 * 6371.007180918475;
            S = i;
            return +l;
          }
          function yb(a2, c2) {
            a2 = a2 | 0;
            c2 = c2 | 0;
            var d2 = 0, f = 0, g2 = 0, h = 0, i = 0, j = 0, k = 0, l = 0;
            i = S;
            S = S + 176 | 0;
            h = i;
            gc(a2, c2, h);
            a2 = b[h >> 2] | 0;
            if ((a2 | 0) <= 1) {
              g2 = 0;
              S = i;
              return +g2;
            }
            c2 = a2 + -1 | 0;
            a2 = 0;
            d2 = 0;
            f = +e[h + 8 >> 3];
            g2 = +e[h + 16 >> 3];
            do {
              a2 = a2 + 1 | 0;
              k = f;
              f = +e[h + 8 + (a2 << 4) >> 3];
              l = +u(+((f - k) * 0.5));
              j = g2;
              g2 = +e[h + 8 + (a2 << 4) + 8 >> 3];
              j = +u(+((g2 - j) * 0.5));
              j = l * l + j * (+t(+k) * +t(+f) * j);
              d2 = d2 + +z(+ +r(+j), + +r(+(1 - j))) * 2;
            } while ((a2 | 0) != (c2 | 0));
            l = d2 * 6371.007180918475 * 1e3;
            S = i;
            return +l;
          }
          function zb(a2, b2) {
            a2 = a2 | 0;
            b2 = b2 | 0;
            b2 = cd(a2 | 0, b2 | 0, 52) | 0;
            G() | 0;
            return b2 & 15 | 0;
          }
          function Ab(a2, b2) {
            a2 = a2 | 0;
            b2 = b2 | 0;
            b2 = cd(a2 | 0, b2 | 0, 45) | 0;
            G() | 0;
            return b2 & 127 | 0;
          }
          function Bb(a2, b2) {
            a2 = a2 | 0;
            b2 = b2 | 0;
            var c2 = 0, d2 = 0, e2 = 0, f = 0, g2 = 0, h = 0;
            if (!(true & (b2 & -16777216 | 0) == 134217728)) {
              b2 = 0;
              return b2 | 0;
            }
            g2 = cd(a2 | 0, b2 | 0, 45) | 0;
            G() | 0;
            g2 = g2 & 127;
            if (g2 >>> 0 > 121) {
              b2 = 0;
              return b2 | 0;
            }
            c2 = cd(a2 | 0, b2 | 0, 52) | 0;
            G() | 0;
            c2 = c2 & 15;
            do {
              if (c2 | 0) {
                e2 = 1;
                d2 = 0;
                while (1) {
                  f = cd(a2 | 0, b2 | 0, (15 - e2 | 0) * 3 | 0) | 0;
                  G() | 0;
                  f = f & 7;
                  if ((f | 0) != 0 & (d2 ^ 1)) {
                    if ((f | 0) == 1 & (la(g2) | 0) != 0) {
                      h = 0;
                      d2 = 13;
                      break;
                    } else {
                      d2 = 1;
                    }
                  }
                  if ((f | 0) == 7) {
                    h = 0;
                    d2 = 13;
                    break;
                  }
                  if (e2 >>> 0 < c2 >>> 0) {
                    e2 = e2 + 1 | 0;
                  } else {
                    d2 = 9;
                    break;
                  }
                }
                if ((d2 | 0) == 9) {
                  if ((c2 | 0) == 15) {
                    h = 1;
                  } else {
                    break;
                  }
                  return h | 0;
                } else if ((d2 | 0) == 13) {
                  return h | 0;
                }
              }
            } while (0);
            while (1) {
              h = cd(a2 | 0, b2 | 0, (14 - c2 | 0) * 3 | 0) | 0;
              G() | 0;
              if (!((h & 7 | 0) == 7 & true)) {
                h = 0;
                d2 = 13;
                break;
              }
              if (c2 >>> 0 < 14) {
                c2 = c2 + 1 | 0;
              } else {
                h = 1;
                d2 = 13;
                break;
              }
            }
            if ((d2 | 0) == 13) {
              return h | 0;
            }
            return 0;
          }
          function Cb(a2, b2, c2) {
            a2 = a2 | 0;
            b2 = b2 | 0;
            c2 = c2 | 0;
            var d2 = 0, e2 = 0;
            d2 = cd(a2 | 0, b2 | 0, 52) | 0;
            G() | 0;
            d2 = d2 & 15;
            if ((d2 | 0) >= (c2 | 0)) {
              if ((d2 | 0) != (c2 | 0)) {
                if (c2 >>> 0 <= 15) {
                  e2 = dd(c2 | 0, 0, 52) | 0;
                  a2 = e2 | a2;
                  b2 = G() | 0 | b2 & -15728641;
                  if ((d2 | 0) > (c2 | 0)) {
                    do {
                      e2 = dd(7, 0, (14 - c2 | 0) * 3 | 0) | 0;
                      c2 = c2 + 1 | 0;
                      a2 = e2 | a2;
                      b2 = G() | 0 | b2;
                    } while ((c2 | 0) < (d2 | 0));
                  }
                } else {
                  b2 = 0;
                  a2 = 0;
                }
              }
            } else {
              b2 = 0;
              a2 = 0;
            }
            F(b2 | 0);
            return a2 | 0;
          }
          function Db(a2, b2, c2) {
            a2 = a2 | 0;
            b2 = b2 | 0;
            c2 = c2 | 0;
            a2 = cd(a2 | 0, b2 | 0, 52) | 0;
            G() | 0;
            a2 = a2 & 15;
            if (!((c2 | 0) < 16 & (a2 | 0) <= (c2 | 0))) {
              c2 = 0;
              return c2 | 0;
            }
            c2 = tc(7, c2 - a2 | 0) | 0;
            return c2 | 0;
          }
          function Eb(a2, c2, d2, e2) {
            a2 = a2 | 0;
            c2 = c2 | 0;
            d2 = d2 | 0;
            e2 = e2 | 0;
            var f = 0, g2 = 0, h = 0, i = 0, j = 0, k = 0, l = 0;
            h = cd(a2 | 0, c2 | 0, 52) | 0;
            G() | 0;
            h = h & 15;
            if (!((d2 | 0) < 16 & (h | 0) <= (d2 | 0))) {
              return;
            }
            if ((h | 0) == (d2 | 0)) {
              d2 = e2;
              b[d2 >> 2] = a2;
              b[d2 + 4 >> 2] = c2;
              return;
            }
            j = tc(7, d2 - h | 0) | 0;
            k = (j | 0) / 7 | 0;
            i = cd(a2 | 0, c2 | 0, 45) | 0;
            G() | 0;
            if (!(la(i & 127) | 0)) {
              g2 = 0;
            } else {
              a: do {
                if (!h) {
                  f = 0;
                } else {
                  g2 = 1;
                  while (1) {
                    f = cd(a2 | 0, c2 | 0, (15 - g2 | 0) * 3 | 0) | 0;
                    G() | 0;
                    f = f & 7;
                    if (f | 0) {
                      break a;
                    }
                    if (g2 >>> 0 < h >>> 0) {
                      g2 = g2 + 1 | 0;
                    } else {
                      f = 0;
                      break;
                    }
                  }
                }
              } while (0);
              g2 = (f | 0) == 0;
            }
            l = dd(h + 1 | 0, 0, 52) | 0;
            f = G() | 0 | c2 & -15728641;
            i = (14 - h | 0) * 3 | 0;
            c2 = dd(7, 0, i | 0) | 0;
            c2 = (l | a2) & ~c2;
            h = f & ~(G() | 0);
            Eb(c2, h, d2, e2);
            f = e2 + (k << 3) | 0;
            if (!g2) {
              l = dd(1, 0, i | 0) | 0;
              Eb(l | c2, G() | 0 | h, d2, f);
              l = f + (k << 3) | 0;
              j = dd(2, 0, i | 0) | 0;
              Eb(j | c2, G() | 0 | h, d2, l);
              l = l + (k << 3) | 0;
              j = dd(3, 0, i | 0) | 0;
              Eb(j | c2, G() | 0 | h, d2, l);
              l = l + (k << 3) | 0;
              j = dd(4, 0, i | 0) | 0;
              Eb(j | c2, G() | 0 | h, d2, l);
              l = l + (k << 3) | 0;
              j = dd(5, 0, i | 0) | 0;
              Eb(j | c2, G() | 0 | h, d2, l);
              j = dd(6, 0, i | 0) | 0;
              Eb(j | c2, G() | 0 | h, d2, l + (k << 3) | 0);
              return;
            }
            g2 = f + (k << 3) | 0;
            if ((j | 0) > 6) {
              j = f + 8 | 0;
              l = (g2 >>> 0 > j >>> 0 ? g2 : j) + -1 + (0 - f) | 0;
              hd(f | 0, 0, l + 8 & -8 | 0) | 0;
              f = j + (l >>> 3 << 3) | 0;
            }
            l = dd(2, 0, i | 0) | 0;
            Eb(l | c2, G() | 0 | h, d2, f);
            l = f + (k << 3) | 0;
            j = dd(3, 0, i | 0) | 0;
            Eb(j | c2, G() | 0 | h, d2, l);
            l = l + (k << 3) | 0;
            j = dd(4, 0, i | 0) | 0;
            Eb(j | c2, G() | 0 | h, d2, l);
            l = l + (k << 3) | 0;
            j = dd(5, 0, i | 0) | 0;
            Eb(j | c2, G() | 0 | h, d2, l);
            j = dd(6, 0, i | 0) | 0;
            Eb(j | c2, G() | 0 | h, d2, l + (k << 3) | 0);
            return;
          }
          function Fb(a2, b2) {
            a2 = a2 | 0;
            b2 = b2 | 0;
            var c2 = 0, d2 = 0, e2 = 0;
            e2 = cd(a2 | 0, b2 | 0, 45) | 0;
            G() | 0;
            if (!(la(e2 & 127) | 0)) {
              e2 = 0;
              return e2 | 0;
            }
            e2 = cd(a2 | 0, b2 | 0, 52) | 0;
            G() | 0;
            e2 = e2 & 15;
            a: do {
              if (!e2) {
                c2 = 0;
              } else {
                d2 = 1;
                while (1) {
                  c2 = cd(a2 | 0, b2 | 0, (15 - d2 | 0) * 3 | 0) | 0;
                  G() | 0;
                  c2 = c2 & 7;
                  if (c2 | 0) {
                    break a;
                  }
                  if (d2 >>> 0 < e2 >>> 0) {
                    d2 = d2 + 1 | 0;
                  } else {
                    c2 = 0;
                    break;
                  }
                }
              }
            } while (0);
            e2 = (c2 | 0) == 0 & 1;
            return e2 | 0;
          }
          function Gb(a2, b2, c2) {
            a2 = a2 | 0;
            b2 = b2 | 0;
            c2 = c2 | 0;
            var d2 = 0, e2 = 0;
            d2 = cd(a2 | 0, b2 | 0, 52) | 0;
            G() | 0;
            d2 = d2 & 15;
            if ((c2 | 0) < 16 & (d2 | 0) <= (c2 | 0)) {
              if ((d2 | 0) != (c2 | 0)) {
                e2 = dd(c2 | 0, 0, 52) | 0;
                a2 = e2 | a2;
                b2 = G() | 0 | b2 & -15728641;
                if ((d2 | 0) < (c2 | 0)) {
                  do {
                    e2 = dd(7, 0, (14 - d2 | 0) * 3 | 0) | 0;
                    d2 = d2 + 1 | 0;
                    a2 = a2 & ~e2;
                    b2 = b2 & ~(G() | 0);
                  } while ((d2 | 0) < (c2 | 0));
                }
              }
            } else {
              b2 = 0;
              a2 = 0;
            }
            F(b2 | 0);
            return a2 | 0;
          }
          function Hb(a2, c2, d2) {
            a2 = a2 | 0;
            c2 = c2 | 0;
            d2 = d2 | 0;
            var e2 = 0, f = 0, g2 = 0, h = 0, i = 0, j = 0, k = 0, l = 0, m = 0, n = 0, o = 0, p2 = 0, q2 = 0, r2 = 0, s2 = 0, t2 = 0, u2 = 0, v2 = 0, w2 = 0, x2 = 0, y2 = 0, z2 = 0;
            if (!d2) {
              y2 = 0;
              return y2 | 0;
            }
            f = a2;
            e2 = b[f >> 2] | 0;
            f = b[f + 4 >> 2] | 0;
            if (true & (f & 15728640 | 0) == 0) {
              if ((d2 | 0) <= 0) {
                y2 = 0;
                return y2 | 0;
              }
              y2 = c2;
              b[y2 >> 2] = e2;
              b[y2 + 4 >> 2] = f;
              if ((d2 | 0) == 1) {
                y2 = 0;
                return y2 | 0;
              }
              e2 = 1;
              do {
                w2 = a2 + (e2 << 3) | 0;
                x2 = b[w2 + 4 >> 2] | 0;
                y2 = c2 + (e2 << 3) | 0;
                b[y2 >> 2] = b[w2 >> 2];
                b[y2 + 4 >> 2] = x2;
                e2 = e2 + 1 | 0;
              } while ((e2 | 0) != (d2 | 0));
              e2 = 0;
              return e2 | 0;
            }
            w2 = d2 << 3;
            x2 = Wc(w2) | 0;
            if (!x2) {
              y2 = -3;
              return y2 | 0;
            }
            gd(x2 | 0, a2 | 0, w2 | 0) | 0;
            v2 = Yc(d2, 8) | 0;
            if (!v2) {
              Xc(x2);
              y2 = -3;
              return y2 | 0;
            }
            e2 = d2;
            a: while (1) {
              h = x2;
              l = b[h >> 2] | 0;
              h = b[h + 4 >> 2] | 0;
              t2 = cd(l | 0, h | 0, 52) | 0;
              G() | 0;
              t2 = t2 & 15;
              u2 = t2 + -1 | 0;
              s2 = (e2 | 0) > 0;
              b: do {
                if (s2) {
                  r2 = ((e2 | 0) < 0) << 31 >> 31;
                  p2 = dd(u2 | 0, 0, 52) | 0;
                  q2 = G() | 0;
                  if (u2 >>> 0 > 15) {
                    f = 0;
                    a2 = l;
                    d2 = h;
                    while (1) {
                      if (!((a2 | 0) == 0 & (d2 | 0) == 0)) {
                        g2 = cd(a2 | 0, d2 | 0, 52) | 0;
                        G() | 0;
                        g2 = g2 & 15;
                        i = (g2 | 0) < (u2 | 0);
                        g2 = (g2 | 0) == (u2 | 0);
                        k = i ? 0 : g2 ? a2 : 0;
                        a2 = i ? 0 : g2 ? d2 : 0;
                        d2 = bd(k | 0, a2 | 0, e2 | 0, r2 | 0) | 0;
                        G() | 0;
                        g2 = v2 + (d2 << 3) | 0;
                        i = g2;
                        j = b[i >> 2] | 0;
                        i = b[i + 4 >> 2] | 0;
                        if ((j | 0) == 0 & (i | 0) == 0) {
                          d2 = k;
                        } else {
                          p2 = 0;
                          o = d2;
                          n = i;
                          d2 = k;
                          while (1) {
                            if ((p2 | 0) > (e2 | 0)) {
                              y2 = 41;
                              break a;
                            }
                            if ((j | 0) == (d2 | 0) & (n & -117440513 | 0) == (a2 | 0)) {
                              k = cd(j | 0, n | 0, 56) | 0;
                              G() | 0;
                              k = k & 7;
                              m = k + 1 | 0;
                              q2 = cd(j | 0, n | 0, 45) | 0;
                              G() | 0;
                              c: do {
                                if (!(la(q2 & 127) | 0)) {
                                  i = 7;
                                } else {
                                  j = cd(j | 0, n | 0, 52) | 0;
                                  G() | 0;
                                  j = j & 15;
                                  if (!j) {
                                    i = 6;
                                    break;
                                  }
                                  i = 1;
                                  while (1) {
                                    q2 = dd(7, 0, (15 - i | 0) * 3 | 0) | 0;
                                    if (!((q2 & d2 | 0) == 0 & ((G() | 0) & a2 | 0) == 0)) {
                                      i = 7;
                                      break c;
                                    }
                                    if (i >>> 0 < j >>> 0) {
                                      i = i + 1 | 0;
                                    } else {
                                      i = 6;
                                      break;
                                    }
                                  }
                                }
                              } while (0);
                              if ((k + 2 | 0) >>> 0 > i >>> 0) {
                                y2 = 51;
                                break a;
                              }
                              q2 = dd(m | 0, 0, 56) | 0;
                              a2 = G() | 0 | a2 & -117440513;
                              i = g2;
                              b[i >> 2] = 0;
                              b[i + 4 >> 2] = 0;
                              i = o;
                              d2 = q2 | d2;
                            } else {
                              i = (o + 1 | 0) % (e2 | 0) | 0;
                            }
                            g2 = v2 + (i << 3) | 0;
                            n = g2;
                            j = b[n >> 2] | 0;
                            n = b[n + 4 >> 2] | 0;
                            if ((j | 0) == 0 & (n | 0) == 0) {
                              break;
                            } else {
                              p2 = p2 + 1 | 0;
                              o = i;
                            }
                          }
                        }
                        q2 = g2;
                        b[q2 >> 2] = d2;
                        b[q2 + 4 >> 2] = a2;
                      }
                      f = f + 1 | 0;
                      if ((f | 0) >= (e2 | 0)) {
                        break b;
                      }
                      d2 = x2 + (f << 3) | 0;
                      a2 = b[d2 >> 2] | 0;
                      d2 = b[d2 + 4 >> 2] | 0;
                    }
                  }
                  f = 0;
                  a2 = l;
                  d2 = h;
                  while (1) {
                    if (!((a2 | 0) == 0 & (d2 | 0) == 0)) {
                      i = cd(a2 | 0, d2 | 0, 52) | 0;
                      G() | 0;
                      i = i & 15;
                      if ((i | 0) >= (u2 | 0)) {
                        if ((i | 0) != (u2 | 0)) {
                          a2 = a2 | p2;
                          d2 = d2 & -15728641 | q2;
                          if (i >>> 0 >= t2 >>> 0) {
                            g2 = u2;
                            do {
                              o = dd(7, 0, (14 - g2 | 0) * 3 | 0) | 0;
                              g2 = g2 + 1 | 0;
                              a2 = o | a2;
                              d2 = G() | 0 | d2;
                            } while (g2 >>> 0 < i >>> 0);
                          }
                        }
                      } else {
                        a2 = 0;
                        d2 = 0;
                      }
                      i = bd(a2 | 0, d2 | 0, e2 | 0, r2 | 0) | 0;
                      G() | 0;
                      g2 = v2 + (i << 3) | 0;
                      j = g2;
                      k = b[j >> 2] | 0;
                      j = b[j + 4 >> 2] | 0;
                      if (!((k | 0) == 0 & (j | 0) == 0)) {
                        o = 0;
                        while (1) {
                          if ((o | 0) > (e2 | 0)) {
                            y2 = 41;
                            break a;
                          }
                          if ((k | 0) == (a2 | 0) & (j & -117440513 | 0) == (d2 | 0)) {
                            m = cd(k | 0, j | 0, 56) | 0;
                            G() | 0;
                            m = m & 7;
                            n = m + 1 | 0;
                            z2 = cd(k | 0, j | 0, 45) | 0;
                            G() | 0;
                            d: do {
                              if (!(la(z2 & 127) | 0)) {
                                j = 7;
                              } else {
                                k = cd(k | 0, j | 0, 52) | 0;
                                G() | 0;
                                k = k & 15;
                                if (!k) {
                                  j = 6;
                                  break;
                                }
                                j = 1;
                                while (1) {
                                  z2 = dd(7, 0, (15 - j | 0) * 3 | 0) | 0;
                                  if (!((z2 & a2 | 0) == 0 & ((G() | 0) & d2 | 0) == 0)) {
                                    j = 7;
                                    break d;
                                  }
                                  if (j >>> 0 < k >>> 0) {
                                    j = j + 1 | 0;
                                  } else {
                                    j = 6;
                                    break;
                                  }
                                }
                              }
                            } while (0);
                            if ((m + 2 | 0) >>> 0 > j >>> 0) {
                              y2 = 51;
                              break a;
                            }
                            z2 = dd(n | 0, 0, 56) | 0;
                            d2 = G() | 0 | d2 & -117440513;
                            n = g2;
                            b[n >> 2] = 0;
                            b[n + 4 >> 2] = 0;
                            a2 = z2 | a2;
                          } else {
                            i = (i + 1 | 0) % (e2 | 0) | 0;
                          }
                          g2 = v2 + (i << 3) | 0;
                          j = g2;
                          k = b[j >> 2] | 0;
                          j = b[j + 4 >> 2] | 0;
                          if ((k | 0) == 0 & (j | 0) == 0) {
                            break;
                          } else {
                            o = o + 1 | 0;
                          }
                        }
                      }
                      z2 = g2;
                      b[z2 >> 2] = a2;
                      b[z2 + 4 >> 2] = d2;
                    }
                    f = f + 1 | 0;
                    if ((f | 0) >= (e2 | 0)) {
                      break b;
                    }
                    d2 = x2 + (f << 3) | 0;
                    a2 = b[d2 >> 2] | 0;
                    d2 = b[d2 + 4 >> 2] | 0;
                  }
                }
              } while (0);
              if ((e2 + 5 | 0) >>> 0 < 11) {
                y2 = 99;
                break;
              }
              q2 = Yc((e2 | 0) / 6 | 0, 8) | 0;
              if (!q2) {
                y2 = 58;
                break;
              }
              e: do {
                if (s2) {
                  o = 0;
                  n = 0;
                  do {
                    i = v2 + (o << 3) | 0;
                    a2 = i;
                    f = b[a2 >> 2] | 0;
                    a2 = b[a2 + 4 >> 2] | 0;
                    if (!((f | 0) == 0 & (a2 | 0) == 0)) {
                      j = cd(f | 0, a2 | 0, 56) | 0;
                      G() | 0;
                      j = j & 7;
                      d2 = j + 1 | 0;
                      k = a2 & -117440513;
                      z2 = cd(f | 0, a2 | 0, 45) | 0;
                      G() | 0;
                      f: do {
                        if (la(z2 & 127) | 0) {
                          m = cd(f | 0, a2 | 0, 52) | 0;
                          G() | 0;
                          m = m & 15;
                          if (m | 0) {
                            g2 = 1;
                            while (1) {
                              z2 = dd(7, 0, (15 - g2 | 0) * 3 | 0) | 0;
                              if (!((f & z2 | 0) == 0 & (k & (G() | 0) | 0) == 0)) {
                                break f;
                              }
                              if (g2 >>> 0 < m >>> 0) {
                                g2 = g2 + 1 | 0;
                              } else {
                                break;
                              }
                            }
                          }
                          a2 = dd(d2 | 0, 0, 56) | 0;
                          f = a2 | f;
                          a2 = G() | 0 | k;
                          d2 = i;
                          b[d2 >> 2] = f;
                          b[d2 + 4 >> 2] = a2;
                          d2 = j + 2 | 0;
                        }
                      } while (0);
                      if ((d2 | 0) == 7) {
                        z2 = q2 + (n << 3) | 0;
                        b[z2 >> 2] = f;
                        b[z2 + 4 >> 2] = a2 & -117440513;
                        n = n + 1 | 0;
                      }
                    }
                    o = o + 1 | 0;
                  } while ((o | 0) != (e2 | 0));
                  if (s2) {
                    p2 = ((e2 | 0) < 0) << 31 >> 31;
                    m = dd(u2 | 0, 0, 52) | 0;
                    o = G() | 0;
                    if (u2 >>> 0 > 15) {
                      a2 = 0;
                      f = 0;
                      while (1) {
                        do {
                          if (!((l | 0) == 0 & (h | 0) == 0)) {
                            j = cd(l | 0, h | 0, 52) | 0;
                            G() | 0;
                            j = j & 15;
                            g2 = (j | 0) < (u2 | 0);
                            j = (j | 0) == (u2 | 0);
                            i = g2 ? 0 : j ? l : 0;
                            j = g2 ? 0 : j ? h : 0;
                            g2 = bd(i | 0, j | 0, e2 | 0, p2 | 0) | 0;
                            G() | 0;
                            d2 = 0;
                            while (1) {
                              if ((d2 | 0) > (e2 | 0)) {
                                y2 = 98;
                                break a;
                              }
                              z2 = v2 + (g2 << 3) | 0;
                              k = b[z2 + 4 >> 2] | 0;
                              if ((k & -117440513 | 0) == (j | 0) ? (b[z2 >> 2] | 0) == (i | 0) : 0) {
                                y2 = 70;
                                break;
                              }
                              g2 = (g2 + 1 | 0) % (e2 | 0) | 0;
                              z2 = v2 + (g2 << 3) | 0;
                              if ((b[z2 >> 2] | 0) == (i | 0) ? (b[z2 + 4 >> 2] | 0) == (j | 0) : 0) {
                                break;
                              } else {
                                d2 = d2 + 1 | 0;
                              }
                            }
                            if ((y2 | 0) == 70 ? (y2 = 0, true & (k & 117440512 | 0) == 100663296) : 0) {
                              break;
                            }
                            z2 = c2 + (f << 3) | 0;
                            b[z2 >> 2] = l;
                            b[z2 + 4 >> 2] = h;
                            f = f + 1 | 0;
                          }
                        } while (0);
                        a2 = a2 + 1 | 0;
                        if ((a2 | 0) >= (e2 | 0)) {
                          e2 = n;
                          break e;
                        }
                        h = x2 + (a2 << 3) | 0;
                        l = b[h >> 2] | 0;
                        h = b[h + 4 >> 2] | 0;
                      }
                    }
                    a2 = 0;
                    f = 0;
                    while (1) {
                      do {
                        if (!((l | 0) == 0 & (h | 0) == 0)) {
                          j = cd(l | 0, h | 0, 52) | 0;
                          G() | 0;
                          j = j & 15;
                          if ((j | 0) >= (u2 | 0)) {
                            if ((j | 0) != (u2 | 0)) {
                              d2 = l | m;
                              g2 = h & -15728641 | o;
                              if (j >>> 0 < t2 >>> 0) {
                                j = g2;
                              } else {
                                i = u2;
                                do {
                                  z2 = dd(7, 0, (14 - i | 0) * 3 | 0) | 0;
                                  i = i + 1 | 0;
                                  d2 = z2 | d2;
                                  g2 = G() | 0 | g2;
                                } while (i >>> 0 < j >>> 0);
                                j = g2;
                              }
                            } else {
                              d2 = l;
                              j = h;
                            }
                          } else {
                            d2 = 0;
                            j = 0;
                          }
                          i = bd(d2 | 0, j | 0, e2 | 0, p2 | 0) | 0;
                          G() | 0;
                          g2 = 0;
                          while (1) {
                            if ((g2 | 0) > (e2 | 0)) {
                              y2 = 98;
                              break a;
                            }
                            z2 = v2 + (i << 3) | 0;
                            k = b[z2 + 4 >> 2] | 0;
                            if ((k & -117440513 | 0) == (j | 0) ? (b[z2 >> 2] | 0) == (d2 | 0) : 0) {
                              y2 = 93;
                              break;
                            }
                            i = (i + 1 | 0) % (e2 | 0) | 0;
                            z2 = v2 + (i << 3) | 0;
                            if ((b[z2 >> 2] | 0) == (d2 | 0) ? (b[z2 + 4 >> 2] | 0) == (j | 0) : 0) {
                              break;
                            } else {
                              g2 = g2 + 1 | 0;
                            }
                          }
                          if ((y2 | 0) == 93 ? (y2 = 0, true & (k & 117440512 | 0) == 100663296) : 0) {
                            break;
                          }
                          z2 = c2 + (f << 3) | 0;
                          b[z2 >> 2] = l;
                          b[z2 + 4 >> 2] = h;
                          f = f + 1 | 0;
                        }
                      } while (0);
                      a2 = a2 + 1 | 0;
                      if ((a2 | 0) >= (e2 | 0)) {
                        e2 = n;
                        break e;
                      }
                      h = x2 + (a2 << 3) | 0;
                      l = b[h >> 2] | 0;
                      h = b[h + 4 >> 2] | 0;
                    }
                  } else {
                    f = 0;
                    e2 = n;
                  }
                } else {
                  f = 0;
                  e2 = 0;
                }
              } while (0);
              hd(v2 | 0, 0, w2 | 0) | 0;
              gd(x2 | 0, q2 | 0, e2 << 3 | 0) | 0;
              Xc(q2);
              if (!e2) {
                break;
              } else {
                c2 = c2 + (f << 3) | 0;
              }
            }
            if ((y2 | 0) == 41) {
              Xc(x2);
              Xc(v2);
              z2 = -1;
              return z2 | 0;
            } else if ((y2 | 0) == 51) {
              Xc(x2);
              Xc(v2);
              z2 = -2;
              return z2 | 0;
            } else if ((y2 | 0) == 58) {
              Xc(x2);
              Xc(v2);
              z2 = -3;
              return z2 | 0;
            } else if ((y2 | 0) == 98) {
              Xc(q2);
              Xc(x2);
              Xc(v2);
              z2 = -1;
              return z2 | 0;
            } else if ((y2 | 0) == 99) {
              gd(c2 | 0, x2 | 0, e2 << 3 | 0) | 0;
            }
            Xc(x2);
            Xc(v2);
            z2 = 0;
            return z2 | 0;
          }
          function Ib(a2, c2, d2, e2, f) {
            a2 = a2 | 0;
            c2 = c2 | 0;
            d2 = d2 | 0;
            e2 = e2 | 0;
            f = f | 0;
            var g2 = 0, h = 0, i = 0, j = 0, k = 0, l = 0;
            if ((c2 | 0) <= 0) {
              f = 0;
              return f | 0;
            }
            if ((f | 0) >= 16) {
              g2 = 0;
              while (1) {
                l = a2 + (g2 << 3) | 0;
                if (!((b[l >> 2] | 0) == 0 & (b[l + 4 >> 2] | 0) == 0)) {
                  g2 = 14;
                  break;
                }
                g2 = g2 + 1 | 0;
                if ((g2 | 0) >= (c2 | 0)) {
                  h = 0;
                  g2 = 16;
                  break;
                }
              }
              if ((g2 | 0) == 14) {
                return ((e2 | 0) > 0 ? -2 : -1) | 0;
              } else if ((g2 | 0) == 16) {
                return h | 0;
              }
            }
            g2 = 0;
            l = 0;
            a: while (1) {
              k = a2 + (l << 3) | 0;
              i = k;
              h = b[i >> 2] | 0;
              i = b[i + 4 >> 2] | 0;
              do {
                if (!((h | 0) == 0 & (i | 0) == 0)) {
                  if ((g2 | 0) >= (e2 | 0)) {
                    h = -1;
                    g2 = 16;
                    break a;
                  }
                  j = cd(h | 0, i | 0, 52) | 0;
                  G() | 0;
                  j = j & 15;
                  if ((j | 0) > (f | 0)) {
                    h = -2;
                    g2 = 16;
                    break a;
                  }
                  if ((j | 0) == (f | 0)) {
                    k = d2 + (g2 << 3) | 0;
                    b[k >> 2] = h;
                    b[k + 4 >> 2] = i;
                    g2 = g2 + 1 | 0;
                    break;
                  }
                  h = (tc(7, f - j | 0) | 0) + g2 | 0;
                  if ((h | 0) > (e2 | 0)) {
                    h = -1;
                    g2 = 16;
                    break a;
                  }
                  Eb(b[k >> 2] | 0, b[k + 4 >> 2] | 0, f, d2 + (g2 << 3) | 0);
                  g2 = h;
                }
              } while (0);
              l = l + 1 | 0;
              if ((l | 0) >= (c2 | 0)) {
                h = 0;
                g2 = 16;
                break;
              }
            }
            if ((g2 | 0) == 16) {
              return h | 0;
            }
            return 0;
          }
          function Jb(a2, c2, d2) {
            a2 = a2 | 0;
            c2 = c2 | 0;
            d2 = d2 | 0;
            var e2 = 0, f = 0, g2 = 0, h = 0;
            if ((c2 | 0) <= 0) {
              d2 = 0;
              return d2 | 0;
            }
            if ((d2 | 0) >= 16) {
              e2 = 0;
              while (1) {
                h = a2 + (e2 << 3) | 0;
                if (!((b[h >> 2] | 0) == 0 & (b[h + 4 >> 2] | 0) == 0)) {
                  e2 = -1;
                  f = 13;
                  break;
                }
                e2 = e2 + 1 | 0;
                if ((e2 | 0) >= (c2 | 0)) {
                  e2 = 0;
                  f = 13;
                  break;
                }
              }
              if ((f | 0) == 13) {
                return e2 | 0;
              }
            }
            e2 = 0;
            h = 0;
            a: while (1) {
              f = a2 + (h << 3) | 0;
              g2 = b[f >> 2] | 0;
              f = b[f + 4 >> 2] | 0;
              do {
                if (!((g2 | 0) == 0 & (f | 0) == 0)) {
                  f = cd(g2 | 0, f | 0, 52) | 0;
                  G() | 0;
                  f = f & 15;
                  if ((f | 0) > (d2 | 0)) {
                    e2 = -1;
                    f = 13;
                    break a;
                  }
                  if ((f | 0) == (d2 | 0)) {
                    e2 = e2 + 1 | 0;
                    break;
                  } else {
                    e2 = (tc(7, d2 - f | 0) | 0) + e2 | 0;
                    break;
                  }
                }
              } while (0);
              h = h + 1 | 0;
              if ((h | 0) >= (c2 | 0)) {
                f = 13;
                break;
              }
            }
            if ((f | 0) == 13) {
              return e2 | 0;
            }
            return 0;
          }
          function Kb(a2, b2) {
            a2 = a2 | 0;
            b2 = b2 | 0;
            b2 = cd(a2 | 0, b2 | 0, 52) | 0;
            G() | 0;
            return b2 & 1 | 0;
          }
          function Lb(a2, b2) {
            a2 = a2 | 0;
            b2 = b2 | 0;
            var c2 = 0, d2 = 0, e2 = 0;
            e2 = cd(a2 | 0, b2 | 0, 52) | 0;
            G() | 0;
            e2 = e2 & 15;
            if (!e2) {
              e2 = 0;
              return e2 | 0;
            }
            d2 = 1;
            while (1) {
              c2 = cd(a2 | 0, b2 | 0, (15 - d2 | 0) * 3 | 0) | 0;
              G() | 0;
              c2 = c2 & 7;
              if (c2 | 0) {
                d2 = 5;
                break;
              }
              if (d2 >>> 0 < e2 >>> 0) {
                d2 = d2 + 1 | 0;
              } else {
                c2 = 0;
                d2 = 5;
                break;
              }
            }
            if ((d2 | 0) == 5) {
              return c2 | 0;
            }
            return 0;
          }
          function Mb(a2, b2) {
            a2 = a2 | 0;
            b2 = b2 | 0;
            var c2 = 0, d2 = 0, e2 = 0, f = 0, g2 = 0, h = 0, i = 0;
            i = cd(a2 | 0, b2 | 0, 52) | 0;
            G() | 0;
            i = i & 15;
            if (!i) {
              h = b2;
              i = a2;
              F(h | 0);
              return i | 0;
            }
            h = 1;
            c2 = 0;
            while (1) {
              f = (15 - h | 0) * 3 | 0;
              d2 = dd(7, 0, f | 0) | 0;
              e2 = G() | 0;
              g2 = cd(a2 | 0, b2 | 0, f | 0) | 0;
              G() | 0;
              f = dd(Pa(g2 & 7) | 0, 0, f | 0) | 0;
              g2 = G() | 0;
              a2 = f | a2 & ~d2;
              b2 = g2 | b2 & ~e2;
              a: do {
                if (!c2) {
                  if (!((f & d2 | 0) == 0 & (g2 & e2 | 0) == 0)) {
                    d2 = cd(a2 | 0, b2 | 0, 52) | 0;
                    G() | 0;
                    d2 = d2 & 15;
                    if (!d2) {
                      c2 = 1;
                    } else {
                      c2 = 1;
                      b: while (1) {
                        g2 = cd(a2 | 0, b2 | 0, (15 - c2 | 0) * 3 | 0) | 0;
                        G() | 0;
                        switch (g2 & 7) {
                          case 1:
                            break b;
                          case 0:
                            break;
                          default: {
                            c2 = 1;
                            break a;
                          }
                        }
                        if (c2 >>> 0 < d2 >>> 0) {
                          c2 = c2 + 1 | 0;
                        } else {
                          c2 = 1;
                          break a;
                        }
                      }
                      c2 = 1;
                      while (1) {
                        g2 = (15 - c2 | 0) * 3 | 0;
                        e2 = cd(a2 | 0, b2 | 0, g2 | 0) | 0;
                        G() | 0;
                        f = dd(7, 0, g2 | 0) | 0;
                        b2 = b2 & ~(G() | 0);
                        g2 = dd(Pa(e2 & 7) | 0, 0, g2 | 0) | 0;
                        a2 = a2 & ~f | g2;
                        b2 = b2 | (G() | 0);
                        if (c2 >>> 0 < d2 >>> 0) {
                          c2 = c2 + 1 | 0;
                        } else {
                          c2 = 1;
                          break;
                        }
                      }
                    }
                  } else {
                    c2 = 0;
                  }
                }
              } while (0);
              if (h >>> 0 < i >>> 0) {
                h = h + 1 | 0;
              } else {
                break;
              }
            }
            F(b2 | 0);
            return a2 | 0;
          }
          function Nb(a2, b2) {
            a2 = a2 | 0;
            b2 = b2 | 0;
            var c2 = 0, d2 = 0, e2 = 0, f = 0, g2 = 0;
            d2 = cd(a2 | 0, b2 | 0, 52) | 0;
            G() | 0;
            d2 = d2 & 15;
            if (!d2) {
              c2 = b2;
              d2 = a2;
              F(c2 | 0);
              return d2 | 0;
            }
            c2 = 1;
            while (1) {
              f = (15 - c2 | 0) * 3 | 0;
              g2 = cd(a2 | 0, b2 | 0, f | 0) | 0;
              G() | 0;
              e2 = dd(7, 0, f | 0) | 0;
              b2 = b2 & ~(G() | 0);
              f = dd(Pa(g2 & 7) | 0, 0, f | 0) | 0;
              a2 = f | a2 & ~e2;
              b2 = G() | 0 | b2;
              if (c2 >>> 0 < d2 >>> 0) {
                c2 = c2 + 1 | 0;
              } else {
                break;
              }
            }
            F(b2 | 0);
            return a2 | 0;
          }
          function Ob(a2, b2) {
            a2 = a2 | 0;
            b2 = b2 | 0;
            var c2 = 0, d2 = 0, e2 = 0, f = 0, g2 = 0, h = 0, i = 0;
            i = cd(a2 | 0, b2 | 0, 52) | 0;
            G() | 0;
            i = i & 15;
            if (!i) {
              h = b2;
              i = a2;
              F(h | 0);
              return i | 0;
            }
            h = 1;
            c2 = 0;
            while (1) {
              f = (15 - h | 0) * 3 | 0;
              d2 = dd(7, 0, f | 0) | 0;
              e2 = G() | 0;
              g2 = cd(a2 | 0, b2 | 0, f | 0) | 0;
              G() | 0;
              f = dd(Qa(g2 & 7) | 0, 0, f | 0) | 0;
              g2 = G() | 0;
              a2 = f | a2 & ~d2;
              b2 = g2 | b2 & ~e2;
              a: do {
                if (!c2) {
                  if (!((f & d2 | 0) == 0 & (g2 & e2 | 0) == 0)) {
                    d2 = cd(a2 | 0, b2 | 0, 52) | 0;
                    G() | 0;
                    d2 = d2 & 15;
                    if (!d2) {
                      c2 = 1;
                    } else {
                      c2 = 1;
                      b: while (1) {
                        g2 = cd(a2 | 0, b2 | 0, (15 - c2 | 0) * 3 | 0) | 0;
                        G() | 0;
                        switch (g2 & 7) {
                          case 1:
                            break b;
                          case 0:
                            break;
                          default: {
                            c2 = 1;
                            break a;
                          }
                        }
                        if (c2 >>> 0 < d2 >>> 0) {
                          c2 = c2 + 1 | 0;
                        } else {
                          c2 = 1;
                          break a;
                        }
                      }
                      c2 = 1;
                      while (1) {
                        e2 = (15 - c2 | 0) * 3 | 0;
                        f = dd(7, 0, e2 | 0) | 0;
                        g2 = b2 & ~(G() | 0);
                        b2 = cd(a2 | 0, b2 | 0, e2 | 0) | 0;
                        G() | 0;
                        b2 = dd(Qa(b2 & 7) | 0, 0, e2 | 0) | 0;
                        a2 = a2 & ~f | b2;
                        b2 = g2 | (G() | 0);
                        if (c2 >>> 0 < d2 >>> 0) {
                          c2 = c2 + 1 | 0;
                        } else {
                          c2 = 1;
                          break;
                        }
                      }
                    }
                  } else {
                    c2 = 0;
                  }
                }
              } while (0);
              if (h >>> 0 < i >>> 0) {
                h = h + 1 | 0;
              } else {
                break;
              }
            }
            F(b2 | 0);
            return a2 | 0;
          }
          function Pb(a2, b2) {
            a2 = a2 | 0;
            b2 = b2 | 0;
            var c2 = 0, d2 = 0, e2 = 0, f = 0, g2 = 0;
            d2 = cd(a2 | 0, b2 | 0, 52) | 0;
            G() | 0;
            d2 = d2 & 15;
            if (!d2) {
              c2 = b2;
              d2 = a2;
              F(c2 | 0);
              return d2 | 0;
            }
            c2 = 1;
            while (1) {
              g2 = (15 - c2 | 0) * 3 | 0;
              f = dd(7, 0, g2 | 0) | 0;
              e2 = b2 & ~(G() | 0);
              b2 = cd(a2 | 0, b2 | 0, g2 | 0) | 0;
              G() | 0;
              b2 = dd(Qa(b2 & 7) | 0, 0, g2 | 0) | 0;
              a2 = b2 | a2 & ~f;
              b2 = G() | 0 | e2;
              if (c2 >>> 0 < d2 >>> 0) {
                c2 = c2 + 1 | 0;
              } else {
                break;
              }
            }
            F(b2 | 0);
            return a2 | 0;
          }
          function Qb(a2, c2) {
            a2 = a2 | 0;
            c2 = c2 | 0;
            var d2 = 0, e2 = 0, f = 0, g2 = 0, h = 0, i = 0, j = 0, k = 0, l = 0;
            j = S;
            S = S + 64 | 0;
            i = j + 40 | 0;
            e2 = j + 24 | 0;
            f = j + 12 | 0;
            g2 = j;
            dd(c2 | 0, 0, 52) | 0;
            d2 = G() | 0 | 134225919;
            if (!c2) {
              if ((b[a2 + 4 >> 2] | 0) > 2) {
                h = 0;
                i = 0;
                F(h | 0);
                S = j;
                return i | 0;
              }
              if ((b[a2 + 8 >> 2] | 0) > 2) {
                h = 0;
                i = 0;
                F(h | 0);
                S = j;
                return i | 0;
              }
              if ((b[a2 + 12 >> 2] | 0) > 2) {
                h = 0;
                i = 0;
                F(h | 0);
                S = j;
                return i | 0;
              }
              dd(na(a2) | 0, 0, 45) | 0;
              h = G() | 0 | d2;
              i = -1;
              F(h | 0);
              S = j;
              return i | 0;
            }
            b[i >> 2] = b[a2 >> 2];
            b[i + 4 >> 2] = b[a2 + 4 >> 2];
            b[i + 8 >> 2] = b[a2 + 8 >> 2];
            b[i + 12 >> 2] = b[a2 + 12 >> 2];
            h = i + 4 | 0;
            if ((c2 | 0) > 0) {
              a2 = -1;
              while (1) {
                b[e2 >> 2] = b[h >> 2];
                b[e2 + 4 >> 2] = b[h + 4 >> 2];
                b[e2 + 8 >> 2] = b[h + 8 >> 2];
                if (!(c2 & 1)) {
                  Ja(h);
                  b[f >> 2] = b[h >> 2];
                  b[f + 4 >> 2] = b[h + 4 >> 2];
                  b[f + 8 >> 2] = b[h + 8 >> 2];
                  La(f);
                } else {
                  Ia(h);
                  b[f >> 2] = b[h >> 2];
                  b[f + 4 >> 2] = b[h + 4 >> 2];
                  b[f + 8 >> 2] = b[h + 8 >> 2];
                  Ka(f);
                }
                Fa(e2, f, g2);
                Ca(g2);
                l = (15 - c2 | 0) * 3 | 0;
                k = dd(7, 0, l | 0) | 0;
                d2 = d2 & ~(G() | 0);
                l = dd(Ha(g2) | 0, 0, l | 0) | 0;
                a2 = l | a2 & ~k;
                d2 = G() | 0 | d2;
                if ((c2 | 0) > 1) {
                  c2 = c2 + -1 | 0;
                } else {
                  break;
                }
              }
            } else {
              a2 = -1;
            }
            a: do {
              if (((b[h >> 2] | 0) <= 2 ? (b[i + 8 >> 2] | 0) <= 2 : 0) ? (b[i + 12 >> 2] | 0) <= 2 : 0) {
                e2 = na(i) | 0;
                c2 = dd(e2 | 0, 0, 45) | 0;
                c2 = c2 | a2;
                a2 = G() | 0 | d2 & -1040385;
                g2 = oa(i) | 0;
                if (!(la(e2) | 0)) {
                  if ((g2 | 0) <= 0) {
                    break;
                  }
                  f = 0;
                  while (1) {
                    e2 = cd(c2 | 0, a2 | 0, 52) | 0;
                    G() | 0;
                    e2 = e2 & 15;
                    if (e2) {
                      d2 = 1;
                      while (1) {
                        l = (15 - d2 | 0) * 3 | 0;
                        i = cd(c2 | 0, a2 | 0, l | 0) | 0;
                        G() | 0;
                        k = dd(7, 0, l | 0) | 0;
                        a2 = a2 & ~(G() | 0);
                        l = dd(Pa(i & 7) | 0, 0, l | 0) | 0;
                        c2 = c2 & ~k | l;
                        a2 = a2 | (G() | 0);
                        if (d2 >>> 0 < e2 >>> 0) {
                          d2 = d2 + 1 | 0;
                        } else {
                          break;
                        }
                      }
                    }
                    f = f + 1 | 0;
                    if ((f | 0) == (g2 | 0)) {
                      break a;
                    }
                  }
                }
                f = cd(c2 | 0, a2 | 0, 52) | 0;
                G() | 0;
                f = f & 15;
                b: do {
                  if (f) {
                    d2 = 1;
                    c: while (1) {
                      l = cd(c2 | 0, a2 | 0, (15 - d2 | 0) * 3 | 0) | 0;
                      G() | 0;
                      switch (l & 7) {
                        case 1:
                          break c;
                        case 0:
                          break;
                        default:
                          break b;
                      }
                      if (d2 >>> 0 < f >>> 0) {
                        d2 = d2 + 1 | 0;
                      } else {
                        break b;
                      }
                    }
                    if (ra(e2, b[i >> 2] | 0) | 0) {
                      d2 = 1;
                      while (1) {
                        i = (15 - d2 | 0) * 3 | 0;
                        k = dd(7, 0, i | 0) | 0;
                        l = a2 & ~(G() | 0);
                        a2 = cd(c2 | 0, a2 | 0, i | 0) | 0;
                        G() | 0;
                        a2 = dd(Qa(a2 & 7) | 0, 0, i | 0) | 0;
                        c2 = c2 & ~k | a2;
                        a2 = l | (G() | 0);
                        if (d2 >>> 0 < f >>> 0) {
                          d2 = d2 + 1 | 0;
                        } else {
                          break;
                        }
                      }
                    } else {
                      d2 = 1;
                      while (1) {
                        l = (15 - d2 | 0) * 3 | 0;
                        i = cd(c2 | 0, a2 | 0, l | 0) | 0;
                        G() | 0;
                        k = dd(7, 0, l | 0) | 0;
                        a2 = a2 & ~(G() | 0);
                        l = dd(Pa(i & 7) | 0, 0, l | 0) | 0;
                        c2 = c2 & ~k | l;
                        a2 = a2 | (G() | 0);
                        if (d2 >>> 0 < f >>> 0) {
                          d2 = d2 + 1 | 0;
                        } else {
                          break;
                        }
                      }
                    }
                  }
                } while (0);
                if ((g2 | 0) > 0) {
                  d2 = 0;
                  do {
                    c2 = Mb(c2, a2) | 0;
                    a2 = G() | 0;
                    d2 = d2 + 1 | 0;
                  } while ((d2 | 0) != (g2 | 0));
                }
              } else {
                c2 = 0;
                a2 = 0;
              }
            } while (0);
            k = a2;
            l = c2;
            F(k | 0);
            S = j;
            return l | 0;
          }
          function Rb(a2) {
            a2 = a2 | 0;
            return (a2 | 0) % 2 | 0 | 0;
          }
          function Sb(a2, c2) {
            a2 = a2 | 0;
            c2 = c2 | 0;
            var d2 = 0, e2 = 0;
            e2 = S;
            S = S + 16 | 0;
            d2 = e2;
            if ((c2 >>> 0 <= 15 ? !((b[a2 + 4 >> 2] & 2146435072 | 0) == 2146435072) : 0) ? !((b[a2 + 8 + 4 >> 2] & 2146435072 | 0) == 2146435072) : 0) {
              Ya(a2, c2, d2);
              c2 = Qb(d2, c2) | 0;
              a2 = G() | 0;
            } else {
              a2 = 0;
              c2 = 0;
            }
            F(a2 | 0);
            S = e2;
            return c2 | 0;
          }
          function Tb(a2, c2, d2) {
            a2 = a2 | 0;
            c2 = c2 | 0;
            d2 = d2 | 0;
            var e2 = 0, f = 0, g2 = 0, h = 0;
            f = d2 + 4 | 0;
            g2 = cd(a2 | 0, c2 | 0, 52) | 0;
            G() | 0;
            g2 = g2 & 15;
            h = cd(a2 | 0, c2 | 0, 45) | 0;
            G() | 0;
            e2 = (g2 | 0) == 0;
            if (!(la(h & 127) | 0)) {
              if (e2) {
                h = 0;
                return h | 0;
              }
              if ((b[f >> 2] | 0) == 0 ? (b[d2 + 8 >> 2] | 0) == 0 : 0) {
                e2 = (b[d2 + 12 >> 2] | 0) != 0 & 1;
              } else {
                e2 = 1;
              }
            } else if (e2) {
              h = 1;
              return h | 0;
            } else {
              e2 = 1;
            }
            d2 = 1;
            while (1) {
              if (!(d2 & 1)) {
                La(f);
              } else {
                Ka(f);
              }
              h = cd(a2 | 0, c2 | 0, (15 - d2 | 0) * 3 | 0) | 0;
              G() | 0;
              Ma(f, h & 7);
              if (d2 >>> 0 < g2 >>> 0) {
                d2 = d2 + 1 | 0;
              } else {
                break;
              }
            }
            return e2 | 0;
          }
          function Ub(a2, c2, d2) {
            a2 = a2 | 0;
            c2 = c2 | 0;
            d2 = d2 | 0;
            var e2 = 0, f = 0, g2 = 0, h = 0, i = 0, j = 0, k = 0, l = 0;
            l = S;
            S = S + 16 | 0;
            j = l;
            k = cd(a2 | 0, c2 | 0, 45) | 0;
            G() | 0;
            k = k & 127;
            a: do {
              if ((la(k) | 0) != 0 ? (g2 = cd(a2 | 0, c2 | 0, 52) | 0, G() | 0, g2 = g2 & 15, (g2 | 0) != 0) : 0) {
                e2 = 1;
                b: while (1) {
                  i = cd(a2 | 0, c2 | 0, (15 - e2 | 0) * 3 | 0) | 0;
                  G() | 0;
                  switch (i & 7) {
                    case 5:
                      break b;
                    case 0:
                      break;
                    default: {
                      e2 = c2;
                      break a;
                    }
                  }
                  if (e2 >>> 0 < g2 >>> 0) {
                    e2 = e2 + 1 | 0;
                  } else {
                    e2 = c2;
                    break a;
                  }
                }
                f = 1;
                e2 = c2;
                while (1) {
                  c2 = (15 - f | 0) * 3 | 0;
                  h = dd(7, 0, c2 | 0) | 0;
                  i = e2 & ~(G() | 0);
                  e2 = cd(a2 | 0, e2 | 0, c2 | 0) | 0;
                  G() | 0;
                  e2 = dd(Qa(e2 & 7) | 0, 0, c2 | 0) | 0;
                  a2 = a2 & ~h | e2;
                  e2 = i | (G() | 0);
                  if (f >>> 0 < g2 >>> 0) {
                    f = f + 1 | 0;
                  } else {
                    break;
                  }
                }
              } else {
                e2 = c2;
              }
            } while (0);
            i = 7728 + (k * 28 | 0) | 0;
            b[d2 >> 2] = b[i >> 2];
            b[d2 + 4 >> 2] = b[i + 4 >> 2];
            b[d2 + 8 >> 2] = b[i + 8 >> 2];
            b[d2 + 12 >> 2] = b[i + 12 >> 2];
            if (!(Tb(a2, e2, d2) | 0)) {
              S = l;
              return;
            }
            h = d2 + 4 | 0;
            b[j >> 2] = b[h >> 2];
            b[j + 4 >> 2] = b[h + 4 >> 2];
            b[j + 8 >> 2] = b[h + 8 >> 2];
            g2 = cd(a2 | 0, e2 | 0, 52) | 0;
            G() | 0;
            i = g2 & 15;
            if (!(g2 & 1)) {
              g2 = i;
            } else {
              La(h);
              g2 = i + 1 | 0;
            }
            if (!(la(k) | 0)) {
              e2 = 0;
            } else {
              c: do {
                if (!i) {
                  e2 = 0;
                } else {
                  c2 = 1;
                  while (1) {
                    f = cd(a2 | 0, e2 | 0, (15 - c2 | 0) * 3 | 0) | 0;
                    G() | 0;
                    f = f & 7;
                    if (f | 0) {
                      e2 = f;
                      break c;
                    }
                    if (c2 >>> 0 < i >>> 0) {
                      c2 = c2 + 1 | 0;
                    } else {
                      e2 = 0;
                      break;
                    }
                  }
                }
              } while (0);
              e2 = (e2 | 0) == 4 & 1;
            }
            if (!(cb(d2, g2, e2, 0) | 0)) {
              if ((g2 | 0) != (i | 0)) {
                b[h >> 2] = b[j >> 2];
                b[h + 4 >> 2] = b[j + 4 >> 2];
                b[h + 8 >> 2] = b[j + 8 >> 2];
              }
            } else {
              if (la(k) | 0) {
                do {
                } while ((cb(d2, g2, 0, 0) | 0) != 0);
              }
              if ((g2 | 0) != (i | 0)) {
                Ja(h);
              }
            }
            S = l;
            return;
          }
          function Vb(a2, b2, c2) {
            a2 = a2 | 0;
            b2 = b2 | 0;
            c2 = c2 | 0;
            var d2 = 0, e2 = 0;
            d2 = S;
            S = S + 16 | 0;
            e2 = d2;
            Ub(a2, b2, e2);
            b2 = cd(a2 | 0, b2 | 0, 52) | 0;
            G() | 0;
            $a(e2, b2 & 15, c2);
            S = d2;
            return;
          }
          function Wb(a2, b2, c2) {
            a2 = a2 | 0;
            b2 = b2 | 0;
            c2 = c2 | 0;
            var d2 = 0, e2 = 0, f = 0, g2 = 0, h = 0;
            g2 = S;
            S = S + 16 | 0;
            f = g2;
            Ub(a2, b2, f);
            d2 = cd(a2 | 0, b2 | 0, 45) | 0;
            G() | 0;
            d2 = (la(d2 & 127) | 0) == 0;
            e2 = cd(a2 | 0, b2 | 0, 52) | 0;
            G() | 0;
            e2 = e2 & 15;
            a: do {
              if (!d2) {
                if (e2 | 0) {
                  d2 = 1;
                  while (1) {
                    h = dd(7, 0, (15 - d2 | 0) * 3 | 0) | 0;
                    if (!((h & a2 | 0) == 0 & ((G() | 0) & b2 | 0) == 0)) {
                      break a;
                    }
                    if (d2 >>> 0 < e2 >>> 0) {
                      d2 = d2 + 1 | 0;
                    } else {
                      break;
                    }
                  }
                }
                ab(f, e2, 0, 5, c2);
                S = g2;
                return;
              }
            } while (0);
            eb(f, e2, 0, 6, c2);
            S = g2;
            return;
          }
          function Xb(a2, b2) {
            a2 = a2 | 0;
            b2 = b2 | 0;
            var c2 = 0, d2 = 0, e2 = 0;
            d2 = cd(a2 | 0, b2 | 0, 45) | 0;
            G() | 0;
            if (!(la(d2 & 127) | 0)) {
              d2 = 2;
              return d2 | 0;
            }
            d2 = cd(a2 | 0, b2 | 0, 52) | 0;
            G() | 0;
            d2 = d2 & 15;
            if (!d2) {
              d2 = 5;
              return d2 | 0;
            }
            c2 = 1;
            while (1) {
              e2 = dd(7, 0, (15 - c2 | 0) * 3 | 0) | 0;
              if (!((e2 & a2 | 0) == 0 & ((G() | 0) & b2 | 0) == 0)) {
                c2 = 2;
                a2 = 6;
                break;
              }
              if (c2 >>> 0 < d2 >>> 0) {
                c2 = c2 + 1 | 0;
              } else {
                c2 = 5;
                a2 = 6;
                break;
              }
            }
            if ((a2 | 0) == 6) {
              return c2 | 0;
            }
            return 0;
          }
          function Yb(a2, c2, d2) {
            a2 = a2 | 0;
            c2 = c2 | 0;
            d2 = d2 | 0;
            var e2 = 0, f = 0, g2 = 0, h = 0, i = 0, j = 0, k = 0, l = 0, m = 0;
            m = S;
            S = S + 128 | 0;
            k = m + 112 | 0;
            g2 = m + 96 | 0;
            l = m;
            f = cd(a2 | 0, c2 | 0, 52) | 0;
            G() | 0;
            i = f & 15;
            b[k >> 2] = i;
            h = cd(a2 | 0, c2 | 0, 45) | 0;
            G() | 0;
            h = h & 127;
            a: do {
              if (la(h) | 0) {
                if (i | 0) {
                  e2 = 1;
                  while (1) {
                    j = dd(7, 0, (15 - e2 | 0) * 3 | 0) | 0;
                    if (!((j & a2 | 0) == 0 & ((G() | 0) & c2 | 0) == 0)) {
                      f = 0;
                      break a;
                    }
                    if (e2 >>> 0 < i >>> 0) {
                      e2 = e2 + 1 | 0;
                    } else {
                      break;
                    }
                  }
                }
                if (!(f & 1)) {
                  j = dd(i + 1 | 0, 0, 52) | 0;
                  l = G() | 0 | c2 & -15728641;
                  k = dd(7, 0, (14 - i | 0) * 3 | 0) | 0;
                  Yb((j | a2) & ~k, l & ~(G() | 0), d2);
                  S = m;
                  return;
                } else {
                  f = 1;
                }
              } else {
                f = 0;
              }
            } while (0);
            Ub(a2, c2, g2);
            if (f) {
              bb(g2, k, l);
              j = 5;
            } else {
              fb(g2, k, l);
              j = 6;
            }
            b: do {
              if (la(h) | 0) {
                if (!i) {
                  e2 = 20;
                } else {
                  e2 = 1;
                  while (1) {
                    h = dd(7, 0, (15 - e2 | 0) * 3 | 0) | 0;
                    if (!((h & a2 | 0) == 0 & ((G() | 0) & c2 | 0) == 0)) {
                      e2 = 8;
                      break b;
                    }
                    if (e2 >>> 0 < i >>> 0) {
                      e2 = e2 + 1 | 0;
                    } else {
                      e2 = 20;
                      break;
                    }
                  }
                }
              } else {
                e2 = 8;
              }
            } while (0);
            hd(d2 | 0, -1, e2 | 0) | 0;
            if (f) {
              f = 0;
              do {
                g2 = l + (f << 4) | 0;
                db(g2, b[k >> 2] | 0) | 0;
                g2 = b[g2 >> 2] | 0;
                e2 = 0;
                while (1) {
                  h = d2 + (e2 << 2) | 0;
                  i = b[h >> 2] | 0;
                  if ((i | 0) == -1 | (i | 0) == (g2 | 0)) {
                    break;
                  } else {
                    e2 = e2 + 1 | 0;
                  }
                }
                b[h >> 2] = g2;
                f = f + 1 | 0;
              } while ((f | 0) != (j | 0));
            } else {
              f = 0;
              do {
                g2 = l + (f << 4) | 0;
                cb(g2, b[k >> 2] | 0, 0, 1) | 0;
                g2 = b[g2 >> 2] | 0;
                e2 = 0;
                while (1) {
                  h = d2 + (e2 << 2) | 0;
                  i = b[h >> 2] | 0;
                  if ((i | 0) == -1 | (i | 0) == (g2 | 0)) {
                    break;
                  } else {
                    e2 = e2 + 1 | 0;
                  }
                }
                b[h >> 2] = g2;
                f = f + 1 | 0;
              } while ((f | 0) != (j | 0));
            }
            S = m;
            return;
          }
          function Zb() {
            return 12;
          }
          function _b(a2, c2) {
            a2 = a2 | 0;
            c2 = c2 | 0;
            var d2 = 0, e2 = 0, f = 0, g2 = 0, h = 0, i = 0, j = 0;
            dd(a2 | 0, 0, 52) | 0;
            i = G() | 0 | 134225919;
            if ((a2 | 0) < 1) {
              e2 = 0;
              d2 = 0;
              do {
                if (la(e2) | 0) {
                  dd(e2 | 0, 0, 45) | 0;
                  h = i | (G() | 0);
                  a2 = c2 + (d2 << 3) | 0;
                  b[a2 >> 2] = -1;
                  b[a2 + 4 >> 2] = h;
                  d2 = d2 + 1 | 0;
                }
                e2 = e2 + 1 | 0;
              } while ((e2 | 0) != 122);
              return;
            }
            h = 0;
            d2 = 0;
            do {
              if (la(h) | 0) {
                dd(h | 0, 0, 45) | 0;
                e2 = 1;
                f = -1;
                g2 = i | (G() | 0);
                while (1) {
                  j = dd(7, 0, (15 - e2 | 0) * 3 | 0) | 0;
                  f = f & ~j;
                  g2 = g2 & ~(G() | 0);
                  if ((e2 | 0) == (a2 | 0)) {
                    break;
                  } else {
                    e2 = e2 + 1 | 0;
                  }
                }
                j = c2 + (d2 << 3) | 0;
                b[j >> 2] = f;
                b[j + 4 >> 2] = g2;
                d2 = d2 + 1 | 0;
              }
              h = h + 1 | 0;
            } while ((h | 0) != 122);
            return;
          }
          function $b(a2, c2, d2, e2) {
            a2 = a2 | 0;
            c2 = c2 | 0;
            d2 = d2 | 0;
            e2 = e2 | 0;
            var f = 0, g2 = 0, h = 0, i = 0, j = 0, k = 0;
            i = S;
            S = S + 64 | 0;
            h = i;
            if ((a2 | 0) == (d2 | 0) & (c2 | 0) == (e2 | 0) | (false | (c2 & 2013265920 | 0) != 134217728 | (false | (e2 & 2013265920 | 0) != 134217728))) {
              h = 0;
              S = i;
              return h | 0;
            }
            f = cd(a2 | 0, c2 | 0, 52) | 0;
            G() | 0;
            f = f & 15;
            g2 = cd(d2 | 0, e2 | 0, 52) | 0;
            G() | 0;
            if ((f | 0) != (g2 & 15 | 0)) {
              h = 0;
              S = i;
              return h | 0;
            }
            g2 = f + -1 | 0;
            if (f >>> 0 > 1 ? (k = Cb(a2, c2, g2) | 0, j = G() | 0, g2 = Cb(d2, e2, g2) | 0, (k | 0) == (g2 | 0) & (j | 0) == (G() | 0)) : 0) {
              g2 = (f ^ 15) * 3 | 0;
              f = cd(a2 | 0, c2 | 0, g2 | 0) | 0;
              G() | 0;
              f = f & 7;
              g2 = cd(d2 | 0, e2 | 0, g2 | 0) | 0;
              G() | 0;
              g2 = g2 & 7;
              if ((f | 0) == 0 | (g2 | 0) == 0) {
                k = 1;
                S = i;
                return k | 0;
              }
              if ((b[21136 + (f << 2) >> 2] | 0) == (g2 | 0)) {
                k = 1;
                S = i;
                return k | 0;
              }
              if ((b[21168 + (f << 2) >> 2] | 0) == (g2 | 0)) {
                k = 1;
                S = i;
                return k | 0;
              }
            }
            f = h;
            g2 = f + 56 | 0;
            do {
              b[f >> 2] = 0;
              f = f + 4 | 0;
            } while ((f | 0) < (g2 | 0));
            $(a2, c2, 1, h);
            k = h;
            if (((((!((b[k >> 2] | 0) == (d2 | 0) ? (b[k + 4 >> 2] | 0) == (e2 | 0) : 0) ? (k = h + 8 | 0, !((b[k >> 2] | 0) == (d2 | 0) ? (b[k + 4 >> 2] | 0) == (e2 | 0) : 0)) : 0) ? (k = h + 16 | 0, !((b[k >> 2] | 0) == (d2 | 0) ? (b[k + 4 >> 2] | 0) == (e2 | 0) : 0)) : 0) ? (k = h + 24 | 0, !((b[k >> 2] | 0) == (d2 | 0) ? (b[k + 4 >> 2] | 0) == (e2 | 0) : 0)) : 0) ? (k = h + 32 | 0, !((b[k >> 2] | 0) == (d2 | 0) ? (b[k + 4 >> 2] | 0) == (e2 | 0) : 0)) : 0) ? (k = h + 40 | 0, !((b[k >> 2] | 0) == (d2 | 0) ? (b[k + 4 >> 2] | 0) == (e2 | 0) : 0)) : 0) {
              f = h + 48 | 0;
              f = ((b[f >> 2] | 0) == (d2 | 0) ? (b[f + 4 >> 2] | 0) == (e2 | 0) : 0) & 1;
            } else {
              f = 1;
            }
            k = f;
            S = i;
            return k | 0;
          }
          function ac(a2, c2, d2, e2) {
            a2 = a2 | 0;
            c2 = c2 | 0;
            d2 = d2 | 0;
            e2 = e2 | 0;
            var f = 0, g2 = 0, h = 0, i = 0, j = 0, k = 0, l = 0;
            k = S;
            S = S + 16 | 0;
            h = k;
            if (!($b(a2, c2, d2, e2) | 0)) {
              i = 0;
              j = 0;
              F(i | 0);
              S = k;
              return j | 0;
            }
            i = c2 & -2130706433;
            f = (Fb(a2, c2) | 0) == 0;
            f = f ? 1 : 2;
            while (1) {
              b[h >> 2] = 0;
              l = da(a2, c2, f, h) | 0;
              g2 = f + 1 | 0;
              if ((l | 0) == (d2 | 0) & (G() | 0) == (e2 | 0)) {
                break;
              }
              if (g2 >>> 0 < 7) {
                f = g2;
              } else {
                f = 0;
                a2 = 0;
                j = 6;
                break;
              }
            }
            if ((j | 0) == 6) {
              F(f | 0);
              S = k;
              return a2 | 0;
            }
            l = dd(f | 0, 0, 56) | 0;
            j = i | (G() | 0) | 268435456;
            l = a2 | l;
            F(j | 0);
            S = k;
            return l | 0;
          }
          function bc(a2, b2) {
            a2 = a2 | 0;
            b2 = b2 | 0;
            var c2 = 0;
            c2 = true & (b2 & 2013265920 | 0) == 268435456;
            F((c2 ? b2 & -2130706433 | 134217728 : 0) | 0);
            return (c2 ? a2 : 0) | 0;
          }
          function cc(a2, c2) {
            a2 = a2 | 0;
            c2 = c2 | 0;
            var d2 = 0, e2 = 0, f = 0;
            e2 = S;
            S = S + 16 | 0;
            d2 = e2;
            if (!(true & (c2 & 2013265920 | 0) == 268435456)) {
              c2 = 0;
              d2 = 0;
              F(c2 | 0);
              S = e2;
              return d2 | 0;
            }
            f = cd(a2 | 0, c2 | 0, 56) | 0;
            G() | 0;
            b[d2 >> 2] = 0;
            d2 = da(a2, c2 & -2130706433 | 134217728, f & 7, d2) | 0;
            c2 = G() | 0;
            F(c2 | 0);
            S = e2;
            return d2 | 0;
          }
          function dc(a2, b2) {
            a2 = a2 | 0;
            b2 = b2 | 0;
            var c2 = 0;
            if (!(true & (b2 & 2013265920 | 0) == 268435456)) {
              c2 = 0;
              return c2 | 0;
            }
            c2 = cd(a2 | 0, b2 | 0, 56) | 0;
            G() | 0;
            switch (c2 & 7) {
              case 0:
              case 7: {
                c2 = 0;
                return c2 | 0;
              }
              default:
            }
            c2 = b2 & -2130706433 | 134217728;
            if (true & (b2 & 117440512 | 0) == 16777216 & (Fb(a2, c2) | 0) != 0) {
              c2 = 0;
              return c2 | 0;
            }
            c2 = Bb(a2, c2) | 0;
            return c2 | 0;
          }
          function ec(a2, c2, d2) {
            a2 = a2 | 0;
            c2 = c2 | 0;
            d2 = d2 | 0;
            var e2 = 0, f = 0, g2 = 0, h = 0, i = 0;
            g2 = S;
            S = S + 16 | 0;
            e2 = g2;
            h = true & (c2 & 2013265920 | 0) == 268435456;
            f = c2 & -2130706433 | 134217728;
            i = d2;
            b[i >> 2] = h ? a2 : 0;
            b[i + 4 >> 2] = h ? f : 0;
            if (h) {
              c2 = cd(a2 | 0, c2 | 0, 56) | 0;
              G() | 0;
              b[e2 >> 2] = 0;
              a2 = da(a2, f, c2 & 7, e2) | 0;
              c2 = G() | 0;
            } else {
              a2 = 0;
              c2 = 0;
            }
            i = d2 + 8 | 0;
            b[i >> 2] = a2;
            b[i + 4 >> 2] = c2;
            S = g2;
            return;
          }
          function fc(a2, c2, d2) {
            a2 = a2 | 0;
            c2 = c2 | 0;
            d2 = d2 | 0;
            var e2 = 0, f = 0;
            f = (Fb(a2, c2) | 0) == 0;
            c2 = c2 & -2130706433;
            e2 = d2;
            b[e2 >> 2] = f ? a2 : 0;
            b[e2 + 4 >> 2] = f ? c2 | 285212672 : 0;
            e2 = d2 + 8 | 0;
            b[e2 >> 2] = a2;
            b[e2 + 4 >> 2] = c2 | 301989888;
            e2 = d2 + 16 | 0;
            b[e2 >> 2] = a2;
            b[e2 + 4 >> 2] = c2 | 318767104;
            e2 = d2 + 24 | 0;
            b[e2 >> 2] = a2;
            b[e2 + 4 >> 2] = c2 | 335544320;
            e2 = d2 + 32 | 0;
            b[e2 >> 2] = a2;
            b[e2 + 4 >> 2] = c2 | 352321536;
            d2 = d2 + 40 | 0;
            b[d2 >> 2] = a2;
            b[d2 + 4 >> 2] = c2 | 369098752;
            return;
          }
          function gc(a2, c2, d2) {
            a2 = a2 | 0;
            c2 = c2 | 0;
            d2 = d2 | 0;
            var e2 = 0, f = 0, g2 = 0, h = 0, i = 0;
            h = S;
            S = S + 16 | 0;
            g2 = h;
            e2 = cd(a2 | 0, c2 | 0, 56) | 0;
            G() | 0;
            i = true & (c2 & 2013265920 | 0) == 268435456;
            f = i ? a2 : 0;
            a2 = i ? c2 & -2130706433 | 134217728 : 0;
            c2 = Lc(f, a2, e2 & 7) | 0;
            if ((c2 | 0) == -1) {
              b[d2 >> 2] = 0;
              S = h;
              return;
            }
            Ub(f, a2, g2);
            e2 = cd(f | 0, a2 | 0, 52) | 0;
            G() | 0;
            e2 = e2 & 15;
            if (!(Fb(f, a2) | 0)) {
              eb(g2, e2, c2, 2, d2);
            } else {
              ab(g2, e2, c2, 2, d2);
            }
            S = h;
            return;
          }
          function hc(a2) {
            a2 = a2 | 0;
            var c2 = 0, d2 = 0, e2 = 0;
            c2 = Yc(1, 12) | 0;
            if (!c2) {
              H(22691, 22646, 49, 22704);
            }
            d2 = a2 + 4 | 0;
            e2 = b[d2 >> 2] | 0;
            if (e2 | 0) {
              e2 = e2 + 8 | 0;
              b[e2 >> 2] = c2;
              b[d2 >> 2] = c2;
              return c2 | 0;
            }
            if (b[a2 >> 2] | 0) {
              H(22721, 22646, 61, 22744);
            }
            e2 = a2;
            b[e2 >> 2] = c2;
            b[d2 >> 2] = c2;
            return c2 | 0;
          }
          function ic(a2, c2) {
            a2 = a2 | 0;
            c2 = c2 | 0;
            var d2 = 0, e2 = 0;
            e2 = Wc(24) | 0;
            if (!e2) {
              H(22758, 22646, 78, 22772);
            }
            b[e2 >> 2] = b[c2 >> 2];
            b[e2 + 4 >> 2] = b[c2 + 4 >> 2];
            b[e2 + 8 >> 2] = b[c2 + 8 >> 2];
            b[e2 + 12 >> 2] = b[c2 + 12 >> 2];
            b[e2 + 16 >> 2] = 0;
            c2 = a2 + 4 | 0;
            d2 = b[c2 >> 2] | 0;
            if (d2 | 0) {
              b[d2 + 16 >> 2] = e2;
              b[c2 >> 2] = e2;
              return e2 | 0;
            }
            if (b[a2 >> 2] | 0) {
              H(22787, 22646, 82, 22772);
            }
            b[a2 >> 2] = e2;
            b[c2 >> 2] = e2;
            return e2 | 0;
          }
          function jc(a2) {
            a2 = a2 | 0;
            var c2 = 0, d2 = 0, e2 = 0, f = 0;
            if (!a2) {
              return;
            }
            e2 = 1;
            while (1) {
              c2 = b[a2 >> 2] | 0;
              if (c2 | 0) {
                do {
                  d2 = b[c2 >> 2] | 0;
                  if (d2 | 0) {
                    do {
                      f = d2;
                      d2 = b[d2 + 16 >> 2] | 0;
                      Xc(f);
                    } while ((d2 | 0) != 0);
                  }
                  f = c2;
                  c2 = b[c2 + 8 >> 2] | 0;
                  Xc(f);
                } while ((c2 | 0) != 0);
              }
              c2 = a2;
              a2 = b[a2 + 8 >> 2] | 0;
              if (!e2) {
                Xc(c2);
              }
              if (!a2) {
                break;
              } else {
                e2 = 0;
              }
            }
            return;
          }
          function kc(a2) {
            a2 = a2 | 0;
            var c2 = 0, d2 = 0, f = 0, g2 = 0, h = 0, i = 0, j = 0, k = 0, l = 0, m = 0, n = 0, o = 0, p2 = 0, r2 = 0, s2 = 0, t2 = 0, u2 = 0, v2 = 0, w2 = 0, x2 = 0, y2 = 0, z2 = 0, A2 = 0, B2 = 0, C2 = 0, D2 = 0, E = 0, F2 = 0, G2 = 0, I2 = 0, J2 = 0, K2 = 0;
            g2 = a2 + 8 | 0;
            if (b[g2 >> 2] | 0) {
              K2 = 1;
              return K2 | 0;
            }
            f = b[a2 >> 2] | 0;
            if (!f) {
              K2 = 0;
              return K2 | 0;
            }
            c2 = f;
            d2 = 0;
            do {
              d2 = d2 + 1 | 0;
              c2 = b[c2 + 8 >> 2] | 0;
            } while ((c2 | 0) != 0);
            if (d2 >>> 0 < 2) {
              K2 = 0;
              return K2 | 0;
            }
            I2 = Wc(d2 << 2) | 0;
            if (!I2) {
              H(22807, 22646, 317, 22826);
            }
            G2 = Wc(d2 << 5) | 0;
            if (!G2) {
              H(22848, 22646, 321, 22826);
            }
            b[a2 >> 2] = 0;
            z2 = a2 + 4 | 0;
            b[z2 >> 2] = 0;
            b[g2 >> 2] = 0;
            d2 = 0;
            F2 = 0;
            y2 = 0;
            n = 0;
            a: while (1) {
              m = b[f >> 2] | 0;
              if (m) {
                h = 0;
                i = m;
                do {
                  k = +e[i + 8 >> 3];
                  c2 = i;
                  i = b[i + 16 >> 2] | 0;
                  l = (i | 0) == 0;
                  g2 = l ? m : i;
                  j = +e[g2 + 8 >> 3];
                  if (+q(+(k - j)) > 3.141592653589793) {
                    K2 = 14;
                    break;
                  }
                  h = h + (j - k) * (+e[c2 >> 3] + +e[g2 >> 3]);
                } while (!l);
                if ((K2 | 0) == 14) {
                  K2 = 0;
                  h = 0;
                  c2 = m;
                  do {
                    x2 = +e[c2 + 8 >> 3];
                    E = c2 + 16 | 0;
                    D2 = b[E >> 2] | 0;
                    D2 = (D2 | 0) == 0 ? m : D2;
                    w2 = +e[D2 + 8 >> 3];
                    h = h + (+e[c2 >> 3] + +e[D2 >> 3]) * ((w2 < 0 ? w2 + 6.283185307179586 : w2) - (x2 < 0 ? x2 + 6.283185307179586 : x2));
                    c2 = b[((c2 | 0) == 0 ? f : E) >> 2] | 0;
                  } while ((c2 | 0) != 0);
                }
                if (h > 0) {
                  b[I2 + (F2 << 2) >> 2] = f;
                  F2 = F2 + 1 | 0;
                  g2 = y2;
                  c2 = n;
                } else {
                  K2 = 19;
                }
              } else {
                K2 = 19;
              }
              if ((K2 | 0) == 19) {
                K2 = 0;
                do {
                  if (!d2) {
                    if (!n) {
                      if (!(b[a2 >> 2] | 0)) {
                        g2 = z2;
                        i = a2;
                        c2 = f;
                        d2 = a2;
                        break;
                      } else {
                        K2 = 27;
                        break a;
                      }
                    } else {
                      g2 = z2;
                      i = n + 8 | 0;
                      c2 = f;
                      d2 = a2;
                      break;
                    }
                  } else {
                    c2 = d2 + 8 | 0;
                    if (b[c2 >> 2] | 0) {
                      K2 = 21;
                      break a;
                    }
                    d2 = Yc(1, 12) | 0;
                    if (!d2) {
                      K2 = 23;
                      break a;
                    }
                    b[c2 >> 2] = d2;
                    g2 = d2 + 4 | 0;
                    i = d2;
                    c2 = n;
                  }
                } while (0);
                b[i >> 2] = f;
                b[g2 >> 2] = f;
                i = G2 + (y2 << 5) | 0;
                l = b[f >> 2] | 0;
                if (l) {
                  m = G2 + (y2 << 5) + 8 | 0;
                  e[m >> 3] = 17976931348623157e292;
                  n = G2 + (y2 << 5) + 24 | 0;
                  e[n >> 3] = 17976931348623157e292;
                  e[i >> 3] = -17976931348623157e292;
                  o = G2 + (y2 << 5) + 16 | 0;
                  e[o >> 3] = -17976931348623157e292;
                  u2 = 17976931348623157e292;
                  v2 = -17976931348623157e292;
                  g2 = 0;
                  p2 = l;
                  k = 17976931348623157e292;
                  s2 = 17976931348623157e292;
                  t2 = -17976931348623157e292;
                  j = -17976931348623157e292;
                  while (1) {
                    h = +e[p2 >> 3];
                    x2 = +e[p2 + 8 >> 3];
                    p2 = b[p2 + 16 >> 2] | 0;
                    r2 = (p2 | 0) == 0;
                    w2 = +e[(r2 ? l : p2) + 8 >> 3];
                    if (h < k) {
                      e[m >> 3] = h;
                      k = h;
                    }
                    if (x2 < s2) {
                      e[n >> 3] = x2;
                      s2 = x2;
                    }
                    if (h > t2) {
                      e[i >> 3] = h;
                    } else {
                      h = t2;
                    }
                    if (x2 > j) {
                      e[o >> 3] = x2;
                      j = x2;
                    }
                    u2 = x2 > 0 & x2 < u2 ? x2 : u2;
                    v2 = x2 < 0 & x2 > v2 ? x2 : v2;
                    g2 = g2 | +q(+(x2 - w2)) > 3.141592653589793;
                    if (r2) {
                      break;
                    } else {
                      t2 = h;
                    }
                  }
                  if (g2) {
                    e[o >> 3] = v2;
                    e[n >> 3] = u2;
                  }
                } else {
                  b[i >> 2] = 0;
                  b[i + 4 >> 2] = 0;
                  b[i + 8 >> 2] = 0;
                  b[i + 12 >> 2] = 0;
                  b[i + 16 >> 2] = 0;
                  b[i + 20 >> 2] = 0;
                  b[i + 24 >> 2] = 0;
                  b[i + 28 >> 2] = 0;
                }
                g2 = y2 + 1 | 0;
              }
              E = f + 8 | 0;
              f = b[E >> 2] | 0;
              b[E >> 2] = 0;
              if (!f) {
                K2 = 45;
                break;
              } else {
                y2 = g2;
                n = c2;
              }
            }
            if ((K2 | 0) == 21) {
              H(22624, 22646, 35, 22658);
            } else if ((K2 | 0) == 23) {
              H(22678, 22646, 37, 22658);
            } else if ((K2 | 0) == 27) {
              H(22721, 22646, 61, 22744);
            } else if ((K2 | 0) == 45) {
              b: do {
                if ((F2 | 0) > 0) {
                  E = (g2 | 0) == 0;
                  C2 = g2 << 2;
                  D2 = (a2 | 0) == 0;
                  B2 = 0;
                  c2 = 0;
                  while (1) {
                    A2 = b[I2 + (B2 << 2) >> 2] | 0;
                    if (!E) {
                      y2 = Wc(C2) | 0;
                      if (!y2) {
                        K2 = 50;
                        break;
                      }
                      z2 = Wc(C2) | 0;
                      if (!z2) {
                        K2 = 52;
                        break;
                      }
                      c: do {
                        if (!D2) {
                          g2 = 0;
                          d2 = 0;
                          i = a2;
                          while (1) {
                            f = G2 + (g2 << 5) | 0;
                            if (lc(b[i >> 2] | 0, f, b[A2 >> 2] | 0) | 0) {
                              b[y2 + (d2 << 2) >> 2] = i;
                              b[z2 + (d2 << 2) >> 2] = f;
                              r2 = d2 + 1 | 0;
                            } else {
                              r2 = d2;
                            }
                            i = b[i + 8 >> 2] | 0;
                            if (!i) {
                              break;
                            } else {
                              g2 = g2 + 1 | 0;
                              d2 = r2;
                            }
                          }
                          if ((r2 | 0) > 0) {
                            f = b[y2 >> 2] | 0;
                            if ((r2 | 0) == 1) {
                              d2 = f;
                            } else {
                              o = 0;
                              p2 = -1;
                              d2 = f;
                              n = f;
                              while (1) {
                                l = b[n >> 2] | 0;
                                f = 0;
                                i = 0;
                                while (1) {
                                  g2 = b[b[y2 + (i << 2) >> 2] >> 2] | 0;
                                  if ((g2 | 0) == (l | 0)) {
                                    m = f;
                                  } else {
                                    m = f + ((lc(g2, b[z2 + (i << 2) >> 2] | 0, b[l >> 2] | 0) | 0) & 1) | 0;
                                  }
                                  i = i + 1 | 0;
                                  if ((i | 0) == (r2 | 0)) {
                                    break;
                                  } else {
                                    f = m;
                                  }
                                }
                                g2 = (m | 0) > (p2 | 0);
                                d2 = g2 ? n : d2;
                                f = o + 1 | 0;
                                if ((f | 0) == (r2 | 0)) {
                                  break c;
                                }
                                o = f;
                                p2 = g2 ? m : p2;
                                n = b[y2 + (f << 2) >> 2] | 0;
                              }
                            }
                          } else {
                            d2 = 0;
                          }
                        } else {
                          d2 = 0;
                        }
                      } while (0);
                      Xc(y2);
                      Xc(z2);
                      if (d2) {
                        g2 = d2 + 4 | 0;
                        f = b[g2 >> 2] | 0;
                        if (!f) {
                          if (b[d2 >> 2] | 0) {
                            K2 = 70;
                            break;
                          }
                        } else {
                          d2 = f + 8 | 0;
                        }
                        b[d2 >> 2] = A2;
                        b[g2 >> 2] = A2;
                      } else {
                        K2 = 73;
                      }
                    } else {
                      K2 = 73;
                    }
                    if ((K2 | 0) == 73) {
                      K2 = 0;
                      c2 = b[A2 >> 2] | 0;
                      if (c2 | 0) {
                        do {
                          z2 = c2;
                          c2 = b[c2 + 16 >> 2] | 0;
                          Xc(z2);
                        } while ((c2 | 0) != 0);
                      }
                      Xc(A2);
                      c2 = 2;
                    }
                    B2 = B2 + 1 | 0;
                    if ((B2 | 0) >= (F2 | 0)) {
                      J2 = c2;
                      break b;
                    }
                  }
                  if ((K2 | 0) == 50) {
                    H(22863, 22646, 249, 22882);
                  } else if ((K2 | 0) == 52) {
                    H(22901, 22646, 252, 22882);
                  } else if ((K2 | 0) == 70) {
                    H(22721, 22646, 61, 22744);
                  }
                } else {
                  J2 = 0;
                }
              } while (0);
              Xc(I2);
              Xc(G2);
              K2 = J2;
              return K2 | 0;
            }
            return 0;
          }
          function lc(a2, c2, d2) {
            a2 = a2 | 0;
            c2 = c2 | 0;
            d2 = d2 | 0;
            var f = 0, g2 = 0, h = 0, i = 0, j = 0, k = 0, l = 0, m = 0;
            if (!(xa(c2, d2) | 0)) {
              a2 = 0;
              return a2 | 0;
            }
            c2 = wa(c2) | 0;
            l = +e[d2 >> 3];
            f = +e[d2 + 8 >> 3];
            f = c2 & f < 0 ? f + 6.283185307179586 : f;
            a2 = b[a2 >> 2] | 0;
            if (!a2) {
              a2 = 0;
              return a2 | 0;
            }
            if (c2) {
              c2 = 0;
              d2 = a2;
              a: while (1) {
                while (1) {
                  i = +e[d2 >> 3];
                  k = +e[d2 + 8 >> 3];
                  d2 = d2 + 16 | 0;
                  m = b[d2 >> 2] | 0;
                  m = (m | 0) == 0 ? a2 : m;
                  h = +e[m >> 3];
                  g2 = +e[m + 8 >> 3];
                  if (i > h) {
                    j = i;
                    i = k;
                  } else {
                    j = h;
                    h = i;
                    i = g2;
                    g2 = k;
                  }
                  if (!(l < h | l > j)) {
                    break;
                  }
                  d2 = b[d2 >> 2] | 0;
                  if (!d2) {
                    d2 = 22;
                    break a;
                  }
                }
                k = g2 < 0 ? g2 + 6.283185307179586 : g2;
                i = i < 0 ? i + 6.283185307179586 : i;
                f = i == f | k == f ? f + -2220446049250313e-31 : f;
                k = k + (l - h) / (j - h) * (i - k);
                if ((k < 0 ? k + 6.283185307179586 : k) > f) {
                  c2 = c2 ^ 1;
                }
                d2 = b[d2 >> 2] | 0;
                if (!d2) {
                  d2 = 22;
                  break;
                }
              }
              if ((d2 | 0) == 22) {
                return c2 | 0;
              }
            } else {
              c2 = 0;
              d2 = a2;
              b: while (1) {
                while (1) {
                  i = +e[d2 >> 3];
                  k = +e[d2 + 8 >> 3];
                  d2 = d2 + 16 | 0;
                  m = b[d2 >> 2] | 0;
                  m = (m | 0) == 0 ? a2 : m;
                  h = +e[m >> 3];
                  g2 = +e[m + 8 >> 3];
                  if (i > h) {
                    j = i;
                    i = k;
                  } else {
                    j = h;
                    h = i;
                    i = g2;
                    g2 = k;
                  }
                  if (!(l < h | l > j)) {
                    break;
                  }
                  d2 = b[d2 >> 2] | 0;
                  if (!d2) {
                    d2 = 22;
                    break b;
                  }
                }
                f = i == f | g2 == f ? f + -2220446049250313e-31 : f;
                if (g2 + (l - h) / (j - h) * (i - g2) > f) {
                  c2 = c2 ^ 1;
                }
                d2 = b[d2 >> 2] | 0;
                if (!d2) {
                  d2 = 22;
                  break;
                }
              }
              if ((d2 | 0) == 22) {
                return c2 | 0;
              }
            }
            return 0;
          }
          function mc(c2, d2, e2, f, g2) {
            c2 = c2 | 0;
            d2 = d2 | 0;
            e2 = e2 | 0;
            f = f | 0;
            g2 = g2 | 0;
            var h = 0, i = 0, j = 0, k = 0, l = 0, m = 0, n = 0, o = 0, p2 = 0, q2 = 0, r2 = 0, s2 = 0, t2 = 0, u2 = 0;
            u2 = S;
            S = S + 32 | 0;
            t2 = u2 + 16 | 0;
            s2 = u2;
            h = cd(c2 | 0, d2 | 0, 52) | 0;
            G() | 0;
            h = h & 15;
            p2 = cd(e2 | 0, f | 0, 52) | 0;
            G() | 0;
            if ((h | 0) != (p2 & 15 | 0)) {
              t2 = 1;
              S = u2;
              return t2 | 0;
            }
            l = cd(c2 | 0, d2 | 0, 45) | 0;
            G() | 0;
            l = l & 127;
            m = cd(e2 | 0, f | 0, 45) | 0;
            G() | 0;
            m = m & 127;
            p2 = (l | 0) != (m | 0);
            if (p2) {
              j = ta(l, m) | 0;
              if ((j | 0) == 7) {
                t2 = 2;
                S = u2;
                return t2 | 0;
              }
              k = ta(m, l) | 0;
              if ((k | 0) == 7) {
                H(22925, 22949, 151, 22959);
              } else {
                q2 = j;
                i = k;
              }
            } else {
              q2 = 0;
              i = 0;
            }
            n = la(l) | 0;
            o = la(m) | 0;
            b[t2 >> 2] = 0;
            b[t2 + 4 >> 2] = 0;
            b[t2 + 8 >> 2] = 0;
            b[t2 + 12 >> 2] = 0;
            do {
              if (!q2) {
                Tb(e2, f, t2) | 0;
                if ((n | 0) != 0 & (o | 0) != 0) {
                  if ((m | 0) != (l | 0)) {
                    H(23077, 22949, 243, 22959);
                  }
                  i = Lb(c2, d2) | 0;
                  h = Lb(e2, f) | 0;
                  if (!(a[22032 + (i * 7 | 0) + h >> 0] | 0)) {
                    i = b[21200 + (i * 28 | 0) + (h << 2) >> 2] | 0;
                    if ((i | 0) > 0) {
                      j = t2 + 4 | 0;
                      h = 0;
                      do {
                        Oa(j);
                        h = h + 1 | 0;
                      } while ((h | 0) != (i | 0));
                      r2 = 50;
                    } else {
                      r2 = 50;
                    }
                  } else {
                    h = 5;
                  }
                } else {
                  r2 = 50;
                }
              } else {
                m = b[4304 + (l * 28 | 0) + (q2 << 2) >> 2] | 0;
                j = (m | 0) > 0;
                if (!o) {
                  if (j) {
                    l = 0;
                    k = e2;
                    j = f;
                    do {
                      k = Pb(k, j) | 0;
                      j = G() | 0;
                      i = Qa(i) | 0;
                      l = l + 1 | 0;
                    } while ((l | 0) != (m | 0));
                    m = i;
                    l = k;
                    k = j;
                  } else {
                    m = i;
                    l = e2;
                    k = f;
                  }
                } else if (j) {
                  l = 0;
                  k = e2;
                  j = f;
                  do {
                    k = Ob(k, j) | 0;
                    j = G() | 0;
                    i = Qa(i) | 0;
                    if ((i | 0) == 1) {
                      i = Qa(1) | 0;
                    }
                    l = l + 1 | 0;
                  } while ((l | 0) != (m | 0));
                  m = i;
                  l = k;
                  k = j;
                } else {
                  m = i;
                  l = e2;
                  k = f;
                }
                Tb(l, k, t2) | 0;
                if (!p2) {
                  H(22972, 22949, 181, 22959);
                }
                j = (n | 0) != 0;
                i = (o | 0) != 0;
                if (j & i) {
                  H(22999, 22949, 182, 22959);
                }
                if (!j) {
                  if (i) {
                    i = Lb(l, k) | 0;
                    if (a[22032 + (i * 7 | 0) + m >> 0] | 0) {
                      h = 4;
                      break;
                    }
                    l = 0;
                    k = b[21200 + (m * 28 | 0) + (i << 2) >> 2] | 0;
                    r2 = 26;
                  } else {
                    i = 0;
                  }
                } else {
                  i = Lb(c2, d2) | 0;
                  if (a[22032 + (i * 7 | 0) + q2 >> 0] | 0) {
                    h = 3;
                    break;
                  }
                  k = b[21200 + (i * 28 | 0) + (q2 << 2) >> 2] | 0;
                  l = k;
                  r2 = 26;
                }
                if ((r2 | 0) == 26) {
                  if ((k | 0) <= -1) {
                    H(23030, 22949, 212, 22959);
                  }
                  if ((l | 0) <= -1) {
                    H(23053, 22949, 213, 22959);
                  }
                  if ((k | 0) > 0) {
                    j = t2 + 4 | 0;
                    i = 0;
                    do {
                      Oa(j);
                      i = i + 1 | 0;
                    } while ((i | 0) != (k | 0));
                    i = l;
                  } else {
                    i = l;
                  }
                }
                b[s2 >> 2] = 0;
                b[s2 + 4 >> 2] = 0;
                b[s2 + 8 >> 2] = 0;
                Ma(s2, q2);
                if (h | 0) {
                  while (1) {
                    if (!(Rb(h) | 0)) {
                      La(s2);
                    } else {
                      Ka(s2);
                    }
                    if ((h | 0) > 1) {
                      h = h + -1 | 0;
                    } else {
                      break;
                    }
                  }
                }
                if ((i | 0) > 0) {
                  h = 0;
                  do {
                    Oa(s2);
                    h = h + 1 | 0;
                  } while ((h | 0) != (i | 0));
                }
                r2 = t2 + 4 | 0;
                Ea(r2, s2, r2);
                Ca(r2);
                r2 = 50;
              }
            } while (0);
            if ((r2 | 0) == 50) {
              h = t2 + 4 | 0;
              b[g2 >> 2] = b[h >> 2];
              b[g2 + 4 >> 2] = b[h + 4 >> 2];
              b[g2 + 8 >> 2] = b[h + 8 >> 2];
              h = 0;
            }
            t2 = h;
            S = u2;
            return t2 | 0;
          }
          function nc(a2, c2, d2, e2) {
            a2 = a2 | 0;
            c2 = c2 | 0;
            d2 = d2 | 0;
            e2 = e2 | 0;
            var f = 0, g2 = 0, h = 0, i = 0, j = 0, k = 0, l = 0, m = 0, n = 0, o = 0, p2 = 0, q2 = 0, r2 = 0, s2 = 0, t2 = 0;
            p2 = S;
            S = S + 48 | 0;
            h = p2 + 36 | 0;
            i = p2 + 24 | 0;
            j = p2 + 12 | 0;
            k = p2;
            g2 = cd(a2 | 0, c2 | 0, 52) | 0;
            G() | 0;
            g2 = g2 & 15;
            n = cd(a2 | 0, c2 | 0, 45) | 0;
            G() | 0;
            n = n & 127;
            l = la(n) | 0;
            dd(g2 | 0, 0, 52) | 0;
            r2 = G() | 0 | 134225919;
            q2 = e2;
            b[q2 >> 2] = -1;
            b[q2 + 4 >> 2] = r2;
            if (!g2) {
              if ((b[d2 >> 2] | 0) > 1) {
                r2 = 1;
                S = p2;
                return r2 | 0;
              }
              if ((b[d2 + 4 >> 2] | 0) > 1) {
                r2 = 1;
                S = p2;
                return r2 | 0;
              }
              if ((b[d2 + 8 >> 2] | 0) > 1) {
                r2 = 1;
                S = p2;
                return r2 | 0;
              }
              f = sa(n, Ha(d2) | 0) | 0;
              if ((f | 0) == 127) {
                r2 = 1;
                S = p2;
                return r2 | 0;
              }
              o = dd(f | 0, 0, 45) | 0;
              q2 = G() | 0;
              n = e2;
              q2 = b[n + 4 >> 2] & -1040385 | q2;
              r2 = e2;
              b[r2 >> 2] = b[n >> 2] | o;
              b[r2 + 4 >> 2] = q2;
              r2 = 0;
              S = p2;
              return r2 | 0;
            }
            b[h >> 2] = b[d2 >> 2];
            b[h + 4 >> 2] = b[d2 + 4 >> 2];
            b[h + 8 >> 2] = b[d2 + 8 >> 2];
            while (1) {
              b[i >> 2] = b[h >> 2];
              b[i + 4 >> 2] = b[h + 4 >> 2];
              b[i + 8 >> 2] = b[h + 8 >> 2];
              if (!(Rb(g2) | 0)) {
                Ja(h);
                b[j >> 2] = b[h >> 2];
                b[j + 4 >> 2] = b[h + 4 >> 2];
                b[j + 8 >> 2] = b[h + 8 >> 2];
                La(j);
              } else {
                Ia(h);
                b[j >> 2] = b[h >> 2];
                b[j + 4 >> 2] = b[h + 4 >> 2];
                b[j + 8 >> 2] = b[h + 8 >> 2];
                Ka(j);
              }
              Fa(i, j, k);
              Ca(k);
              q2 = e2;
              s2 = b[q2 >> 2] | 0;
              q2 = b[q2 + 4 >> 2] | 0;
              t2 = (15 - g2 | 0) * 3 | 0;
              d2 = dd(7, 0, t2 | 0) | 0;
              q2 = q2 & ~(G() | 0);
              t2 = dd(Ha(k) | 0, 0, t2 | 0) | 0;
              q2 = G() | 0 | q2;
              r2 = e2;
              b[r2 >> 2] = t2 | s2 & ~d2;
              b[r2 + 4 >> 2] = q2;
              if ((g2 | 0) > 1) {
                g2 = g2 + -1 | 0;
              } else {
                break;
              }
            }
            a: do {
              if (((b[h >> 2] | 0) <= 1 ? (b[h + 4 >> 2] | 0) <= 1 : 0) ? (b[h + 8 >> 2] | 0) <= 1 : 0) {
                g2 = Ha(h) | 0;
                i = sa(n, g2) | 0;
                if ((i | 0) == 127) {
                  k = 0;
                } else {
                  k = la(i) | 0;
                }
                b: do {
                  if (!g2) {
                    if ((l | 0) != 0 & (k | 0) != 0) {
                      t2 = Lb(a2, c2) | 0;
                      g2 = e2;
                      g2 = 21408 + (t2 * 28 | 0) + ((Lb(b[g2 >> 2] | 0, b[g2 + 4 >> 2] | 0) | 0) << 2) | 0;
                      g2 = b[g2 >> 2] | 0;
                      if ((g2 | 0) <= -1) {
                        H(23201, 22949, 433, 23134);
                      }
                      if (!g2) {
                        f = i;
                        g2 = 55;
                      } else {
                        h = e2;
                        f = 0;
                        d2 = b[h >> 2] | 0;
                        h = b[h + 4 >> 2] | 0;
                        do {
                          d2 = Nb(d2, h) | 0;
                          h = G() | 0;
                          t2 = e2;
                          b[t2 >> 2] = d2;
                          b[t2 + 4 >> 2] = h;
                          f = f + 1 | 0;
                        } while ((f | 0) < (g2 | 0));
                        f = i;
                        g2 = 54;
                      }
                    } else {
                      f = i;
                      g2 = 54;
                    }
                  } else {
                    if (l) {
                      h = 21408 + ((Lb(a2, c2) | 0) * 28 | 0) + (g2 << 2) | 0;
                      h = b[h >> 2] | 0;
                      if ((h | 0) > 0) {
                        d2 = 0;
                        do {
                          g2 = Pa(g2) | 0;
                          d2 = d2 + 1 | 0;
                        } while ((d2 | 0) != (h | 0));
                      }
                      if ((g2 | 0) == 1) {
                        f = 3;
                        break a;
                      }
                      d2 = sa(n, g2) | 0;
                      if ((d2 | 0) == 127) {
                        H(23104, 22949, 376, 23134);
                      }
                      if (!(la(d2) | 0)) {
                        o = h;
                        m = g2;
                        f = d2;
                      } else {
                        H(23147, 22949, 377, 23134);
                      }
                    } else {
                      o = 0;
                      m = g2;
                      f = i;
                    }
                    j = b[4304 + (n * 28 | 0) + (m << 2) >> 2] | 0;
                    if ((j | 0) <= -1) {
                      H(23178, 22949, 384, 23134);
                    }
                    if (!k) {
                      if ((o | 0) <= -1) {
                        H(23030, 22949, 417, 23134);
                      }
                      if (o | 0) {
                        h = e2;
                        g2 = 0;
                        d2 = b[h >> 2] | 0;
                        h = b[h + 4 >> 2] | 0;
                        do {
                          d2 = Nb(d2, h) | 0;
                          h = G() | 0;
                          t2 = e2;
                          b[t2 >> 2] = d2;
                          b[t2 + 4 >> 2] = h;
                          g2 = g2 + 1 | 0;
                        } while ((g2 | 0) < (o | 0));
                      }
                      if ((j | 0) <= 0) {
                        g2 = 54;
                        break;
                      }
                      h = e2;
                      g2 = 0;
                      d2 = b[h >> 2] | 0;
                      h = b[h + 4 >> 2] | 0;
                      while (1) {
                        d2 = Nb(d2, h) | 0;
                        h = G() | 0;
                        t2 = e2;
                        b[t2 >> 2] = d2;
                        b[t2 + 4 >> 2] = h;
                        g2 = g2 + 1 | 0;
                        if ((g2 | 0) == (j | 0)) {
                          g2 = 54;
                          break b;
                        }
                      }
                    }
                    i = ta(f, n) | 0;
                    if ((i | 0) == 7) {
                      H(22925, 22949, 393, 23134);
                    }
                    g2 = e2;
                    d2 = b[g2 >> 2] | 0;
                    g2 = b[g2 + 4 >> 2] | 0;
                    if ((j | 0) > 0) {
                      h = 0;
                      do {
                        d2 = Nb(d2, g2) | 0;
                        g2 = G() | 0;
                        t2 = e2;
                        b[t2 >> 2] = d2;
                        b[t2 + 4 >> 2] = g2;
                        h = h + 1 | 0;
                      } while ((h | 0) != (j | 0));
                    }
                    d2 = Lb(d2, g2) | 0;
                    t2 = ma(f) | 0;
                    d2 = b[(t2 ? 21824 : 21616) + (i * 28 | 0) + (d2 << 2) >> 2] | 0;
                    if ((d2 | 0) <= -1) {
                      H(23030, 22949, 412, 23134);
                    }
                    if (!d2) {
                      g2 = 54;
                    } else {
                      i = e2;
                      g2 = 0;
                      h = b[i >> 2] | 0;
                      i = b[i + 4 >> 2] | 0;
                      do {
                        h = Mb(h, i) | 0;
                        i = G() | 0;
                        t2 = e2;
                        b[t2 >> 2] = h;
                        b[t2 + 4 >> 2] = i;
                        g2 = g2 + 1 | 0;
                      } while ((g2 | 0) < (d2 | 0));
                      g2 = 54;
                    }
                  }
                } while (0);
                if ((g2 | 0) == 54) {
                  if (k) {
                    g2 = 55;
                  }
                }
                if ((g2 | 0) == 55) {
                  t2 = e2;
                  if ((Lb(b[t2 >> 2] | 0, b[t2 + 4 >> 2] | 0) | 0) == 1) {
                    f = 4;
                    break;
                  }
                }
                t2 = e2;
                r2 = b[t2 >> 2] | 0;
                t2 = b[t2 + 4 >> 2] & -1040385;
                s2 = dd(f | 0, 0, 45) | 0;
                t2 = t2 | (G() | 0);
                f = e2;
                b[f >> 2] = r2 | s2;
                b[f + 4 >> 2] = t2;
                f = 0;
              } else {
                f = 2;
              }
            } while (0);
            t2 = f;
            S = p2;
            return t2 | 0;
          }
          function oc(a2, b2, c2, d2, e2) {
            a2 = a2 | 0;
            b2 = b2 | 0;
            c2 = c2 | 0;
            d2 = d2 | 0;
            e2 = e2 | 0;
            var f = 0, g2 = 0;
            g2 = S;
            S = S + 16 | 0;
            f = g2;
            a2 = mc(a2, b2, c2, d2, f) | 0;
            if (!a2) {
              Ua(f, e2);
              a2 = 0;
            }
            S = g2;
            return a2 | 0;
          }
          function pc(a2, b2, c2, d2) {
            a2 = a2 | 0;
            b2 = b2 | 0;
            c2 = c2 | 0;
            d2 = d2 | 0;
            var e2 = 0, f = 0;
            e2 = S;
            S = S + 16 | 0;
            f = e2;
            Va(c2, f);
            d2 = nc(a2, b2, f, d2) | 0;
            S = e2;
            return d2 | 0;
          }
          function qc(a2, b2, c2, d2) {
            a2 = a2 | 0;
            b2 = b2 | 0;
            c2 = c2 | 0;
            d2 = d2 | 0;
            var e2 = 0, f = 0, g2 = 0;
            g2 = S;
            S = S + 32 | 0;
            e2 = g2 + 12 | 0;
            f = g2;
            if ((mc(a2, b2, a2, b2, e2) | 0) == 0 ? (mc(a2, b2, c2, d2, f) | 0) == 0 : 0) {
              a2 = Ta(e2, f) | 0;
            } else {
              a2 = -1;
            }
            S = g2;
            return a2 | 0;
          }
          function rc(a2, b2, c2, d2) {
            a2 = a2 | 0;
            b2 = b2 | 0;
            c2 = c2 | 0;
            d2 = d2 | 0;
            var e2 = 0, f = 0, g2 = 0;
            g2 = S;
            S = S + 32 | 0;
            e2 = g2 + 12 | 0;
            f = g2;
            if ((mc(a2, b2, a2, b2, e2) | 0) == 0 ? (mc(a2, b2, c2, d2, f) | 0) == 0 : 0) {
              a2 = Ta(e2, f) | 0;
            } else {
              a2 = -1;
            }
            S = g2;
            return (a2 >>> 31 ^ 1) + a2 | 0;
          }
          function sc(a2, c2, d2, e2, f) {
            a2 = a2 | 0;
            c2 = c2 | 0;
            d2 = d2 | 0;
            e2 = e2 | 0;
            f = f | 0;
            var g2 = 0, h = 0, i = 0, j = 0, k = 0, l = 0, m = 0, n = 0, o = 0, p2 = 0, r2 = 0, s2 = 0, t2 = 0, u2 = 0, v2 = 0, w2 = 0, x2 = 0, y2 = 0;
            x2 = S;
            S = S + 48 | 0;
            g2 = x2 + 24 | 0;
            h = x2 + 12 | 0;
            w2 = x2;
            if ((mc(a2, c2, a2, c2, g2) | 0) == 0 ? (mc(a2, c2, d2, e2, h) | 0) == 0 : 0) {
              v2 = Ta(g2, h) | 0;
              if ((v2 | 0) < 0) {
                w2 = v2;
                S = x2;
                return w2 | 0;
              }
              b[g2 >> 2] = 0;
              b[g2 + 4 >> 2] = 0;
              b[g2 + 8 >> 2] = 0;
              b[h >> 2] = 0;
              b[h + 4 >> 2] = 0;
              b[h + 8 >> 2] = 0;
              mc(a2, c2, a2, c2, g2) | 0;
              mc(a2, c2, d2, e2, h) | 0;
              Wa(g2);
              Wa(h);
              if (!v2) {
                e2 = g2 + 4 | 0;
                n = g2 + 8 | 0;
                s2 = e2;
                t2 = n;
                u2 = g2;
                d2 = b[g2 >> 2] | 0;
                e2 = b[e2 >> 2] | 0;
                g2 = b[n >> 2] | 0;
                p2 = 0;
                r2 = 0;
                o = 0;
              } else {
                l = b[g2 >> 2] | 0;
                o = +(v2 | 0);
                s2 = g2 + 4 | 0;
                m = b[s2 >> 2] | 0;
                t2 = g2 + 8 | 0;
                n = b[t2 >> 2] | 0;
                u2 = g2;
                d2 = l;
                e2 = m;
                g2 = n;
                p2 = +((b[h >> 2] | 0) - l | 0) / o;
                r2 = +((b[h + 4 >> 2] | 0) - m | 0) / o;
                o = +((b[h + 8 >> 2] | 0) - n | 0) / o;
              }
              b[w2 >> 2] = d2;
              n = w2 + 4 | 0;
              b[n >> 2] = e2;
              m = w2 + 8 | 0;
              b[m >> 2] = g2;
              l = 0;
              while (1) {
                j = +(l | 0);
                y2 = p2 * j + +(d2 | 0);
                i = r2 * j + +(b[s2 >> 2] | 0);
                j = o * j + +(b[t2 >> 2] | 0);
                e2 = ~~+fd(+y2);
                h = ~~+fd(+i);
                d2 = ~~+fd(+j);
                y2 = +q(+(+(e2 | 0) - y2));
                i = +q(+(+(h | 0) - i));
                j = +q(+(+(d2 | 0) - j));
                do {
                  if (!(y2 > i & y2 > j)) {
                    k = 0 - e2 | 0;
                    if (i > j) {
                      g2 = k - d2 | 0;
                      break;
                    } else {
                      g2 = h;
                      d2 = k - h | 0;
                      break;
                    }
                  } else {
                    e2 = 0 - (h + d2) | 0;
                    g2 = h;
                  }
                } while (0);
                b[w2 >> 2] = e2;
                b[n >> 2] = g2;
                b[m >> 2] = d2;
                Xa(w2);
                nc(a2, c2, w2, f + (l << 3) | 0) | 0;
                if ((l | 0) == (v2 | 0)) {
                  break;
                }
                l = l + 1 | 0;
                d2 = b[u2 >> 2] | 0;
              }
              w2 = 0;
              S = x2;
              return w2 | 0;
            }
            w2 = -1;
            S = x2;
            return w2 | 0;
          }
          function tc(a2, b2) {
            a2 = a2 | 0;
            b2 = b2 | 0;
            var c2 = 0;
            if (!b2) {
              c2 = 1;
              return c2 | 0;
            }
            c2 = a2;
            a2 = 1;
            do {
              a2 = B((b2 & 1 | 0) == 0 ? 1 : c2, a2) | 0;
              b2 = b2 >> 1;
              c2 = B(c2, c2) | 0;
            } while ((b2 | 0) != 0);
            return a2 | 0;
          }
          function uc(a2, c2, d2) {
            a2 = a2 | 0;
            c2 = c2 | 0;
            d2 = d2 | 0;
            var f = 0, g2 = 0, h = 0, i = 0, j = 0, k = 0, l = 0, m = 0, n = 0, o = 0;
            if (!(xa(c2, d2) | 0)) {
              n = 0;
              return n | 0;
            }
            c2 = wa(c2) | 0;
            o = +e[d2 >> 3];
            f = +e[d2 + 8 >> 3];
            f = c2 & f < 0 ? f + 6.283185307179586 : f;
            n = b[a2 >> 2] | 0;
            if ((n | 0) <= 0) {
              n = 0;
              return n | 0;
            }
            m = b[a2 + 4 >> 2] | 0;
            if (c2) {
              c2 = 0;
              d2 = -1;
              a2 = 0;
              a: while (1) {
                l = a2;
                while (1) {
                  i = +e[m + (l << 4) >> 3];
                  k = +e[m + (l << 4) + 8 >> 3];
                  a2 = (d2 + 2 | 0) % (n | 0) | 0;
                  h = +e[m + (a2 << 4) >> 3];
                  g2 = +e[m + (a2 << 4) + 8 >> 3];
                  if (i > h) {
                    j = i;
                    i = k;
                  } else {
                    j = h;
                    h = i;
                    i = g2;
                    g2 = k;
                  }
                  if (!(o < h | o > j)) {
                    break;
                  }
                  d2 = l + 1 | 0;
                  if ((d2 | 0) < (n | 0)) {
                    a2 = l;
                    l = d2;
                    d2 = a2;
                  } else {
                    d2 = 22;
                    break a;
                  }
                }
                k = g2 < 0 ? g2 + 6.283185307179586 : g2;
                i = i < 0 ? i + 6.283185307179586 : i;
                f = i == f | k == f ? f + -2220446049250313e-31 : f;
                k = k + (o - h) / (j - h) * (i - k);
                if ((k < 0 ? k + 6.283185307179586 : k) > f) {
                  c2 = c2 ^ 1;
                }
                a2 = l + 1 | 0;
                if ((a2 | 0) >= (n | 0)) {
                  d2 = 22;
                  break;
                } else {
                  d2 = l;
                }
              }
              if ((d2 | 0) == 22) {
                return c2 | 0;
              }
            } else {
              c2 = 0;
              d2 = -1;
              a2 = 0;
              b: while (1) {
                l = a2;
                while (1) {
                  i = +e[m + (l << 4) >> 3];
                  k = +e[m + (l << 4) + 8 >> 3];
                  a2 = (d2 + 2 | 0) % (n | 0) | 0;
                  h = +e[m + (a2 << 4) >> 3];
                  g2 = +e[m + (a2 << 4) + 8 >> 3];
                  if (i > h) {
                    j = i;
                    i = k;
                  } else {
                    j = h;
                    h = i;
                    i = g2;
                    g2 = k;
                  }
                  if (!(o < h | o > j)) {
                    break;
                  }
                  d2 = l + 1 | 0;
                  if ((d2 | 0) < (n | 0)) {
                    a2 = l;
                    l = d2;
                    d2 = a2;
                  } else {
                    d2 = 22;
                    break b;
                  }
                }
                f = i == f | g2 == f ? f + -2220446049250313e-31 : f;
                if (g2 + (o - h) / (j - h) * (i - g2) > f) {
                  c2 = c2 ^ 1;
                }
                a2 = l + 1 | 0;
                if ((a2 | 0) >= (n | 0)) {
                  d2 = 22;
                  break;
                } else {
                  d2 = l;
                }
              }
              if ((d2 | 0) == 22) {
                return c2 | 0;
              }
            }
            return 0;
          }
          function vc(a2, c2) {
            a2 = a2 | 0;
            c2 = c2 | 0;
            var d2 = 0, f = 0, g2 = 0, h = 0, i = 0, j = 0, k = 0, l = 0, m = 0, n = 0, o = 0, p2 = 0, r2 = 0, s2 = 0, t2 = 0, u2 = 0, v2 = 0;
            r2 = b[a2 >> 2] | 0;
            if (!r2) {
              b[c2 >> 2] = 0;
              b[c2 + 4 >> 2] = 0;
              b[c2 + 8 >> 2] = 0;
              b[c2 + 12 >> 2] = 0;
              b[c2 + 16 >> 2] = 0;
              b[c2 + 20 >> 2] = 0;
              b[c2 + 24 >> 2] = 0;
              b[c2 + 28 >> 2] = 0;
              return;
            }
            s2 = c2 + 8 | 0;
            e[s2 >> 3] = 17976931348623157e292;
            t2 = c2 + 24 | 0;
            e[t2 >> 3] = 17976931348623157e292;
            e[c2 >> 3] = -17976931348623157e292;
            u2 = c2 + 16 | 0;
            e[u2 >> 3] = -17976931348623157e292;
            if ((r2 | 0) <= 0) {
              return;
            }
            o = b[a2 + 4 >> 2] | 0;
            l = 17976931348623157e292;
            m = -17976931348623157e292;
            n = 0;
            a2 = -1;
            h = 17976931348623157e292;
            i = 17976931348623157e292;
            k = -17976931348623157e292;
            f = -17976931348623157e292;
            p2 = 0;
            while (1) {
              d2 = +e[o + (p2 << 4) >> 3];
              j = +e[o + (p2 << 4) + 8 >> 3];
              a2 = a2 + 2 | 0;
              g2 = +e[o + (((a2 | 0) == (r2 | 0) ? 0 : a2) << 4) + 8 >> 3];
              if (d2 < h) {
                e[s2 >> 3] = d2;
                h = d2;
              }
              if (j < i) {
                e[t2 >> 3] = j;
                i = j;
              }
              if (d2 > k) {
                e[c2 >> 3] = d2;
              } else {
                d2 = k;
              }
              if (j > f) {
                e[u2 >> 3] = j;
                f = j;
              }
              l = j > 0 & j < l ? j : l;
              m = j < 0 & j > m ? j : m;
              n = n | +q(+(j - g2)) > 3.141592653589793;
              a2 = p2 + 1 | 0;
              if ((a2 | 0) == (r2 | 0)) {
                break;
              } else {
                v2 = p2;
                k = d2;
                p2 = a2;
                a2 = v2;
              }
            }
            if (!n) {
              return;
            }
            e[u2 >> 3] = m;
            e[t2 >> 3] = l;
            return;
          }
          function wc(a2, c2) {
            a2 = a2 | 0;
            c2 = c2 | 0;
            var d2 = 0, f = 0, g2 = 0, h = 0, i = 0, j = 0, k = 0, l = 0, m = 0, n = 0, o = 0, p2 = 0, r2 = 0, s2 = 0, t2 = 0, u2 = 0, v2 = 0, w2 = 0, x2 = 0, y2 = 0, z2 = 0, A2 = 0;
            r2 = b[a2 >> 2] | 0;
            if (r2) {
              s2 = c2 + 8 | 0;
              e[s2 >> 3] = 17976931348623157e292;
              t2 = c2 + 24 | 0;
              e[t2 >> 3] = 17976931348623157e292;
              e[c2 >> 3] = -17976931348623157e292;
              u2 = c2 + 16 | 0;
              e[u2 >> 3] = -17976931348623157e292;
              if ((r2 | 0) > 0) {
                g2 = b[a2 + 4 >> 2] | 0;
                o = 17976931348623157e292;
                p2 = -17976931348623157e292;
                f = 0;
                d2 = -1;
                k = 17976931348623157e292;
                l = 17976931348623157e292;
                n = -17976931348623157e292;
                i = -17976931348623157e292;
                v2 = 0;
                while (1) {
                  h = +e[g2 + (v2 << 4) >> 3];
                  m = +e[g2 + (v2 << 4) + 8 >> 3];
                  z2 = d2 + 2 | 0;
                  j = +e[g2 + (((z2 | 0) == (r2 | 0) ? 0 : z2) << 4) + 8 >> 3];
                  if (h < k) {
                    e[s2 >> 3] = h;
                    k = h;
                  }
                  if (m < l) {
                    e[t2 >> 3] = m;
                    l = m;
                  }
                  if (h > n) {
                    e[c2 >> 3] = h;
                  } else {
                    h = n;
                  }
                  if (m > i) {
                    e[u2 >> 3] = m;
                    i = m;
                  }
                  o = m > 0 & m < o ? m : o;
                  p2 = m < 0 & m > p2 ? m : p2;
                  f = f | +q(+(m - j)) > 3.141592653589793;
                  d2 = v2 + 1 | 0;
                  if ((d2 | 0) == (r2 | 0)) {
                    break;
                  } else {
                    z2 = v2;
                    n = h;
                    v2 = d2;
                    d2 = z2;
                  }
                }
                if (f) {
                  e[u2 >> 3] = p2;
                  e[t2 >> 3] = o;
                }
              }
            } else {
              b[c2 >> 2] = 0;
              b[c2 + 4 >> 2] = 0;
              b[c2 + 8 >> 2] = 0;
              b[c2 + 12 >> 2] = 0;
              b[c2 + 16 >> 2] = 0;
              b[c2 + 20 >> 2] = 0;
              b[c2 + 24 >> 2] = 0;
              b[c2 + 28 >> 2] = 0;
            }
            z2 = a2 + 8 | 0;
            d2 = b[z2 >> 2] | 0;
            if ((d2 | 0) <= 0) {
              return;
            }
            y2 = a2 + 12 | 0;
            x2 = 0;
            do {
              g2 = b[y2 >> 2] | 0;
              f = x2;
              x2 = x2 + 1 | 0;
              t2 = c2 + (x2 << 5) | 0;
              u2 = b[g2 + (f << 3) >> 2] | 0;
              if (u2) {
                v2 = c2 + (x2 << 5) + 8 | 0;
                e[v2 >> 3] = 17976931348623157e292;
                a2 = c2 + (x2 << 5) + 24 | 0;
                e[a2 >> 3] = 17976931348623157e292;
                e[t2 >> 3] = -17976931348623157e292;
                w2 = c2 + (x2 << 5) + 16 | 0;
                e[w2 >> 3] = -17976931348623157e292;
                if ((u2 | 0) > 0) {
                  r2 = b[g2 + (f << 3) + 4 >> 2] | 0;
                  o = 17976931348623157e292;
                  p2 = -17976931348623157e292;
                  g2 = 0;
                  f = -1;
                  s2 = 0;
                  k = 17976931348623157e292;
                  l = 17976931348623157e292;
                  m = -17976931348623157e292;
                  i = -17976931348623157e292;
                  while (1) {
                    h = +e[r2 + (s2 << 4) >> 3];
                    n = +e[r2 + (s2 << 4) + 8 >> 3];
                    f = f + 2 | 0;
                    j = +e[r2 + (((f | 0) == (u2 | 0) ? 0 : f) << 4) + 8 >> 3];
                    if (h < k) {
                      e[v2 >> 3] = h;
                      k = h;
                    }
                    if (n < l) {
                      e[a2 >> 3] = n;
                      l = n;
                    }
                    if (h > m) {
                      e[t2 >> 3] = h;
                    } else {
                      h = m;
                    }
                    if (n > i) {
                      e[w2 >> 3] = n;
                      i = n;
                    }
                    o = n > 0 & n < o ? n : o;
                    p2 = n < 0 & n > p2 ? n : p2;
                    g2 = g2 | +q(+(n - j)) > 3.141592653589793;
                    f = s2 + 1 | 0;
                    if ((f | 0) == (u2 | 0)) {
                      break;
                    } else {
                      A2 = s2;
                      s2 = f;
                      m = h;
                      f = A2;
                    }
                  }
                  if (g2) {
                    e[w2 >> 3] = p2;
                    e[a2 >> 3] = o;
                  }
                }
              } else {
                b[t2 >> 2] = 0;
                b[t2 + 4 >> 2] = 0;
                b[t2 + 8 >> 2] = 0;
                b[t2 + 12 >> 2] = 0;
                b[t2 + 16 >> 2] = 0;
                b[t2 + 20 >> 2] = 0;
                b[t2 + 24 >> 2] = 0;
                b[t2 + 28 >> 2] = 0;
                d2 = b[z2 >> 2] | 0;
              }
            } while ((x2 | 0) < (d2 | 0));
            return;
          }
          function xc(a2, c2, d2) {
            a2 = a2 | 0;
            c2 = c2 | 0;
            d2 = d2 | 0;
            var e2 = 0, f = 0, g2 = 0;
            if (!(uc(a2, c2, d2) | 0)) {
              f = 0;
              return f | 0;
            }
            f = a2 + 8 | 0;
            if ((b[f >> 2] | 0) <= 0) {
              f = 1;
              return f | 0;
            }
            e2 = a2 + 12 | 0;
            a2 = 0;
            while (1) {
              g2 = a2;
              a2 = a2 + 1 | 0;
              if (uc((b[e2 >> 2] | 0) + (g2 << 3) | 0, c2 + (a2 << 5) | 0, d2) | 0) {
                a2 = 0;
                e2 = 6;
                break;
              }
              if ((a2 | 0) >= (b[f >> 2] | 0)) {
                a2 = 1;
                e2 = 6;
                break;
              }
            }
            if ((e2 | 0) == 6) {
              return a2 | 0;
            }
            return 0;
          }
          function yc() {
            return 8;
          }
          function zc() {
            return 16;
          }
          function Ac() {
            return 168;
          }
          function Bc() {
            return 8;
          }
          function Cc() {
            return 16;
          }
          function Dc() {
            return 12;
          }
          function Ec() {
            return 8;
          }
          function Fc(a2) {
            a2 = a2 | 0;
            var b2 = 0, c2 = 0;
            c2 = +e[a2 >> 3];
            b2 = +e[a2 + 8 >> 3];
            return + +r(+(c2 * c2 + b2 * b2));
          }
          function Gc(a2, b2, c2, d2, f) {
            a2 = a2 | 0;
            b2 = b2 | 0;
            c2 = c2 | 0;
            d2 = d2 | 0;
            f = f | 0;
            var g2 = 0, h = 0, i = 0, j = 0, k = 0, l = 0, m = 0, n = 0;
            k = +e[a2 >> 3];
            j = +e[b2 >> 3] - k;
            i = +e[a2 + 8 >> 3];
            h = +e[b2 + 8 >> 3] - i;
            m = +e[c2 >> 3];
            g2 = +e[d2 >> 3] - m;
            n = +e[c2 + 8 >> 3];
            l = +e[d2 + 8 >> 3] - n;
            g2 = (g2 * (i - n) - (k - m) * l) / (j * l - h * g2);
            e[f >> 3] = k + j * g2;
            e[f + 8 >> 3] = i + h * g2;
            return;
          }
          function Hc(a2, b2) {
            a2 = a2 | 0;
            b2 = b2 | 0;
            if (!(+e[a2 >> 3] == +e[b2 >> 3])) {
              b2 = 0;
              return b2 | 0;
            }
            b2 = +e[a2 + 8 >> 3] == +e[b2 + 8 >> 3];
            return b2 | 0;
          }
          function Ic(a2, b2) {
            a2 = a2 | 0;
            b2 = b2 | 0;
            var c2 = 0, d2 = 0, f = 0;
            f = +e[a2 >> 3] - +e[b2 >> 3];
            d2 = +e[a2 + 8 >> 3] - +e[b2 + 8 >> 3];
            c2 = +e[a2 + 16 >> 3] - +e[b2 + 16 >> 3];
            return +(f * f + d2 * d2 + c2 * c2);
          }
          function Jc(a2, b2) {
            a2 = a2 | 0;
            b2 = b2 | 0;
            var c2 = 0, d2 = 0, f = 0;
            c2 = +e[a2 >> 3];
            d2 = +t(+c2);
            c2 = +u(+c2);
            e[b2 + 16 >> 3] = c2;
            c2 = +e[a2 + 8 >> 3];
            f = d2 * +t(+c2);
            e[b2 >> 3] = f;
            c2 = d2 * +u(+c2);
            e[b2 + 8 >> 3] = c2;
            return;
          }
          function Kc(a2, c2) {
            a2 = a2 | 0;
            c2 = c2 | 0;
            var d2 = 0, e2 = 0, f = 0, g2 = 0, h = 0, i = 0, j = 0, k = 0;
            k = S;
            S = S + 32 | 0;
            f = k + 16 | 0;
            g2 = k;
            Ub(a2, c2, f);
            h = Ab(a2, c2) | 0;
            j = Lb(a2, c2) | 0;
            pa(h, g2);
            c2 = qa(h, b[f >> 2] | 0) | 0;
            if (!(la(h) | 0)) {
              j = c2;
              S = k;
              return j | 0;
            }
            do {
              switch (h | 0) {
                case 4: {
                  a2 = 0;
                  d2 = 14;
                  break;
                }
                case 14: {
                  a2 = 1;
                  d2 = 14;
                  break;
                }
                case 24: {
                  a2 = 2;
                  d2 = 14;
                  break;
                }
                case 38: {
                  a2 = 3;
                  d2 = 14;
                  break;
                }
                case 49: {
                  a2 = 4;
                  d2 = 14;
                  break;
                }
                case 58: {
                  a2 = 5;
                  d2 = 14;
                  break;
                }
                case 63: {
                  a2 = 6;
                  d2 = 14;
                  break;
                }
                case 72: {
                  a2 = 7;
                  d2 = 14;
                  break;
                }
                case 83: {
                  a2 = 8;
                  d2 = 14;
                  break;
                }
                case 97: {
                  a2 = 9;
                  d2 = 14;
                  break;
                }
                case 107: {
                  a2 = 10;
                  d2 = 14;
                  break;
                }
                case 117: {
                  a2 = 11;
                  d2 = 14;
                  break;
                }
                default: {
                  i = 0;
                  e2 = 0;
                }
              }
            } while (0);
            if ((d2 | 0) == 14) {
              i = b[22096 + (a2 * 24 | 0) + 8 >> 2] | 0;
              e2 = b[22096 + (a2 * 24 | 0) + 16 >> 2] | 0;
            }
            a2 = b[f >> 2] | 0;
            if ((a2 | 0) != (b[g2 >> 2] | 0)) {
              h = ma(h) | 0;
              a2 = b[f >> 2] | 0;
              if (h | (a2 | 0) == (e2 | 0)) {
                c2 = (c2 + 1 | 0) % 6 | 0;
              }
            }
            if ((j | 0) == 3 & (a2 | 0) == (e2 | 0)) {
              j = (c2 + 5 | 0) % 6 | 0;
              S = k;
              return j | 0;
            }
            if (!((j | 0) == 5 & (a2 | 0) == (i | 0))) {
              j = c2;
              S = k;
              return j | 0;
            }
            j = (c2 + 1 | 0) % 6 | 0;
            S = k;
            return j | 0;
          }
          function Lc(a2, c2, d2) {
            a2 = a2 | 0;
            c2 = c2 | 0;
            d2 = d2 | 0;
            var e2 = 0, f = 0;
            e2 = Fb(a2, c2) | 0;
            if ((d2 + -1 | 0) >>> 0 > 5) {
              d2 = -1;
              return d2 | 0;
            }
            f = (e2 | 0) != 0;
            if ((d2 | 0) == 1 & f) {
              d2 = -1;
              return d2 | 0;
            }
            e2 = Kc(a2, c2) | 0;
            if (f) {
              d2 = (5 - e2 + (b[22384 + (d2 << 2) >> 2] | 0) | 0) % 5 | 0;
              return d2 | 0;
            } else {
              d2 = (6 - e2 + (b[22416 + (d2 << 2) >> 2] | 0) | 0) % 6 | 0;
              return d2 | 0;
            }
            return 0;
          }
          function Mc(a2, c2, d2) {
            a2 = a2 | 0;
            c2 = c2 | 0;
            d2 = d2 | 0;
            var e2 = 0;
            if ((c2 | 0) > 0) {
              e2 = Yc(c2, 4) | 0;
              b[a2 >> 2] = e2;
              if (!e2) {
                H(23230, 23253, 40, 23267);
              }
            } else {
              b[a2 >> 2] = 0;
            }
            b[a2 + 4 >> 2] = c2;
            b[a2 + 8 >> 2] = 0;
            b[a2 + 12 >> 2] = d2;
            return;
          }
          function Nc(a2) {
            a2 = a2 | 0;
            var c2 = 0, d2 = 0, f = 0, g2 = 0, h = 0, i = 0, j = 0;
            g2 = a2 + 4 | 0;
            h = a2 + 12 | 0;
            i = a2 + 8 | 0;
            a: while (1) {
              d2 = b[g2 >> 2] | 0;
              c2 = 0;
              while (1) {
                if ((c2 | 0) >= (d2 | 0)) {
                  break a;
                }
                f = b[a2 >> 2] | 0;
                j = b[f + (c2 << 2) >> 2] | 0;
                if (!j) {
                  c2 = c2 + 1 | 0;
                } else {
                  break;
                }
              }
              c2 = f + (~~(+q(+(+s(10, + +(15 - (b[h >> 2] | 0) | 0)) * (+e[j >> 3] + +e[j + 8 >> 3]))) % +(d2 | 0)) >>> 0 << 2) | 0;
              d2 = b[c2 >> 2] | 0;
              b: do {
                if (d2 | 0) {
                  f = j + 32 | 0;
                  if ((d2 | 0) == (j | 0)) {
                    b[c2 >> 2] = b[f >> 2];
                  } else {
                    d2 = d2 + 32 | 0;
                    c2 = b[d2 >> 2] | 0;
                    if (!c2) {
                      break;
                    }
                    while (1) {
                      if ((c2 | 0) == (j | 0)) {
                        break;
                      }
                      d2 = c2 + 32 | 0;
                      c2 = b[d2 >> 2] | 0;
                      if (!c2) {
                        break b;
                      }
                    }
                    b[d2 >> 2] = b[f >> 2];
                  }
                  Xc(j);
                  b[i >> 2] = (b[i >> 2] | 0) + -1;
                }
              } while (0);
            }
            Xc(b[a2 >> 2] | 0);
            return;
          }
          function Oc(a2) {
            a2 = a2 | 0;
            var c2 = 0, d2 = 0, e2 = 0;
            e2 = b[a2 + 4 >> 2] | 0;
            d2 = 0;
            while (1) {
              if ((d2 | 0) >= (e2 | 0)) {
                c2 = 0;
                d2 = 4;
                break;
              }
              c2 = b[(b[a2 >> 2] | 0) + (d2 << 2) >> 2] | 0;
              if (!c2) {
                d2 = d2 + 1 | 0;
              } else {
                d2 = 4;
                break;
              }
            }
            if ((d2 | 0) == 4) {
              return c2 | 0;
            }
            return 0;
          }
          function Pc(a2, c2) {
            a2 = a2 | 0;
            c2 = c2 | 0;
            var d2 = 0, f = 0, g2 = 0, h = 0;
            d2 = ~~(+q(+(+s(10, + +(15 - (b[a2 + 12 >> 2] | 0) | 0)) * (+e[c2 >> 3] + +e[c2 + 8 >> 3]))) % +(b[a2 + 4 >> 2] | 0)) >>> 0;
            d2 = (b[a2 >> 2] | 0) + (d2 << 2) | 0;
            f = b[d2 >> 2] | 0;
            if (!f) {
              h = 1;
              return h | 0;
            }
            h = c2 + 32 | 0;
            do {
              if ((f | 0) != (c2 | 0)) {
                d2 = b[f + 32 >> 2] | 0;
                if (!d2) {
                  h = 1;
                  return h | 0;
                }
                g2 = d2;
                while (1) {
                  if ((g2 | 0) == (c2 | 0)) {
                    g2 = 8;
                    break;
                  }
                  d2 = b[g2 + 32 >> 2] | 0;
                  if (!d2) {
                    d2 = 1;
                    g2 = 10;
                    break;
                  } else {
                    f = g2;
                    g2 = d2;
                  }
                }
                if ((g2 | 0) == 8) {
                  b[f + 32 >> 2] = b[h >> 2];
                  break;
                } else if ((g2 | 0) == 10) {
                  return d2 | 0;
                }
              } else {
                b[d2 >> 2] = b[h >> 2];
              }
            } while (0);
            Xc(c2);
            h = a2 + 8 | 0;
            b[h >> 2] = (b[h >> 2] | 0) + -1;
            h = 0;
            return h | 0;
          }
          function Qc(a2, c2, d2) {
            a2 = a2 | 0;
            c2 = c2 | 0;
            d2 = d2 | 0;
            var f = 0, g2 = 0, h = 0, i = 0;
            h = Wc(40) | 0;
            if (!h) {
              H(23283, 23253, 98, 23296);
            }
            b[h >> 2] = b[c2 >> 2];
            b[h + 4 >> 2] = b[c2 + 4 >> 2];
            b[h + 8 >> 2] = b[c2 + 8 >> 2];
            b[h + 12 >> 2] = b[c2 + 12 >> 2];
            g2 = h + 16 | 0;
            b[g2 >> 2] = b[d2 >> 2];
            b[g2 + 4 >> 2] = b[d2 + 4 >> 2];
            b[g2 + 8 >> 2] = b[d2 + 8 >> 2];
            b[g2 + 12 >> 2] = b[d2 + 12 >> 2];
            b[h + 32 >> 2] = 0;
            g2 = ~~(+q(+(+s(10, + +(15 - (b[a2 + 12 >> 2] | 0) | 0)) * (+e[c2 >> 3] + +e[c2 + 8 >> 3]))) % +(b[a2 + 4 >> 2] | 0)) >>> 0;
            g2 = (b[a2 >> 2] | 0) + (g2 << 2) | 0;
            f = b[g2 >> 2] | 0;
            do {
              if (!f) {
                b[g2 >> 2] = h;
              } else {
                while (1) {
                  if (hb(f, c2) | 0 ? hb(f + 16 | 0, d2) | 0 : 0) {
                    break;
                  }
                  g2 = b[f + 32 >> 2] | 0;
                  f = (g2 | 0) == 0 ? f : g2;
                  if (!(b[f + 32 >> 2] | 0)) {
                    i = 10;
                    break;
                  }
                }
                if ((i | 0) == 10) {
                  b[f + 32 >> 2] = h;
                  break;
                }
                Xc(h);
                i = f;
                return i | 0;
              }
            } while (0);
            i = a2 + 8 | 0;
            b[i >> 2] = (b[i >> 2] | 0) + 1;
            i = h;
            return i | 0;
          }
          function Rc(a2, c2, d2) {
            a2 = a2 | 0;
            c2 = c2 | 0;
            d2 = d2 | 0;
            var f = 0, g2 = 0;
            g2 = ~~(+q(+(+s(10, + +(15 - (b[a2 + 12 >> 2] | 0) | 0)) * (+e[c2 >> 3] + +e[c2 + 8 >> 3]))) % +(b[a2 + 4 >> 2] | 0)) >>> 0;
            g2 = b[(b[a2 >> 2] | 0) + (g2 << 2) >> 2] | 0;
            if (!g2) {
              d2 = 0;
              return d2 | 0;
            }
            if (!d2) {
              a2 = g2;
              while (1) {
                if (hb(a2, c2) | 0) {
                  f = 10;
                  break;
                }
                a2 = b[a2 + 32 >> 2] | 0;
                if (!a2) {
                  a2 = 0;
                  f = 10;
                  break;
                }
              }
              if ((f | 0) == 10) {
                return a2 | 0;
              }
            }
            a2 = g2;
            while (1) {
              if (hb(a2, c2) | 0 ? hb(a2 + 16 | 0, d2) | 0 : 0) {
                f = 10;
                break;
              }
              a2 = b[a2 + 32 >> 2] | 0;
              if (!a2) {
                a2 = 0;
                f = 10;
                break;
              }
            }
            if ((f | 0) == 10) {
              return a2 | 0;
            }
            return 0;
          }
          function Sc(a2, c2) {
            a2 = a2 | 0;
            c2 = c2 | 0;
            var d2 = 0;
            d2 = ~~(+q(+(+s(10, + +(15 - (b[a2 + 12 >> 2] | 0) | 0)) * (+e[c2 >> 3] + +e[c2 + 8 >> 3]))) % +(b[a2 + 4 >> 2] | 0)) >>> 0;
            a2 = b[(b[a2 >> 2] | 0) + (d2 << 2) >> 2] | 0;
            if (!a2) {
              d2 = 0;
              return d2 | 0;
            }
            while (1) {
              if (hb(a2, c2) | 0) {
                c2 = 5;
                break;
              }
              a2 = b[a2 + 32 >> 2] | 0;
              if (!a2) {
                a2 = 0;
                c2 = 5;
                break;
              }
            }
            if ((c2 | 0) == 5) {
              return a2 | 0;
            }
            return 0;
          }
          function Tc() {
            return 23312;
          }
          function Uc(a2) {
            a2 = +a2;
            return + +id(+a2);
          }
          function Vc(a2) {
            a2 = +a2;
            return ~~+Uc(a2) | 0;
          }
          function Wc(a2) {
            a2 = a2 | 0;
            var c2 = 0, d2 = 0, e2 = 0, f = 0, g2 = 0, h = 0, i = 0, j = 0, k = 0, l = 0, m = 0, n = 0, o = 0, p2 = 0, q2 = 0, r2 = 0, s2 = 0, t2 = 0, u2 = 0, v2 = 0, w2 = 0;
            w2 = S;
            S = S + 16 | 0;
            n = w2;
            do {
              if (a2 >>> 0 < 245) {
                k = a2 >>> 0 < 11 ? 16 : a2 + 11 & -8;
                a2 = k >>> 3;
                m = b[5829] | 0;
                d2 = m >>> a2;
                if (d2 & 3 | 0) {
                  c2 = (d2 & 1 ^ 1) + a2 | 0;
                  a2 = 23356 + (c2 << 1 << 2) | 0;
                  d2 = a2 + 8 | 0;
                  e2 = b[d2 >> 2] | 0;
                  f = e2 + 8 | 0;
                  g2 = b[f >> 2] | 0;
                  if ((g2 | 0) == (a2 | 0)) {
                    b[5829] = m & ~(1 << c2);
                  } else {
                    b[g2 + 12 >> 2] = a2;
                    b[d2 >> 2] = g2;
                  }
                  v2 = c2 << 3;
                  b[e2 + 4 >> 2] = v2 | 3;
                  v2 = e2 + v2 + 4 | 0;
                  b[v2 >> 2] = b[v2 >> 2] | 1;
                  v2 = f;
                  S = w2;
                  return v2 | 0;
                }
                l = b[5831] | 0;
                if (k >>> 0 > l >>> 0) {
                  if (d2 | 0) {
                    c2 = 2 << a2;
                    c2 = d2 << a2 & (c2 | 0 - c2);
                    c2 = (c2 & 0 - c2) + -1 | 0;
                    i = c2 >>> 12 & 16;
                    c2 = c2 >>> i;
                    d2 = c2 >>> 5 & 8;
                    c2 = c2 >>> d2;
                    g2 = c2 >>> 2 & 4;
                    c2 = c2 >>> g2;
                    a2 = c2 >>> 1 & 2;
                    c2 = c2 >>> a2;
                    e2 = c2 >>> 1 & 1;
                    e2 = (d2 | i | g2 | a2 | e2) + (c2 >>> e2) | 0;
                    c2 = 23356 + (e2 << 1 << 2) | 0;
                    a2 = c2 + 8 | 0;
                    g2 = b[a2 >> 2] | 0;
                    i = g2 + 8 | 0;
                    d2 = b[i >> 2] | 0;
                    if ((d2 | 0) == (c2 | 0)) {
                      a2 = m & ~(1 << e2);
                      b[5829] = a2;
                    } else {
                      b[d2 + 12 >> 2] = c2;
                      b[a2 >> 2] = d2;
                      a2 = m;
                    }
                    v2 = e2 << 3;
                    h = v2 - k | 0;
                    b[g2 + 4 >> 2] = k | 3;
                    f = g2 + k | 0;
                    b[f + 4 >> 2] = h | 1;
                    b[g2 + v2 >> 2] = h;
                    if (l | 0) {
                      e2 = b[5834] | 0;
                      c2 = l >>> 3;
                      d2 = 23356 + (c2 << 1 << 2) | 0;
                      c2 = 1 << c2;
                      if (!(a2 & c2)) {
                        b[5829] = a2 | c2;
                        c2 = d2;
                        a2 = d2 + 8 | 0;
                      } else {
                        a2 = d2 + 8 | 0;
                        c2 = b[a2 >> 2] | 0;
                      }
                      b[a2 >> 2] = e2;
                      b[c2 + 12 >> 2] = e2;
                      b[e2 + 8 >> 2] = c2;
                      b[e2 + 12 >> 2] = d2;
                    }
                    b[5831] = h;
                    b[5834] = f;
                    v2 = i;
                    S = w2;
                    return v2 | 0;
                  }
                  g2 = b[5830] | 0;
                  if (g2) {
                    d2 = (g2 & 0 - g2) + -1 | 0;
                    f = d2 >>> 12 & 16;
                    d2 = d2 >>> f;
                    e2 = d2 >>> 5 & 8;
                    d2 = d2 >>> e2;
                    h = d2 >>> 2 & 4;
                    d2 = d2 >>> h;
                    i = d2 >>> 1 & 2;
                    d2 = d2 >>> i;
                    j = d2 >>> 1 & 1;
                    j = b[23620 + ((e2 | f | h | i | j) + (d2 >>> j) << 2) >> 2] | 0;
                    d2 = j;
                    i = j;
                    j = (b[j + 4 >> 2] & -8) - k | 0;
                    while (1) {
                      a2 = b[d2 + 16 >> 2] | 0;
                      if (!a2) {
                        a2 = b[d2 + 20 >> 2] | 0;
                        if (!a2) {
                          break;
                        }
                      }
                      h = (b[a2 + 4 >> 2] & -8) - k | 0;
                      f = h >>> 0 < j >>> 0;
                      d2 = a2;
                      i = f ? a2 : i;
                      j = f ? h : j;
                    }
                    h = i + k | 0;
                    if (h >>> 0 > i >>> 0) {
                      f = b[i + 24 >> 2] | 0;
                      c2 = b[i + 12 >> 2] | 0;
                      do {
                        if ((c2 | 0) == (i | 0)) {
                          a2 = i + 20 | 0;
                          c2 = b[a2 >> 2] | 0;
                          if (!c2) {
                            a2 = i + 16 | 0;
                            c2 = b[a2 >> 2] | 0;
                            if (!c2) {
                              d2 = 0;
                              break;
                            }
                          }
                          while (1) {
                            e2 = c2 + 20 | 0;
                            d2 = b[e2 >> 2] | 0;
                            if (!d2) {
                              e2 = c2 + 16 | 0;
                              d2 = b[e2 >> 2] | 0;
                              if (!d2) {
                                break;
                              } else {
                                c2 = d2;
                                a2 = e2;
                              }
                            } else {
                              c2 = d2;
                              a2 = e2;
                            }
                          }
                          b[a2 >> 2] = 0;
                          d2 = c2;
                        } else {
                          d2 = b[i + 8 >> 2] | 0;
                          b[d2 + 12 >> 2] = c2;
                          b[c2 + 8 >> 2] = d2;
                          d2 = c2;
                        }
                      } while (0);
                      do {
                        if (f | 0) {
                          c2 = b[i + 28 >> 2] | 0;
                          a2 = 23620 + (c2 << 2) | 0;
                          if ((i | 0) == (b[a2 >> 2] | 0)) {
                            b[a2 >> 2] = d2;
                            if (!d2) {
                              b[5830] = g2 & ~(1 << c2);
                              break;
                            }
                          } else {
                            v2 = f + 16 | 0;
                            b[((b[v2 >> 2] | 0) == (i | 0) ? v2 : f + 20 | 0) >> 2] = d2;
                            if (!d2) {
                              break;
                            }
                          }
                          b[d2 + 24 >> 2] = f;
                          c2 = b[i + 16 >> 2] | 0;
                          if (c2 | 0) {
                            b[d2 + 16 >> 2] = c2;
                            b[c2 + 24 >> 2] = d2;
                          }
                          c2 = b[i + 20 >> 2] | 0;
                          if (c2 | 0) {
                            b[d2 + 20 >> 2] = c2;
                            b[c2 + 24 >> 2] = d2;
                          }
                        }
                      } while (0);
                      if (j >>> 0 < 16) {
                        v2 = j + k | 0;
                        b[i + 4 >> 2] = v2 | 3;
                        v2 = i + v2 + 4 | 0;
                        b[v2 >> 2] = b[v2 >> 2] | 1;
                      } else {
                        b[i + 4 >> 2] = k | 3;
                        b[h + 4 >> 2] = j | 1;
                        b[h + j >> 2] = j;
                        if (l | 0) {
                          e2 = b[5834] | 0;
                          c2 = l >>> 3;
                          d2 = 23356 + (c2 << 1 << 2) | 0;
                          c2 = 1 << c2;
                          if (!(c2 & m)) {
                            b[5829] = c2 | m;
                            c2 = d2;
                            a2 = d2 + 8 | 0;
                          } else {
                            a2 = d2 + 8 | 0;
                            c2 = b[a2 >> 2] | 0;
                          }
                          b[a2 >> 2] = e2;
                          b[c2 + 12 >> 2] = e2;
                          b[e2 + 8 >> 2] = c2;
                          b[e2 + 12 >> 2] = d2;
                        }
                        b[5831] = j;
                        b[5834] = h;
                      }
                      v2 = i + 8 | 0;
                      S = w2;
                      return v2 | 0;
                    } else {
                      m = k;
                    }
                  } else {
                    m = k;
                  }
                } else {
                  m = k;
                }
              } else if (a2 >>> 0 <= 4294967231) {
                a2 = a2 + 11 | 0;
                k = a2 & -8;
                e2 = b[5830] | 0;
                if (e2) {
                  f = 0 - k | 0;
                  a2 = a2 >>> 8;
                  if (a2) {
                    if (k >>> 0 > 16777215) {
                      j = 31;
                    } else {
                      m = (a2 + 1048320 | 0) >>> 16 & 8;
                      q2 = a2 << m;
                      i = (q2 + 520192 | 0) >>> 16 & 4;
                      q2 = q2 << i;
                      j = (q2 + 245760 | 0) >>> 16 & 2;
                      j = 14 - (i | m | j) + (q2 << j >>> 15) | 0;
                      j = k >>> (j + 7 | 0) & 1 | j << 1;
                    }
                  } else {
                    j = 0;
                  }
                  d2 = b[23620 + (j << 2) >> 2] | 0;
                  a: do {
                    if (!d2) {
                      d2 = 0;
                      a2 = 0;
                      q2 = 61;
                    } else {
                      a2 = 0;
                      i = k << ((j | 0) == 31 ? 0 : 25 - (j >>> 1) | 0);
                      g2 = 0;
                      while (1) {
                        h = (b[d2 + 4 >> 2] & -8) - k | 0;
                        if (h >>> 0 < f >>> 0) {
                          if (!h) {
                            a2 = d2;
                            f = 0;
                            q2 = 65;
                            break a;
                          } else {
                            a2 = d2;
                            f = h;
                          }
                        }
                        q2 = b[d2 + 20 >> 2] | 0;
                        d2 = b[d2 + 16 + (i >>> 31 << 2) >> 2] | 0;
                        g2 = (q2 | 0) == 0 | (q2 | 0) == (d2 | 0) ? g2 : q2;
                        if (!d2) {
                          d2 = g2;
                          q2 = 61;
                          break;
                        } else {
                          i = i << 1;
                        }
                      }
                    }
                  } while (0);
                  if ((q2 | 0) == 61) {
                    if ((d2 | 0) == 0 & (a2 | 0) == 0) {
                      a2 = 2 << j;
                      a2 = (a2 | 0 - a2) & e2;
                      if (!a2) {
                        m = k;
                        break;
                      }
                      m = (a2 & 0 - a2) + -1 | 0;
                      h = m >>> 12 & 16;
                      m = m >>> h;
                      g2 = m >>> 5 & 8;
                      m = m >>> g2;
                      i = m >>> 2 & 4;
                      m = m >>> i;
                      j = m >>> 1 & 2;
                      m = m >>> j;
                      d2 = m >>> 1 & 1;
                      a2 = 0;
                      d2 = b[23620 + ((g2 | h | i | j | d2) + (m >>> d2) << 2) >> 2] | 0;
                    }
                    if (!d2) {
                      i = a2;
                      h = f;
                    } else {
                      q2 = 65;
                    }
                  }
                  if ((q2 | 0) == 65) {
                    g2 = d2;
                    while (1) {
                      m = (b[g2 + 4 >> 2] & -8) - k | 0;
                      d2 = m >>> 0 < f >>> 0;
                      f = d2 ? m : f;
                      a2 = d2 ? g2 : a2;
                      d2 = b[g2 + 16 >> 2] | 0;
                      if (!d2) {
                        d2 = b[g2 + 20 >> 2] | 0;
                      }
                      if (!d2) {
                        i = a2;
                        h = f;
                        break;
                      } else {
                        g2 = d2;
                      }
                    }
                  }
                  if (((i | 0) != 0 ? h >>> 0 < ((b[5831] | 0) - k | 0) >>> 0 : 0) ? (l = i + k | 0, l >>> 0 > i >>> 0) : 0) {
                    g2 = b[i + 24 >> 2] | 0;
                    c2 = b[i + 12 >> 2] | 0;
                    do {
                      if ((c2 | 0) == (i | 0)) {
                        a2 = i + 20 | 0;
                        c2 = b[a2 >> 2] | 0;
                        if (!c2) {
                          a2 = i + 16 | 0;
                          c2 = b[a2 >> 2] | 0;
                          if (!c2) {
                            c2 = 0;
                            break;
                          }
                        }
                        while (1) {
                          f = c2 + 20 | 0;
                          d2 = b[f >> 2] | 0;
                          if (!d2) {
                            f = c2 + 16 | 0;
                            d2 = b[f >> 2] | 0;
                            if (!d2) {
                              break;
                            } else {
                              c2 = d2;
                              a2 = f;
                            }
                          } else {
                            c2 = d2;
                            a2 = f;
                          }
                        }
                        b[a2 >> 2] = 0;
                      } else {
                        v2 = b[i + 8 >> 2] | 0;
                        b[v2 + 12 >> 2] = c2;
                        b[c2 + 8 >> 2] = v2;
                      }
                    } while (0);
                    do {
                      if (g2) {
                        a2 = b[i + 28 >> 2] | 0;
                        d2 = 23620 + (a2 << 2) | 0;
                        if ((i | 0) == (b[d2 >> 2] | 0)) {
                          b[d2 >> 2] = c2;
                          if (!c2) {
                            e2 = e2 & ~(1 << a2);
                            b[5830] = e2;
                            break;
                          }
                        } else {
                          v2 = g2 + 16 | 0;
                          b[((b[v2 >> 2] | 0) == (i | 0) ? v2 : g2 + 20 | 0) >> 2] = c2;
                          if (!c2) {
                            break;
                          }
                        }
                        b[c2 + 24 >> 2] = g2;
                        a2 = b[i + 16 >> 2] | 0;
                        if (a2 | 0) {
                          b[c2 + 16 >> 2] = a2;
                          b[a2 + 24 >> 2] = c2;
                        }
                        a2 = b[i + 20 >> 2] | 0;
                        if (a2) {
                          b[c2 + 20 >> 2] = a2;
                          b[a2 + 24 >> 2] = c2;
                        }
                      }
                    } while (0);
                    b: do {
                      if (h >>> 0 < 16) {
                        v2 = h + k | 0;
                        b[i + 4 >> 2] = v2 | 3;
                        v2 = i + v2 + 4 | 0;
                        b[v2 >> 2] = b[v2 >> 2] | 1;
                      } else {
                        b[i + 4 >> 2] = k | 3;
                        b[l + 4 >> 2] = h | 1;
                        b[l + h >> 2] = h;
                        c2 = h >>> 3;
                        if (h >>> 0 < 256) {
                          d2 = 23356 + (c2 << 1 << 2) | 0;
                          a2 = b[5829] | 0;
                          c2 = 1 << c2;
                          if (!(a2 & c2)) {
                            b[5829] = a2 | c2;
                            c2 = d2;
                            a2 = d2 + 8 | 0;
                          } else {
                            a2 = d2 + 8 | 0;
                            c2 = b[a2 >> 2] | 0;
                          }
                          b[a2 >> 2] = l;
                          b[c2 + 12 >> 2] = l;
                          b[l + 8 >> 2] = c2;
                          b[l + 12 >> 2] = d2;
                          break;
                        }
                        c2 = h >>> 8;
                        if (c2) {
                          if (h >>> 0 > 16777215) {
                            d2 = 31;
                          } else {
                            u2 = (c2 + 1048320 | 0) >>> 16 & 8;
                            v2 = c2 << u2;
                            t2 = (v2 + 520192 | 0) >>> 16 & 4;
                            v2 = v2 << t2;
                            d2 = (v2 + 245760 | 0) >>> 16 & 2;
                            d2 = 14 - (t2 | u2 | d2) + (v2 << d2 >>> 15) | 0;
                            d2 = h >>> (d2 + 7 | 0) & 1 | d2 << 1;
                          }
                        } else {
                          d2 = 0;
                        }
                        c2 = 23620 + (d2 << 2) | 0;
                        b[l + 28 >> 2] = d2;
                        a2 = l + 16 | 0;
                        b[a2 + 4 >> 2] = 0;
                        b[a2 >> 2] = 0;
                        a2 = 1 << d2;
                        if (!(e2 & a2)) {
                          b[5830] = e2 | a2;
                          b[c2 >> 2] = l;
                          b[l + 24 >> 2] = c2;
                          b[l + 12 >> 2] = l;
                          b[l + 8 >> 2] = l;
                          break;
                        }
                        c2 = b[c2 >> 2] | 0;
                        c: do {
                          if ((b[c2 + 4 >> 2] & -8 | 0) != (h | 0)) {
                            e2 = h << ((d2 | 0) == 31 ? 0 : 25 - (d2 >>> 1) | 0);
                            while (1) {
                              d2 = c2 + 16 + (e2 >>> 31 << 2) | 0;
                              a2 = b[d2 >> 2] | 0;
                              if (!a2) {
                                break;
                              }
                              if ((b[a2 + 4 >> 2] & -8 | 0) == (h | 0)) {
                                c2 = a2;
                                break c;
                              } else {
                                e2 = e2 << 1;
                                c2 = a2;
                              }
                            }
                            b[d2 >> 2] = l;
                            b[l + 24 >> 2] = c2;
                            b[l + 12 >> 2] = l;
                            b[l + 8 >> 2] = l;
                            break b;
                          }
                        } while (0);
                        u2 = c2 + 8 | 0;
                        v2 = b[u2 >> 2] | 0;
                        b[v2 + 12 >> 2] = l;
                        b[u2 >> 2] = l;
                        b[l + 8 >> 2] = v2;
                        b[l + 12 >> 2] = c2;
                        b[l + 24 >> 2] = 0;
                      }
                    } while (0);
                    v2 = i + 8 | 0;
                    S = w2;
                    return v2 | 0;
                  } else {
                    m = k;
                  }
                } else {
                  m = k;
                }
              } else {
                m = -1;
              }
            } while (0);
            d2 = b[5831] | 0;
            if (d2 >>> 0 >= m >>> 0) {
              c2 = d2 - m | 0;
              a2 = b[5834] | 0;
              if (c2 >>> 0 > 15) {
                v2 = a2 + m | 0;
                b[5834] = v2;
                b[5831] = c2;
                b[v2 + 4 >> 2] = c2 | 1;
                b[a2 + d2 >> 2] = c2;
                b[a2 + 4 >> 2] = m | 3;
              } else {
                b[5831] = 0;
                b[5834] = 0;
                b[a2 + 4 >> 2] = d2 | 3;
                v2 = a2 + d2 + 4 | 0;
                b[v2 >> 2] = b[v2 >> 2] | 1;
              }
              v2 = a2 + 8 | 0;
              S = w2;
              return v2 | 0;
            }
            h = b[5832] | 0;
            if (h >>> 0 > m >>> 0) {
              t2 = h - m | 0;
              b[5832] = t2;
              v2 = b[5835] | 0;
              u2 = v2 + m | 0;
              b[5835] = u2;
              b[u2 + 4 >> 2] = t2 | 1;
              b[v2 + 4 >> 2] = m | 3;
              v2 = v2 + 8 | 0;
              S = w2;
              return v2 | 0;
            }
            if (!(b[5947] | 0)) {
              b[5949] = 4096;
              b[5948] = 4096;
              b[5950] = -1;
              b[5951] = -1;
              b[5952] = 0;
              b[5940] = 0;
              b[5947] = n & -16 ^ 1431655768;
              a2 = 4096;
            } else {
              a2 = b[5949] | 0;
            }
            i = m + 48 | 0;
            j = m + 47 | 0;
            g2 = a2 + j | 0;
            f = 0 - a2 | 0;
            k = g2 & f;
            if (k >>> 0 <= m >>> 0) {
              v2 = 0;
              S = w2;
              return v2 | 0;
            }
            a2 = b[5939] | 0;
            if (a2 | 0 ? (l = b[5937] | 0, n = l + k | 0, n >>> 0 <= l >>> 0 | n >>> 0 > a2 >>> 0) : 0) {
              v2 = 0;
              S = w2;
              return v2 | 0;
            }
            d: do {
              if (!(b[5940] & 4)) {
                d2 = b[5835] | 0;
                e: do {
                  if (d2) {
                    e2 = 23764;
                    while (1) {
                      n = b[e2 >> 2] | 0;
                      if (n >>> 0 <= d2 >>> 0 ? (n + (b[e2 + 4 >> 2] | 0) | 0) >>> 0 > d2 >>> 0 : 0) {
                        break;
                      }
                      a2 = b[e2 + 8 >> 2] | 0;
                      if (!a2) {
                        q2 = 128;
                        break e;
                      } else {
                        e2 = a2;
                      }
                    }
                    c2 = g2 - h & f;
                    if (c2 >>> 0 < 2147483647) {
                      a2 = jd(c2 | 0) | 0;
                      if ((a2 | 0) == ((b[e2 >> 2] | 0) + (b[e2 + 4 >> 2] | 0) | 0)) {
                        if ((a2 | 0) != (-1 | 0)) {
                          h = c2;
                          g2 = a2;
                          q2 = 145;
                          break d;
                        }
                      } else {
                        e2 = a2;
                        q2 = 136;
                      }
                    } else {
                      c2 = 0;
                    }
                  } else {
                    q2 = 128;
                  }
                } while (0);
                do {
                  if ((q2 | 0) == 128) {
                    d2 = jd(0) | 0;
                    if ((d2 | 0) != (-1 | 0) ? (c2 = d2, o = b[5948] | 0, p2 = o + -1 | 0, c2 = ((p2 & c2 | 0) == 0 ? 0 : (p2 + c2 & 0 - o) - c2 | 0) + k | 0, o = b[5937] | 0, p2 = c2 + o | 0, c2 >>> 0 > m >>> 0 & c2 >>> 0 < 2147483647) : 0) {
                      n = b[5939] | 0;
                      if (n | 0 ? p2 >>> 0 <= o >>> 0 | p2 >>> 0 > n >>> 0 : 0) {
                        c2 = 0;
                        break;
                      }
                      a2 = jd(c2 | 0) | 0;
                      if ((a2 | 0) == (d2 | 0)) {
                        h = c2;
                        g2 = d2;
                        q2 = 145;
                        break d;
                      } else {
                        e2 = a2;
                        q2 = 136;
                      }
                    } else {
                      c2 = 0;
                    }
                  }
                } while (0);
                do {
                  if ((q2 | 0) == 136) {
                    d2 = 0 - c2 | 0;
                    if (!(i >>> 0 > c2 >>> 0 & (c2 >>> 0 < 2147483647 & (e2 | 0) != (-1 | 0)))) {
                      if ((e2 | 0) == (-1 | 0)) {
                        c2 = 0;
                        break;
                      } else {
                        h = c2;
                        g2 = e2;
                        q2 = 145;
                        break d;
                      }
                    }
                    a2 = b[5949] | 0;
                    a2 = j - c2 + a2 & 0 - a2;
                    if (a2 >>> 0 >= 2147483647) {
                      h = c2;
                      g2 = e2;
                      q2 = 145;
                      break d;
                    }
                    if ((jd(a2 | 0) | 0) == (-1 | 0)) {
                      jd(d2 | 0) | 0;
                      c2 = 0;
                      break;
                    } else {
                      h = a2 + c2 | 0;
                      g2 = e2;
                      q2 = 145;
                      break d;
                    }
                  }
                } while (0);
                b[5940] = b[5940] | 4;
                q2 = 143;
              } else {
                c2 = 0;
                q2 = 143;
              }
            } while (0);
            if (((q2 | 0) == 143 ? k >>> 0 < 2147483647 : 0) ? (t2 = jd(k | 0) | 0, p2 = jd(0) | 0, r2 = p2 - t2 | 0, s2 = r2 >>> 0 > (m + 40 | 0) >>> 0, !((t2 | 0) == (-1 | 0) | s2 ^ 1 | t2 >>> 0 < p2 >>> 0 & ((t2 | 0) != (-1 | 0) & (p2 | 0) != (-1 | 0)) ^ 1)) : 0) {
              h = s2 ? r2 : c2;
              g2 = t2;
              q2 = 145;
            }
            if ((q2 | 0) == 145) {
              c2 = (b[5937] | 0) + h | 0;
              b[5937] = c2;
              if (c2 >>> 0 > (b[5938] | 0) >>> 0) {
                b[5938] = c2;
              }
              j = b[5835] | 0;
              f: do {
                if (j) {
                  c2 = 23764;
                  while (1) {
                    a2 = b[c2 >> 2] | 0;
                    d2 = b[c2 + 4 >> 2] | 0;
                    if ((g2 | 0) == (a2 + d2 | 0)) {
                      q2 = 154;
                      break;
                    }
                    e2 = b[c2 + 8 >> 2] | 0;
                    if (!e2) {
                      break;
                    } else {
                      c2 = e2;
                    }
                  }
                  if (((q2 | 0) == 154 ? (u2 = c2 + 4 | 0, (b[c2 + 12 >> 2] & 8 | 0) == 0) : 0) ? g2 >>> 0 > j >>> 0 & a2 >>> 0 <= j >>> 0 : 0) {
                    b[u2 >> 2] = d2 + h;
                    v2 = (b[5832] | 0) + h | 0;
                    t2 = j + 8 | 0;
                    t2 = (t2 & 7 | 0) == 0 ? 0 : 0 - t2 & 7;
                    u2 = j + t2 | 0;
                    t2 = v2 - t2 | 0;
                    b[5835] = u2;
                    b[5832] = t2;
                    b[u2 + 4 >> 2] = t2 | 1;
                    b[j + v2 + 4 >> 2] = 40;
                    b[5836] = b[5951];
                    break;
                  }
                  if (g2 >>> 0 < (b[5833] | 0) >>> 0) {
                    b[5833] = g2;
                  }
                  d2 = g2 + h | 0;
                  c2 = 23764;
                  while (1) {
                    if ((b[c2 >> 2] | 0) == (d2 | 0)) {
                      q2 = 162;
                      break;
                    }
                    a2 = b[c2 + 8 >> 2] | 0;
                    if (!a2) {
                      break;
                    } else {
                      c2 = a2;
                    }
                  }
                  if ((q2 | 0) == 162 ? (b[c2 + 12 >> 2] & 8 | 0) == 0 : 0) {
                    b[c2 >> 2] = g2;
                    l = c2 + 4 | 0;
                    b[l >> 2] = (b[l >> 2] | 0) + h;
                    l = g2 + 8 | 0;
                    l = g2 + ((l & 7 | 0) == 0 ? 0 : 0 - l & 7) | 0;
                    c2 = d2 + 8 | 0;
                    c2 = d2 + ((c2 & 7 | 0) == 0 ? 0 : 0 - c2 & 7) | 0;
                    k = l + m | 0;
                    i = c2 - l - m | 0;
                    b[l + 4 >> 2] = m | 3;
                    g: do {
                      if ((j | 0) == (c2 | 0)) {
                        v2 = (b[5832] | 0) + i | 0;
                        b[5832] = v2;
                        b[5835] = k;
                        b[k + 4 >> 2] = v2 | 1;
                      } else {
                        if ((b[5834] | 0) == (c2 | 0)) {
                          v2 = (b[5831] | 0) + i | 0;
                          b[5831] = v2;
                          b[5834] = k;
                          b[k + 4 >> 2] = v2 | 1;
                          b[k + v2 >> 2] = v2;
                          break;
                        }
                        a2 = b[c2 + 4 >> 2] | 0;
                        if ((a2 & 3 | 0) == 1) {
                          h = a2 & -8;
                          e2 = a2 >>> 3;
                          h: do {
                            if (a2 >>> 0 < 256) {
                              a2 = b[c2 + 8 >> 2] | 0;
                              d2 = b[c2 + 12 >> 2] | 0;
                              if ((d2 | 0) == (a2 | 0)) {
                                b[5829] = b[5829] & ~(1 << e2);
                                break;
                              } else {
                                b[a2 + 12 >> 2] = d2;
                                b[d2 + 8 >> 2] = a2;
                                break;
                              }
                            } else {
                              g2 = b[c2 + 24 >> 2] | 0;
                              a2 = b[c2 + 12 >> 2] | 0;
                              do {
                                if ((a2 | 0) == (c2 | 0)) {
                                  d2 = c2 + 16 | 0;
                                  e2 = d2 + 4 | 0;
                                  a2 = b[e2 >> 2] | 0;
                                  if (!a2) {
                                    a2 = b[d2 >> 2] | 0;
                                    if (!a2) {
                                      a2 = 0;
                                      break;
                                    }
                                  } else {
                                    d2 = e2;
                                  }
                                  while (1) {
                                    f = a2 + 20 | 0;
                                    e2 = b[f >> 2] | 0;
                                    if (!e2) {
                                      f = a2 + 16 | 0;
                                      e2 = b[f >> 2] | 0;
                                      if (!e2) {
                                        break;
                                      } else {
                                        a2 = e2;
                                        d2 = f;
                                      }
                                    } else {
                                      a2 = e2;
                                      d2 = f;
                                    }
                                  }
                                  b[d2 >> 2] = 0;
                                } else {
                                  v2 = b[c2 + 8 >> 2] | 0;
                                  b[v2 + 12 >> 2] = a2;
                                  b[a2 + 8 >> 2] = v2;
                                }
                              } while (0);
                              if (!g2) {
                                break;
                              }
                              d2 = b[c2 + 28 >> 2] | 0;
                              e2 = 23620 + (d2 << 2) | 0;
                              do {
                                if ((b[e2 >> 2] | 0) != (c2 | 0)) {
                                  v2 = g2 + 16 | 0;
                                  b[((b[v2 >> 2] | 0) == (c2 | 0) ? v2 : g2 + 20 | 0) >> 2] = a2;
                                  if (!a2) {
                                    break h;
                                  }
                                } else {
                                  b[e2 >> 2] = a2;
                                  if (a2 | 0) {
                                    break;
                                  }
                                  b[5830] = b[5830] & ~(1 << d2);
                                  break h;
                                }
                              } while (0);
                              b[a2 + 24 >> 2] = g2;
                              d2 = c2 + 16 | 0;
                              e2 = b[d2 >> 2] | 0;
                              if (e2 | 0) {
                                b[a2 + 16 >> 2] = e2;
                                b[e2 + 24 >> 2] = a2;
                              }
                              d2 = b[d2 + 4 >> 2] | 0;
                              if (!d2) {
                                break;
                              }
                              b[a2 + 20 >> 2] = d2;
                              b[d2 + 24 >> 2] = a2;
                            }
                          } while (0);
                          c2 = c2 + h | 0;
                          f = h + i | 0;
                        } else {
                          f = i;
                        }
                        c2 = c2 + 4 | 0;
                        b[c2 >> 2] = b[c2 >> 2] & -2;
                        b[k + 4 >> 2] = f | 1;
                        b[k + f >> 2] = f;
                        c2 = f >>> 3;
                        if (f >>> 0 < 256) {
                          d2 = 23356 + (c2 << 1 << 2) | 0;
                          a2 = b[5829] | 0;
                          c2 = 1 << c2;
                          if (!(a2 & c2)) {
                            b[5829] = a2 | c2;
                            c2 = d2;
                            a2 = d2 + 8 | 0;
                          } else {
                            a2 = d2 + 8 | 0;
                            c2 = b[a2 >> 2] | 0;
                          }
                          b[a2 >> 2] = k;
                          b[c2 + 12 >> 2] = k;
                          b[k + 8 >> 2] = c2;
                          b[k + 12 >> 2] = d2;
                          break;
                        }
                        c2 = f >>> 8;
                        do {
                          if (!c2) {
                            e2 = 0;
                          } else {
                            if (f >>> 0 > 16777215) {
                              e2 = 31;
                              break;
                            }
                            u2 = (c2 + 1048320 | 0) >>> 16 & 8;
                            v2 = c2 << u2;
                            t2 = (v2 + 520192 | 0) >>> 16 & 4;
                            v2 = v2 << t2;
                            e2 = (v2 + 245760 | 0) >>> 16 & 2;
                            e2 = 14 - (t2 | u2 | e2) + (v2 << e2 >>> 15) | 0;
                            e2 = f >>> (e2 + 7 | 0) & 1 | e2 << 1;
                          }
                        } while (0);
                        c2 = 23620 + (e2 << 2) | 0;
                        b[k + 28 >> 2] = e2;
                        a2 = k + 16 | 0;
                        b[a2 + 4 >> 2] = 0;
                        b[a2 >> 2] = 0;
                        a2 = b[5830] | 0;
                        d2 = 1 << e2;
                        if (!(a2 & d2)) {
                          b[5830] = a2 | d2;
                          b[c2 >> 2] = k;
                          b[k + 24 >> 2] = c2;
                          b[k + 12 >> 2] = k;
                          b[k + 8 >> 2] = k;
                          break;
                        }
                        c2 = b[c2 >> 2] | 0;
                        i: do {
                          if ((b[c2 + 4 >> 2] & -8 | 0) != (f | 0)) {
                            e2 = f << ((e2 | 0) == 31 ? 0 : 25 - (e2 >>> 1) | 0);
                            while (1) {
                              d2 = c2 + 16 + (e2 >>> 31 << 2) | 0;
                              a2 = b[d2 >> 2] | 0;
                              if (!a2) {
                                break;
                              }
                              if ((b[a2 + 4 >> 2] & -8 | 0) == (f | 0)) {
                                c2 = a2;
                                break i;
                              } else {
                                e2 = e2 << 1;
                                c2 = a2;
                              }
                            }
                            b[d2 >> 2] = k;
                            b[k + 24 >> 2] = c2;
                            b[k + 12 >> 2] = k;
                            b[k + 8 >> 2] = k;
                            break g;
                          }
                        } while (0);
                        u2 = c2 + 8 | 0;
                        v2 = b[u2 >> 2] | 0;
                        b[v2 + 12 >> 2] = k;
                        b[u2 >> 2] = k;
                        b[k + 8 >> 2] = v2;
                        b[k + 12 >> 2] = c2;
                        b[k + 24 >> 2] = 0;
                      }
                    } while (0);
                    v2 = l + 8 | 0;
                    S = w2;
                    return v2 | 0;
                  }
                  c2 = 23764;
                  while (1) {
                    a2 = b[c2 >> 2] | 0;
                    if (a2 >>> 0 <= j >>> 0 ? (v2 = a2 + (b[c2 + 4 >> 2] | 0) | 0, v2 >>> 0 > j >>> 0) : 0) {
                      break;
                    }
                    c2 = b[c2 + 8 >> 2] | 0;
                  }
                  f = v2 + -47 | 0;
                  a2 = f + 8 | 0;
                  a2 = f + ((a2 & 7 | 0) == 0 ? 0 : 0 - a2 & 7) | 0;
                  f = j + 16 | 0;
                  a2 = a2 >>> 0 < f >>> 0 ? j : a2;
                  c2 = a2 + 8 | 0;
                  d2 = h + -40 | 0;
                  t2 = g2 + 8 | 0;
                  t2 = (t2 & 7 | 0) == 0 ? 0 : 0 - t2 & 7;
                  u2 = g2 + t2 | 0;
                  t2 = d2 - t2 | 0;
                  b[5835] = u2;
                  b[5832] = t2;
                  b[u2 + 4 >> 2] = t2 | 1;
                  b[g2 + d2 + 4 >> 2] = 40;
                  b[5836] = b[5951];
                  d2 = a2 + 4 | 0;
                  b[d2 >> 2] = 27;
                  b[c2 >> 2] = b[5941];
                  b[c2 + 4 >> 2] = b[5942];
                  b[c2 + 8 >> 2] = b[5943];
                  b[c2 + 12 >> 2] = b[5944];
                  b[5941] = g2;
                  b[5942] = h;
                  b[5944] = 0;
                  b[5943] = c2;
                  c2 = a2 + 24 | 0;
                  do {
                    u2 = c2;
                    c2 = c2 + 4 | 0;
                    b[c2 >> 2] = 7;
                  } while ((u2 + 8 | 0) >>> 0 < v2 >>> 0);
                  if ((a2 | 0) != (j | 0)) {
                    g2 = a2 - j | 0;
                    b[d2 >> 2] = b[d2 >> 2] & -2;
                    b[j + 4 >> 2] = g2 | 1;
                    b[a2 >> 2] = g2;
                    c2 = g2 >>> 3;
                    if (g2 >>> 0 < 256) {
                      d2 = 23356 + (c2 << 1 << 2) | 0;
                      a2 = b[5829] | 0;
                      c2 = 1 << c2;
                      if (!(a2 & c2)) {
                        b[5829] = a2 | c2;
                        c2 = d2;
                        a2 = d2 + 8 | 0;
                      } else {
                        a2 = d2 + 8 | 0;
                        c2 = b[a2 >> 2] | 0;
                      }
                      b[a2 >> 2] = j;
                      b[c2 + 12 >> 2] = j;
                      b[j + 8 >> 2] = c2;
                      b[j + 12 >> 2] = d2;
                      break;
                    }
                    c2 = g2 >>> 8;
                    if (c2) {
                      if (g2 >>> 0 > 16777215) {
                        e2 = 31;
                      } else {
                        u2 = (c2 + 1048320 | 0) >>> 16 & 8;
                        v2 = c2 << u2;
                        t2 = (v2 + 520192 | 0) >>> 16 & 4;
                        v2 = v2 << t2;
                        e2 = (v2 + 245760 | 0) >>> 16 & 2;
                        e2 = 14 - (t2 | u2 | e2) + (v2 << e2 >>> 15) | 0;
                        e2 = g2 >>> (e2 + 7 | 0) & 1 | e2 << 1;
                      }
                    } else {
                      e2 = 0;
                    }
                    d2 = 23620 + (e2 << 2) | 0;
                    b[j + 28 >> 2] = e2;
                    b[j + 20 >> 2] = 0;
                    b[f >> 2] = 0;
                    c2 = b[5830] | 0;
                    a2 = 1 << e2;
                    if (!(c2 & a2)) {
                      b[5830] = c2 | a2;
                      b[d2 >> 2] = j;
                      b[j + 24 >> 2] = d2;
                      b[j + 12 >> 2] = j;
                      b[j + 8 >> 2] = j;
                      break;
                    }
                    c2 = b[d2 >> 2] | 0;
                    j: do {
                      if ((b[c2 + 4 >> 2] & -8 | 0) != (g2 | 0)) {
                        e2 = g2 << ((e2 | 0) == 31 ? 0 : 25 - (e2 >>> 1) | 0);
                        while (1) {
                          d2 = c2 + 16 + (e2 >>> 31 << 2) | 0;
                          a2 = b[d2 >> 2] | 0;
                          if (!a2) {
                            break;
                          }
                          if ((b[a2 + 4 >> 2] & -8 | 0) == (g2 | 0)) {
                            c2 = a2;
                            break j;
                          } else {
                            e2 = e2 << 1;
                            c2 = a2;
                          }
                        }
                        b[d2 >> 2] = j;
                        b[j + 24 >> 2] = c2;
                        b[j + 12 >> 2] = j;
                        b[j + 8 >> 2] = j;
                        break f;
                      }
                    } while (0);
                    u2 = c2 + 8 | 0;
                    v2 = b[u2 >> 2] | 0;
                    b[v2 + 12 >> 2] = j;
                    b[u2 >> 2] = j;
                    b[j + 8 >> 2] = v2;
                    b[j + 12 >> 2] = c2;
                    b[j + 24 >> 2] = 0;
                  }
                } else {
                  v2 = b[5833] | 0;
                  if ((v2 | 0) == 0 | g2 >>> 0 < v2 >>> 0) {
                    b[5833] = g2;
                  }
                  b[5941] = g2;
                  b[5942] = h;
                  b[5944] = 0;
                  b[5838] = b[5947];
                  b[5837] = -1;
                  b[5842] = 23356;
                  b[5841] = 23356;
                  b[5844] = 23364;
                  b[5843] = 23364;
                  b[5846] = 23372;
                  b[5845] = 23372;
                  b[5848] = 23380;
                  b[5847] = 23380;
                  b[5850] = 23388;
                  b[5849] = 23388;
                  b[5852] = 23396;
                  b[5851] = 23396;
                  b[5854] = 23404;
                  b[5853] = 23404;
                  b[5856] = 23412;
                  b[5855] = 23412;
                  b[5858] = 23420;
                  b[5857] = 23420;
                  b[5860] = 23428;
                  b[5859] = 23428;
                  b[5862] = 23436;
                  b[5861] = 23436;
                  b[5864] = 23444;
                  b[5863] = 23444;
                  b[5866] = 23452;
                  b[5865] = 23452;
                  b[5868] = 23460;
                  b[5867] = 23460;
                  b[5870] = 23468;
                  b[5869] = 23468;
                  b[5872] = 23476;
                  b[5871] = 23476;
                  b[5874] = 23484;
                  b[5873] = 23484;
                  b[5876] = 23492;
                  b[5875] = 23492;
                  b[5878] = 23500;
                  b[5877] = 23500;
                  b[5880] = 23508;
                  b[5879] = 23508;
                  b[5882] = 23516;
                  b[5881] = 23516;
                  b[5884] = 23524;
                  b[5883] = 23524;
                  b[5886] = 23532;
                  b[5885] = 23532;
                  b[5888] = 23540;
                  b[5887] = 23540;
                  b[5890] = 23548;
                  b[5889] = 23548;
                  b[5892] = 23556;
                  b[5891] = 23556;
                  b[5894] = 23564;
                  b[5893] = 23564;
                  b[5896] = 23572;
                  b[5895] = 23572;
                  b[5898] = 23580;
                  b[5897] = 23580;
                  b[5900] = 23588;
                  b[5899] = 23588;
                  b[5902] = 23596;
                  b[5901] = 23596;
                  b[5904] = 23604;
                  b[5903] = 23604;
                  v2 = h + -40 | 0;
                  t2 = g2 + 8 | 0;
                  t2 = (t2 & 7 | 0) == 0 ? 0 : 0 - t2 & 7;
                  u2 = g2 + t2 | 0;
                  t2 = v2 - t2 | 0;
                  b[5835] = u2;
                  b[5832] = t2;
                  b[u2 + 4 >> 2] = t2 | 1;
                  b[g2 + v2 + 4 >> 2] = 40;
                  b[5836] = b[5951];
                }
              } while (0);
              c2 = b[5832] | 0;
              if (c2 >>> 0 > m >>> 0) {
                t2 = c2 - m | 0;
                b[5832] = t2;
                v2 = b[5835] | 0;
                u2 = v2 + m | 0;
                b[5835] = u2;
                b[u2 + 4 >> 2] = t2 | 1;
                b[v2 + 4 >> 2] = m | 3;
                v2 = v2 + 8 | 0;
                S = w2;
                return v2 | 0;
              }
            }
            v2 = Tc() | 0;
            b[v2 >> 2] = 12;
            v2 = 0;
            S = w2;
            return v2 | 0;
          }
          function Xc(a2) {
            a2 = a2 | 0;
            var c2 = 0, d2 = 0, e2 = 0, f = 0, g2 = 0, h = 0, i = 0, j = 0;
            if (!a2) {
              return;
            }
            d2 = a2 + -8 | 0;
            f = b[5833] | 0;
            a2 = b[a2 + -4 >> 2] | 0;
            c2 = a2 & -8;
            j = d2 + c2 | 0;
            do {
              if (!(a2 & 1)) {
                e2 = b[d2 >> 2] | 0;
                if (!(a2 & 3)) {
                  return;
                }
                h = d2 + (0 - e2) | 0;
                g2 = e2 + c2 | 0;
                if (h >>> 0 < f >>> 0) {
                  return;
                }
                if ((b[5834] | 0) == (h | 0)) {
                  a2 = j + 4 | 0;
                  c2 = b[a2 >> 2] | 0;
                  if ((c2 & 3 | 0) != 3) {
                    i = h;
                    c2 = g2;
                    break;
                  }
                  b[5831] = g2;
                  b[a2 >> 2] = c2 & -2;
                  b[h + 4 >> 2] = g2 | 1;
                  b[h + g2 >> 2] = g2;
                  return;
                }
                d2 = e2 >>> 3;
                if (e2 >>> 0 < 256) {
                  a2 = b[h + 8 >> 2] | 0;
                  c2 = b[h + 12 >> 2] | 0;
                  if ((c2 | 0) == (a2 | 0)) {
                    b[5829] = b[5829] & ~(1 << d2);
                    i = h;
                    c2 = g2;
                    break;
                  } else {
                    b[a2 + 12 >> 2] = c2;
                    b[c2 + 8 >> 2] = a2;
                    i = h;
                    c2 = g2;
                    break;
                  }
                }
                f = b[h + 24 >> 2] | 0;
                a2 = b[h + 12 >> 2] | 0;
                do {
                  if ((a2 | 0) == (h | 0)) {
                    c2 = h + 16 | 0;
                    d2 = c2 + 4 | 0;
                    a2 = b[d2 >> 2] | 0;
                    if (!a2) {
                      a2 = b[c2 >> 2] | 0;
                      if (!a2) {
                        a2 = 0;
                        break;
                      }
                    } else {
                      c2 = d2;
                    }
                    while (1) {
                      e2 = a2 + 20 | 0;
                      d2 = b[e2 >> 2] | 0;
                      if (!d2) {
                        e2 = a2 + 16 | 0;
                        d2 = b[e2 >> 2] | 0;
                        if (!d2) {
                          break;
                        } else {
                          a2 = d2;
                          c2 = e2;
                        }
                      } else {
                        a2 = d2;
                        c2 = e2;
                      }
                    }
                    b[c2 >> 2] = 0;
                  } else {
                    i = b[h + 8 >> 2] | 0;
                    b[i + 12 >> 2] = a2;
                    b[a2 + 8 >> 2] = i;
                  }
                } while (0);
                if (f) {
                  c2 = b[h + 28 >> 2] | 0;
                  d2 = 23620 + (c2 << 2) | 0;
                  if ((b[d2 >> 2] | 0) == (h | 0)) {
                    b[d2 >> 2] = a2;
                    if (!a2) {
                      b[5830] = b[5830] & ~(1 << c2);
                      i = h;
                      c2 = g2;
                      break;
                    }
                  } else {
                    i = f + 16 | 0;
                    b[((b[i >> 2] | 0) == (h | 0) ? i : f + 20 | 0) >> 2] = a2;
                    if (!a2) {
                      i = h;
                      c2 = g2;
                      break;
                    }
                  }
                  b[a2 + 24 >> 2] = f;
                  c2 = h + 16 | 0;
                  d2 = b[c2 >> 2] | 0;
                  if (d2 | 0) {
                    b[a2 + 16 >> 2] = d2;
                    b[d2 + 24 >> 2] = a2;
                  }
                  c2 = b[c2 + 4 >> 2] | 0;
                  if (c2) {
                    b[a2 + 20 >> 2] = c2;
                    b[c2 + 24 >> 2] = a2;
                    i = h;
                    c2 = g2;
                  } else {
                    i = h;
                    c2 = g2;
                  }
                } else {
                  i = h;
                  c2 = g2;
                }
              } else {
                i = d2;
                h = d2;
              }
            } while (0);
            if (h >>> 0 >= j >>> 0) {
              return;
            }
            a2 = j + 4 | 0;
            e2 = b[a2 >> 2] | 0;
            if (!(e2 & 1)) {
              return;
            }
            if (!(e2 & 2)) {
              if ((b[5835] | 0) == (j | 0)) {
                j = (b[5832] | 0) + c2 | 0;
                b[5832] = j;
                b[5835] = i;
                b[i + 4 >> 2] = j | 1;
                if ((i | 0) != (b[5834] | 0)) {
                  return;
                }
                b[5834] = 0;
                b[5831] = 0;
                return;
              }
              if ((b[5834] | 0) == (j | 0)) {
                j = (b[5831] | 0) + c2 | 0;
                b[5831] = j;
                b[5834] = h;
                b[i + 4 >> 2] = j | 1;
                b[h + j >> 2] = j;
                return;
              }
              f = (e2 & -8) + c2 | 0;
              d2 = e2 >>> 3;
              do {
                if (e2 >>> 0 < 256) {
                  c2 = b[j + 8 >> 2] | 0;
                  a2 = b[j + 12 >> 2] | 0;
                  if ((a2 | 0) == (c2 | 0)) {
                    b[5829] = b[5829] & ~(1 << d2);
                    break;
                  } else {
                    b[c2 + 12 >> 2] = a2;
                    b[a2 + 8 >> 2] = c2;
                    break;
                  }
                } else {
                  g2 = b[j + 24 >> 2] | 0;
                  a2 = b[j + 12 >> 2] | 0;
                  do {
                    if ((a2 | 0) == (j | 0)) {
                      c2 = j + 16 | 0;
                      d2 = c2 + 4 | 0;
                      a2 = b[d2 >> 2] | 0;
                      if (!a2) {
                        a2 = b[c2 >> 2] | 0;
                        if (!a2) {
                          d2 = 0;
                          break;
                        }
                      } else {
                        c2 = d2;
                      }
                      while (1) {
                        e2 = a2 + 20 | 0;
                        d2 = b[e2 >> 2] | 0;
                        if (!d2) {
                          e2 = a2 + 16 | 0;
                          d2 = b[e2 >> 2] | 0;
                          if (!d2) {
                            break;
                          } else {
                            a2 = d2;
                            c2 = e2;
                          }
                        } else {
                          a2 = d2;
                          c2 = e2;
                        }
                      }
                      b[c2 >> 2] = 0;
                      d2 = a2;
                    } else {
                      d2 = b[j + 8 >> 2] | 0;
                      b[d2 + 12 >> 2] = a2;
                      b[a2 + 8 >> 2] = d2;
                      d2 = a2;
                    }
                  } while (0);
                  if (g2 | 0) {
                    a2 = b[j + 28 >> 2] | 0;
                    c2 = 23620 + (a2 << 2) | 0;
                    if ((b[c2 >> 2] | 0) == (j | 0)) {
                      b[c2 >> 2] = d2;
                      if (!d2) {
                        b[5830] = b[5830] & ~(1 << a2);
                        break;
                      }
                    } else {
                      e2 = g2 + 16 | 0;
                      b[((b[e2 >> 2] | 0) == (j | 0) ? e2 : g2 + 20 | 0) >> 2] = d2;
                      if (!d2) {
                        break;
                      }
                    }
                    b[d2 + 24 >> 2] = g2;
                    a2 = j + 16 | 0;
                    c2 = b[a2 >> 2] | 0;
                    if (c2 | 0) {
                      b[d2 + 16 >> 2] = c2;
                      b[c2 + 24 >> 2] = d2;
                    }
                    a2 = b[a2 + 4 >> 2] | 0;
                    if (a2 | 0) {
                      b[d2 + 20 >> 2] = a2;
                      b[a2 + 24 >> 2] = d2;
                    }
                  }
                }
              } while (0);
              b[i + 4 >> 2] = f | 1;
              b[h + f >> 2] = f;
              if ((i | 0) == (b[5834] | 0)) {
                b[5831] = f;
                return;
              }
            } else {
              b[a2 >> 2] = e2 & -2;
              b[i + 4 >> 2] = c2 | 1;
              b[h + c2 >> 2] = c2;
              f = c2;
            }
            a2 = f >>> 3;
            if (f >>> 0 < 256) {
              d2 = 23356 + (a2 << 1 << 2) | 0;
              c2 = b[5829] | 0;
              a2 = 1 << a2;
              if (!(c2 & a2)) {
                b[5829] = c2 | a2;
                a2 = d2;
                c2 = d2 + 8 | 0;
              } else {
                c2 = d2 + 8 | 0;
                a2 = b[c2 >> 2] | 0;
              }
              b[c2 >> 2] = i;
              b[a2 + 12 >> 2] = i;
              b[i + 8 >> 2] = a2;
              b[i + 12 >> 2] = d2;
              return;
            }
            a2 = f >>> 8;
            if (a2) {
              if (f >>> 0 > 16777215) {
                e2 = 31;
              } else {
                h = (a2 + 1048320 | 0) >>> 16 & 8;
                j = a2 << h;
                g2 = (j + 520192 | 0) >>> 16 & 4;
                j = j << g2;
                e2 = (j + 245760 | 0) >>> 16 & 2;
                e2 = 14 - (g2 | h | e2) + (j << e2 >>> 15) | 0;
                e2 = f >>> (e2 + 7 | 0) & 1 | e2 << 1;
              }
            } else {
              e2 = 0;
            }
            a2 = 23620 + (e2 << 2) | 0;
            b[i + 28 >> 2] = e2;
            b[i + 20 >> 2] = 0;
            b[i + 16 >> 2] = 0;
            c2 = b[5830] | 0;
            d2 = 1 << e2;
            a: do {
              if (!(c2 & d2)) {
                b[5830] = c2 | d2;
                b[a2 >> 2] = i;
                b[i + 24 >> 2] = a2;
                b[i + 12 >> 2] = i;
                b[i + 8 >> 2] = i;
              } else {
                a2 = b[a2 >> 2] | 0;
                b: do {
                  if ((b[a2 + 4 >> 2] & -8 | 0) != (f | 0)) {
                    e2 = f << ((e2 | 0) == 31 ? 0 : 25 - (e2 >>> 1) | 0);
                    while (1) {
                      d2 = a2 + 16 + (e2 >>> 31 << 2) | 0;
                      c2 = b[d2 >> 2] | 0;
                      if (!c2) {
                        break;
                      }
                      if ((b[c2 + 4 >> 2] & -8 | 0) == (f | 0)) {
                        a2 = c2;
                        break b;
                      } else {
                        e2 = e2 << 1;
                        a2 = c2;
                      }
                    }
                    b[d2 >> 2] = i;
                    b[i + 24 >> 2] = a2;
                    b[i + 12 >> 2] = i;
                    b[i + 8 >> 2] = i;
                    break a;
                  }
                } while (0);
                h = a2 + 8 | 0;
                j = b[h >> 2] | 0;
                b[j + 12 >> 2] = i;
                b[h >> 2] = i;
                b[i + 8 >> 2] = j;
                b[i + 12 >> 2] = a2;
                b[i + 24 >> 2] = 0;
              }
            } while (0);
            j = (b[5837] | 0) + -1 | 0;
            b[5837] = j;
            if (j | 0) {
              return;
            }
            a2 = 23772;
            while (1) {
              a2 = b[a2 >> 2] | 0;
              if (!a2) {
                break;
              } else {
                a2 = a2 + 8 | 0;
              }
            }
            b[5837] = -1;
            return;
          }
          function Yc(a2, c2) {
            a2 = a2 | 0;
            c2 = c2 | 0;
            var d2 = 0;
            if (a2) {
              d2 = B(c2, a2) | 0;
              if ((c2 | a2) >>> 0 > 65535) {
                d2 = ((d2 >>> 0) / (a2 >>> 0) | 0 | 0) == (c2 | 0) ? d2 : -1;
              }
            } else {
              d2 = 0;
            }
            a2 = Wc(d2) | 0;
            if (!a2) {
              return a2 | 0;
            }
            if (!(b[a2 + -4 >> 2] & 3)) {
              return a2 | 0;
            }
            hd(a2 | 0, 0, d2 | 0) | 0;
            return a2 | 0;
          }
          function Zc(a2, b2, c2, d2) {
            a2 = a2 | 0;
            b2 = b2 | 0;
            c2 = c2 | 0;
            d2 = d2 | 0;
            c2 = a2 + c2 >>> 0;
            return (F(b2 + d2 + (c2 >>> 0 < a2 >>> 0 | 0) >>> 0 | 0), c2 | 0) | 0;
          }
          function _c(a2, b2, c2, d2) {
            a2 = a2 | 0;
            b2 = b2 | 0;
            c2 = c2 | 0;
            d2 = d2 | 0;
            d2 = b2 - d2 - (c2 >>> 0 > a2 >>> 0 | 0) >>> 0;
            return (F(d2 | 0), a2 - c2 >>> 0 | 0) | 0;
          }
          function $c(a2) {
            a2 = a2 | 0;
            return (a2 ? 31 - (D(a2 ^ a2 - 1) | 0) | 0 : 32) | 0;
          }
          function ad(a2, c2, d2, e2, f) {
            a2 = a2 | 0;
            c2 = c2 | 0;
            d2 = d2 | 0;
            e2 = e2 | 0;
            f = f | 0;
            var g2 = 0, h = 0, i = 0, j = 0, k = 0, l = 0, m = 0, n = 0, o = 0, p2 = 0;
            l = a2;
            j = c2;
            k = j;
            h = d2;
            n = e2;
            i = n;
            if (!k) {
              g2 = (f | 0) != 0;
              if (!i) {
                if (g2) {
                  b[f >> 2] = (l >>> 0) % (h >>> 0);
                  b[f + 4 >> 2] = 0;
                }
                n = 0;
                f = (l >>> 0) / (h >>> 0) >>> 0;
                return (F(n | 0), f) | 0;
              } else {
                if (!g2) {
                  n = 0;
                  f = 0;
                  return (F(n | 0), f) | 0;
                }
                b[f >> 2] = a2 | 0;
                b[f + 4 >> 2] = c2 & 0;
                n = 0;
                f = 0;
                return (F(n | 0), f) | 0;
              }
            }
            g2 = (i | 0) == 0;
            do {
              if (h) {
                if (!g2) {
                  g2 = (D(i | 0) | 0) - (D(k | 0) | 0) | 0;
                  if (g2 >>> 0 <= 31) {
                    m = g2 + 1 | 0;
                    i = 31 - g2 | 0;
                    c2 = g2 - 31 >> 31;
                    h = m;
                    a2 = l >>> (m >>> 0) & c2 | k << i;
                    c2 = k >>> (m >>> 0) & c2;
                    g2 = 0;
                    i = l << i;
                    break;
                  }
                  if (!f) {
                    n = 0;
                    f = 0;
                    return (F(n | 0), f) | 0;
                  }
                  b[f >> 2] = a2 | 0;
                  b[f + 4 >> 2] = j | c2 & 0;
                  n = 0;
                  f = 0;
                  return (F(n | 0), f) | 0;
                }
                g2 = h - 1 | 0;
                if (g2 & h | 0) {
                  i = (D(h | 0) | 0) + 33 - (D(k | 0) | 0) | 0;
                  p2 = 64 - i | 0;
                  m = 32 - i | 0;
                  j = m >> 31;
                  o = i - 32 | 0;
                  c2 = o >> 31;
                  h = i;
                  a2 = m - 1 >> 31 & k >>> (o >>> 0) | (k << m | l >>> (i >>> 0)) & c2;
                  c2 = c2 & k >>> (i >>> 0);
                  g2 = l << p2 & j;
                  i = (k << p2 | l >>> (o >>> 0)) & j | l << m & i - 33 >> 31;
                  break;
                }
                if (f | 0) {
                  b[f >> 2] = g2 & l;
                  b[f + 4 >> 2] = 0;
                }
                if ((h | 0) == 1) {
                  o = j | c2 & 0;
                  p2 = a2 | 0 | 0;
                  return (F(o | 0), p2) | 0;
                } else {
                  p2 = $c(h | 0) | 0;
                  o = k >>> (p2 >>> 0) | 0;
                  p2 = k << 32 - p2 | l >>> (p2 >>> 0) | 0;
                  return (F(o | 0), p2) | 0;
                }
              } else {
                if (g2) {
                  if (f | 0) {
                    b[f >> 2] = (k >>> 0) % (h >>> 0);
                    b[f + 4 >> 2] = 0;
                  }
                  o = 0;
                  p2 = (k >>> 0) / (h >>> 0) >>> 0;
                  return (F(o | 0), p2) | 0;
                }
                if (!l) {
                  if (f | 0) {
                    b[f >> 2] = 0;
                    b[f + 4 >> 2] = (k >>> 0) % (i >>> 0);
                  }
                  o = 0;
                  p2 = (k >>> 0) / (i >>> 0) >>> 0;
                  return (F(o | 0), p2) | 0;
                }
                g2 = i - 1 | 0;
                if (!(g2 & i)) {
                  if (f | 0) {
                    b[f >> 2] = a2 | 0;
                    b[f + 4 >> 2] = g2 & k | c2 & 0;
                  }
                  o = 0;
                  p2 = k >>> (($c(i | 0) | 0) >>> 0);
                  return (F(o | 0), p2) | 0;
                }
                g2 = (D(i | 0) | 0) - (D(k | 0) | 0) | 0;
                if (g2 >>> 0 <= 30) {
                  c2 = g2 + 1 | 0;
                  i = 31 - g2 | 0;
                  h = c2;
                  a2 = k << i | l >>> (c2 >>> 0);
                  c2 = k >>> (c2 >>> 0);
                  g2 = 0;
                  i = l << i;
                  break;
                }
                if (!f) {
                  o = 0;
                  p2 = 0;
                  return (F(o | 0), p2) | 0;
                }
                b[f >> 2] = a2 | 0;
                b[f + 4 >> 2] = j | c2 & 0;
                o = 0;
                p2 = 0;
                return (F(o | 0), p2) | 0;
              }
            } while (0);
            if (!h) {
              k = i;
              j = 0;
              i = 0;
            } else {
              m = d2 | 0 | 0;
              l = n | e2 & 0;
              k = Zc(m | 0, l | 0, -1, -1) | 0;
              d2 = G() | 0;
              j = i;
              i = 0;
              do {
                e2 = j;
                j = g2 >>> 31 | j << 1;
                g2 = i | g2 << 1;
                e2 = a2 << 1 | e2 >>> 31 | 0;
                n = a2 >>> 31 | c2 << 1 | 0;
                _c(k | 0, d2 | 0, e2 | 0, n | 0) | 0;
                p2 = G() | 0;
                o = p2 >> 31 | ((p2 | 0) < 0 ? -1 : 0) << 1;
                i = o & 1;
                a2 = _c(e2 | 0, n | 0, o & m | 0, (((p2 | 0) < 0 ? -1 : 0) >> 31 | ((p2 | 0) < 0 ? -1 : 0) << 1) & l | 0) | 0;
                c2 = G() | 0;
                h = h - 1 | 0;
              } while ((h | 0) != 0);
              k = j;
              j = 0;
            }
            h = 0;
            if (f | 0) {
              b[f >> 2] = a2;
              b[f + 4 >> 2] = c2;
            }
            o = (g2 | 0) >>> 31 | (k | h) << 1 | (h << 1 | g2 >>> 31) & 0 | j;
            p2 = (g2 << 1 | 0 >>> 31) & -2 | i;
            return (F(o | 0), p2) | 0;
          }
          function bd(a2, c2, d2, e2) {
            a2 = a2 | 0;
            c2 = c2 | 0;
            d2 = d2 | 0;
            e2 = e2 | 0;
            var f = 0, g2 = 0;
            g2 = S;
            S = S + 16 | 0;
            f = g2 | 0;
            ad(a2, c2, d2, e2, f) | 0;
            S = g2;
            return (F(b[f + 4 >> 2] | 0), b[f >> 2] | 0) | 0;
          }
          function cd(a2, b2, c2) {
            a2 = a2 | 0;
            b2 = b2 | 0;
            c2 = c2 | 0;
            if ((c2 | 0) < 32) {
              F(b2 >>> c2 | 0);
              return a2 >>> c2 | (b2 & (1 << c2) - 1) << 32 - c2;
            }
            F(0);
            return b2 >>> c2 - 32 | 0;
          }
          function dd(a2, b2, c2) {
            a2 = a2 | 0;
            b2 = b2 | 0;
            c2 = c2 | 0;
            if ((c2 | 0) < 32) {
              F(b2 << c2 | (a2 & (1 << c2) - 1 << 32 - c2) >>> 32 - c2 | 0);
              return a2 << c2;
            }
            F(a2 << c2 - 32 | 0);
            return 0;
          }
          function ed(a2, b2) {
            a2 = +a2;
            b2 = +b2;
            if (a2 != a2) {
              return +b2;
            }
            if (b2 != b2) {
              return +a2;
            }
            return +C(+a2, +b2);
          }
          function fd(a2) {
            a2 = +a2;
            return a2 >= 0 ? +p(a2 + 0.5) : +A(a2 - 0.5);
          }
          function gd(c2, d2, e2) {
            c2 = c2 | 0;
            d2 = d2 | 0;
            e2 = e2 | 0;
            var f = 0, g2 = 0, h = 0;
            if ((e2 | 0) >= 8192) {
              K(c2 | 0, d2 | 0, e2 | 0) | 0;
              return c2 | 0;
            }
            h = c2 | 0;
            g2 = c2 + e2 | 0;
            if ((c2 & 3) == (d2 & 3)) {
              while (c2 & 3) {
                if (!e2) {
                  return h | 0;
                }
                a[c2 >> 0] = a[d2 >> 0] | 0;
                c2 = c2 + 1 | 0;
                d2 = d2 + 1 | 0;
                e2 = e2 - 1 | 0;
              }
              e2 = g2 & -4 | 0;
              f = e2 - 64 | 0;
              while ((c2 | 0) <= (f | 0)) {
                b[c2 >> 2] = b[d2 >> 2];
                b[c2 + 4 >> 2] = b[d2 + 4 >> 2];
                b[c2 + 8 >> 2] = b[d2 + 8 >> 2];
                b[c2 + 12 >> 2] = b[d2 + 12 >> 2];
                b[c2 + 16 >> 2] = b[d2 + 16 >> 2];
                b[c2 + 20 >> 2] = b[d2 + 20 >> 2];
                b[c2 + 24 >> 2] = b[d2 + 24 >> 2];
                b[c2 + 28 >> 2] = b[d2 + 28 >> 2];
                b[c2 + 32 >> 2] = b[d2 + 32 >> 2];
                b[c2 + 36 >> 2] = b[d2 + 36 >> 2];
                b[c2 + 40 >> 2] = b[d2 + 40 >> 2];
                b[c2 + 44 >> 2] = b[d2 + 44 >> 2];
                b[c2 + 48 >> 2] = b[d2 + 48 >> 2];
                b[c2 + 52 >> 2] = b[d2 + 52 >> 2];
                b[c2 + 56 >> 2] = b[d2 + 56 >> 2];
                b[c2 + 60 >> 2] = b[d2 + 60 >> 2];
                c2 = c2 + 64 | 0;
                d2 = d2 + 64 | 0;
              }
              while ((c2 | 0) < (e2 | 0)) {
                b[c2 >> 2] = b[d2 >> 2];
                c2 = c2 + 4 | 0;
                d2 = d2 + 4 | 0;
              }
            } else {
              e2 = g2 - 4 | 0;
              while ((c2 | 0) < (e2 | 0)) {
                a[c2 >> 0] = a[d2 >> 0] | 0;
                a[c2 + 1 >> 0] = a[d2 + 1 >> 0] | 0;
                a[c2 + 2 >> 0] = a[d2 + 2 >> 0] | 0;
                a[c2 + 3 >> 0] = a[d2 + 3 >> 0] | 0;
                c2 = c2 + 4 | 0;
                d2 = d2 + 4 | 0;
              }
            }
            while ((c2 | 0) < (g2 | 0)) {
              a[c2 >> 0] = a[d2 >> 0] | 0;
              c2 = c2 + 1 | 0;
              d2 = d2 + 1 | 0;
            }
            return h | 0;
          }
          function hd(c2, d2, e2) {
            c2 = c2 | 0;
            d2 = d2 | 0;
            e2 = e2 | 0;
            var f = 0, g2 = 0, h = 0, i = 0;
            h = c2 + e2 | 0;
            d2 = d2 & 255;
            if ((e2 | 0) >= 67) {
              while (c2 & 3) {
                a[c2 >> 0] = d2;
                c2 = c2 + 1 | 0;
              }
              f = h & -4 | 0;
              i = d2 | d2 << 8 | d2 << 16 | d2 << 24;
              g2 = f - 64 | 0;
              while ((c2 | 0) <= (g2 | 0)) {
                b[c2 >> 2] = i;
                b[c2 + 4 >> 2] = i;
                b[c2 + 8 >> 2] = i;
                b[c2 + 12 >> 2] = i;
                b[c2 + 16 >> 2] = i;
                b[c2 + 20 >> 2] = i;
                b[c2 + 24 >> 2] = i;
                b[c2 + 28 >> 2] = i;
                b[c2 + 32 >> 2] = i;
                b[c2 + 36 >> 2] = i;
                b[c2 + 40 >> 2] = i;
                b[c2 + 44 >> 2] = i;
                b[c2 + 48 >> 2] = i;
                b[c2 + 52 >> 2] = i;
                b[c2 + 56 >> 2] = i;
                b[c2 + 60 >> 2] = i;
                c2 = c2 + 64 | 0;
              }
              while ((c2 | 0) < (f | 0)) {
                b[c2 >> 2] = i;
                c2 = c2 + 4 | 0;
              }
            }
            while ((c2 | 0) < (h | 0)) {
              a[c2 >> 0] = d2;
              c2 = c2 + 1 | 0;
            }
            return h - e2 | 0;
          }
          function id(a2) {
            a2 = +a2;
            return a2 >= 0 ? +p(a2 + 0.5) : +A(a2 - 0.5);
          }
          function jd(a2) {
            a2 = a2 | 0;
            var c2 = 0, d2 = 0, e2 = 0;
            e2 = J() | 0;
            d2 = b[g >> 2] | 0;
            c2 = d2 + a2 | 0;
            if ((a2 | 0) > 0 & (c2 | 0) < (d2 | 0) | (c2 | 0) < 0) {
              M(c2 | 0) | 0;
              I(12);
              return -1;
            }
            if ((c2 | 0) > (e2 | 0)) {
              if (!(L(c2 | 0) | 0)) {
                I(12);
                return -1;
              }
            }
            b[g >> 2] = c2;
            return d2 | 0;
          }
          return {
            ___uremdi3: bd,
            _bitshift64Lshr: cd,
            _bitshift64Shl: dd,
            _calloc: Yc,
            _cellAreaKm2: ub,
            _cellAreaM2: vb,
            _cellAreaRads2: tb,
            _compact: Hb,
            _destroyLinkedPolygon: jc,
            _edgeLengthKm: pb,
            _edgeLengthM: qb,
            _emscripten_replace_memory: V,
            _exactEdgeLengthKm: xb,
            _exactEdgeLengthM: yb,
            _exactEdgeLengthRads: wb,
            _experimentalH3ToLocalIj: oc,
            _experimentalLocalIjToH3: pc,
            _free: Xc,
            _geoToH3: Sb,
            _getDestinationH3IndexFromUnidirectionalEdge: cc,
            _getH3IndexesFromUnidirectionalEdge: ec,
            _getH3UnidirectionalEdge: ac,
            _getH3UnidirectionalEdgeBoundary: gc,
            _getH3UnidirectionalEdgesFromHexagon: fc,
            _getOriginH3IndexFromUnidirectionalEdge: bc,
            _getPentagonIndexes: _b,
            _getRes0Indexes: va,
            _h3Distance: qc,
            _h3GetBaseCell: Ab,
            _h3GetFaces: Yb,
            _h3GetResolution: zb,
            _h3IndexesAreNeighbors: $b,
            _h3IsPentagon: Fb,
            _h3IsResClassIII: Kb,
            _h3IsValid: Bb,
            _h3Line: sc,
            _h3LineSize: rc,
            _h3SetToLinkedGeo: ka,
            _h3ToCenterChild: Gb,
            _h3ToChildren: Eb,
            _h3ToGeo: Vb,
            _h3ToGeoBoundary: Wb,
            _h3ToParent: Cb,
            _h3UnidirectionalEdgeIsValid: dc,
            _hexAreaKm2: nb,
            _hexAreaM2: ob,
            _hexRing: ea,
            _i64Subtract: _c,
            _kRing: $,
            _kRingDistances: aa,
            _llvm_minnum_f64: ed,
            _llvm_round_f64: fd,
            _malloc: Wc,
            _maxFaceCount: Xb,
            _maxH3ToChildrenSize: Db,
            _maxKringSize: _,
            _maxPolyfillSize: fa,
            _maxUncompactSize: Jb,
            _memcpy: gd,
            _memset: hd,
            _numHexagons: rb,
            _pentagonIndexCount: Zb,
            _pointDistKm: jb,
            _pointDistM: kb,
            _pointDistRads: ib,
            _polyfill: ga,
            _res0IndexCount: ua,
            _round: id,
            _sbrk: jd,
            _sizeOfCoordIJ: Ec,
            _sizeOfGeoBoundary: Ac,
            _sizeOfGeoCoord: zc,
            _sizeOfGeoPolygon: Cc,
            _sizeOfGeofence: Bc,
            _sizeOfH3Index: yc,
            _sizeOfLinkedGeoPolygon: Dc,
            _uncompact: Ib,
            establishStackSpace: Z,
            stackAlloc: W,
            stackRestore: Y,
            stackSave: X
          };
        }(
          // EMSCRIPTEN_END_ASM
          asmGlobalArg,
          asmLibraryArg,
          buffer
        )
      );
      var ___uremdi3 = Module["___uremdi3"] = asm["___uremdi3"];
      var _bitshift64Lshr = Module["_bitshift64Lshr"] = asm["_bitshift64Lshr"];
      var _bitshift64Shl = Module["_bitshift64Shl"] = asm["_bitshift64Shl"];
      var _calloc = Module["_calloc"] = asm["_calloc"];
      var _cellAreaKm2 = Module["_cellAreaKm2"] = asm["_cellAreaKm2"];
      var _cellAreaM2 = Module["_cellAreaM2"] = asm["_cellAreaM2"];
      var _cellAreaRads2 = Module["_cellAreaRads2"] = asm["_cellAreaRads2"];
      var _compact = Module["_compact"] = asm["_compact"];
      var _destroyLinkedPolygon = Module["_destroyLinkedPolygon"] = asm["_destroyLinkedPolygon"];
      var _edgeLengthKm = Module["_edgeLengthKm"] = asm["_edgeLengthKm"];
      var _edgeLengthM = Module["_edgeLengthM"] = asm["_edgeLengthM"];
      var _emscripten_replace_memory = Module["_emscripten_replace_memory"] = asm["_emscripten_replace_memory"];
      var _exactEdgeLengthKm = Module["_exactEdgeLengthKm"] = asm["_exactEdgeLengthKm"];
      var _exactEdgeLengthM = Module["_exactEdgeLengthM"] = asm["_exactEdgeLengthM"];
      var _exactEdgeLengthRads = Module["_exactEdgeLengthRads"] = asm["_exactEdgeLengthRads"];
      var _experimentalH3ToLocalIj = Module["_experimentalH3ToLocalIj"] = asm["_experimentalH3ToLocalIj"];
      var _experimentalLocalIjToH3 = Module["_experimentalLocalIjToH3"] = asm["_experimentalLocalIjToH3"];
      var _free = Module["_free"] = asm["_free"];
      var _geoToH3 = Module["_geoToH3"] = asm["_geoToH3"];
      var _getDestinationH3IndexFromUnidirectionalEdge = Module["_getDestinationH3IndexFromUnidirectionalEdge"] = asm["_getDestinationH3IndexFromUnidirectionalEdge"];
      var _getH3IndexesFromUnidirectionalEdge = Module["_getH3IndexesFromUnidirectionalEdge"] = asm["_getH3IndexesFromUnidirectionalEdge"];
      var _getH3UnidirectionalEdge = Module["_getH3UnidirectionalEdge"] = asm["_getH3UnidirectionalEdge"];
      var _getH3UnidirectionalEdgeBoundary = Module["_getH3UnidirectionalEdgeBoundary"] = asm["_getH3UnidirectionalEdgeBoundary"];
      var _getH3UnidirectionalEdgesFromHexagon = Module["_getH3UnidirectionalEdgesFromHexagon"] = asm["_getH3UnidirectionalEdgesFromHexagon"];
      var _getOriginH3IndexFromUnidirectionalEdge = Module["_getOriginH3IndexFromUnidirectionalEdge"] = asm["_getOriginH3IndexFromUnidirectionalEdge"];
      var _getPentagonIndexes = Module["_getPentagonIndexes"] = asm["_getPentagonIndexes"];
      var _getRes0Indexes = Module["_getRes0Indexes"] = asm["_getRes0Indexes"];
      var _h3Distance = Module["_h3Distance"] = asm["_h3Distance"];
      var _h3GetBaseCell = Module["_h3GetBaseCell"] = asm["_h3GetBaseCell"];
      var _h3GetFaces = Module["_h3GetFaces"] = asm["_h3GetFaces"];
      var _h3GetResolution = Module["_h3GetResolution"] = asm["_h3GetResolution"];
      var _h3IndexesAreNeighbors = Module["_h3IndexesAreNeighbors"] = asm["_h3IndexesAreNeighbors"];
      var _h3IsPentagon = Module["_h3IsPentagon"] = asm["_h3IsPentagon"];
      var _h3IsResClassIII = Module["_h3IsResClassIII"] = asm["_h3IsResClassIII"];
      var _h3IsValid = Module["_h3IsValid"] = asm["_h3IsValid"];
      var _h3Line = Module["_h3Line"] = asm["_h3Line"];
      var _h3LineSize = Module["_h3LineSize"] = asm["_h3LineSize"];
      var _h3SetToLinkedGeo = Module["_h3SetToLinkedGeo"] = asm["_h3SetToLinkedGeo"];
      var _h3ToCenterChild = Module["_h3ToCenterChild"] = asm["_h3ToCenterChild"];
      var _h3ToChildren = Module["_h3ToChildren"] = asm["_h3ToChildren"];
      var _h3ToGeo = Module["_h3ToGeo"] = asm["_h3ToGeo"];
      var _h3ToGeoBoundary = Module["_h3ToGeoBoundary"] = asm["_h3ToGeoBoundary"];
      var _h3ToParent = Module["_h3ToParent"] = asm["_h3ToParent"];
      var _h3UnidirectionalEdgeIsValid = Module["_h3UnidirectionalEdgeIsValid"] = asm["_h3UnidirectionalEdgeIsValid"];
      var _hexAreaKm2 = Module["_hexAreaKm2"] = asm["_hexAreaKm2"];
      var _hexAreaM2 = Module["_hexAreaM2"] = asm["_hexAreaM2"];
      var _hexRing = Module["_hexRing"] = asm["_hexRing"];
      var _i64Subtract = Module["_i64Subtract"] = asm["_i64Subtract"];
      var _kRing = Module["_kRing"] = asm["_kRing"];
      var _kRingDistances = Module["_kRingDistances"] = asm["_kRingDistances"];
      var _llvm_minnum_f64 = Module["_llvm_minnum_f64"] = asm["_llvm_minnum_f64"];
      var _llvm_round_f64 = Module["_llvm_round_f64"] = asm["_llvm_round_f64"];
      var _malloc = Module["_malloc"] = asm["_malloc"];
      var _maxFaceCount = Module["_maxFaceCount"] = asm["_maxFaceCount"];
      var _maxH3ToChildrenSize = Module["_maxH3ToChildrenSize"] = asm["_maxH3ToChildrenSize"];
      var _maxKringSize = Module["_maxKringSize"] = asm["_maxKringSize"];
      var _maxPolyfillSize = Module["_maxPolyfillSize"] = asm["_maxPolyfillSize"];
      var _maxUncompactSize = Module["_maxUncompactSize"] = asm["_maxUncompactSize"];
      var _memcpy = Module["_memcpy"] = asm["_memcpy"];
      var _memset = Module["_memset"] = asm["_memset"];
      var _numHexagons = Module["_numHexagons"] = asm["_numHexagons"];
      var _pentagonIndexCount = Module["_pentagonIndexCount"] = asm["_pentagonIndexCount"];
      var _pointDistKm = Module["_pointDistKm"] = asm["_pointDistKm"];
      var _pointDistM = Module["_pointDistM"] = asm["_pointDistM"];
      var _pointDistRads = Module["_pointDistRads"] = asm["_pointDistRads"];
      var _polyfill = Module["_polyfill"] = asm["_polyfill"];
      var _res0IndexCount = Module["_res0IndexCount"] = asm["_res0IndexCount"];
      var _round = Module["_round"] = asm["_round"];
      var _sbrk = Module["_sbrk"] = asm["_sbrk"];
      var _sizeOfCoordIJ = Module["_sizeOfCoordIJ"] = asm["_sizeOfCoordIJ"];
      var _sizeOfGeoBoundary = Module["_sizeOfGeoBoundary"] = asm["_sizeOfGeoBoundary"];
      var _sizeOfGeoCoord = Module["_sizeOfGeoCoord"] = asm["_sizeOfGeoCoord"];
      var _sizeOfGeoPolygon = Module["_sizeOfGeoPolygon"] = asm["_sizeOfGeoPolygon"];
      var _sizeOfGeofence = Module["_sizeOfGeofence"] = asm["_sizeOfGeofence"];
      var _sizeOfH3Index = Module["_sizeOfH3Index"] = asm["_sizeOfH3Index"];
      var _sizeOfLinkedGeoPolygon = Module["_sizeOfLinkedGeoPolygon"] = asm["_sizeOfLinkedGeoPolygon"];
      var _uncompact = Module["_uncompact"] = asm["_uncompact"];
      var establishStackSpace = Module["establishStackSpace"] = asm["establishStackSpace"];
      var stackAlloc = Module["stackAlloc"] = asm["stackAlloc"];
      var stackRestore = Module["stackRestore"] = asm["stackRestore"];
      var stackSave = Module["stackSave"] = asm["stackSave"];
      Module["asm"] = asm;
      Module["cwrap"] = cwrap;
      Module["setValue"] = setValue;
      Module["getValue"] = getValue;
      Module["getTempRet0"] = getTempRet0;
      if (memoryInitializer) {
        if (!isDataURI(memoryInitializer)) {
          memoryInitializer = locateFile(memoryInitializer);
        }
        {
          addRunDependency("memory initializer");
          var applyMemoryInitializer = function(data) {
            if (data.byteLength) {
              data = new Uint8Array(data);
            }
            HEAPU8.set(data, GLOBAL_BASE);
            if (Module["memoryInitializerRequest"]) {
              delete Module["memoryInitializerRequest"].response;
            }
            removeRunDependency("memory initializer");
          };
          var doBrowserLoad = function() {
            readAsync(memoryInitializer, applyMemoryInitializer, function() {
              throw "could not load memory initializer " + memoryInitializer;
            });
          };
          var memoryInitializerBytes = tryParseAsDataURI(memoryInitializer);
          if (memoryInitializerBytes) {
            applyMemoryInitializer(memoryInitializerBytes.buffer);
          } else if (Module["memoryInitializerRequest"]) {
            var useRequest = function() {
              var request = Module["memoryInitializerRequest"];
              var response = request.response;
              if (request.status !== 200 && request.status !== 0) {
                var data = tryParseAsDataURI(Module["memoryInitializerRequestURL"]);
                if (data) {
                  response = data.buffer;
                } else {
                  console.warn("a problem seems to have happened with Module.memoryInitializerRequest, status: " + request.status + ", retrying " + memoryInitializer);
                  doBrowserLoad();
                  return;
                }
              }
              applyMemoryInitializer(response);
            };
            if (Module["memoryInitializerRequest"].response) {
              setTimeout(useRequest, 0);
            } else {
              Module["memoryInitializerRequest"].addEventListener("load", useRequest);
            }
          } else {
            doBrowserLoad();
          }
        }
      }
      var calledRun;
      dependenciesFulfilled = function runCaller() {
        if (!calledRun) {
          run();
        }
        if (!calledRun) {
          dependenciesFulfilled = runCaller;
        }
      };
      function run(args) {
        args = args || arguments_;
        if (runDependencies > 0) {
          return;
        }
        preRun();
        if (runDependencies > 0) {
          return;
        }
        function doRun() {
          if (calledRun) {
            return;
          }
          calledRun = true;
          if (ABORT) {
            return;
          }
          initRuntime();
          preMain();
          if (Module["onRuntimeInitialized"]) {
            Module["onRuntimeInitialized"]();
          }
          postRun();
        }
        if (Module["setStatus"]) {
          Module["setStatus"]("Running...");
          setTimeout(function() {
            setTimeout(function() {
              Module["setStatus"]("");
            }, 1);
            doRun();
          }, 1);
        } else {
          doRun();
        }
      }
      Module["run"] = run;
      function abort(what) {
        if (Module["onAbort"]) {
          Module["onAbort"](what);
        }
        what += "";
        out(what);
        err(what);
        ABORT = true;
        throw "abort(" + what + "). Build with -s ASSERTIONS=1 for more info.";
      }
      Module["abort"] = abort;
      if (Module["preInit"]) {
        if (typeof Module["preInit"] == "function") {
          Module["preInit"] = [Module["preInit"]];
        }
        while (Module["preInit"].length > 0) {
          Module["preInit"].pop()();
        }
      }
      run();
      return libh32;
    }(typeof libh3 === "object" ? libh3 : {});
    NUMBER = "number";
    BOOLEAN = NUMBER;
    H3_LOWER = NUMBER;
    H3_UPPER = NUMBER;
    RESOLUTION = NUMBER;
    POINTER = NUMBER;
    BINDINGS = [
      // The size functions are inserted via build/sizes.h
      ["sizeOfH3Index", NUMBER],
      ["sizeOfGeoCoord", NUMBER],
      ["sizeOfGeoBoundary", NUMBER],
      ["sizeOfGeoPolygon", NUMBER],
      ["sizeOfGeofence", NUMBER],
      ["sizeOfLinkedGeoPolygon", NUMBER],
      ["sizeOfCoordIJ", NUMBER],
      // The remaining functions are defined in the core lib in h3Api.h
      ["h3IsValid", BOOLEAN, [H3_LOWER, H3_UPPER]],
      ["geoToH3", H3_LOWER, [NUMBER, NUMBER, RESOLUTION]],
      ["h3ToGeo", null, [H3_LOWER, H3_UPPER, POINTER]],
      ["h3ToGeoBoundary", null, [H3_LOWER, H3_UPPER, POINTER]],
      ["maxKringSize", NUMBER, [NUMBER]],
      ["kRing", null, [H3_LOWER, H3_UPPER, NUMBER, POINTER]],
      ["kRingDistances", null, [H3_LOWER, H3_UPPER, NUMBER, POINTER, POINTER]],
      ["hexRing", null, [H3_LOWER, H3_UPPER, NUMBER, POINTER]],
      ["maxPolyfillSize", NUMBER, [POINTER, RESOLUTION]],
      ["polyfill", null, [POINTER, RESOLUTION, POINTER]],
      ["h3SetToLinkedGeo", null, [POINTER, NUMBER, POINTER]],
      ["destroyLinkedPolygon", null, [POINTER]],
      ["compact", NUMBER, [POINTER, POINTER, NUMBER]],
      ["uncompact", NUMBER, [POINTER, NUMBER, POINTER, NUMBER, RESOLUTION]],
      ["maxUncompactSize", NUMBER, [POINTER, NUMBER, RESOLUTION]],
      ["h3IsPentagon", BOOLEAN, [H3_LOWER, H3_UPPER]],
      ["h3IsResClassIII", BOOLEAN, [H3_LOWER, H3_UPPER]],
      ["h3GetBaseCell", NUMBER, [H3_LOWER, H3_UPPER]],
      ["h3GetResolution", NUMBER, [H3_LOWER, H3_UPPER]],
      ["maxFaceCount", NUMBER, [H3_LOWER, H3_UPPER]],
      ["h3GetFaces", null, [H3_LOWER, H3_UPPER, POINTER]],
      ["h3ToParent", H3_LOWER, [H3_LOWER, H3_UPPER, RESOLUTION]],
      ["h3ToChildren", null, [H3_LOWER, H3_UPPER, RESOLUTION, POINTER]],
      ["h3ToCenterChild", H3_LOWER, [H3_LOWER, H3_UPPER, RESOLUTION]],
      ["maxH3ToChildrenSize", NUMBER, [H3_LOWER, H3_UPPER, RESOLUTION]],
      ["h3IndexesAreNeighbors", BOOLEAN, [H3_LOWER, H3_UPPER, H3_LOWER, H3_UPPER]],
      ["getH3UnidirectionalEdge", H3_LOWER, [H3_LOWER, H3_UPPER, H3_LOWER, H3_UPPER]],
      ["getOriginH3IndexFromUnidirectionalEdge", H3_LOWER, [H3_LOWER, H3_UPPER]],
      ["getDestinationH3IndexFromUnidirectionalEdge", H3_LOWER, [H3_LOWER, H3_UPPER]],
      ["h3UnidirectionalEdgeIsValid", BOOLEAN, [H3_LOWER, H3_UPPER]],
      ["getH3IndexesFromUnidirectionalEdge", null, [H3_LOWER, H3_UPPER, POINTER]],
      ["getH3UnidirectionalEdgesFromHexagon", null, [H3_LOWER, H3_UPPER, POINTER]],
      ["getH3UnidirectionalEdgeBoundary", null, [H3_LOWER, H3_UPPER, POINTER]],
      ["h3Distance", NUMBER, [H3_LOWER, H3_UPPER, H3_LOWER, H3_UPPER]],
      ["h3Line", NUMBER, [H3_LOWER, H3_UPPER, H3_LOWER, H3_UPPER, POINTER]],
      ["h3LineSize", NUMBER, [H3_LOWER, H3_UPPER, H3_LOWER, H3_UPPER]],
      ["experimentalH3ToLocalIj", NUMBER, [H3_LOWER, H3_UPPER, H3_LOWER, H3_UPPER, POINTER]],
      ["experimentalLocalIjToH3", NUMBER, [H3_LOWER, H3_UPPER, POINTER, POINTER]],
      ["hexAreaM2", NUMBER, [RESOLUTION]],
      ["hexAreaKm2", NUMBER, [RESOLUTION]],
      ["edgeLengthM", NUMBER, [RESOLUTION]],
      ["edgeLengthKm", NUMBER, [RESOLUTION]],
      ["pointDistM", NUMBER, [POINTER, POINTER]],
      ["pointDistKm", NUMBER, [POINTER, POINTER]],
      ["pointDistRads", NUMBER, [POINTER, POINTER]],
      ["cellAreaM2", NUMBER, [H3_LOWER, H3_UPPER]],
      ["cellAreaKm2", NUMBER, [H3_LOWER, H3_UPPER]],
      ["cellAreaRads2", NUMBER, [H3_LOWER, H3_UPPER]],
      ["exactEdgeLengthM", NUMBER, [H3_LOWER, H3_UPPER]],
      ["exactEdgeLengthKm", NUMBER, [H3_LOWER, H3_UPPER]],
      ["exactEdgeLengthRads", NUMBER, [H3_LOWER, H3_UPPER]],
      ["numHexagons", NUMBER, [RESOLUTION]],
      ["getRes0Indexes", null, [POINTER]],
      ["res0IndexCount", NUMBER],
      ["getPentagonIndexes", null, [NUMBER, POINTER]],
      ["pentagonIndexCount", NUMBER]
    ];
    H3 = {};
    BINDINGS.forEach(function bind(def) {
      H3[def[0]] = libh3.cwrap.apply(libh3, def);
    });
    BASE_16 = 16;
    SZ_INT = 4;
    SZ_PTR = 4;
    SZ_DBL = 8;
    SZ_H3INDEX = H3.sizeOfH3Index();
    SZ_GEOCOORD = H3.sizeOfGeoCoord();
    SZ_GEOBOUNDARY = H3.sizeOfGeoBoundary();
    SZ_GEOPOLYGON = H3.sizeOfGeoPolygon();
    SZ_GEOFENCE = H3.sizeOfGeofence();
    SZ_LINKED_GEOPOLYGON = H3.sizeOfLinkedGeoPolygon();
    SZ_COORDIJ = H3.sizeOfCoordIJ();
    UNITS = {
      m: "m",
      m2: "m2",
      km: "km",
      km2: "km2",
      rads: "rads",
      rads2: "rads2"
    };
    INVALID_HEXIDECIMAL_CHAR = /[^0-9a-fA-F]/;
  }
});

// node_modules/geojson2h3/dist/src/geojson2h3.js
var require_geojson2h3 = __commonJS({
  "node_modules/geojson2h3/dist/src/geojson2h3.js"(exports, module) {
    var h3 = (init_h3_js_es(), __toCommonJS(h3_js_es_exports));
    var FEATURE = "Feature";
    var FEATURE_COLLECTION = "FeatureCollection";
    var POLYGON = "Polygon";
    var MULTI_POLYGON = "MultiPolygon";
    function flatten(arrays) {
      var out = null;
      for (var i = 0; i < arrays.length; i++) {
        if (out !== null) {
          for (var j = 0; j < arrays[i].length; j++) {
            out.push(arrays[i][j]);
          }
        } else {
          out = arrays[i];
        }
      }
      return Array.from(new Set(out));
    }
    function centroid(polygon) {
      var lngSum = 0;
      var latSum = 0;
      var count = 0;
      var loop = polygon[0];
      for (var i = 0; i < loop.length; i++) {
        lngSum += loop[i][0];
        latSum += loop[i][1];
        count++;
      }
      return [lngSum / count, latSum / count];
    }
    function featureCollectionToH3Set(featureCollection, resolution) {
      var features = featureCollection.features;
      if (!features) {
        throw new Error("No features found");
      }
      return flatten(features.map(function(feature) {
        return featureToH3Set(feature, resolution);
      }));
    }
    function featureToH3Set(feature, resolution, options) {
      if (options === void 0) options = {};
      var type = feature.type;
      var geometry = feature.geometry;
      var geometryType = geometry && geometry.type;
      if (type === FEATURE_COLLECTION) {
        return featureCollectionToH3Set(feature, resolution);
      }
      if (type !== FEATURE) {
        throw new Error("Unhandled type: " + type);
      }
      if (geometryType !== POLYGON && geometryType !== MULTI_POLYGON) {
        throw new Error("Unhandled geometry type: " + geometryType);
      }
      var polygons = geometryType === POLYGON ? [geometry.coordinates] : geometry.coordinates;
      return flatten(
        polygons.map(function(polygon) {
          var result = h3.polyfill(polygon, resolution, true);
          if (result.length || !options.ensureOutput) {
            return result;
          }
          var ref = centroid(polygon);
          var lng = ref[0];
          var lat = ref[1];
          return [h3.geoToH3(lat, lng, resolution)];
        })
      );
    }
    function h3ToFeature(h3Index, properties) {
      if (properties === void 0) properties = {};
      var coordinates = [h3.h3ToGeoBoundary(h3Index, true)];
      return {
        type: FEATURE,
        id: h3Index,
        properties,
        geometry: {
          type: POLYGON,
          coordinates
        }
      };
    }
    function h3SetToFeature(hexagons, properties) {
      if (properties === void 0) properties = {};
      var polygons = h3.h3SetToMultiPolygon(hexagons, true);
      var isMultiPolygon = polygons.length > 1;
      var type = isMultiPolygon ? MULTI_POLYGON : POLYGON;
      var coordinates = isMultiPolygon ? polygons : polygons[0] || [];
      return {
        type: FEATURE,
        properties,
        geometry: {
          type,
          coordinates
        }
      };
    }
    function h3SetToMultiPolygonFeature(hexagons, properties) {
      if (properties === void 0) properties = {};
      var coordinates = hexagons.map(
        function(h3Index) {
          return [h3.h3ToGeoBoundary(h3Index, { geoJson: true })];
        }
      );
      return {
        type: FEATURE,
        properties,
        geometry: {
          type: MULTI_POLYGON,
          coordinates
        }
      };
    }
    function h3SetToFeatureCollection(hexagons, getProperties) {
      var features = [];
      for (var i = 0; i < hexagons.length; i++) {
        var h3Index = hexagons[i];
        var properties = getProperties ? getProperties(h3Index) : {};
        features.push(h3ToFeature(h3Index, properties));
      }
      return {
        type: FEATURE_COLLECTION,
        features
      };
    }
    module.exports = {
      featureToH3Set,
      h3ToFeature,
      h3SetToFeature,
      h3SetToMultiPolygonFeature,
      h3SetToFeatureCollection
    };
  }
});

// node_modules/geojson2h3/index.js
var require_geojson2h32 = __commonJS({
  "node_modules/geojson2h3/index.js"(exports, module) {
    module.exports = require_geojson2h3();
  }
});
export default require_geojson2h32();
//# sourceMappingURL=geojson2h3.js.map
