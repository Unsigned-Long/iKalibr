/*M///////////////////////////////////////////////////////////////////////////////////////
 //
 //  IMPORTANT: READ BEFORE DOWNLOADING, COPYING, INSTALLING OR USING.
 //
 //  By downloading, copying, installing or using the software you agree to this license.
 //  If you do not agree to this license, do not download, install,
 //  copy or use the software.
 //
 //
 //                           License Agreement
 //                For Open Source Computer Vision Library
 //
 // Copyright (C) 2000-2008, Intel Corporation, all rights reserved.
 // Copyright (C) 2009, Willow Garage Inc., all rights reserved.
 // Third party copyrights are property of their respective owners.
 //
 // Redistribution and use in source and binary forms, with or without modification,
 // are permitted provided that the following conditions are met:
 //
 //   * Redistribution's of source code must retain the above copyright notice,
 //     this list of conditions and the following disclaimer.
 //
 //   * Redistribution's in binary form must reproduce the above copyright notice,
 //     this list of conditions and the following disclaimer in the documentation
 //     and/or other materials provided with the distribution.
 //
 //   * The name of the copyright holders may not be used to endorse or promote products
 //     derived from this software without specific prior written permission.
 //
 // This software is provided by the copyright holders and contributors "as is" and
 // any express or implied warranties, including, but not limited to, the implied
 // warranties of merchantability and fitness for a particular purpose are disclaimed.
 // In no event shall the Intel Corporation or contributors be liable for any direct,
 // indirect, incidental, special, exemplary, or consequential damages
 // (including, but not limited to, procurement of substitute goods or services;
 // loss of use, data, or profits; or business interruption) however caused
 // and on any theory of liability, whether in contract, strict liability,
 // or tort (including negligence or otherwise) arising in any way out of
 // the use of this software, even if advised of the possibility of such damage.
 //
 //M*/

// #include "precomp.hpp"
#include "util/circlesgrid.h"
#include <limits>
#include "opencv2/calib3d.hpp"

#include <iostream>
#include <opencv2/imgproc.hpp>

// Requires CMake flag: DEBUG_opencv_calib3d=ON
// #define DEBUG_CIRCLES

#ifdef DEBUG_CIRCLES
    #include <iostream>
    #include "opencv2/opencv_modules.hpp"
    #ifdef HAVE_OPENCV_HIGHGUI
        #include "opencv2/highgui.hpp"
    #else
        #undef DEBUG_CIRCLES
    #endif
#endif

namespace ns_cv_helper {

using namespace cv;

#ifdef DEBUG_CIRCLES
void drawPoints(const std::vector<Point2f> &points,
                Mat &outImage,
                int radius = 2,
                Scalar color = Scalar::all(255),
                int thickness = -1) {
    for (size_t i = 0; i < points.size(); i++) {
        circle(outImage, points[i], radius, color, thickness);
    }
}
#endif

Graph::Graph(size_t n) {
    for (size_t i = 0; i < n; i++) {
        addVertex(i);
    }
}

bool Graph::doesVertexExist(size_t id) const { return (vertices.find(id) != vertices.end()); }

void Graph::addVertex(size_t id) {
    CV_Assert(!doesVertexExist(id));

    vertices.emplace(id, Vertex());
}

void Graph::addEdge(size_t id1, size_t id2) {
    CV_Assert(doesVertexExist(id1));
    CV_Assert(doesVertexExist(id2));

    vertices[id1].neighbors.insert(id2);
    vertices[id2].neighbors.insert(id1);
}

void Graph::removeEdge(size_t id1, size_t id2) {
    CV_Assert(doesVertexExist(id1));
    CV_Assert(doesVertexExist(id2));

    vertices[id1].neighbors.erase(id2);
    vertices[id2].neighbors.erase(id1);
}

bool Graph::areVerticesAdjacent(size_t id1, size_t id2) const {
    Vertices::const_iterator it = vertices.find(id1);
    CV_Assert(it != vertices.end());
    const Neighbors &neighbors = it->second.neighbors;
    return neighbors.find(id2) != neighbors.end();
}

size_t Graph::getVerticesCount() const { return vertices.size(); }

size_t Graph::getDegree(size_t id) const {
    Vertices::const_iterator it = vertices.find(id);
    CV_Assert(it != vertices.end());
    return it->second.neighbors.size();
}

void Graph::floydWarshall(cv::Mat &distanceMatrix, int infinity) const {
    const int edgeWeight = 1;

    const int n = (int)getVerticesCount();
    distanceMatrix.create(n, n, CV_32SC1);
    distanceMatrix.setTo(infinity);
    for (Vertices::const_iterator it1 = vertices.begin(); it1 != vertices.end(); ++it1) {
        distanceMatrix.at<int>((int)it1->first, (int)it1->first) = 0;
        for (Neighbors::const_iterator it2 = it1->second.neighbors.begin();
             it2 != it1->second.neighbors.end(); ++it2) {
            CV_Assert(it1->first != *it2);
            distanceMatrix.at<int>((int)it1->first, (int)*it2) = edgeWeight;
        }
    }

    for (Vertices::const_iterator it1 = vertices.begin(); it1 != vertices.end(); ++it1) {
        for (Vertices::const_iterator it2 = vertices.begin(); it2 != vertices.end(); ++it2) {
            for (Vertices::const_iterator it3 = vertices.begin(); it3 != vertices.end(); ++it3) {
                int i1 = (int)it1->first, i2 = (int)it2->first, i3 = (int)it3->first;
                int val1 = distanceMatrix.at<int>(i2, i3);
                int val2;
                if (distanceMatrix.at<int>(i2, i1) == infinity ||
                    distanceMatrix.at<int>(i1, i3) == infinity)
                    val2 = val1;
                else {
                    val2 = distanceMatrix.at<int>(i2, i1) + distanceMatrix.at<int>(i1, i3);
                }
                distanceMatrix.at<int>(i2, i3) = (val1 == infinity) ? val2 : std::min(val1, val2);
            }
        }
    }
}

const Graph::Neighbors &Graph::getNeighbors(size_t id) const {
    Vertices::const_iterator it = vertices.find(id);
    CV_Assert(it != vertices.end());
    return it->second.neighbors;
}

CirclesGridFinder::Segment::Segment(cv::Point2f _s, cv::Point2f _e)
    : s(_s),
      e(_e) {}

void computeShortestPath(Mat &predecessorMatrix, int v1, int v2, std::vector<int> &path);
void computePredecessorMatrix(const Mat &dm, int verticesCount, Mat &predecessorMatrix);

CirclesGridFinderParameters::CirclesGridFinderParameters() {
    minDensity = 10;
    densityNeighborhoodSize = Size2f(16, 16);
    minDistanceToAddKeypoint = 20;
    kmeansAttempts = 100;
    convexHullFactor = 1.1f;
    keypointScale = 1;

    minGraphConfidence = 9;
    vertexGain = 1;
    vertexPenalty = -0.6f;
    edgeGain = 1;
    edgePenalty = -0.6f;
    existingVertexGain = 10000;

    minRNGEdgeSwitchDist = 5.f;
    gridType = SYMMETRIC_GRID;

    squareSize = 1.0f;
    maxRectifiedDistance = squareSize / 2.0f;
}

CirclesGridFinder::CirclesGridFinder(Size _patternSize,
                                     const std::vector<Point2f> &testKeypoints,
                                     const CirclesGridFinderParameters &_parameters)
    : patternSize(static_cast<size_t>(_patternSize.width),
                  static_cast<size_t>(_patternSize.height)) {
    CV_Assert(_patternSize.height >= 0 && _patternSize.width >= 0);

    keypoints = testKeypoints;
    parameters = _parameters;
    largeHoles = 0;
    smallHoles = 0;
}

bool CirclesGridFinder::findHoles() {
    switch (parameters.gridType) {
        case CirclesGridFinderParameters::SYMMETRIC_GRID: {
            std::vector<Point2f> vectors, filteredVectors, basis;
            Graph rng(0);
            computeRNG(rng, vectors);
            filterOutliersByDensity(vectors, filteredVectors);
            std::vector<Graph> basisGraphs;
            findBasis(filteredVectors, basis, basisGraphs);
            findMCS(basis, basisGraphs);
            break;
        }

        case CirclesGridFinderParameters::ASYMMETRIC_GRID: {
            std::vector<Point2f> vectors, tmpVectors, filteredVectors, basis;
            Graph rng(0);
            computeRNG(rng, tmpVectors);
            rng2gridGraph(rng, vectors);
            filterOutliersByDensity(vectors, filteredVectors);
            std::vector<Graph> basisGraphs;
            findBasis(filteredVectors, basis, basisGraphs);
            findMCS(basis, basisGraphs);
            eraseUsedGraph(basisGraphs);
            holes2 = holes;
            holes.clear();
            findMCS(basis, basisGraphs);
            break;
        }

        default:
            CV_Error(Error::StsBadArg, "Unknown pattern type");
    }
    return (isDetectionCorrect());
    // CV_Error( 0, "Detection is not correct" );
}

void CirclesGridFinder::rng2gridGraph(Graph &rng, std::vector<cv::Point2f> &vectors) const {
    for (size_t i = 0; i < rng.getVerticesCount(); i++) {
        Graph::Neighbors neighbors1 = rng.getNeighbors(i);
        for (Graph::Neighbors::iterator it1 = neighbors1.begin(); it1 != neighbors1.end(); ++it1) {
            Graph::Neighbors neighbors2 = rng.getNeighbors(*it1);
            for (Graph::Neighbors::iterator it2 = neighbors2.begin(); it2 != neighbors2.end();
                 ++it2) {
                if (i < *it2) {
                    Point2f vec1 = keypoints[i] - keypoints[*it1];
                    Point2f vec2 = keypoints[*it1] - keypoints[*it2];
                    if (norm(vec1 - vec2) < parameters.minRNGEdgeSwitchDist ||
                        norm(vec1 + vec2) < parameters.minRNGEdgeSwitchDist)
                        continue;

                    vectors.push_back(keypoints[i] - keypoints[*it2]);
                    vectors.push_back(keypoints[*it2] - keypoints[i]);
                }
            }
        }
    }
}

void CirclesGridFinder::eraseUsedGraph(std::vector<Graph> &basisGraphs) const {
    for (size_t i = 0; i < holes.size(); i++) {
        for (size_t j = 0; j < holes[i].size(); j++) {
            for (size_t k = 0; k < basisGraphs.size(); k++) {
                if (i != holes.size() - 1 &&
                    basisGraphs[k].areVerticesAdjacent(holes[i][j], holes[i + 1][j])) {
                    basisGraphs[k].removeEdge(holes[i][j], holes[i + 1][j]);
                }

                if (j != holes[i].size() - 1 &&
                    basisGraphs[k].areVerticesAdjacent(holes[i][j], holes[i][j + 1])) {
                    basisGraphs[k].removeEdge(holes[i][j], holes[i][j + 1]);
                }
            }
        }
    }
}

bool CirclesGridFinder::isDetectionCorrect() {
    switch (parameters.gridType) {
        case CirclesGridFinderParameters::SYMMETRIC_GRID: {
            if (holes.size() != patternSize.height) return false;

            std::set<size_t> vertices;
            for (size_t i = 0; i < holes.size(); i++) {
                if (holes[i].size() != patternSize.width) return false;

                for (size_t j = 0; j < holes[i].size(); j++) {
                    vertices.insert(holes[i][j]);
                }
            }

            return vertices.size() == patternSize.area();
        }

        case CirclesGridFinderParameters::ASYMMETRIC_GRID: {
            if (holes.size() < holes2.size() || holes[0].size() < holes2[0].size()) {
                largeHoles = &holes2;
                smallHoles = &holes;
            } else {
                largeHoles = &holes;
                smallHoles = &holes2;
            }

            size_t largeWidth = patternSize.width;
            size_t largeHeight = (size_t)ceil(patternSize.height / 2.);
            size_t smallWidth = patternSize.width;
            size_t smallHeight = (size_t)floor(patternSize.height / 2.);

            size_t sw = smallWidth, sh = smallHeight, lw = largeWidth, lh = largeHeight;
            if (largeHoles->size() != largeHeight) {
                std::swap(lh, lw);
            }
            if (smallHoles->size() != smallHeight) {
                std::swap(sh, sw);
            }

            if (largeHoles->size() != lh || smallHoles->size() != sh) {
                return false;
            }

            std::set<size_t> vertices;
            for (size_t i = 0; i < largeHoles->size(); i++) {
                if (largeHoles->at(i).size() != lw) {
                    return false;
                }

                for (size_t j = 0; j < largeHoles->at(i).size(); j++) {
                    vertices.insert(largeHoles->at(i)[j]);
                }

                if (i < smallHoles->size()) {
                    if (smallHoles->at(i).size() != sw) {
                        return false;
                    }

                    for (size_t j = 0; j < smallHoles->at(i).size(); j++) {
                        vertices.insert(smallHoles->at(i)[j]);
                    }
                }
            }
            return (vertices.size() == largeHeight * largeWidth + smallHeight * smallWidth);
        }
    }
    CV_Error(Error::StsBadArg, "Unknown pattern type");
}

void CirclesGridFinder::findMCS(const std::vector<Point2f> &basis,
                                std::vector<Graph> &basisGraphs) {
    holes.clear();
    Path longestPath;
    size_t bestGraphIdx = findLongestPath(basisGraphs, longestPath);
    std::vector<size_t> holesRow = longestPath.vertices;

    while (holesRow.size() > std::max(patternSize.width, patternSize.height)) {
        holesRow.pop_back();
        holesRow.erase(holesRow.begin());
    }

    if (bestGraphIdx == 0) {
        holes.push_back(holesRow);
        size_t w = holes[0].size();
        size_t h = holes.size();

        // parameters.minGraphConfidence = holes[0].size() * parameters.vertexGain +
        // (holes[0].size() - 1) * parameters.edgeGain; parameters.minGraphConfidence =
        // holes[0].size() * parameters.vertexGain + (holes[0].size() / 2) * parameters.edgeGain;
        // parameters.minGraphConfidence = holes[0].size() * parameters.existingVertexGain +
        // (holes[0].size() / 2) * parameters.edgeGain;
        parameters.minGraphConfidence = holes[0].size() * parameters.existingVertexGain;
        for (size_t i = h; i < patternSize.height; i++) {
            addHolesByGraph(basisGraphs, true, basis[1]);
        }

        // parameters.minGraphConfidence = holes.size() * parameters.existingVertexGain +
        // (holes.size() / 2) * parameters.edgeGain;
        parameters.minGraphConfidence = holes.size() * parameters.existingVertexGain;

        for (size_t i = w; i < patternSize.width; i++) {
            addHolesByGraph(basisGraphs, false, basis[0]);
        }
    } else {
        holes.resize(holesRow.size());
        for (size_t i = 0; i < holesRow.size(); i++) holes[i].push_back(holesRow[i]);

        size_t w = holes[0].size();
        size_t h = holes.size();

        parameters.minGraphConfidence = holes.size() * parameters.existingVertexGain;
        for (size_t i = w; i < patternSize.width; i++) {
            addHolesByGraph(basisGraphs, false, basis[0]);
        }

        parameters.minGraphConfidence = holes[0].size() * parameters.existingVertexGain;
        for (size_t i = h; i < patternSize.height; i++) {
            addHolesByGraph(basisGraphs, true, basis[1]);
        }
    }
}

Mat CirclesGridFinder::rectifyGrid(Size detectedGridSize,
                                   const std::vector<Point2f> &centers,
                                   const std::vector<Point2f> &keypoints,
                                   std::vector<Point2f> &warpedKeypoints) {
    CV_Assert(!centers.empty());
    const float edgeLength = 30;
    const Point2f offset(150, 150);

    std::vector<Point2f> dstPoints;
    bool isClockwiseBefore = getDirection(centers[0], centers[detectedGridSize.width - 1],
                                          centers[centers.size() - 1]) < 0;

    int iStart = isClockwiseBefore ? 0 : detectedGridSize.height - 1;
    int iEnd = isClockwiseBefore ? detectedGridSize.height : -1;
    int iStep = isClockwiseBefore ? 1 : -1;
    for (int i = iStart; i != iEnd; i += iStep) {
        for (int j = 0; j < detectedGridSize.width; j++) {
            dstPoints.push_back(offset + Point2f(edgeLength * j, edgeLength * i));
        }
    }

    Mat H = findHomography(centers, dstPoints, RANSAC);
    // Mat H = findHomography(corners, dstPoints);

    if (H.empty()) {
        H = Mat::zeros(3, 3, CV_64FC1);
        warpedKeypoints.clear();
        return H;
    }

    std::vector<Point2f> srcKeypoints;
    for (size_t i = 0; i < keypoints.size(); i++) {
        srcKeypoints.push_back(keypoints[i]);
    }

    Mat dstKeypointsMat;
    transform(srcKeypoints, dstKeypointsMat, H);
    std::vector<Point2f> dstKeypoints;
    convertPointsFromHomogeneous(dstKeypointsMat, dstKeypoints);

    warpedKeypoints.clear();
    for (auto &pt : dstKeypoints) {
        warpedKeypoints.emplace_back(std::move(pt));
    }

    return H;
}

size_t CirclesGridFinder::findNearestKeypoint(Point2f pt) const {
    size_t bestIdx = 0;
    double minDist = std::numeric_limits<double>::max();
    for (size_t i = 0; i < keypoints.size(); i++) {
        double dist = norm(pt - keypoints[i]);
        if (dist < minDist) {
            minDist = dist;
            bestIdx = i;
        }
    }
    return bestIdx;
}

void CirclesGridFinder::addPoint(Point2f pt, std::vector<size_t> &points) {
    size_t ptIdx = findNearestKeypoint(pt);
    if (norm(keypoints[ptIdx] - pt) > parameters.minDistanceToAddKeypoint) {
        Point2f kpt = Point2f(pt);
        keypoints.push_back(kpt);
        points.push_back(keypoints.size() - 1);
    } else {
        points.push_back(ptIdx);
    }
}

void CirclesGridFinder::findCandidateLine(std::vector<size_t> &line,
                                          size_t seedLineIdx,
                                          bool addRow,
                                          Point2f basisVec,
                                          std::vector<size_t> &seeds) {
    line.clear();
    seeds.clear();

    if (addRow) {
        for (size_t i = 0; i < holes[seedLineIdx].size(); i++) {
            Point2f pt = keypoints[holes[seedLineIdx][i]] + basisVec;
            addPoint(pt, line);
            seeds.push_back(holes[seedLineIdx][i]);
        }
    } else {
        for (size_t i = 0; i < holes.size(); i++) {
            Point2f pt = keypoints[holes[i][seedLineIdx]] + basisVec;
            addPoint(pt, line);
            seeds.push_back(holes[i][seedLineIdx]);
        }
    }

    CV_Assert(line.size() == seeds.size());
}

void CirclesGridFinder::findCandidateHoles(std::vector<size_t> &above,
                                           std::vector<size_t> &below,
                                           bool addRow,
                                           Point2f basisVec,
                                           std::vector<size_t> &aboveSeeds,
                                           std::vector<size_t> &belowSeeds) {
    above.clear();
    below.clear();
    aboveSeeds.clear();
    belowSeeds.clear();

    findCandidateLine(above, 0, addRow, -basisVec, aboveSeeds);
    size_t lastIdx = addRow ? holes.size() - 1 : holes[0].size() - 1;
    findCandidateLine(below, lastIdx, addRow, basisVec, belowSeeds);

    CV_Assert(below.size() == above.size());
    CV_Assert(belowSeeds.size() == aboveSeeds.size());
    CV_Assert(below.size() == belowSeeds.size());
}

bool CirclesGridFinder::areCentersNew(const std::vector<size_t> &newCenters,
                                      const std::vector<std::vector<size_t>> &holes) {
    for (size_t i = 0; i < newCenters.size(); i++) {
        for (size_t j = 0; j < holes.size(); j++) {
            if (holes[j].end() != std::find(holes[j].begin(), holes[j].end(), newCenters[i])) {
                return false;
            }
        }
    }

    return true;
}

void CirclesGridFinder::insertWinner(float aboveConfidence,
                                     float belowConfidence,
                                     float minConfidence,
                                     bool addRow,
                                     const std::vector<size_t> &above,
                                     const std::vector<size_t> &below,
                                     std::vector<std::vector<size_t>> &holes) {
    if (aboveConfidence < minConfidence && belowConfidence < minConfidence) return;

    if (addRow) {
        if (aboveConfidence >= belowConfidence) {
            if (!areCentersNew(above, holes)) CV_Error(0, "Centers are not new");

            holes.insert(holes.begin(), above);
        } else {
            if (!areCentersNew(below, holes)) CV_Error(0, "Centers are not new");

            holes.insert(holes.end(), below);
        }
    } else {
        if (aboveConfidence >= belowConfidence) {
            if (!areCentersNew(above, holes)) CV_Error(0, "Centers are not new");

            for (size_t i = 0; i < holes.size(); i++) {
                holes[i].insert(holes[i].begin(), above[i]);
            }
        } else {
            if (!areCentersNew(below, holes)) CV_Error(0, "Centers are not new");

            for (size_t i = 0; i < holes.size(); i++) {
                holes[i].insert(holes[i].end(), below[i]);
            }
        }
    }
}

float CirclesGridFinder::computeGraphConfidence(const std::vector<Graph> &basisGraphs,
                                                bool addRow,
                                                const std::vector<size_t> &points,
                                                const std::vector<size_t> &seeds) {
    CV_Assert(points.size() == seeds.size());
    float confidence = 0;
    const size_t vCount = basisGraphs[0].getVerticesCount();
    CV_Assert(basisGraphs[0].getVerticesCount() == basisGraphs[1].getVerticesCount());

    for (size_t i = 0; i < seeds.size(); i++) {
        if (seeds[i] < vCount && points[i] < vCount) {
            if (!basisGraphs[addRow].areVerticesAdjacent(seeds[i], points[i])) {
                confidence += parameters.vertexPenalty;
            } else {
                confidence += parameters.vertexGain;
            }
        }

        if (points[i] < vCount) {
            confidence += parameters.existingVertexGain;
        }
    }

    for (size_t i = 1; i < points.size(); i++) {
        if (points[i - 1] < vCount && points[i] < vCount) {
            if (!basisGraphs[!addRow].areVerticesAdjacent(points[i - 1], points[i])) {
                confidence += parameters.edgePenalty;
            } else {
                confidence += parameters.edgeGain;
            }
        }
    }
    return confidence;
}

void CirclesGridFinder::addHolesByGraph(const std::vector<Graph> &basisGraphs,
                                        bool addRow,
                                        Point2f basisVec) {
    std::vector<size_t> above, below, aboveSeeds, belowSeeds;
    findCandidateHoles(above, below, addRow, basisVec, aboveSeeds, belowSeeds);
    float aboveConfidence = computeGraphConfidence(basisGraphs, addRow, above, aboveSeeds);
    float belowConfidence = computeGraphConfidence(basisGraphs, addRow, below, belowSeeds);

    insertWinner(aboveConfidence, belowConfidence, parameters.minGraphConfidence, addRow, above,
                 below, holes);
}

void CirclesGridFinder::filterOutliersByDensity(const std::vector<Point2f> &samples,
                                                std::vector<Point2f> &filteredSamples) {
    if (samples.empty()) CV_Error(0, "samples is empty");

    filteredSamples.clear();

    for (size_t i = 0; i < samples.size(); i++) {
        Rect_<float> rect(samples[i] - Point2f(parameters.densityNeighborhoodSize) * 0.5,
                          parameters.densityNeighborhoodSize);
        int neighborsCount = 0;
        for (size_t j = 0; j < samples.size(); j++) {
            if (rect.contains(samples[j])) neighborsCount++;
        }
        if (neighborsCount >= parameters.minDensity) filteredSamples.push_back(samples[i]);
    }

    if (filteredSamples.empty()) CV_Error(0, "filteredSamples is empty");
}

void CirclesGridFinder::findBasis(const std::vector<Point2f> &samples,
                                  std::vector<Point2f> &basis,
                                  std::vector<Graph> &basisGraphs) {
    basis.clear();
    Mat bestLabels;
    TermCriteria termCriteria;
    Mat centers;
    const int clustersCount = 4;
    kmeans(Mat(samples).reshape(1, 0), clustersCount, bestLabels, termCriteria,
           parameters.kmeansAttempts, KMEANS_RANDOM_CENTERS, centers);
    CV_Assert(centers.type() == CV_32FC1);

    std::vector<int> basisIndices;
    // TODO: only remove duplicate
    for (int i = 0; i < clustersCount; i++) {
        int maxIdx = (fabs(centers.at<float>(i, 0)) < fabs(centers.at<float>(i, 1)));
        if (centers.at<float>(i, maxIdx) > 0) {
            Point2f vec(centers.at<float>(i, 0), centers.at<float>(i, 1));
            basis.push_back(vec);
            basisIndices.push_back(i);
        }
    }
    if (basis.size() != 2) CV_Error(0, "Basis size is not 2");

    if (basis[1].x > basis[0].x) {
        std::swap(basis[0], basis[1]);
        std::swap(basisIndices[0], basisIndices[1]);
    }

    const float minBasisDif = 2;
    if (norm(basis[0] - basis[1]) < minBasisDif) CV_Error(0, "degenerate basis");

    std::vector<std::vector<Point2f>> clusters(2), hulls(2);
    for (int k = 0; k < (int)samples.size(); k++) {
        int label = bestLabels.at<int>(k, 0);
        int idx = -1;
        if (label == basisIndices[0]) idx = 0;
        if (label == basisIndices[1]) idx = 1;
        if (idx >= 0) {
            clusters[idx].push_back(basis[idx] +
                                    parameters.convexHullFactor * (samples[k] - basis[idx]));
        }
    }
    for (size_t i = 0; i < basis.size(); i++) {
        convexHull(clusters[i], hulls[i]);
    }

    basisGraphs.resize(basis.size(), Graph(keypoints.size()));
    for (size_t i = 0; i < keypoints.size(); i++) {
        for (size_t j = 0; j < keypoints.size(); j++) {
            if (i == j) continue;

            Point2f vec = keypoints[i] - keypoints[j];

            for (size_t k = 0; k < hulls.size(); k++) {
                if (pointPolygonTest(hulls[k], vec, false) >= 0) {
                    basisGraphs[k].addEdge(i, j);
                }
            }
        }
    }
    if (basisGraphs.size() != 2) CV_Error(0, "Number of basis graphs is not 2");
}

void CirclesGridFinder::computeRNG(Graph &rng,
                                   std::vector<cv::Point2f> &vectors,
                                   Mat *drawImage) const {
    rng = Graph(keypoints.size());
    vectors.clear();

    // TODO: use more fast algorithm instead of naive N^3
    for (size_t i = 0; i < keypoints.size(); i++) {
        for (size_t j = 0; j < keypoints.size(); j++) {
            if (i == j) continue;

            Point2f vec = keypoints[i] - keypoints[j];
            double dist = norm(vec);

            bool isNeighbors = true;
            for (size_t k = 0; k < keypoints.size(); k++) {
                if (k == i || k == j) continue;

                double dist1 = norm(keypoints[i] - keypoints[k]);
                double dist2 = norm(keypoints[j] - keypoints[k]);
                if (dist1 < dist && dist2 < dist) {
                    isNeighbors = false;
                    break;
                }
            }

            if (isNeighbors) {
                rng.addEdge(i, j);
                vectors.push_back(keypoints[i] - keypoints[j]);
                if (drawImage != 0) {
                    line(*drawImage, keypoints[i], keypoints[j], Scalar(255, 0, 0), 2);
                    circle(*drawImage, keypoints[i], 3, Scalar(0, 0, 255), -1);
                    circle(*drawImage, keypoints[j], 3, Scalar(0, 0, 255), -1);
                }
            }
        }
    }
}

void computePredecessorMatrix(const Mat &dm, int verticesCount, Mat &predecessorMatrix) {
    CV_Assert(dm.type() == CV_32SC1);
    predecessorMatrix.create(verticesCount, verticesCount, CV_32SC1);
    predecessorMatrix = -1;
    for (int i = 0; i < predecessorMatrix.rows; i++) {
        for (int j = 0; j < predecessorMatrix.cols; j++) {
            int dist = dm.at<int>(i, j);
            for (int k = 0; k < verticesCount; k++) {
                if (dm.at<int>(i, k) == dist - 1 && dm.at<int>(k, j) == 1) {
                    predecessorMatrix.at<int>(i, j) = k;
                    break;
                }
            }
        }
    }
}

static void computeShortestPath(Mat &predecessorMatrix,
                                size_t v1,
                                size_t v2,
                                std::vector<size_t> &path) {
    if (predecessorMatrix.at<int>((int)v1, (int)v2) < 0) {
        path.push_back(v1);
        return;
    }

    computeShortestPath(predecessorMatrix, v1, predecessorMatrix.at<int>((int)v1, (int)v2), path);
    path.push_back(v2);
}

size_t CirclesGridFinder::findLongestPath(std::vector<Graph> &basisGraphs, Path &bestPath) {
    std::vector<Path> longestPaths(1);
    std::vector<int> confidences;

    size_t bestGraphIdx = 0;
    const int infinity = -1;
    for (size_t graphIdx = 0; graphIdx < basisGraphs.size(); graphIdx++) {
        const Graph &g = basisGraphs[graphIdx];
        Mat distanceMatrix;
        g.floydWarshall(distanceMatrix, infinity);
        Mat predecessorMatrix;
        computePredecessorMatrix(distanceMatrix, (int)g.getVerticesCount(), predecessorMatrix);

        double maxVal;
        Point maxLoc;
        minMaxLoc(distanceMatrix, 0, &maxVal, 0, &maxLoc);

        if (maxVal > longestPaths[0].length) {
            longestPaths.clear();
            confidences.clear();
            bestGraphIdx = graphIdx;
        }
        if (longestPaths.empty() ||
            (maxVal == longestPaths[0].length && graphIdx == bestGraphIdx)) {
            Path path = Path(maxLoc.x, maxLoc.y, cvRound(maxVal));
            CV_Assert(maxLoc.x >= 0 && maxLoc.y >= 0);
            size_t id1 = static_cast<size_t>(maxLoc.x);
            size_t id2 = static_cast<size_t>(maxLoc.y);
            computeShortestPath(predecessorMatrix, id1, id2, path.vertices);
            longestPaths.push_back(path);

            int conf = 0;
            for (int v2 = 0; v2 < (int)path.vertices.size(); v2++) {
                conf += (int)basisGraphs[1 - (int)graphIdx].getDegree(v2);
            }
            confidences.push_back(conf);
        }
    }
    // if( bestGraphIdx != 0 )
    // CV_Error( 0, "" );

    int maxConf = -1;
    int bestPathIdx = -1;
    for (int i = 0; i < (int)confidences.size(); i++) {
        if (confidences[i] > maxConf) {
            maxConf = confidences[i];
            bestPathIdx = i;
        }
    }

    // int bestPathIdx = rand() % longestPaths.size();
    bestPath = longestPaths.at(bestPathIdx);
    bool needReverse =
        (bestGraphIdx == 0 &&
         keypoints[bestPath.lastVertex].x < keypoints[bestPath.firstVertex].x) ||
        (bestGraphIdx == 1 && keypoints[bestPath.lastVertex].y < keypoints[bestPath.firstVertex].y);
    if (needReverse) {
        std::swap(bestPath.lastVertex, bestPath.firstVertex);
        std::reverse(bestPath.vertices.begin(), bestPath.vertices.end());
    }
    return bestGraphIdx;
}

void CirclesGridFinder::drawBasis(const std::vector<Point2f> &basis,
                                  Point2f origin,
                                  Mat &drawImg) const {
    for (size_t i = 0; i < basis.size(); i++) {
        Point2f pt(basis[i]);
        line(drawImg, origin, origin + pt, Scalar(0, (double)(i * 255), 0), 2);
    }
}

void CirclesGridFinder::drawBasisGraphs(const std::vector<Graph> &basisGraphs,
                                        Mat &drawImage,
                                        bool drawEdges,
                                        bool drawVertices) const {
    // const int vertexRadius = 1;
    const int vertexRadius = 3;
    const Scalar vertexColor = Scalar(0, 0, 255);
    const int vertexThickness = -1;

    const Scalar edgeColor = Scalar(255, 0, 0);
    // const int edgeThickness = 1;
    const int edgeThickness = 2;

    if (drawEdges) {
        for (size_t i = 0; i < basisGraphs.size(); i++) {
            for (size_t v1 = 0; v1 < basisGraphs[i].getVerticesCount(); v1++) {
                for (size_t v2 = 0; v2 < basisGraphs[i].getVerticesCount(); v2++) {
                    if (basisGraphs[i].areVerticesAdjacent(v1, v2)) {
                        line(drawImage, keypoints[v1], keypoints[v2], edgeColor, edgeThickness);
                    }
                }
            }
        }
    }
    if (drawVertices) {
        for (size_t v = 0; v < basisGraphs[0].getVerticesCount(); v++) {
            circle(drawImage, keypoints[v], vertexRadius, vertexColor, vertexThickness);
        }
    }
}

void CirclesGridFinder::drawHoles(const Mat &srcImage, Mat &drawImage) const {
    // const int holeRadius = 4;
    // const int holeRadius = 2;
    // const int holeThickness = 1;
    const int holeRadius = 3;
    const int holeThickness = -1;

    // const Scalar holeColor = Scalar(0, 0, 255);
    const Scalar holeColor = Scalar(0, 255, 0);

    if (srcImage.channels() == 1)
        cvtColor(srcImage, drawImage, COLOR_GRAY2RGB);
    else
        srcImage.copyTo(drawImage);

    for (size_t i = 0; i < holes.size(); i++) {
        for (size_t j = 0; j < holes[i].size(); j++) {
            if (j != holes[i].size() - 1)
                line(drawImage, keypoints[holes[i][j]], keypoints[holes[i][j + 1]],
                     Scalar(255, 0, 0), 2);
            if (i != holes.size() - 1)
                line(drawImage, keypoints[holes[i][j]], keypoints[holes[i + 1][j]],
                     Scalar(255, 0, 0), 2);

            circle(drawImage, keypoints[holes[i][j]], holeRadius, holeColor, holeThickness);
        }
    }
}

Size CirclesGridFinder::getDetectedGridSize() const {
    if (holes.size() == 0) return Size(0, 0);

    return Size((int)holes[0].size(), (int)holes.size());
}

void CirclesGridFinder::getHoles(std::vector<Point2f> &outHoles) const {
    outHoles.clear();

    for (size_t i = 0; i < holes.size(); i++) {
        for (size_t j = 0; j < holes[i].size(); j++) {
            outHoles.push_back(keypoints[holes[i][j]]);
        }
    }
}

static bool areIndicesCorrect(Point pos, std::vector<std::vector<size_t>> *points) {
    if (pos.y < 0 || pos.x < 0) return false;
    return (static_cast<size_t>(pos.y) < points->size() &&
            static_cast<size_t>(pos.x) < points->at(pos.y).size());
}

void CirclesGridFinder::getAsymmetricHoles(std::vector<cv::Point2f> &outHoles) const {
    outHoles.clear();

    std::vector<Point> largeCornerIndices, smallCornerIndices;
    std::vector<Point> firstSteps, secondSteps;
    size_t cornerIdx =
        getFirstCorner(largeCornerIndices, smallCornerIndices, firstSteps, secondSteps);
    CV_Assert(largeHoles != 0 && smallHoles != 0);

    Point srcLargePos = largeCornerIndices[cornerIdx];
    Point srcSmallPos = smallCornerIndices[cornerIdx];

    while (areIndicesCorrect(srcLargePos, largeHoles) ||
           areIndicesCorrect(srcSmallPos, smallHoles)) {
        Point largePos = srcLargePos;
        while (areIndicesCorrect(largePos, largeHoles)) {
            outHoles.push_back(keypoints[largeHoles->at(largePos.y)[largePos.x]]);
            largePos += firstSteps[cornerIdx];
        }
        srcLargePos += secondSteps[cornerIdx];

        Point smallPos = srcSmallPos;
        while (areIndicesCorrect(smallPos, smallHoles)) {
            outHoles.push_back(keypoints[smallHoles->at(smallPos.y)[smallPos.x]]);
            smallPos += firstSteps[cornerIdx];
        }
        srcSmallPos += secondSteps[cornerIdx];
    }
}

double CirclesGridFinder::getDirection(Point2f p1, Point2f p2, Point2f p3) {
    Point2f a = p3 - p1;
    Point2f b = p2 - p1;
    return a.x * b.y - a.y * b.x;
}

bool CirclesGridFinder::areSegmentsIntersecting(Segment seg1, Segment seg2) {
    bool doesStraddle1 =
        (getDirection(seg2.s, seg2.e, seg1.s) * getDirection(seg2.s, seg2.e, seg1.e)) < 0;
    bool doesStraddle2 =
        (getDirection(seg1.s, seg1.e, seg2.s) * getDirection(seg1.s, seg1.e, seg2.e)) < 0;
    return doesStraddle1 && doesStraddle2;

    /*
     Point2f t1 = e1-s1;
     Point2f n1(t1.y, -t1.x);
     double c1 = -n1.ddot(s1);

     Point2f t2 = e2-s2;
     Point2f n2(t2.y, -t2.x);
     double c2 = -n2.ddot(s2);

     bool seg1 = ((n1.ddot(s2) + c1) * (n1.ddot(e2) + c1)) <= 0;
     bool seg1 = ((n2.ddot(s1) + c2) * (n2.ddot(e1) + c2)) <= 0;

     return seg1 && seg2;
     */
}

void CirclesGridFinder::getCornerSegments(const std::vector<std::vector<size_t>> &points,
                                          std::vector<std::vector<Segment>> &segments,
                                          std::vector<Point> &cornerIndices,
                                          std::vector<Point> &firstSteps,
                                          std::vector<Point> &secondSteps) const {
    segments.clear();
    cornerIndices.clear();
    firstSteps.clear();
    secondSteps.clear();
    int h = (int)points.size();
    int w = (int)points[0].size();
    CV_Assert(h >= 2 && w >= 2);

    // all 8 segments with one end in a corner
    std::vector<Segment> corner;
    corner.emplace_back(keypoints[points[1][0]], keypoints[points[0][0]]);
    corner.emplace_back(keypoints[points[0][0]], keypoints[points[0][1]]);
    segments.push_back(corner);
    cornerIndices.emplace_back(0, 0);
    firstSteps.emplace_back(1, 0);
    secondSteps.emplace_back(0, 1);
    corner.clear();

    corner.emplace_back(keypoints[points[0][w - 2]], keypoints[points[0][w - 1]]);
    corner.emplace_back(keypoints[points[0][w - 1]], keypoints[points[1][w - 1]]);
    segments.push_back(corner);
    cornerIndices.emplace_back(w - 1, 0);
    firstSteps.emplace_back(0, 1);
    secondSteps.emplace_back(-1, 0);
    corner.clear();

    corner.emplace_back(keypoints[points[h - 2][w - 1]], keypoints[points[h - 1][w - 1]]);
    corner.emplace_back(keypoints[points[h - 1][w - 1]], keypoints[points[h - 1][w - 2]]);
    segments.push_back(corner);
    cornerIndices.emplace_back(w - 1, h - 1);
    firstSteps.emplace_back(-1, 0);
    secondSteps.emplace_back(0, -1);
    corner.clear();

    corner.emplace_back(keypoints[points[h - 1][1]], keypoints[points[h - 1][0]]);
    corner.emplace_back(keypoints[points[h - 1][0]], keypoints[points[h - 2][0]]);
    cornerIndices.emplace_back(0, h - 1);
    firstSteps.emplace_back(0, -1);
    secondSteps.emplace_back(1, 0);
    segments.push_back(corner);
    corner.clear();

    // y axis is inverted in computer vision so we check < 0
    bool isClockwise = getDirection(keypoints[points[0][0]], keypoints[points[0][w - 1]],
                                    keypoints[points[h - 1][w - 1]]) < 0;
    if (!isClockwise) {
#ifdef DEBUG_CIRCLES
        std::cout << "Corners are counterclockwise" << std::endl;
#endif
        std::reverse(segments.begin(), segments.end());
        std::reverse(cornerIndices.begin(), cornerIndices.end());
        std::reverse(firstSteps.begin(), firstSteps.end());
        std::reverse(secondSteps.begin(), secondSteps.end());
        std::swap(firstSteps, secondSteps);
    }
}

bool CirclesGridFinder::doesIntersectionExist(const std::vector<Segment> &corner,
                                              const std::vector<std::vector<Segment>> &segments) {
    for (size_t i = 0; i < corner.size(); i++) {
        for (size_t j = 0; j < segments.size(); j++) {
            for (size_t k = 0; k < segments[j].size(); k++) {
                if (areSegmentsIntersecting(corner[i], segments[j][k])) return true;
            }
        }
    }

    return false;
}

size_t CirclesGridFinder::getFirstCorner(std::vector<Point> &largeCornerIndices,
                                         std::vector<Point> &smallCornerIndices,
                                         std::vector<Point> &firstSteps,
                                         std::vector<Point> &secondSteps) const {
    std::vector<std::vector<Segment>> largeSegments;
    std::vector<std::vector<Segment>> smallSegments;

    getCornerSegments(*largeHoles, largeSegments, largeCornerIndices, firstSteps, secondSteps);
    getCornerSegments(*smallHoles, smallSegments, smallCornerIndices, firstSteps, secondSteps);

    const size_t cornersCount = 4;
    CV_Assert(largeSegments.size() == cornersCount);

    bool isInsider[cornersCount];

    for (size_t i = 0; i < cornersCount; i++) {
        isInsider[i] = doesIntersectionExist(largeSegments[i], smallSegments);
    }

    int cornerIdx = 0;
    bool waitOutsider = true;

    for (size_t i = 0; i < cornersCount * 2; ++i) {
        if (waitOutsider) {
            if (!isInsider[(cornerIdx + 1) % cornersCount]) waitOutsider = false;
        } else {
            if (isInsider[(cornerIdx + 1) % cornersCount]) return cornerIdx;
        }

        cornerIdx = (cornerIdx + 1) % cornersCount;
    }

    CV_Error(Error::StsNoConv, "isInsider array has the same values");
}

bool FindCirclesGrid(cv::InputArray _image,
                     cv::Size patternSize,
                     cv::OutputArray _centers,
                     int flags,
                     const cv::Ptr<cv::FeatureDetector> &blobDetector,
                     const ns_cv_helper::CirclesGridFinderParameters &parameters_) {
    // using namespace cv;
    // CV_INSTRUMENT_REGION();

    ns_cv_helper::CirclesGridFinderParameters parameters =
        parameters_;  // parameters.gridType is amended below

    bool isAsymmetricGrid = (flags & cv::CALIB_CB_ASYMMETRIC_GRID) ? true : false;
    bool isSymmetricGrid = (flags & cv::CALIB_CB_SYMMETRIC_GRID) ? true : false;
    CV_Assert(isAsymmetricGrid ^ isSymmetricGrid);

    std::vector<cv::Point2f> centers;

    std::vector<cv::Point2f> points;
    if (blobDetector) {
        std::vector<cv::KeyPoint> keypoints;
        blobDetector->detect(_image, keypoints);
        for (size_t i = 0; i < keypoints.size(); i++) {
            points.push_back(keypoints[i].pt);
        }
    } else {
        CV_CheckTypeEQ(_image.type(), CV_32FC2,
                       "blobDetector must be provided or image must contains Point2f array "
                       "(std::vector<Point2f>) with candidates");
        _image.copyTo(points);
    }

    if (flags & cv::CALIB_CB_ASYMMETRIC_GRID)
        parameters.gridType = ns_cv_helper::CirclesGridFinderParameters::ASYMMETRIC_GRID;
    if (flags & cv::CALIB_CB_SYMMETRIC_GRID)
        parameters.gridType = ns_cv_helper::CirclesGridFinderParameters::SYMMETRIC_GRID;

    bool isValid = false;
    const int attempts = 2;
    const size_t minHomographyPoints = 4;
    cv::Mat H;
    for (int i = 0; i < attempts; i++) {
        centers.clear();
        ns_cv_helper::CirclesGridFinder boxFinder(patternSize, points, parameters);
        try {
            bool isFound = boxFinder.findHoles();
            if (isFound) {
                switch (parameters.gridType) {
                    case ns_cv_helper::CirclesGridFinderParameters::SYMMETRIC_GRID:
                        boxFinder.getHoles(centers);
                        break;
                    case ns_cv_helper::CirclesGridFinderParameters::ASYMMETRIC_GRID:
                        boxFinder.getAsymmetricHoles(centers);
                        break;
                    default:
                        CV_Error(cv::Error::StsBadArg, "Unknown pattern type");
                }

                isValid = true;
                break;  // done, return result
            }
        } catch (const cv::Exception &e) {
            CV_UNUSED(e);
            // CV_LOG_DEBUG(NULL, "findCirclesGrid2: attempt=" << i << ": " << e.what());
            // nothing, next attempt
        }

        boxFinder.getHoles(centers);
        if (i != attempts - 1) {
            if (centers.size() < minHomographyPoints) break;
            H = ns_cv_helper::CirclesGridFinder::rectifyGrid(boxFinder.getDetectedGridSize(),
                                                             centers, points, points);
        }
    }

    if (!centers.empty() && !H.empty())  // undone rectification
    {
        cv::Mat orgPointsMat;
        transform(centers, orgPointsMat, H.inv());
        convertPointsFromHomogeneous(orgPointsMat, centers);
    }
    cv::Mat(centers).copyTo(_centers);
    return isValid;
}
}  // namespace ns_cv_helper