
#include "curve2d.hpp"

#include <utility>
#include <map>
#include <cmath>

#include "transform2d.hpp"
#include "propagate2d.hpp"
#include "plotanalysis.hpp"

int numberOfNeighbors(const Image2D<char>& mask_original, Int2 loc){
    auto neighbors = 0;
    auto w = mask_original.width();
    auto h = mask_original.height();

    for(int xd = -1; xd <= 1; xd++) {
        for(int yd = -1; yd <= 1; yd++) {
            if(xd == 0 && yd == 0) continue;
            Int2 neighbour { loc[0] + xd, loc[1] + yd };
            if(neighbour[0] < 0 || neighbour[0] >= w) continue;
            if(neighbour[1] < 0 || neighbour[1] >= h) continue;
            if (!mask_original(neighbour[0], neighbour[1])) continue;
            neighbors++;
        }
    }

    return neighbors;
}


void tryLongestPathFrom(const Image2D<char>& mask_original, Int2 start, float& currentMaxLength, std::vector<Int2>& currentMaxPath, int pixPerSegment) {
    auto w = mask_original.width();
    auto h = mask_original.height();

    Image2D<char> mask = mask_original;
    Image2D<Int2> movement {w, h};
    movement.fill({0, 0});
    std::multimap<float, std::pair<Int2, Int2>> queue;
    queue.insert({0.0f, {0, start}});

    float maxDist = -1;
    Int2 maxPos = start;

    //breadth first search
    while(!queue.empty()) {
        float dist;
        Int2 mov, pos;
        {
            auto begin_iterator = queue.begin();
            dist = begin_iterator->first;
            mov = begin_iterator->second.first;
            pos =  begin_iterator->second.second;
            queue.erase(begin_iterator);
        }
        if(!mask(pos[0], pos[1])) continue;

        mask(pos[0], pos[1]) = 0;
        movement(pos[0], pos[1]) = mov;

        if(dist > maxDist) {
            maxDist = dist;
            maxPos = pos;
        }

        for(int xd = -1; xd <= 1; xd++) {
            for(int yd = -1; yd <= 1; yd++) {
                if(xd == 0 && yd == 0) continue;
                Int2 neighbour { pos[0] + xd, pos[1] + yd };
                if(neighbour[0] < 0 || neighbour[0] >= w) continue;
                if(neighbour[1] < 0 || neighbour[1] >= h) continue;
                //if(!mask(neighbour[0], neighbour[1])) continue;

                queue.insert({dist + 1, {{xd, yd}, neighbour}});
            }
        }
    }

    if(maxDist > currentMaxLength) {
        currentMaxLength = maxDist;
        currentMaxPath.resize(0);
        auto pos = maxPos;
        currentMaxPath.push_back(pos);
        int i = 0;

        while(movement(pos[0], pos[1]) != Int2({0, 0})) {
            pos -= movement(pos[0], pos[1]);
            i++;
            i %= pixPerSegment;
            if(i == 0) {
                currentMaxPath.push_back(pos);
            }
        }
    }
}

Float2 curveDirectionVector(const std::vector<Int2>& curve, int index) {
    Float2 result;
    if(index == 0) {
        result = curve[1] - curve[0];
    } else if(index == curve.size() -1) {
        result = curve[index] - curve[index-1];
    } else {
        Float2 a = curve[index+1] - curve[index];
        a /= sqrt(vectorLengthSquared(a));
        Float2 b = curve[index] - curve[index-1];
        b /= sqrt(vectorLengthSquared(b));
        result = a + b;
    }
    result /= sqrt(vectorLengthSquared(result));
    return result;
}

float curveLocalLength(const std::vector<Int2>& curve, int index) {

    if(index == 0) {
        return sqrt(vectorLengthSquared(curve[1] - curve[0]));
    } else if(index == curve.size() - 1) {
        return sqrt(vectorLengthSquared(curve[index] - curve[index - 1]));
    } else {
        return 0.5 * sqrt(vectorLengthSquared(curve[index + 1] - curve[index]))
            + 0.5 * sqrt(vectorLengthSquared(curve[index] - curve[index - 1]));
    }
}


std::vector<AugmentedPoint> augmentCurve(const std::vector<Int2>& curve, const Image2D<unsigned>& plot, int smoothCount) {
    std::vector<AugmentedPoint> result;
    result.resize(curve.size());

    auto width = plot.width();
    auto height = plot.height();

    unsigned total = 0;
    for(int x = 0; x < width; x++) {
        for(int y = 0; y < height; y++) {
            total += plot(x, y);
        }
    }

    Image2D<int> pointsImage { width, height };
    Image2D<float> pointsDistance { width, height };
    Image2D<char> dummyMask { width, height };

    pointsImage.fill(-1);
    pointsDistance.fill(INFINITY);
    dummyMask.fill(1);

    for(int i = 0; i < curve.size(); i++) {
        pointsImage(curve[i][0], curve[i][1]) = i;
        pointsDistance(curve[i][0], curve[i][1]) = 0;
    }

    auto origins = propagate2D(dummyMask, pointsDistance);
    auto labels = indexImage(pointsImage, origins);

    auto regions = regionStats(plot, labels, curve.size());

    for(int i = 0; i < curve.size(); i++) {
        auto direction = curveDirectionVector(curve, i);
        auto ortho = Float2({direction[1], -direction[0]});

        RegionStat smoothedStat;
        smoothedStat.pointStat.count = 0;
        smoothedStat.productStat.count = 0;
        for(int j = std::max(i - smoothCount, 0); j <= i + smoothCount && j < regions.size(); j++) {
            smoothedStat = smoothedStat + regions[j];
        }

        AugmentedPoint ap;
        ap.point = Float2({ curve[i][0] * 1.0f, curve[i][1] * 1.0f });
        ap.direction = direction;
        ap.amount = regions[i].pointStat.count * 1.0f / total;
        ap.amountPerLength = ap.amount / curveLocalLength(curve, i);

        if(!smoothedStat.empty()) {
            ap.orthoVariance = smoothedStat.covarianceMatrix().varianceAlong(ortho);
        }

        result[i] = ap;
    }

    return result;
}

Float2 couplingDistanceAugmented_nvec(Float2 vec, Float2 orig, float scale, float arc) {
    vec -= orig;
    vec /= scale;
    auto s = std::sin(-arc);
    auto c = std::cos(-arc);
    return Float2({c * vec[0] - s * vec[1], s * vec[0] + c * vec[1]});
}

float couplingDistanceAugmentedSqured(const std::vector<AugmentedPoint>& curve1, const std::vector<AugmentedPoint>& curve2, float vf, float af, bool enableMirror) {
    Image2D<float> table { curve1.size(), curve2.size() };

    if ( curve1.empty() || curve2.empty())
        return 999999999999;

    //normalization
    Float2 orig1 = curve1[0].point;
    Float2 orig2_array[2] = {curve2[0].point, curve2[curve2.size() - 1].point};
    Float2 way1 = curve1[curve1.size() - 1].point - curve1[0].point;
    Float2 way2 = curve2[curve2.size() - 1].point - curve2[0].point;
    float scale1 = sqrt(vectorLengthSquared(way1));
    float scale2 = sqrt(vectorLengthSquared(way2));
    float arc1 = std::atan2(way1[1], way1[0]);
    float arc2_base = std::atan2(way2[1], way2[0]);
    float arc2_array[2] = {arc2_base, arc2_base - M_PI};

    float result;

    for(int mirror = 0; mirror <= (enableMirror ? 1 : 0); mirror++) {
        for(int rotation = 0; rotation <= 1; rotation++) {
            Float2 orig2 = orig2_array[rotation];
            float arc2 = arc2_array[rotation];
            for(int x = 0; x < curve1.size(); x++) {
                for(int y = 0; y < curve2.size(); y++) {
                    int yi = y;
                    if(rotation != 0) {
                        yi = curve2.size() - 1 - yi;
                    }
                    auto p1 = couplingDistanceAugmented_nvec(curve1[x].point, orig1, scale1, arc1);
                    auto p2 = couplingDistanceAugmented_nvec(curve2[yi].point, orig2, scale2, arc2);
                    if(mirror) {
                        p2[1] *= -1;
                    }
                    auto distSquared = vectorLengthSquared(p1 - p2);
                    auto v1 = curve1[x].orthoVariance / scale1;
                    auto v2 = curve2[yi].orthoVariance / scale2;
                    auto varianceDistanceSquared = (v1 - v2) * (v1 - v2);
                    auto a1 = curve1[x].amountPerLength * scale1;
                    auto a2 = curve2[yi].amountPerLength * scale2;
                    auto amountDistanceSquared = (a1 - a2) * (a1 - a2);
                    auto value = vf * varianceDistanceSquared + af * amountDistanceSquared + (1.0f - vf - af) * distSquared;

                    if(x == 0 && y == 0) {
                        table(x, y) = value;
                    } else if(x == 0) {
                        table(x, y) = std::max(value, table(x, y-1));
                    } else if(y == 0) {
                        table(x, y) = std::max(value, table(x-1, y));
                    } else {
                        table(x, y) = std::max(value, std::min(table(x-1, y), std::min(table(x, y-1), table(x-1, y-1))) );
                    }
                }
            }
            if(rotation == 0 && mirror == 0) {
                result =  table(curve1.size() - 1, curve2.size() - 1);
            } else {
                result =  std::min(result, table(curve1.size() - 1, curve2.size() - 1));
            }
        }
    }
    return result;
}

float couplingDistanceAugmented(const std::vector<AugmentedPoint>& curve1, const std::vector<AugmentedPoint>& curve2, float vf, float af, bool enableMirror) {
    return sqrt(couplingDistanceAugmentedSqured(curve1, curve2, vf, af, enableMirror));
}

float sharpenssOfAugmentedCurve(const std::vector<AugmentedPoint>& curve) {
    Float2 way = curve[0].point - curve[curve.size() - 1].point;
    float scale = std::sqrt(vectorLengthSquared(way));
    Stat<float> stat;
    for(int i = 0; i < curve.size(); i++) {
        stat.put(curve[i].orthoVariance / scale);
    }
    return stat.average;
}
