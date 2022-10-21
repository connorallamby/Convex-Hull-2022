import math
import sys
from typing import List
from typing import Tuple

EPSILON = sys.float_info.epsilon
Point = Tuple[int, int]

def y_intercept(p1: Point, p2: Point, x: float) -> float:
    """
    Given two points, p1 and p2, an x coordinate from a vertical line,
    compute and return the the y-intercept of the line segment p1->p2
    with the vertical line passing through x.
    """
    x1, y1 = p1
    x2, y2 = p2
    '''if x2 == x1:
        slope = 0
    else:
    '''
    slope = (y2 - y1) / (x2 - x1)
    return y1 + (x - x1) * slope


def triangle_area(a: Point, b: Point, c: Point) -> float:
    """
    Given three points a,b,c,
    computes and returns the area defined by the triangle a,b,c.
    Note that this area will be negative if a,b,c represents a clockwise sequence,
    positive if it is counter-clockwise,
    and zero if the points are collinear.
    """
    ax, ay = a
    bx, by = b
    cx, cy = c
    return ((cx - bx) * (by - ay) - (bx - ax) * (cy - by)) / 2


def is_clockwise(a: Point, b: Point, c: Point) -> bool:
    """
    Given three points a,b,c,
    returns True if and only if a,b,c represents a clockwise sequence
    (subject to floating-point precision)
    """
    return triangle_area(a, b, c) < -EPSILON


def is_counter_clockwise(a: Point, b: Point, c: Point) -> bool:
    """
    Given three points a,b,c,
    returns True if and only if a,b,c represents a counter-clockwise sequence
    (subject to floating-point precision)
    """
    return triangle_area(a, b, c) > EPSILON


def collinear(a: Point, b: Point, c: Point) -> bool:
    """
    Given three points a,b,c,
    returns True if and only if a,b,c are collinear
    (subject to floating-point precision)
    """
    return abs(triangle_area(a, b, c)) <= EPSILON


def clockwise_sort(points: List[Point]):
    """
    Given a list of points, sorts those points in clockwise order about their centroid.
    Note: this function modifies its argument.
    """
    # get mean x coord, mean y coord
    x_mean = sum(p[0] for p in points) / len(points)
    y_mean = sum(p[1] for p in points) / len(points)

    def angle(point: Point):
        return (math.atan2(point[1] - y_mean, point[0] - x_mean) + 2 * math.pi) % (2 * math.pi)

    points.sort(key=angle)
    return


def base_case_hull(points: List[Point]) -> List[Point]:
    """ Base case of the recursive algorithm.
    """
    # TODO: You need to implement this function.
    # empty set - will be filled with points that should remain on the hull
    hull = set()
    length = len(points)
    # The first two loops will iterate through the point list and make line segments between i-j
    for i in range(length - 1):
        for j in range(i + 1, length):
            # create clockwise and counterclockwise counters
            cwCount = 0
            ccCount = 0
            # Loop will test all remaing permutations of segment j-k to prove the current
            # i-j can remain on the hull
            for k in range(length):
                if k != i and k != j:
                    # use the isclockiwise function to determine the type of turn
                    if is_counter_clockwise(points[i], points[j], points[k]):
                        ccCount += 1
                    else:
                        cwCount += 1
            # if all turns made are of the same type i-j can remain on the hull
            if ccCount == 0 or cwCount == 0:
                hull.add(points[i])
                hull.add(points[j])
    # load the set information in output list
    points = [p for p in hull]
    clockwise_sort(points)

    return points


def find_median(points: List[Point]):
    length = len(points)
    if length % 2 == 0:
        return int(length / 2)
    else:
        return math.ceil(length / 2)


def compute_hull(points: List[Point]) -> List[Point]:
    """
    Given a list of points, computes the convex hull around those points
    and returns only the points that are on the hull.
    """
    # TODO: Implement a correct computation of the convex hull
    #  using the divide-and-conquer algorithm
    # TODO: Document your Initialization, Maintenance and Termination invariants.
    allColinear = True
    # check if all points are colinear
    for i in range(len(points)):
        if points[0][1] != points[i][1]:
            allColinear = False
            break

    # case for all colinear points
    if allColinear:
        # find leftmost and right most points and that is the hull
        finalRight, rightMost = findRightMost(points)
        finalLeft, leftMost = findLeftMost(points)

        leftPoint = points[finalLeft]
        rightPoint = points[finalRight]
        points = [leftPoint, rightPoint]
        return points

    # base case
    elif len(points) <= 7:
        clockwise_sort(points)
        points = base_case_hull(points)
        return points

    else:
        # Initialization Invariant: Points are sorted in ascending order and split into a left
        # and right list
        points = sorted(points, key=lambda k: [k[0], k[1]])
        median = find_median(points)
        leftList = []
        rightList = []
        # split the list of points but ensure coordinates with the same x values end up on the same half
        for i in range(len(points)):
            if i > median and points[i][0] != points[median][0]:
                rightList.append(points[i])
            else:
                leftList.append(points[i])

        # Recursive call with each half of the previously too large set of points
        leftList = compute_hull(leftList)
        rightList = compute_hull(rightList)

        clockwise_sort(leftList)
        clockwise_sort(rightList)

        # get the upper and lower tanget points indices
        upperTanLeft, upperTanRight = upperTangent(leftList, rightList)
        lowerTanLeft, lowerTanRight = lowerTangent(leftList, rightList)

        points = merge(lowerTanLeft, lowerTanRight, upperTanLeft, upperTanRight, leftList, rightList)

        clockwise_sort(points)
        return points


def findRightMost(inputList: List[Point]):
    # index where rightmost is found
    currentLeft = 0
    # x value of righmost
    rightMost = inputList[0][0]
    for i in range(len(inputList)):
        if inputList[i][0] >= rightMost:
            rightMost = inputList[i][0]
            currentLeft = i

    return currentLeft, rightMost

def findLeftMost(inputList: List[Point]):
    # index where leftmost is found
    currentRight = 0
    # x value of leftmose
    leftMost = inputList[0][0]
    for i in range(len(inputList)):
        if inputList[i][0] <= leftMost:
            leftMost = inputList[i][0]
            currentRight = i

    return currentRight, leftMost


# Maintenance Invariant - Each iteration through, 4 tangent points are found, 
# upper left, upper right, lower left, lower right
def upperTangent(leftList: List[Point], rightList: List[Point]):
# find a suitable x value for the vertical line - will store it in l
    # index where rightmost is found on the leftList, x value of rightMost
    currentLeft, rightMost = findRightMost(leftList)

    # index where leftmost is found on the rightList, x value of leftMost
    currentRight, leftMost = findLeftMost(rightList)

    # l is the middle of the left list's most right point and the right list's left most point
    l = (leftMost + rightMost) / 2

    currentIntercept = y_intercept(leftList[currentLeft], rightList[currentRight], l)

    # increment the left and right as long as the y intercepts are lower than currentIntercept
    while y_intercept(leftList[currentLeft], rightList[(currentRight + 1) % len(rightList)], l) < currentIntercept or y_intercept(leftList[(currentLeft - 1) % len(leftList)], rightList[currentRight], l) < currentIntercept:
        # if the new intercept with current right + 1 is lower then increment the current right else increment left
        if y_intercept(leftList[currentLeft], rightList[(currentRight + 1) % len(rightList)],l) < currentIntercept:
            currentIntercept = y_intercept(leftList[currentLeft], rightList[(currentRight + 1) % len(rightList)], l)
            currentRight = (currentRight + 1) % len(rightList)
        else:
            currentIntercept = y_intercept(leftList[(currentLeft - 1) % len(leftList)], rightList[currentRight], l)
            currentLeft = (currentLeft - 1) % len(leftList)

    #return the indecies of the tangents
    return currentLeft, currentRight

def lowerTangent(leftList: List[Point], rightList: List[Point]):
# find a sutible x value for the vertical line - will store it in l
    # index where rightmost is found on the leftList, x value of rightMost
    currentLeft, rightMost = findRightMost(leftList)

    # index where leftmost is found on the rightList, x value of leftMost
    currentRight, leftMost = findLeftMost(rightList)

    # l is the middle of the left list's most right point and the right list's left most point
    l = (leftMost + rightMost) / 2

    currentIntercept = y_intercept(leftList[currentLeft], rightList[currentRight], l)

    # increment the left and right as long as the y intercepts are higher than currentIntercept
    while y_intercept(leftList[currentLeft], rightList[(currentRight - 1) % len(rightList)], l) > currentIntercept or y_intercept(leftList[(currentLeft + 1) % len(leftList)], rightList[currentRight], l) > currentIntercept:

        # if the new intercept with current right + 1 is greater then increment the current right else increment left
        if y_intercept(leftList[currentLeft], rightList[(currentRight - 1) % len(rightList)], l) > currentIntercept:
            currentIntercept = y_intercept(leftList[currentLeft], rightList[(currentRight - 1) % len(rightList)], l)
            currentRight = (currentRight - 1) % len(rightList)
        else:
            currentIntercept = y_intercept(leftList[(currentLeft + 1) % len(leftList)], rightList[currentRight], l)
            currentLeft = (currentLeft + 1) % len(leftList)

    return currentLeft, currentRight

# Given a left convex hull and right convex hull, returns one combined hull
# will accept the index of lower left, lower right, upper left, upper right tangents as well as left and right list
def merge(lowerLeftTan: int, lowerRightTan: int, upperLeftTan: int, upperRightTan: int, leftList: List[Point],
          rightList: List[Point]) -> List[Point]:
    '''Delete points that are in between upper and lower tangent points that are not on the hull,
    then combine the two lists into one list'''
    # Termination Invariant - Two hulls are merged and the points that are inside the convex hull 
    # are removed from the points list
    clockwise_sort(leftList)
    clockwise_sort(rightList)

    # pop points that are in between upperLeftTan and lowerLeftTan
    current_point_l_index = (upperLeftTan + 1) % len(leftList)
    lowerLeftPoint = leftList[lowerLeftTan]
    while current_point_l_index != leftList.index(lowerLeftPoint):
        leftList.pop(current_point_l_index)
        current_point_l_index = current_point_l_index % len(leftList)

    # pop points that are in between upperRightTan and lowerRightTan
    current_point_r_index = (upperRightTan - 1) % len(rightList)
    lowerRightPoint = rightList[lowerRightTan]
    while current_point_r_index != rightList.index(lowerRightPoint):
        rightList.pop(current_point_r_index)
        current_point_r_index = (current_point_r_index - 1) % len(rightList)

    # add the completed lists together to create one combined hull
    points = leftList + rightList
    clockwise_sort(points)

    return points