//
// Copyright (c) 2015 Mahyar Khayatkhoei
// Copyright (c) 2009-2014 Shawn Singh, Glen Berseth, Mubbasir Kapadia, Petros Faloutsos, Glenn Reinman
// See license.txt for complete license.
//

#include <algorithm>
#include <vector>
#include <util/Geometry.h>
#include <util/Curve.h>
#include <util/Color.h>
#include <util/DrawLib.h>
#include "Globals.h"

using namespace Util;

Curve::Curve(const CurvePoint& startPoint, int curveType) : type(curveType)
{
	controlPoints.push_back(startPoint);
}

Curve::Curve(const std::vector<CurvePoint>& inputPoints, int curveType) : type(curveType)
{
	controlPoints = inputPoints;
	sortControlPoints();
}

// Add one control point to the vector controlPoints
void Curve::addControlPoint(const CurvePoint& inputPoint)
{
	controlPoints.push_back(inputPoint);
	sortControlPoints();
}

// Add a vector of control points to the vector controlPoints
void Curve::addControlPoints(const std::vector<CurvePoint>& inputPoints)
{
	for (int i = 0; i < inputPoints.size(); i++)
		controlPoints.push_back(inputPoints[i]);
	sortControlPoints();
}

// Draw the curve shape on screen, usign window as step size (bigger window: less accurate shape)
void Curve::drawCurve(Color curveColor, float curveThickness, int window)
{	
	Point & a,b;
#ifdef ENABLE_GUI
	for ( int t = 0; t < controlPoints.end().time; t = t + window)
	{
		Curve::calculatePoint(Point& a ,t);
		Curve::calculatePoint(Point& b ,t+window);
		Drawlib::drawLine(const Point & a, const Point & b, const Color curveColor, float curveThickness);
	}
	//================DELETE THIS PART AND THEN START CODING===================
	//static bool flag = false;
	//if (!flag)
	//{
	//	std::cerr << "ERROR>>>>Member function drawCurve is not implemented!" << std::endl;
	//	flag = true;
	//}
	//=========================================================================

	// Robustness: make sure there is at least two control point: start and end points

	// Move on the curve from t=0 to t=finalPoint, using window as step size, and linearly interpolate the curve points

	return;
#endif
}

// Sort controlPoints vector in ascending order: min-first
void Curve::sortControlPoints()
{
	//binary search and insert sort maybe better,lambda for c++ 11
	std::sort(controlPoints.begin(), controlPoints.end(), [](const CurvePoint &left, const CurvePoint &right) {return left.time < right.time; });

	//================DELETE THIS PART AND THEN START CODING===================
	/*static bool flag = false;
	if (!flag)
	{
	std::cerr << "ERROR>>>>Member function sortControlPoints is not implemented!" << std::endl;
	flag = true;
	}*/
	//=========================================================================

	return;
}

// Calculate the position on curve corresponding to the given time, outputPoint is the resulting position
bool Curve::calculatePoint(Point& outputPoint, float time)
{
	// Robustness: make sure there is at least two control point: start and end points
	if (!checkRobust())
		return false;

	// Define temporary parameters for calculation
	unsigned int nextPoint;
	float normalTime, intervalTime;

	// Find the current interval in time, supposing that controlPoints is sorted (sorting is done whenever control points are added)
	if (!findTimeInterval(nextPoint, time))
		return false;

	// Calculate position at t = time on curve
	if (type == hermiteCurve)
	{
		outputPoint = useHermiteCurve(nextPoint, time);
	}
	else if (type == catmullCurve)
	{
		outputPoint = useCatmullCurve(nextPoint, time);
	}

	// Return
	return true;
}

// Check Roboustness
bool Curve::checkRobust()
{
	if (controlPoints.size() < 2)
		return false;


	//================DELETE THIS PART AND THEN START CODING===================
	/*static bool flag = false;
	if (!flag)
	{
	std::cerr << "ERROR>>>>Member function checkRobust is not implemented!" << std::endl;
	flag = true;
	}*/
	//=========================================================================


	return true;
}

// Find the current time interval (i.e. index of the next control point to follow according to current time)
bool Curve::findTimeInterval(unsigned int& nextPoint, float time)
{

	for (int i = 0; i < controlPoints.size(); i++)
	{
		if (controlPoints[i - 1].time < time && controlPoints[i].time > time)
		{
			nextPoint = i;
			return true;
		}
	}


	//================DELETE THIS PART AND THEN START CODING===================

	//=========================================================================


	return false;
}

// Implement Hermite curve
Point Curve::useHermiteCurve(const unsigned int nextPoint, const float time)
{
	Point newPosition;
	float normalTime, intervalTime;



	CurvePoint &P1 = controlPoints[nextPoint - 1];
	CurvePoint &P2 = controlPoints[nextPoint];
	float s = (time - P1.time) / (P2.time - P1.time);
	float h1 = 2 * pow(s, 3) - 3 * pow(s, 2) + 1;
	float h2 = -2 * pow(s, 3) + 3 * pow(s, 2);
	float h3 = pow(s, 3) - 2 * pow(s, 2) + s;
	float h4 = pow(s, 3) - pow(s, 2);

	newPosition = h1*P1.position + h2*P2.position + h3*P1.tangent*(P2.time-P1.time) + h4*P2.tangent*(P2.time-P1.time);



	//================DELETE THIS PART AND THEN START CODING===================
	/*	static bool flag = false;
	if (!flag)
	{
	std::cerr << "ERROR>>>>Member function useHermiteCurve is not implemented!" << std::endl;
	flag = true;
	}*/
	//=========================================================================


	// Calculate time interval, and normal time required for later curve calculations

	// Calculate position at t = time on Hermite curve

	// Return result
	return newPosition;
}

// Implement Catmull-Rom curve
Point Curve::useCatmullCurve(const unsigned int nextPoint, const float time)
{
	Point newPosition;
	float normalTime, intervalTime;


	CurvePoint &P1 = controlPoints[nextPoint - 1];
	CurvePoint &P2 = controlPoints[nextPoint];

	float s = (time - P1.time) / (P2.time - P1.time);
	float h1 = 2 * pow(s, 3) - 3 * pow(s, 2) + 1;
	float h2 = -2 * pow(s, 3) + 3 * pow(s, 2);
	float h3 = pow(s, 3) - 2 * pow(s, 2) + s;
	float h4 = pow(s, 3) - pow(s, 2);
	Vector T1, T2;

	if (controlPoints.size() == 2)
	{
		T1 = (P2.position - P1.position)/(P2.time - P1.time);
		T2 = T1;
	}
	else if (nextPoint == 1)
	{	
		CurvePoint &P3 = controlPoints[nextPoint + 1];
		T1 = (P2.position - P1.position)/(P2.time - P1.time)*(P3.time - P1.time)/(P3.time - P2.time) - (P3.position - P1.position)/(P3.time - P1.time)*(P2.time-P1.time)/(P3.time - P2.time);
		T2 = (P3.position - P2.position)/(P3.time - P2.time)*(P2.time - P1.time)/(P3.time - P1.time) + (P2.position - P1.position)/(P2.time - P1.time)*(P3.time-P2.time)/(P3.time - P1.time);
	}

	else if (nextPoint == controlPoints.size() - 1)
	{
		
		CurvePoint &P0 = controlPoints[nextPoint - 2];
		T1 = (P2.position - P1.position)/(P2.time - P1.time)*(P1.time - P0.time)/(P2.time - P0.time) + (P1.position - P0.position)/(P1.time - P0.time)*(P2.time-P1.time)/(P2.time - P0.time);
		T2 = (P2.position - P1.position)/(P2.time - P1.time)*(P2.time - P0.time)/(P1.time - P0.time) - (P2.position - P0.position)/(P2.time - P0.time)*(P1.time-P0.time)/(P1.time - P0.time); 
	}
	else {
		CurvePoint &P0 = controlPoints[nextPoint - 2];
		CurvePoint &P3 = controlPoints[nextPoint + 1];
		T1 = (P2.position - P1.position)/(P2.time - P1.time)*(P1.time - P0.time)/(P2.time - P0.time) + (P1.position - P0.position)/(P1.time - P0.time)*(P2.time-P1.time)/(P2.time - P0.time);
		T2 = (P3.position - P2.position)/(P3.time - P2.time)*(P2.time - P1.time)/(P3.time - P1.time) + (P2.position - P1.position)/(P2.time - P1.time)*(P3.time-P2.time)/(P3.time - P1.time);
	}

	newPosition = h1*P1.position + h2*P2.position + h3*T1*(P2.time - P1.time) + h4*T2*(P2.time - P1.time);






	//================DELETE THIS PART AND THEN START CODING===================
	/*	static bool flag = false;
	if (!flag)
	{
	std::cerr << "ERROR>>>>Member function useCatmullCurve is not implemented!" << std::endl;
	flag = true;
	}*/
	//=========================================================================


	// Calculate time interval, and normal time required for later curve calculations

	// Calculate position at t = time on Catmull-Rom curve

	// Return result
	return newPosition;
}
