/*******************************************************************************
 File:				TrackBall.h

 Author: 			Gaspard Petit (gaspardpetit@gmail.com)
					Modified by Yajie.

// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
*******************************************************************************/
#ifndef CAMERA_TRACKBALL_H
#define CAMERA_TRACKBALL_H

//==============================================================================
//	TrackBall.h
//
//	This class implements a 3D trackball that can be used to convert mouse
//	drags on a 2D plane into 3D rotations.  
//
//	How to use:
//
//	//	w, h are current viewport dimensions
//	TrackBall tBall(radius, Point2(xCenter, yCenter));
//
//	tBall.beginTracking(clickPoint);
//	...
//	Quatf qrot = tBall.track(dragPoint);
//
//	//	here, qrot is a quaternion containing the rotation, to convert to
//	//	angles as a rotation matrix:
//
//	myRot = qrot*myRot;	//	initially, myRot would be Quatf(1,0,0,0)
//	Matrix4f rot = myRot.as4x4Matrix();
//	
//==============================================================================

#include <oglplus/all.hpp>

#include <iostream>
#include <algorithm>
#include <QSettings>
#include <memory>
#include "glUtils.h"

using std::cout;
using std::endl;

using namespace oglplus;
class TrackBall
{
public:
	TrackBall()
		: mCenter(Vec2f(0.0f, 0.0f))
		, mRadius(1.0f)
		, is_tracking(false)
		, rotate_mat(Mat4f()), translate_mat(Mat4f()), cam_mat(Mat4f())
	  {}
	TrackBall(GLuint _w, GLuint _h)
		: mCenter(_w/2, _h/2)
		, mRadius(std::min(_w, _h)/2)
		, is_tracking(false)
		, rotate_mat(Mat4f()), translate_mat(Mat4f()), cam_mat(Mat4f())
	{}
	TrackBall(const TrackBall &rhs)
		: mCenter(rhs.mCenter)
		, mRadius(rhs.mRadius)
		, is_tracking(false)
		, rotate_mat(Mat4f()), translate_mat(Mat4f()), cam_mat(Mat4f())
	{}

	TrackBall& operator = (const TrackBall &rhs)
	{
		mCenter = rhs.mCenter;
		mRadius = rhs.mRadius;
		return *this;
	}

	/// save/load current/previous trackball info to/from the given qsetting ptr
	bool saveToSetting(std::shared_ptr<QSettings>& _qsetting)
	{
		QVariantList mat4x4;
		for (int i = 0; i < 16; ++i)
			mat4x4 << cam_mat.Data()[ i ];
		_qsetting->setValue("TrackBall/cam_mat", mat4x4);
		mat4x4.clear();
		for (int i = 0; i < 16; ++i)
			mat4x4 << rotate_mat.Data()[ i ];
		_qsetting->setValue("TrackBall/rotate_mat", mat4x4);
		mat4x4.clear();
		for (int i = 0; i < 16; ++i)
			mat4x4 << translate_mat.Data()[ i ];
		_qsetting->setValue("TrackBall/translate_mat", mat4x4);

		return true;
	}
	bool loadFromSetting(std::shared_ptr<QSettings>& _qsetting)
	{
		auto temp = _qsetting->value("TrackBall/cam_mat").toList();
		assert(temp.size() == 16);

		for (int i = 0; i < temp.size(); ++i)
		{
			cam_mat.Set( i / 4, i % 4, temp[i].toFloat() );
		}
		temp = _qsetting->value("TrackBall/rotate_mat").toList();
		for (int i = 0; i < temp.size(); ++i)
		{
			rotate_mat.Set( i / 4, i % 4, temp[i].toFloat() );
		}
		temp = _qsetting->value("TrackBall/translate_mat").toList();
		for (int i = 0; i < temp.size(); ++i)
		{
			translate_mat.Set( i / 4, i % 4, temp[i].toFloat() );
		}

		return true;
	}

	/// reset to initial camera position
	void reset()
	{
		this->rotate_mat = Mat4f();
		this->translate_mat = Mat4f();
	}

	// need to set/update viewport info whenever gl window reshapes
	void setViewport(int _w, int _h)
	{
		mCenter = Vec2f(_w/2, _h/2);
		mRadius = std::min(_w, _h)/2;
	}

	// set up initial camera matrix
	void setCamera(Vec3f _eye, Vec3f _lookat, Vec3f _up)
	{
		cam_mat = CameraMatrix<float>::LookingAt(_eye, _lookat, _up);
	}

	void stopTracking()
	{
		is_tracking = false;
		is_panning = false;
		is_zooming = false;
	}

	void beginRotating(const Vec2f& pt)
	{
		/*cout << "begin rotating" << endl;*/
		is_tracking = true;

		mLastPos = pt;
		mLastSphPos = Vec3f(mLastPos.x(), mLastPos.y(), mRadius);
	}
	void updateRotating(const Vec2f& pt)
	{
		/*cout << "rotating... " <<pt.x()<<","<<pt.y()<< endl;*/
		Vec2f currPos = pt;
		Vec3f currSphPos = Vec3f(currPos.x(), currPos.y(), mRadius);
		Vec3f centerSph = Vec3f(mCenter.x(), mCenter.y(), 0.0f);
		Vec3f currVec = currSphPos-centerSph; currVec.Normalize();
		Vec3f lastVec = mLastSphPos-centerSph; lastVec.Normalize();

		oglplus::Angle<float> angle(ArcCos( 
			std::max( std::min(Vec3f::DotProduct(currVec, lastVec), 0.9999f ), -0.9999f)
			) );
		Vec3f axis = crossProduct(lastVec, currVec);
		axis.Normalize();

		mLastPos = currPos;
		mLastSphPos = currSphPos;

		rotate_mat = ModelMatrixf::RotationA(axis, angle) * rotate_mat;
	}

	void beginPanning(const Vec2f& pt, const float _unit)
	{
		is_panning = true;

		mLastPos = pt;
		this->mTransUnit = _unit;
	}void updatePanning(const Vec2f& pt)
	{
		Vec2f dxy = pt - mLastPos;
		mLastPos = pt;
		translate_mat = ModelMatrixf::Translation(
			dxy.x()*mTransUnit, 
			dxy.y()*mTransUnit,
			0.0f) * translate_mat;
	}
	
	void beginZooming(const Vec2f& pt, const float _unit)
	{
		is_zooming = true;

		mLastPos = pt;
		this->mTransUnit = _unit;
	}
	void updateZooming(const Vec2f& pt)
	{
		Vec2f dxy = pt - mLastPos;
		mLastPos = pt;
		translate_mat = ModelMatrixf::Translation(
			0.0f, 
			0.0f,
			-dxy.y()*mTransUnit) * translate_mat;
	}
	void update(const Vec2f& pt)
	{
		if (is_tracking)
			updateRotating(pt);
		if (is_panning)
			updatePanning(pt);
		if (is_zooming)
			updateZooming(pt);
	}
	
	// the rotation part of the track ball
	Mat4f getModelMatrix()
	{
		return rotate_mat;
	}
	// the translate part of the track ball
	Mat4f getViewMatrix()
	{
		return translate_mat;
	}
	Mat4f getCamMatrix()
	{
		return cam_mat;
	}
	Mat4f getModelViewMatrix()
	{
		return translate_mat * cam_mat * rotate_mat;
	}
	
	GLuint h()
	{
		return mCenter.y()*2;
	}
public:
	float mTransUnit;
	Vec2f mLastPos;
	Vec3f mLastSphPos;
	Vec2f mCenter;
	float mRadius;

	Mat4f cam_mat;
	Mat4f rotate_mat; // in cam space
	Mat4f translate_mat; // in cam space
/*
private:
	Vec3f proj*/

private:
	bool is_tracking;
	bool is_panning;
	bool is_zooming;
};



#endif // include guard