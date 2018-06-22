// Renderer for spheres
// 
// Copyright (C) 2018 Yajie Yan <danielyan86129@hotmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.

/************************************************************************/
/* Draw sphere. Currently only supports one sphere at a time.           */
/************************************************************************/

#ifndef SPHERE_DRAWER_H
#define SPHERE_DRAWER_H

#include "all_gl.h"
#include <QMouseEvent>
#include <algorithm>
#include <oglplus/shapes/sphere.hpp>
#include "drawable.h"
#include "cameraTrackBall.h"
#include "allTypes.h"
using namespace oglplus;

void createSphereGeometry(float _r, int _subdiv_d, vector<GLfloat>& _pos, vector<GLfloat>& _nrml, vector<GLuint>& _indices);

class SphereDrawer : public Drawable
{
public:
	// by default create a unit sphere centered at origin.
	// can be moved around: setCenter(), setRadius()
	SphereDrawer(/*shared_ptr<oglplus::Program> _prog,*/ shared_ptr<TrackBall> _track_ball);
	~SphereDrawer();

	void reset();
	void setScale(const float _s);

	// set the center	
	void setCenter(const TriPoint& _p);
	void setCenter(const vector<TriPoint>& _pts, float _r);
	void setRadius(float _r);
	// setup a uniform color for all spheres
	void setColor(const TriColor& _c);
	void setColor(float* _color_data, int _n);
	void setPerVertColor(float* _color_data, int _n);

	void reshape(GLuint _w, GLuint _h);
	void render(double _time);

	void MousePress(QMouseEvent* _evnt);
	void MouseMove(GLuint _x, GLuint _y, GLuint _w, GLuint _h);
	void MouseRelease(QMouseEvent* _evnt);

private:
	// current gl context
	oglplus::Context gl;

	// shader related
	shared_ptr<oglplus::Program> m_sphereProg; // draw-points shader program

	// camera settings
	std::shared_ptr<TrackBall> m_trackball;
	// gl-data
	VertexArray m_sphereVAO;
	Buffer m_sphereIBO;
	Buffer m_sphereVertVBO;
	Buffer m_sphereNrmlVBO;
	Buffer m_sphereColorVBO;
	ModelMatrix<float> m_modelMat;

	// model to draw ()

	// flags & related parameters
	bool m_readyToDraw;
	unsigned m_nFaceToDraw;
	TriPoint m_center;
	float m_radius;

	// geometry scale info
	float m_s;
};

#endif