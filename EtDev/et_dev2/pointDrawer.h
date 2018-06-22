// Renderer for points
// 
// Copyright (C) 2018 Yajie Yan <danielyan86129@hotmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef POINT_DRAWER_H
#define POINT_DRAWER_H

#include "all_gl.h"

#include <QMouseEvent>
#include <algorithm>

#include "allTypes.h"
#include "cameraTrackBall.h"
#include "drawable.h"

class PointDrawer : public Drawable
{
public:
	PointDrawer(shared_ptr<oglplus::Program> _prog, shared_ptr<TrackBall> _track_ball);
	~PointDrawer();

	void reset();
	void setScale(const float _s);
	/// prepare for future rendering the given mesh object
	void setPoints(const vector<TriPoint>& _pts);
	void setPointSize(float size);
	void setPerVertColor(float* _color_data, int _n);

	void reshape(GLuint _w, GLuint _h);
	void render(double _time);

	void clear();

	void MousePress(QMouseEvent* _evnt);

	void MouseMove(GLuint _x, GLuint _y, GLuint _w, GLuint _h);
	void MouseRelease(QMouseEvent* _evnt);

private:
	// current gl context
	oglplus::Context gl;

	// shader related
	shared_ptr<oglplus::Program> m_drawPointsProg; // draw-points shader program

	// camera settings
	std::shared_ptr<TrackBall> track_ball;
	// gl-data
	oglplus::VertexArray m_pointsVAO;
	oglplus::Buffer m_ptsVertVBO;
	oglplus::Buffer m_ptsColorVBO;

	// model to draw ()

	// flags & related parameters
	bool m_readyToDraw;
	unsigned m_nPoints;
	float m_drawPointSize;

	// geometry scale info
	float s;
};

#endif