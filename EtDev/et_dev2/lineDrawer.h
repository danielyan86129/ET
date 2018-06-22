// Renderer for lines
// 
// Copyright (C) 2018 Yajie Yan <danielyan86129@hotmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef LINE_DRAWER_H
#define LINE_DRAWER_H

#include "all_gl.h"

#include <QMouseEvent>
#include <algorithm>
#include <cstdint>

#include "allTypes.h"
#include "cameraTrackBall.h"
#include "drawable.h"

using std::shared_ptr;

class LineDrawer : public Drawable
{
public:
	/// _prog: the glsl program this drawer will use. 
	/// _track_ball: the trackball camera this drawer will use to position geometry.
	LineDrawer(shared_ptr<oglplus::Program> _prog, shared_ptr<TrackBall> _track_ball);
	~LineDrawer();

	// reset the line drawer
	void reset();
	// set the geometry related scale. 
	// useful when the scale of the geometry is needed within this drawer
	void setScale(const float _s);
	/// the coordinates of edges' vertices
	void setPoints(const vector<TriPoint>& _pts);
	/// the indices of points each line connects
	void setLines(
		const vector<TriEdge>& _edges);
	/// set the given subset of lines to draw. 
	/// only work correctly when all lines data reside in GPU already
	void setLines(const vector<unsigned>& _edge_indices);
	/// specify line-adjacency info as connectivity
	void setLinesAdjacency(
		const vector<unsigned>& _edge_w_adjacency);
	/// color of the vts of lines
	void setPerVertColor(float* _color_data, int _n);
	/// set/get the transparency enabling status
	void setTransparencyEnabled(const bool _use_alpha);
	bool getTransparencyEnabled();
	/// saliency data for vert (i.e. 2 for each edge)
	void setPerVertSaliency(float* _scalars, int _n);

	/// set the offset value used for resolving z-fighting
	void setZFightOffset(float _r);
	void reshape(GLuint _w, GLuint _h);
	void render(double _time);
	/// render lines as cylinder imposters. 
	/// Note: only call it when current drawer is setup with line-adj data
	/// and the geometry proxy shader program!
	void renderAsProxy(double _time);

	void clear();

	void MousePress(QMouseEvent* _evnt);
	void MouseMove(GLuint _x, GLuint _y, GLuint _w, GLuint _h);
	void MouseRelease(QMouseEvent* _evnt);

	const Mat4f& getCurProjMat() const;
	const Vec2f getNearFar() const;
	std::shared_ptr<TrackBall> getCamera();

private:
	// current gl context
	oglplus::Context gl;
	// shader related
	std::shared_ptr<oglplus::Program> m_linesProg; // draw-points shader program
	
	// camera settings
	std::shared_ptr<TrackBall> track_ball;

	// vts & line connectivity data to draw
	oglplus::VertexArray m_linesVAO;
	oglplus::Buffer m_vtsVBO;
	oglplus::Buffer m_colorVBO;
	oglplus::Buffer m_saliencyVBO;
	oglplus::Buffer m_linesIBO;

	// flags & related parameters/status
	bool m_readyToDraw;
	unsigned m_nPoints;
	unsigned m_nLines;
	float m_lineWidth;
	GLuint m_useAlpha;

	// keep track of the projection matrix
	CameraMatrix<float> m_perspectiveProjMat;
	float m_zNear, m_zFar;
	float m_zfightOffset;

	// the scale of the geometry
	float s;

	// the type of the line connectivity stored in m_linesIBO
	oglplus::PrimitiveType m_lineType;
	unsigned m_elemCount;
};


#endif