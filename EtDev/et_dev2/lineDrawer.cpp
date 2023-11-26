// Renderer for lines
// 
// Copyright (C) 2018 Yajie Yan <danielyan86129@hotmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "lineDrawer.h"

LineDrawer::LineDrawer(
	std::shared_ptr<oglplus::Program> _prog, 
	std::shared_ptr<TrackBall> _track_ball
	)
	: track_ball(_track_ball)
	, m_readyToDraw(false)
	, s(1.0f)
	, m_lineWidth(1.0f)
	, m_zNear(1.0f)
	, m_zFar(10.0f)
	, m_zfightOffset(0.0f)
	, m_useAlpha(0)
{
	using namespace oglplus;
	try {
		//cout << "LineDrawer constructing... " << endl; // debug
		m_linesProg = _prog;

		// enable point size
		gl.Enable(Capability::Blend);
		//cout << "LineDrawer constructed." << endl; // 
	}
	catch(oglplus::ProgramBuildError& pbe)
	{
		std::cerr <<
			"Program build error (in " <<
			pbe.GLSymbol() << ", " <<
			pbe.ClassName() << " '" <<
			pbe.ObjectDescription() << "'): " <<
			pbe.what() << std::endl <<
			pbe.Log() << std::endl;
		pbe.Cleanup();
	}
	catch(oglplus::Error& err)
	{
		std::cerr <<
			"GL error (in " << err.GLSymbol() << ", " <<
			err.ClassName() << ": '" <<
			err.ObjectDescription() << "'): " <<
			err.what() <<
			" [" << err.File() << ":" << err.Line() << "] ";
		std::cerr << std::endl;
		err.Cleanup();
	}
	catch(const std::exception& se)
	{
		std::cerr <<
			"General error: " <<
			se.what() << std::endl;
	}
}

LineDrawer::~LineDrawer()
{

}

void LineDrawer::reset()
{
	m_elemCount = 0;
	m_zfightOffset = 0.0f;

	m_vtsVBO.Bind(Buffer::Target::Array);
	m_vtsVBO.InvalidateData();
	m_vtsVBO.Bind(Buffer::Target::Array);
	m_colorVBO.InvalidateData();
	m_vtsVBO.Bind(Buffer::Target::Array);
	m_saliencyVBO.InvalidateData();
	m_vtsVBO.Bind(Buffer::Target::ElementArray);
	m_linesIBO.InvalidateData();
}

void LineDrawer::setScale(const float _s)
{
	this->s = _s;
}

void LineDrawer::setPoints(const vector<TriPoint>& _pts)
{
	m_nPoints = _pts.size();
	try {
		//cout << "LineDrawer::setPoints() setting up points... " << endl;
		// make current VAO that following affects
		m_linesProg->Use();
		m_linesVAO.Bind();

		// then bind them
		//cout << "allocating mem... ";
		// allocate data for points visualization
		GLfloat * pts_vert_data = new GLfloat[_pts.size() * 3];
		//cout << "Done! " << endl;
		
		for (unsigned i = 0; i < _pts.size(); ++i)
		{
			TriPoint p = _pts[i];
			pts_vert_data[3 * i + 0] = p[0];
			pts_vert_data[3 * i + 1] = p[1];
			pts_vert_data[3 * i + 2] = p[2];
		}

		// setup inputs for line drawer shader program
		m_linesProg->Use();
		m_linesVAO.Bind();
		//m_linesIBO.InvalidateData();

		// cout << "uploading line vts data (in Bytes): " << _pts.size()*3*sizeof(GLfloat) << endl;
		m_vtsVBO.Bind(Buffer::Target::Array);
		//m_vtsVBO.InvalidateData();
		Buffer::Data(Buffer::Target::Array, _pts.size()*3, pts_vert_data);
		VertexAttribArray vert_attr(*m_linesProg, "Position");
		vert_attr.Setup<GLfloat>(3);
		vert_attr.Enable();

		m_vtsVBO.Unbind(Buffer::Target::Array);
		m_linesVAO.Unbind();

		//cout << "freeing mem... ";
		delete [] pts_vert_data;
		//cout << "Done!" << endl;

		m_readyToDraw = true;

		//cout << "LineDrawer::setPoints() done." << endl; // debug
	}
	catch(oglplus::Error& err)
	{
		std::cerr <<
			"GL error (in " << err.GLSymbol() << ", " <<
			err.ClassName() << ": '" <<
			err.ObjectDescription() << "'): " <<
			err.what() <<
			" [" << err.File() << ":" << err.Line() << "] ";
		std::cerr << std::endl;
		err.Cleanup();
	}
	catch(const std::exception& se)
	{
		std::cerr <<
			"General error: " <<
			se.what() << std::endl;
	}
}

void LineDrawer::setPerVertColor(float* _color_data, int _n)
{
	try {
		m_linesVAO.Bind();
		m_colorVBO.Bind(Buffer::Target::Array);
		//m_colorVBO.InvalidateData();
		Buffer::Data(Buffer::Target::Array, _n*3, _color_data);

		m_linesProg->Use();
		VertexAttribArray attr_color(*m_linesProg, "Color");
		attr_color.Setup<GLfloat>(3);
		attr_color.Enable();
		m_colorVBO.Unbind(Buffer::Target::Array);
		m_linesVAO.Unbind();
	}
	catch(oglplus::Error& err)
	{
		std::cerr <<
			"GL error (in " << err.GLSymbol() << ", " <<
			err.ClassName() << ": '" <<
			err.ObjectDescription() << "'): " <<
			err.what() <<
			" [" << err.File() << ":" << err.Line() << "] ";
		std::cerr << std::endl;
		err.Cleanup();
	}
	catch(const std::exception& se)
	{
		std::cerr <<
			"General error: " <<
			se.what() << std::endl;
	}
}

void LineDrawer::setTransparencyEnabled(const bool _use_alpha)
{
	m_useAlpha = _use_alpha ? 1u : 0u;
}

bool LineDrawer::getTransparencyEnabled()
{
	return m_useAlpha == 1u ? true : false;
}

void LineDrawer::setPerVertSaliency(float* _scalars, int _n)
{
	try
	{
		m_linesProg->Use();

		m_linesVAO.Bind();
		m_saliencyVBO.Bind(Buffer::Target::Array);

		Buffer::Data(Buffer::Target::Array, _n, _scalars);
		VertexAttribArray saliency(*m_linesProg, "Saliency");
		saliency.Setup<GLfloat>(1);
		saliency.Enable();

		m_saliencyVBO.Unbind(Buffer::Target::Array);
		m_linesVAO.Unbind();
	}
	catch(oglplus::Error& err)
	{
		std::cerr <<
			"GL error (in " << err.GLSymbol() << ", " <<
			err.ClassName() << ": '" <<
			err.ObjectDescription() << "'): " <<
			err.what() <<
			" [" << err.File() << ":" << err.Line() << "] ";
		std::cerr << std::endl;
		err.Cleanup();
	}
	catch(const std::exception& se)
	{
		std::cerr <<
			"General error: " <<
			se.what() << std::endl;
	}
}

void LineDrawer::setLines(
	const vector<TriEdge>& _edges)
{
	m_nLines = _edges.size();
	try {
		// cout << "LineDrawer::setLines() setting up lines... " << endl; // debug

		m_linesProg->Use();
		m_linesVAO.Bind();

		GLuint* line_indices = new GLuint[m_nLines * 2];
		for (unsigned i = 0; i < m_nLines; ++i)
		{
			TriEdge e = _edges[i];
			line_indices[2*i+0] = e[0];
			line_indices[2*i+1] = e[1];
		}

		// upload data
		// cout << "uploading line connectivity data (in Bytes): " << m_nLines*2*sizeof(GLuint) << endl;
		m_linesIBO.Bind(Buffer::Target::ElementArray);
		//m_linesIBO.InvalidateData();
		Buffer::Data(Buffer::Target::ElementArray, m_nLines*2, line_indices);

		delete [] line_indices;

		m_linesIBO.Unbind(Buffer::Target::ElementArray);
		m_linesVAO.Unbind();

		// cout << "LineDrawer::setLines() done." << endl; // debug
	}
	catch(oglplus::Error& err)
	{
		std::cerr <<
			"GL error (in " << err.GLSymbol() << ", " <<
			err.ClassName() << ": '" <<
			err.ObjectDescription() << "'): " <<
			err.what() <<
			" [" << err.File() << ":" << err.Line() << "] ";
		std::cerr << std::endl;
		err.Cleanup();
	}
	catch(const std::exception& se)
	{
		std::cerr <<
			"General error: " <<
			se.what() << std::endl;
	}

	m_lineType = oglplus::PrimitiveType::Lines;
	m_elemCount = m_nLines * 2;
}

void LineDrawer::setLines(const vector<unsigned>& _edge_indices)
{
	this->m_nLines = _edge_indices.size();
	try {
		// cout << "LineDrawer::setLines() setting up connectivity of lines... " << endl; // debug

		m_linesProg->Use();
		m_linesVAO.Bind();

		GLuint* line_indices = new GLuint[m_nLines * 2];
		for (unsigned i = 0; i < m_nLines; ++i)
		{
			unsigned ei = _edge_indices[i];
			line_indices[2*i+0] = ei*2 + 0;
			line_indices[2*i+1] = ei*2 + 1;
		}

		// upload data
		// cout << "uploading line connectivity data (in Bytes): " << m_nLines*2*sizeof(GLuint) << endl;
		m_linesIBO.Bind(Buffer::Target::ElementArray);
		//m_linesIBO.InvalidateData();
		Buffer::Data(Buffer::Target::ElementArray, m_nLines*2, line_indices);

		delete [] line_indices;

		m_linesIBO.Unbind(Buffer::Target::ElementArray);
		m_linesVAO.Unbind();

		// cout << "LineDrawer::setLines() done." << endl; // debug
	}
	catch(oglplus::Error& err)
	{
		std::cerr <<
			"GL error (in " << err.GLSymbol() << ", " <<
			err.ClassName() << ": '" <<
			err.ObjectDescription() << "'): " <<
			err.what() <<
			" [" << err.File() << ":" << err.Line() << "] ";
		std::cerr << std::endl;
		err.Cleanup();
	}
	catch(const std::exception& se)
	{
		std::cerr <<
			"General error: " <<
			se.what() << std::endl;
	}

	this->m_lineType = PrimitiveType::Lines;
	this->m_elemCount = m_nLines * 2;
}

void LineDrawer::setLinesAdjacency(const vector<unsigned>& _edge_w_adjacency)
{
	try
	{
		cout << "LineDrawer::setLinesAdjacency() ... " << endl; // debug

		m_linesProg->Use();
		m_linesVAO.Bind();

		m_linesIBO.Bind(Buffer::Target::ElementArray);
		Buffer::Data(Buffer::Target::ElementArray, _edge_w_adjacency, BufferUsage::DynamicDraw);
		
		m_linesIBO.Unbind(Buffer::Target::ElementArray);
		m_linesVAO.Unbind();

		this->m_nLines = _edge_w_adjacency.size() / 4;
		this->m_lineType = PrimitiveType::LinesAdjacency;
		this->m_elemCount = _edge_w_adjacency.size();

		cout << "done." << endl;
	}
	catch(oglplus::Error& err)
	{
		std::cerr <<
			"GL error (in " << err.GLSymbol() << ", " <<
			err.ClassName() << ": '" <<
			err.ObjectDescription() << "'): " <<
			err.what() <<
			" [" << err.File() << ":" << err.Line() << "] ";
		std::cerr << std::endl;
		err.Cleanup();
	}
	catch(const std::exception& se)
	{
		std::cerr <<
			"General error: " <<
			se.what() << std::endl;
	}
}

void LineDrawer::setZFightOffset(float _r)
{
	m_zfightOffset = _r*s;
}

void LineDrawer::reshape(GLuint _w, GLuint _h)
{
	track_ball->setViewport(_w, _h);

	using namespace oglplus;

	gl.Viewport(_w, _h);

	m_linesProg->Use();
	m_zNear = (0.01f)*s + m_zfightOffset;
	m_zFar = 1000.0f*s;
	m_perspectiveProjMat = CameraMatrix<float>::PerspectiveX(
		Degrees(60), double(_w)/_h, m_zNear, m_zFar
		);
	LazyUniform<Mat4f> proj_mat(*m_linesProg, "ProjMat");
	proj_mat.Set(m_perspectiveProjMat);

	//cout << "LineDrawer::reshape() done." << endl;
}


CameraMatrix<float> LineDrawer::getPerspectiveProjMat() const
{
	return m_perspectiveProjMat;
}
void LineDrawer::render(double _time)
{
	if (!m_readyToDraw)
		return;

	using namespace oglplus;

	try {
		//LazyUniform<Vec3f>* light_pos;

		//draw lines
		m_linesProg->Use();
		LazyUniform<Mat4f> cam_mat(*m_linesProg, "CamMat");
		LazyUniform<Mat4f> model_mat(*m_linesProg, "ModelMat");
		model_mat.Set(track_ball->getModelMatrix());

		// apply tiny offset to avoid z-fighting between lines and other geometry, e.g. faces.
		//auto offset = ModelMatrixf::Translation(0.0f, 0.0f, -0.01f*s);
		cam_mat.Set(/*offset * */track_ball->getViewMatrix() * track_ball->getCamMatrix());
		LazyUniform<GLuint> use_alpha(*m_linesProg, "UseAlpha");
		use_alpha.Set(m_useAlpha);

		m_linesVAO.Bind();
		m_linesIBO.Bind(Buffer::Target::ElementArray);
		/*gl.DrawArrays(PrimitiveType::Points, 0, m_nPoints * 3);*/
		gl.DrawElements(m_lineType, m_elemCount, DataType::UnsignedInt);
		m_linesIBO.Unbind(Buffer::Target::ElementArray);
		m_linesVAO.Unbind();

		//cout << "m_linesProg in Render." << endl;
	}
	catch(oglplus::Error& err)
	{
		std::cerr <<
			"GL error (in " << err.GLSymbol() << ", " <<
			err.ClassName() << ": '" <<
			err.ObjectDescription() << "'): " <<
			err.what() <<
			" [" << err.File() << ":" << err.Line() << "] ";
		std::cerr << std::endl;
		err.Cleanup();
	}
	catch(const std::exception& se)
	{
		std::cerr <<
			"General error: " <<
			se.what() << std::endl;
	}
}

void LineDrawer::renderAsProxy(double _time)
{
	if (!m_readyToDraw)
		return;
	
	//cout << "render lines as proxy ... " << endl;

	using namespace oglplus;

	try {
		//LazyUniform<Vec3f>* light_pos;

		//draw lines
		m_linesProg->Use();
		LazyUniform<Mat4f> cam_mat(*m_linesProg, "CamMat");
		cam_mat.Set(track_ball->getViewMatrix()*track_ball->getCamMatrix());
		LazyUniform<Mat4f> model_mat(*m_linesProg, "ModelMat");
		model_mat.Set(track_ball->getModelMatrix());

		m_linesVAO.Bind();
		m_linesIBO.Bind(Buffer::Target::ElementArray);
		/*gl.DrawArrays(PrimitiveType::Points, 0, m_nPoints * 3);*/
		gl.DrawElements(/*PrimitiveType::LinesAdjacency*/PrimitiveType::Lines, m_elemCount, DataType::UnsignedInt);
		m_linesIBO.Unbind(Buffer::Target::ElementArray);
		m_linesVAO.Unbind();

		//cout << "m_linesProg in Render." << endl;
	}
	catch(oglplus::Error& err)
	{
		std::cerr <<
			"GL error (in " << err.GLSymbol() << ", " <<
			err.ClassName() << ": '" <<
			err.ObjectDescription() << "'): " <<
			err.what() <<
			" [" << err.File() << ":" << err.Line() << "] ";
		std::cerr << std::endl;
		err.Cleanup();
	}
	catch(const std::exception& se)
	{
		std::cerr <<
			"General error: " <<
			se.what() << std::endl;
	}
}

void LineDrawer::clear()
{
	m_vtsVBO.InvalidateData();
	m_linesIBO.InvalidateData();
}

void LineDrawer::MousePress(QMouseEvent* _evnt)
{
	//cout << "mouse pressed: " << _evnt->x() << "," << track_ball->h()-_evnt->y() << endl;//debug

	if (_evnt->button() == Qt::LeftButton)
		track_ball->beginRotating(Vec2f(_evnt->x(), track_ball->h() - _evnt->y()));
	else if (_evnt->button() == Qt::RightButton)
		track_ball->beginPanning(Vec2f(_evnt->x(), track_ball->h() - _evnt->y()), 0.003f * s);
	else if (_evnt->button() == Qt::MiddleButton)
		track_ball->beginZooming(Vec2f(_evnt->x(), track_ball->h() - _evnt->y()), 0.003f * s);
}

void LineDrawer::MouseMove(GLuint _x, GLuint _y, GLuint _w, GLuint _h)
{
	//cout << "mouse moving: " << _x << "," << _y << endl;//debug

	track_ball->update(Vec2f(_x, _y));
}

void LineDrawer::MouseRelease(QMouseEvent* _evnt)
{
	// TODO
	track_ball->stopTracking();
}

const Mat4f& LineDrawer::getCurProjMat() const
{
	return m_perspectiveProjMat;
}

const Vec2f LineDrawer::getNearFar() const
{
	return Vec2f(m_zNear, m_zFar);
}

std::shared_ptr<TrackBall> LineDrawer::getCamera()
{
	return track_ball;
}