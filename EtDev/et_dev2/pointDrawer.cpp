// Renderer for points
// 
// Copyright (C) 2018 Yajie Yan <danielyan86129@hotmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "pointDrawer.h"

PointDrawer::PointDrawer(
	std::shared_ptr<oglplus::Program> _prog, 
	std::shared_ptr<TrackBall> _track_ball
	)
	: track_ball(_track_ball)
	, m_readyToDraw(false)
	, s(1.0f)
{
	//cout << "PointDrawer constructing... " << endl; // debug

	using namespace oglplus;
	try {
		/*m_drawPointsVS.Source(
			loadShaderSource("./shaders/DrawPoints.vert.glsl")
			).Compile();
		m_drawPointsFS.Source(
			loadShaderSource("./shaders/DrawPoints.frag.glsl")
			).Compile();
		m_drawPointsProg.AttachShaders(Group<Shader>(
			m_drawPointsVS, m_drawPointsFS
			)).Link();*/
		m_drawPointsProg = _prog;

		Program::UseNone();

		// enable point size
		// gl.Enable(Capability::ProgramPointSize);
		glEnable(GL_PROGRAM_POINT_SIZE);
		glEnable(GL_POINT_SMOOTH);
		setPointSize(1.0f);
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

	//cout << "PointDrawer constructed." << endl; // 
}

PointDrawer::~PointDrawer()
{

}

void PointDrawer::reset()
{
	m_nPoints = 0;
}

void PointDrawer::setScale(const float _s)
{
	this->s = _s;
}

void PointDrawer::setPoints(const vector<TriPoint>& _pts)
{
	m_nPoints = _pts.size();
	try {
		//cout << "PointDrawer::setPoints() setting up points... " << endl;
		// make current VAO that following affects
		m_drawPointsProg->Use();
		m_pointsVAO.Bind();

		// then bind them
		//cout << "allocating mem... " << endl;
		// allocate data for points visualization
		GLfloat * pts_vert_data = new GLfloat[_pts.size() * 3];
		GLfloat * pts_color_data = new GLfloat[_pts.size() * 3];
		//cout << "allocating mem Done! " << endl;

		// data for points visualization
		//memcpy(pts_vert_data, verts_data, m->vertices.size()*3*sizeof(GLfloat));
		for (unsigned i = 0; i < _pts.size(); ++i)
		{
			TriPoint p = _pts[i];
			pts_vert_data[3 * i + 0] = p[0];
			pts_vert_data[3 * i + 1] = p[1];
			pts_vert_data[3 * i + 2] = p[2];

			// default color
			TriColor color(0.0f, 0.0f, 1.0f);
			pts_color_data[3*i+0] = color[0];
			pts_color_data[3*i+1] = color[1];
			pts_color_data[3*i+2] = color[2];
		}

		// setup inputs for draw-points shader program
		m_drawPointsProg->Use();
		m_pointsVAO.Bind();

		m_ptsVertVBO.Bind(Buffer::Target::Array);
		Buffer::Data(Buffer::Target::Array, _pts.size()*3, pts_vert_data);
		VertexAttribArray vert_attr3(*m_drawPointsProg, "Position");
		vert_attr3.Setup<GLfloat>(3);
		vert_attr3.Enable();

		m_ptsColorVBO.Bind(Buffer::Target::Array);
		Buffer::Data(Buffer::Target::Array, _pts.size()*3, pts_color_data);
		VertexAttribArray attr_color3(*m_drawPointsProg, "Color");
		attr_color3.Setup<GLfloat>(3);
		attr_color3.Enable();

		m_ptsVertVBO.Unbind(Buffer::Target::Array);
		m_ptsColorVBO.Unbind(Buffer::Target::Array);
		m_pointsVAO.Unbind();

		//cout << "freeing mem..." << endl;

		delete [] pts_vert_data;
		delete [] pts_color_data;
		//cout << "freeing mem Done!" << endl;

		// setup camera
		//initCameraPosAndLens();
		//track_ball->reset();

		m_readyToDraw = true;

		//cout << "PointDrawer::setPoints() done." << endl; // debug
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

void PointDrawer::setPointSize(float size)
{
	this->m_drawPointSize = size;
}

void PointDrawer::setPerVertColor(float* _color_data, int _n_vts)
{
	m_pointsVAO.Bind();
	m_ptsColorVBO.Bind(Buffer::Target::Array);
	Buffer::Data(Buffer::Target::Array, _n_vts*3, _color_data);

	m_drawPointsProg->Use();
	VertexAttribArray attr_color1(*m_drawPointsProg, "Color");
	attr_color1.Setup<GLfloat>(3);
	attr_color1.Enable();
	m_ptsColorVBO.Unbind(Buffer::Target::Array);
	m_pointsVAO.Unbind();
}

void PointDrawer::setShaderSaliency(float* _scalars)
{
	try
	{
		LazyUniform<GLfloat>* saliency;
		m_drawPointsProg->Use();
		saliency = new LazyUniform<GLfloat>(*m_drawPointsProg, "Saliency");
		saliency->Set(*_scalars);
	}
	catch (oglplus::Error& err)
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
	catch (const std::exception& se)
	{
		std::cerr <<
			"General error: " <<
			se.what() << std::endl;
	}
}

void PointDrawer::reshape(GLuint _w, GLuint _h)
{
	track_ball->setViewport(_w, _h);

	using namespace oglplus;

	gl.Viewport(_w, _h);

	m_drawPointsProg->Use();
	LazyUniform<Mat4f> proj_mat(*m_drawPointsProg, "ProjMat");
	proj_mat.Set(
		CameraMatrix<float>::PerspectiveX(
		Degrees(60), double(_w)/_h, (0.01f+0.00001f)*s, 1000.0f*s
		)
		);

	//cout << "m_drawPointsProg in reshape()." << endl;

}

void PointDrawer::render(double _time)
{
	if (!m_readyToDraw)
		return;
	using namespace oglplus;
	try {
		
		LazyUniform<Mat4f>* cam_mat;
		LazyUniform<Mat4f>* model_mat;
		//LazyUniform<Vec3f>* light_pos;
		//draw points
		m_drawPointsProg->Use();
		cam_mat = new LazyUniform<Mat4f>(*m_drawPointsProg, "CamMat");
		model_mat = new LazyUniform<Mat4f>(*m_drawPointsProg, "ModelMat");
		model_mat->Set(track_ball->getModelMatrix());
		cam_mat->Set(track_ball->getViewMatrix()*track_ball->getCamMatrix());

		LazyUniform<GLfloat> pt_size(*m_drawPointsProg, "PointSize");
		pt_size.Set(this->m_drawPointSize);

		m_pointsVAO.Bind();
		gl.DrawArrays(PrimitiveType::Points, 0, m_nPoints);
		m_pointsVAO.Unbind();

		delete cam_mat;
		delete model_mat;
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

void PointDrawer::clear()
{
	m_ptsVertVBO.InvalidateData();
}

void PointDrawer::MousePress(QMouseEvent* _evnt)
{
	//cout << "mouse pressed: " << _evnt->x() << "," << track_ball->h()-_evnt->y() << endl;//debug

	if (_evnt->button() == Qt::LeftButton)
		track_ball->beginRotating(Vec2f(_evnt->x(), track_ball->h() - _evnt->y()));
	else if (_evnt->button() == Qt::RightButton)
		track_ball->beginPanning(Vec2f(_evnt->x(), track_ball->h() - _evnt->y()), 0.003f * s);
	else if (_evnt->button() == Qt::MiddleButton)
		track_ball->beginZooming(Vec2f(_evnt->x(), track_ball->h() - _evnt->y()), 0.003f * s);
}

void PointDrawer::MouseMove(GLuint _x, GLuint _y, GLuint _w, GLuint _h)
{
	//cout << "mouse moving: " << _x << "," << _y << endl;//debug

	track_ball->update(Vec2f(_x, _y));
}

void PointDrawer::MouseRelease(QMouseEvent* _evnt)
{
	// TODO
	track_ball->stopTracking();
}