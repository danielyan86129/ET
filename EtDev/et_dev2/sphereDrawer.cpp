#include <iostream>
#include <vector>

#include "sphereDrawer.h"
#include "glArea.h"

#include <oglplus/shapes/sphere.hpp>

#define X .525731112119133606 
#define Z .850650808352039932

// Icosahedron vertex and indices
static GLfloat vdata[12][3] = {    
	{-X, 0.0, Z}, {X, 0.0, Z}, {-X, 0.0, -Z}, {X, 0.0, -Z},    
	{0.0, Z, X}, {0.0, Z, -X}, {0.0, -Z, X}, {0.0, -Z, -X},    
	{Z, X, 0.0}, {-Z, X, 0.0}, {Z, -X, 0.0}, {-Z, -X, 0.0} 
};
static GLuint tindices[20][3] = { 
	{0,4,1}, {0,9,4}, {9,5,4}, {4,5,8}, {4,8,1},    
	{8,10,1}, {8,3,10}, {5,3,8}, {5,2,3}, {2,7,3},    
	{7,10,3}, {7,6,10}, {7,11,6}, {11,0,6}, {0,1,6}, 
	{6,1,10}, {9,0,11}, {9,11,2}, {9,2,5}, {7,2,11} 
};

void createSphereGeometry(float _r, int _subdiv_d, vector<GLfloat>& _pos, vector<GLfloat>& _nrml, vector<GLuint>& _indices)
{
	for (size_t i = 0; i < 12; ++i)
	{
		_pos.push_back(vdata[i][0]);
		_pos.push_back(vdata[i][1]);
		_pos.push_back(vdata[i][2]);

		_nrml.push_back(vdata[i][0]);
		_nrml.push_back(vdata[i][1]);
		_nrml.push_back(vdata[i][2]);
	}
	
	int d_cnt = _subdiv_d;
	//1. subdivide the default icosahedron to the depth _subdiv_d
	queue<std::pair<TriFace, int> > Q;
	for (size_t i = 0; i < 20; ++i)
	{
		Q.push( make_pair(TriFace(tindices[i]), d_cnt) );
	}
	while (!Q.empty())
	{
		auto pr = Q.front();
		auto f = pr.first;
		auto d_cnt = pr.second;
		Q.pop();

		if (d_cnt == 0)
		{
			_indices.push_back(f[0]);
			_indices.push_back(f[1]);
			_indices.push_back(f[2]);
			continue;
		}

		// subdivide this face f
		auto v0 = TriPoint(vdata[f[0]]);
		auto v1 = TriPoint(vdata[f[1]]);
		auto v2 = TriPoint(vdata[f[2]]);
		int v0i = f[0], v1i = f[1], v2i = f[2];
		int v01i, v12i, v02i;
		auto v01 = trimesh::mix(v0, v1, 0.5f);
		trimesh::normalize(v01);
		auto v12 = trimesh::mix(v1, v2, 0.5f);
		trimesh::normalize(v12);
		auto v02 = trimesh::mix(v0, v2, 0.5f);
		trimesh::normalize(v02);

		// add new vertices to pos list & normal list
		_nrml.push_back(v01[0]);
		_nrml.push_back(v01[1]);
		_nrml.push_back(v01[2]);
		_pos.push_back(v01[0]);
		_pos.push_back(v01[1]);
		_pos.push_back(v01[2]);
		v01i = (_pos.size() / 3) - 1;
		_nrml.push_back(v12[0]);
		_nrml.push_back(v12[1]);
		_nrml.push_back(v12[2]);
		_pos.push_back(v12[0]);
		_pos.push_back(v12[1]);
		_pos.push_back(v12[2]);
		v12i = (_pos.size() / 3) - 1;
		_nrml.push_back(v02[0]);
		_nrml.push_back(v02[1]);
		_nrml.push_back(v02[2]);
		_pos.push_back(v02[0]);
		_pos.push_back(v02[1]);
		_pos.push_back(v02[2]);
		v02i = (_pos.size() / 3) - 1;

		// add subdivided faces to indices list and Q
		Q.push(make_pair(TriFace(v0i, v01i, v02i), d_cnt-1));
		Q.push(make_pair(TriFace(v1i, v01i, v12i), d_cnt-1));
		Q.push(make_pair(TriFace(v2i, v02i, v12i), d_cnt-1));
		Q.push(make_pair(TriFace(v01i, v12i, v02i), d_cnt-1));

	} // end of while()

	// 2. scale the vertices 
	for (auto it = _pos.begin(); it != _pos.end(); ++it)
		*it *= _r;
}

SphereDrawer::SphereDrawer(
	/*shared_ptr<oglplus::Program> _prog, */
	shared_ptr<TrackBall> _track_ball)
	: m_trackball (_track_ball), 
	m_readyToDraw(false), 
	m_s(1.0f)
{
	//cout << "SphereDrawer constructing... " <<endl;

	using namespace oglplus;
	try
	{
		m_sphereProg = GLArea::glslProgram(
			"./shaders/DiffuseOnly2.vert.glsl",
			"./shaders/DiffuseOnly2.frag.glsl"
			);

		/*setCenter(TriPoint(0.0f, 0.0f, 0.0f));*/
		setRadius(1.0f);
		m_readyToDraw = true;
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

	//cout << "SphereDrawer constructed." <<endl;
}

SphereDrawer::~SphereDrawer()
{

}

void SphereDrawer::reset()
{
	m_readyToDraw = false;
	//m_nFaceToDraw = 0;
	// reset the model matrix
	m_modelMat = oglplus::ModelMatrix<GLfloat>();
	m_radius = 1.0f;
	m_center = TriPoint();
}

void SphereDrawer::setScale(const float _s)
{
	m_s = _s;
}

// update model matrix with new center
void SphereDrawer::setCenter(const TriPoint& _p)
{
	m_center = _p;
	m_modelMat = 
		ModelMatrix<GLfloat>::Translation(m_center[0], m_center[1], m_center[2]) * 
		ModelMatrix<GLfloat>::Scale(m_radius, m_radius, m_radius);
	m_readyToDraw = true;
}

void SphereDrawer::setCenter(const vector<TriPoint>& _pts, float _r)
{
	// construct sphere & upload data to GPU
	vector<GLfloat> pos_data;
	vector<GLuint> indices_data;
	vector<GLfloat> nrml_data;
	vector<GLfloat> one_sphere_pos;
	vector<GLuint> one_sphere_indices;
	vector<GLfloat> one_sphere_nrml;
	createSphereGeometry(_r, 1, one_sphere_pos, one_sphere_nrml, one_sphere_indices);

	// bind VAO
	m_sphereVAO.Bind();

	// setup verts VBO, normal VBO, and index IBO
	
	{
		//GLuint n_per_vert = make_sphere.Positions(one_sphere);
		GLuint n_per_vert = 3;
		
		for (size_t i = 0; i < _pts.size(); ++i)
		{
			const auto& p = _pts[i];
			for (size_t j = 0; j < one_sphere_pos.size(); j += n_per_vert)
			{
				pos_data.push_back( one_sphere_pos[j+0] + p[0] );
				pos_data.push_back( one_sphere_pos[j+1] + p[1] );
				pos_data.push_back( one_sphere_pos[j+2] + p[2] );
			}
			nrml_data.insert(nrml_data.end(), one_sphere_nrml.begin(), one_sphere_nrml.end());
			for (auto it = one_sphere_indices.begin(); it != one_sphere_indices.end(); ++it)
			{
				auto index = *it + i * one_sphere_pos.size() / n_per_vert;
				indices_data.push_back(index);
			}
		}

		m_sphereVertVBO.Bind(Buffer::Target::Array);
		Buffer::Data(Buffer::Target::Array, pos_data);
		VertexAttribArray vert_attr(*m_sphereProg, "Position");
		vert_attr.Setup<GLfloat>(n_per_vert);
		vert_attr.Enable();
		pos_data.clear();
		m_sphereVertVBO.Unbind(Buffer::Target::Array);

		m_sphereNrmlVBO.Bind(Buffer::Target::Array);
		Buffer::Data(Buffer::Target::Array, nrml_data);
		VertexAttribArray nrml_attr(*m_sphereProg, "Normal");
		nrml_attr.Setup<GLfloat>(n_per_vert);
		nrml_attr.Enable();
		nrml_data.clear();
		m_sphereNrmlVBO.Unbind(Buffer::Target::Array);

		m_sphereIBO.Bind(Buffer::Target::ElementArray);
		Buffer::Data(Buffer::Target::ElementArray, indices_data);
		m_nFaceToDraw = indices_data.size() / 3;
		indices_data.clear();
		m_sphereIBO.Unbind(Buffer::Target::ElementArray);
	}

	// setup color (uniform color for spheres)
	m_sphereProg->Use();
	LazyUniform<oglplus::Vec3f> uniform_color(*m_sphereProg, "Color");
	uniform_color.Set(oglplus::Vec3f(1.0f, 0.0f, 0.0f));
	Program::UseNone();

	// unbind VAO
	m_sphereVAO.Unbind();
	m_readyToDraw = true;

	//cout << "setCenter(pts) done." << endl;
}

// effectively update the model view matrix by scaling by _r
void SphereDrawer::setRadius(float _r)
{
	m_radius = _r;
	m_modelMat = 
		/*ModelMatrix<GLfloat>::Translation(m_center[0], m_center[1], m_center[2]) * */
		ModelMatrix<GLfloat>::Scale(m_radius, m_radius, m_radius);
	m_readyToDraw = true;
}

void SphereDrawer::setColor(const TriColor& _c)
{
	LazyUniform<oglplus::Vec3f> uniform_color(*m_sphereProg, "Color");
	uniform_color.Set(oglplus::Vec3f(1.0f, 0.0f, 0.0f));
}

void SphereDrawer::setColor(float* _color_data, int _n)
{
	// TODO: support specifying different color for different spheres
}

void SphereDrawer::setPerVertColor(float* _color_data, int _n)
{
}

void SphereDrawer::reshape(GLuint _w, GLuint _h)
{
	m_trackball->setViewport(_w, _h);

	gl.Viewport(_w, _h);

	m_sphereProg->Use();
	LazyUniform<Mat4f> proj_mat(*m_sphereProg, "ProjMat");
	proj_mat.Set(
		CameraMatrix<float>::PerspectiveX(
		Degrees(60), double(_w)/_h, 0.01f*m_s, 1000.0f*m_s
		)
		);
}

// currently only render one sphere. 
// TODO: add support to multiple spheres
void SphereDrawer::render(double _time)
{
	if (!m_readyToDraw)
		return;

	try {
		LazyUniform<Mat4f>* cam_mat;
		LazyUniform<Mat4f>* model_mat;
		//LazyUniform<Vec3f>* light_pos;

		//draw points
		m_sphereProg->Use();
		cam_mat = new LazyUniform<Mat4f>(*m_sphereProg, "CamMat");
		model_mat = new LazyUniform<Mat4f>(*m_sphereProg, "ModelMat");
		model_mat->Set(m_modelMat);
		cam_mat->Set(m_trackball->getViewMatrix()*m_trackball->getCamMatrix()*m_trackball->getModelMatrix());
		/*model_mat->Set(ModelMatrix<GLfloat>());
		cam_mat->Set(m_trackball->getViewMatrix()*m_trackball->getCamMatrix());*/
		LazyUniform<oglplus::Vec3f> light_pos(*m_sphereProg, "LightPosition");
		light_pos.Set(oglplus::Vec3f(0.0f, 0.0f, 50.0f*this->m_s));

		m_sphereVAO.Bind();
		m_sphereIBO.Bind(Buffer::Target::ElementArray);
		//gl.DrawArrays(PrimitiveType::Triangles, 0, m_nFaceToDraw * 3);
		gl.DrawElements(PrimitiveType::Triangles, m_nFaceToDraw * 3, DataType::UnsignedInt);
		m_sphereIBO.Unbind(Buffer::Target::ElementArray);
		m_sphereVAO.Unbind();

		/*m_sphereVAO.Bind();
		m_sphereInstr.Draw(m_sphereIndicesData);
		m_sphereVAO.Unbind();*/

		delete cam_mat;
		delete model_mat;

		//cout << "SphereDrawer::render() called. "/*<<m_nFaceToDraw<<" faces drawn."*/<<endl;
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

void SphereDrawer::MousePress(QMouseEvent* _evnt)
{
	//cout << "mouse pressed: " << _evnt->x() << "," << track_ball->h()-_evnt->y() << endl;//debug

	if (_evnt->button() == Qt::LeftButton)
		m_trackball->beginRotating(Vec2f(_evnt->x(), m_trackball->h() - _evnt->y()));
	else if (_evnt->button() == Qt::RightButton)
		m_trackball->beginPanning(Vec2f(_evnt->x(), m_trackball->h() - _evnt->y()), 0.003f * m_s);
	else if (_evnt->button() == Qt::MiddleButton)
		m_trackball->beginZooming(Vec2f(_evnt->x(), m_trackball->h() - _evnt->y()), 0.003f * m_s);
}

void SphereDrawer::MouseMove(GLuint _x, GLuint _y, GLuint _w, GLuint _h)
{
	//cout << "mouse moving: " << _x << "," << _y << endl;//debug

	m_trackball->update(Vec2f(_x, _y));
}

void SphereDrawer::MouseRelease(QMouseEvent* _evnt)
{
	// TODO
	m_trackball->stopTracking();
}