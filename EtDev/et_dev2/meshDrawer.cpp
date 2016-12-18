#include "meshDrawer.h"

MeshDrawer::MeshDrawer(
	GLuint _w, GLuint _h, 
	std::shared_ptr<oglplus::Program> _simple_prog, 
	std::shared_ptr<oglplus::Program> _diffEdge_prog,
	std::shared_ptr<TrackBall> _track_ball)
{
	//cout << "MeshDrawer constructing... " << endl; // debug
	reset();

	using namespace oglplus;

	track_ball = 
		_track_ball ? _track_ball : std::shared_ptr<TrackBall> (new TrackBall(_w, _h));

	// Vertex shader
	try {
		m_simpleProg = _simple_prog;

		m_diffuseEdgeProg = _diffEdge_prog;

		// set uniforms related to diffuseEdgeProg
		m_diffuseEdgeProg->Use();
		LazyUniform<GLfloat> edge_width(*m_diffuseEdgeProg, "EdgeWidth");
		edge_width.Set(0.5);
		LazyUniform<Vec3f> edge_color(*m_diffuseEdgeProg, "EdgeColor");
		edge_color.Set(Vec3f(0.0f, 0.0f, 0.0f));

		Program::UseNone();

		// initially using simple shader program
		setDrawEdge(false);

		// enable point size
		gl.Enable(Capability::ProgramPointSize);
		// enable lighting
		setLightingEnabled(true);

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

	//cout << "MeshDrawer constructed." << endl; // debug
}

MeshDrawer::MeshDrawer(
	GLuint _w, GLuint _h, 
	std::shared_ptr<oglplus::Program> _simple_prog, 
	std::shared_ptr<oglplus::Program> _diffEdge_prog,
	std::shared_ptr<oglplus::Program> _svp_prog,
	std::shared_ptr<TrackBall> _track_ball)
	: m_readyToDraw(false), m_renderMode(PER_VERT), m_use_alpha(0)
{
	//cout << "MeshDrawer constructing... " << endl; // debug

	using namespace oglplus;

	track_ball = 
		_track_ball ? _track_ball : std::shared_ptr<TrackBall> (new TrackBall(_w, _h));

	// Vertex shader
	try {
		m_simpleProg = _simple_prog;

		m_diffuseEdgeProg = _diffEdge_prog;
		
		// set uniforms related to diffuseEdgeProg
		m_diffuseEdgeProg->Use();
		LazyUniform<GLfloat> edge_width(*m_diffuseEdgeProg, "EdgeWidth");
		edge_width.Set(0.5);
		LazyUniform<Vec3f> edge_color(*m_diffuseEdgeProg, "EdgeColor");
		edge_color.Set(Vec3f(0.0f, 0.0f, 0.0f));

		Program::UseNone();

		// initially using simple shader program
		setDrawEdge(false);

		// enable point size
		gl.Enable(Capability::ProgramPointSize);
		// enable lighting
		setLightingEnabled(true);

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

	//cout << "MeshDrawer constructed." << endl; // debug
}

MeshDrawer::~MeshDrawer()
{
}

void MeshDrawer::reset()
{
	m_nFaces2Draw = 0;
	/*m_meshVAO;
	m_meshDupVAO;*/
	m_vertsVBO.Bind(Buffer::Target::Array);
	//Buffer::Data(Buffer::Target::Array, _vts.size()*3, nullptr);
	m_vertsVBO.InvalidateData();
	m_vertsDupVBO.Bind(Buffer::Target::Array);
	m_vertsDupVBO.InvalidateData();
	m_normVBO.Bind(Buffer::Target::Array);
	m_normVBO.InvalidateData();
	m_normDupVBO.Bind(Buffer::Target::Array);
	m_normDupVBO.InvalidateData();
	m_colorVBO.Bind(Buffer::Target::Array);
	m_colorVBO.InvalidateData();
	m_colorDupVBO.Bind(Buffer::Target::Array);
	m_colorDupVBO.InvalidateData();
	m_saliencyVBO.Bind(Buffer::Target::Array);
	m_saliencyVBO.InvalidateData();
	m_saliencyDupVBO.Bind(Buffer::Target::Array);
	m_saliencyDupVBO.InvalidateData();
	m_indices.Bind(Buffer::Target::ElementArray);
	m_indices.InvalidateData();
	m_indicesDup.Bind(Buffer::Target::ElementArray);
	m_indicesDup.InvalidateData();	
}

bool MeshDrawer::getLightingEnabled()
{
	return m_lightingEnabled == 1u ? true : false;
}

void MeshDrawer::setLightingEnabled(bool _enable)
{
	m_lightingEnabled = _enable ? 1u : 0u;
}

void MeshDrawer::setTransparencyEnabled(const bool _use_alpha)
{
	m_use_alpha = _use_alpha ? 1 : 0;
}

bool MeshDrawer::getTransparencyEnabled()
{
	return m_use_alpha == 1u ? true : false;
}

void MeshDrawer::setScale(const float _s)
{
	this->s = _s;
}

void MeshDrawer::initCameraPosAndLens(void)
{
	using namespace oglplus;

	//Vec3f target = Vec3f(0.0f, 0.0f, 0.0f);
	m->need_bbox();
	trimesh::point c = m->bbox.center();
	float s = m->bbox.size().max();
	Vec3f target = Vec3f(c[0], c[1], c[2]);
	Vec3f eye = Vec3f(target.x(), target.y(), target.z()+1.5f*s);
	track_ball->reset();
	track_ball->setCamera(eye, target, Vec3f(0.0f, 1.0f, 0.0f));

	// very likely to change in TriangleDrawTest::Reshape()
	/*LazyUniform<Mat4f> proj_mat(m_simpleProg, "ProjMat");
	proj_mat.Set(CameraMatrix<float>::PerspectiveX(Degrees(45), 1.0f, 0.2f, 20.0f));*/
}

MeshDrawer::RenderMode MeshDrawer::getRenderMode()
{
	return m_renderMode;
}

void MeshDrawer::setRenderMode(RenderMode _mode)
{
	m_renderMode = _mode;
}

void MeshDrawer::setMeshToDraw(std::shared_ptr<TriMesh> _m)
{
	try {
		this->m = _m;
		this->m_nFaces2Draw = m->faces.size();

		// set up rendering states
		setTransparencyEnabled(false);

		// copy vertices to m_vertsVBO, normals to m_normVBO, connectivity to element_array
		// then bind them
		// allocate data for mesh visualization
		GLfloat * verts_data = new GLfloat[m->vertices.size() * 3];
		GLfloat * norm_data = new GLfloat[m->vertices.size() * 3];
		GLfloat * color_data = new GLfloat[m->vertices.size() * 3];
		GLfloat * saliency_data = new GLfloat[m->vertices.size()];
		GLuint * indices_data = new GLuint[m->faces.size() * 3];

		for (unsigned i = 0; i < m->vertices.size(); i ++)
		{
			trimesh::point p = m->vertices[i];
			verts_data[3 * i + 0] = p[0];
			verts_data[3 * i + 1] = p[1];
			verts_data[3 * i + 2] = p[2];

			trimesh::vec3 n = trimesh::normalize(m->normals[i]);
			norm_data[3*i+0] = n[0];
			norm_data[3*i+1] = n[1];
			norm_data[3*i+2] = n[2];

			TriColor color;
			color = TriColor(1.0f, 0.2f, 0.2f); // TODO: allow customize this mesh color
			color_data[3*i+0] = color[0];
			color_data[3*i+1] = color[1];
			color_data[3*i+2] = color[2];

			saliency_data[i] = 1.0f;
		}

		for (unsigned i = 0; i < m->faces.size(); i ++)
		{
			indices_data[3 * i + 0] = m->faces[i][0];
			indices_data[3 * i + 1] = m->faces[i][1];
			indices_data[3 * i + 2] = m->faces[i][2];
		}

		setPoints(m->vertices);
		setPerVertColor(color_data, m->vertices.size());
		setPerVertNormal(norm_data, m->vertices.size());
		setPerVertSaliency(saliency_data, m->vertices.size());
		m_indices.Bind(Buffer::Target::ElementArray);
		Buffer::Data(Buffer::Target::ElementArray, m->faces.size()*3, indices_data);

		delete [] verts_data;
		delete [] norm_data;
		delete [] color_data;
		delete [] saliency_data;
		delete [] indices_data;

		//cout << "set per vert attributes done." <<endl;

 		verts_data = new GLfloat[m->faces.size() * 3 * 3];
		norm_data = new GLfloat[m->faces.size() * 3];
		color_data = new GLfloat[m->faces.size() * 3];
		saliency_data = new GLfloat[m->faces.size()];
		indices_data = new GLuint[m->faces.size() * 3];
		
		for (unsigned fi = 0; fi < m->faces.size(); fi ++)
		{
			TriFace f = m->faces[fi];

			for (unsigned i = 0; i < 3; ++i)
			{
				trimesh::point p = m->vertices[f[i]];
				verts_data[(3*fi + i)*3 + 0] = p[0];
				verts_data[(3*fi + i)*3 + 1] = p[1];
				verts_data[(3*fi + i)*3 + 2] = p[2];
			}

			trimesh::vec3 n = trimesh::normalize(trimesh::trinorm(
				m->vertices[f[0]], m->vertices[f[1]], m->vertices[f[2]]
				));
			norm_data[3*fi + 0] = n[0];
			norm_data[3*fi + 1] = n[1];
			norm_data[3*fi + 2] = n[2];

			TriColor color;
			color = TriColor(1.0f, 0.2f, 0.2f); // TODO: allow customize this mesh color
			color_data[3*fi + 0] = color[0];
			color_data[3*fi + 1] = color[1];
			color_data[3*fi + 2] = color[2];

			saliency_data[fi] = 0.3f;
		}

		for (unsigned fi = 0; fi < m->faces.size(); fi ++)
		{
			indices_data[3 * fi + 0] = fi*3+0;
			indices_data[3 * fi + 1] = fi*3+1;
			indices_data[3 * fi + 2] = fi*3+2;
		}

		//cout << "set per face attributes ..." <<endl;
		setPointsPerFace(verts_data, m->faces.size()*3*3);
		setPerFaceColor(color_data, m->faces.size());
		setPerFaceNormal(norm_data, m->faces.size());
		setPerFaceSaliency(saliency_data, m->faces.size());
		m_indicesDup.Bind(Buffer::Target::ElementArray);
		m_indicesDup.InvalidateData();
		Buffer::Data(Buffer::Target::ElementArray, m->faces.size()*3, indices_data);

		delete [] verts_data;
		delete [] norm_data;
		delete [] color_data;
		delete [] saliency_data;
		delete [] indices_data;
		//cout << "set per face attributes done." <<endl;

		// setup camera
		initCameraPosAndLens();
		track_ball->reset();

		m_readyToDraw = true;

		//cout << "meshDrawer::setMeshToDraw() done." << endl; // debug
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

void MeshDrawer::setPoints(const vector<TriPoint>& _vts)
{
	unsigned m_nPoints = _vts.size();
	try {
		//cout << "MeshDrawer::setPoints() setting up vertices... " << endl;
		// allocate data for vertices coordinates
		GLfloat * vert_data = new GLfloat[_vts.size() * 3];

		for (unsigned i = 0; i < _vts.size(); ++i)
		{
			TriPoint p = _vts[i];
			vert_data[3 * i + 0] = p[0];
			vert_data[3 * i + 1] = p[1];
			vert_data[3 * i + 2] = p[2];
		}

		// setup inputs for line drawer shader program
		m_simpleProg->Use();
		m_meshVAO.Bind();

		m_vertsVBO.Bind(Buffer::Target::Array);
		Buffer::Data(Buffer::Target::Array, _vts.size()*3, vert_data);
		VertexAttribArray vert_attr(*m_simpleProg, "Position");
		vert_attr.Setup<GLfloat>(3);
		vert_attr.Enable();

		m_vertsVBO.Unbind(Buffer::Target::Array);

		delete [] vert_data;

		m_readyToDraw = true;

		//cout << "MeshDrawer::setPoints() done." << endl; // debug
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

void MeshDrawer::setFaces(const vector<unsigned>& _face_indices)
{
	this->m_nFaces2Draw = std::min(_face_indices.size(), m->faces.size());
	//cout << "# faces to draw: " << this->m_nFaces2Draw << endl; // debug
	GLuint* indices_data = new GLuint[m_nFaces2Draw * 3];

	if (m_renderMode == PER_FACE)
	{
		unsigned fi = 0;
		for (unsigned i = 0; i < m_nFaces2Draw; ++i)
		{
			fi = _face_indices[i];
			indices_data[i*3 + 0] = fi*3 + 0;
			indices_data[i*3 + 1] = fi*3 + 1;
			indices_data[i*3 + 2] = fi*3 + 2;
		}
	}
	else
	{
		TriFace f;
		for (unsigned i = 0; i < m_nFaces2Draw; ++i)
		{
			f = m->faces[_face_indices[i]];
			indices_data[i*3 + 0] = f[0];
			indices_data[i*3 + 1] = f[1];
			indices_data[i*3 + 2] = f[2];
		}
	}
	//cout << "faces' indices set up!" << endl; // debug

	try
	{
		if (m_renderMode == PER_FACE)
		{
			m_meshDupVAO.Bind();
			m_indicesDup.Bind(Buffer::Target::ElementArray);
			Buffer::Data(Buffer::Target::ElementArray, m_nFaces2Draw*3, indices_data);
			m_indicesDup.Unbind(Buffer::Target::ElementArray);
			m_meshDupVAO.Unbind();
		}
		else
		{
			m_meshVAO.Bind();
			m_indices.Bind(Buffer::Target::ElementArray);
			Buffer::Data(Buffer::Target::ElementArray, m_nFaces2Draw*3, indices_data);
			m_indices.Unbind(Buffer::Target::ElementArray);
			m_meshVAO.Unbind();
		}
		//cout << "indices uploaded!" << endl; // debug
		delete [] indices_data;
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

void MeshDrawer::setFaces(const vector<unsigned>& _face_indices, int _faces_must_be_already_uploaded)
{
	this->m_nFaces2Draw = _face_indices.size();
	//cout << "# faces to draw: " << this->m_nFaces2Draw << endl; // debug
	GLuint* indices_data = new GLuint[m_nFaces2Draw * 3];

	unsigned fi = 0;
	for (unsigned i = 0; i < m_nFaces2Draw; ++i)
	{
		fi = _face_indices[i];
		indices_data[i*3 + 0] = fi*3 + 0;
		indices_data[i*3 + 1] = fi*3 + 1;
		indices_data[i*3 + 2] = fi*3 + 2;
	}
	//cout << "faces' indices set up!" << endl; // debug

	try
	{
		if (m_renderMode == PER_FACE)
		{
			m_meshDupVAO.Bind();
			m_indicesDup.Bind(Buffer::Target::ElementArray);
			Buffer::Data(Buffer::Target::ElementArray, m_nFaces2Draw*3, indices_data);
			m_indicesDup.Unbind(Buffer::Target::ElementArray);
			m_meshDupVAO.Unbind();
		}
		else
		{
			m_meshVAO.Bind();
			m_indices.Bind(Buffer::Target::ElementArray);
			Buffer::Data(Buffer::Target::ElementArray, m_nFaces2Draw*3, indices_data);
			m_indices.Unbind(Buffer::Target::ElementArray);
			m_meshVAO.Unbind();
		}
		
		//cout << "indices uploaded!" << endl; // debug
		delete [] indices_data;
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

void MeshDrawer::setFaces(const vector<TriFace>& _faces, bool _use_per_face_attrib)
{
	this->m_nFaces2Draw = _faces.size();
	GLuint* indices_data = new GLuint[3*m_nFaces2Draw];

	// per face drawing
	for (unsigned i = 0; i < m_nFaces2Draw; ++i)
	{
		auto& f = _faces[i];
		indices_data[i*3 + 0] = f[0];
		indices_data[i*3 + 1] = f[1];
		indices_data[i*3 + 2] = f[2];
	}

	try
	{
		if (_use_per_face_attrib)
		{
			m_meshDupVAO.Bind();
			m_indicesDup.Bind(Buffer::Target::ElementArray);
			Buffer::Data(Buffer::Target::ElementArray, m_nFaces2Draw*3, indices_data);
			m_indicesDup.Unbind(Buffer::Target::ElementArray);
			m_meshDupVAO.Unbind();
		}
		else
		{
			m_meshVAO.Bind();
			m_indices.Bind(Buffer::Target::ElementArray);
			Buffer::Data(Buffer::Target::ElementArray, m_nFaces2Draw*3, indices_data);
			m_indices.Unbind(Buffer::Target::ElementArray);
			m_meshVAO.Unbind();
		}
		//cout << "indices uploaded!" << endl; // debug
		delete [] indices_data;
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

void MeshDrawer::setPerVertColor(float* _color_data, int _n)
{
	try{
		if (m_renderMode == PER_VERT)
		{
			m_meshVAO.Bind();
			m_colorVBO.Bind(Buffer::Target::Array);
			m_colorVBO.InvalidateData();
			Buffer::Data(Buffer::Target::Array, _n*3, _color_data);

			m_simpleProg->Use();
			VertexAttribArray attr_color1(*m_simpleProg, "Color");
			attr_color1.Setup<GLfloat>(3);
			attr_color1.Enable();
			m_colorVBO.Unbind(Buffer::Target::Array);
			m_meshVAO.Unbind();
		}
		else // attributes are dup.ed at shared vertices
		{
			m_meshDupVAO.Bind();
			m_colorDupVBO.Bind(Buffer::Target::Array);
			m_colorDupVBO.InvalidateData();
			Buffer::Data(Buffer::Target::Array, _n*3, _color_data);

			m_simpleProg->Use();
			VertexAttribArray attr_color1(*m_simpleProg, "Color");
			attr_color1.Setup<GLfloat>(3);
			attr_color1.Enable();
			m_colorDupVBO.Unbind(Buffer::Target::Array);
			m_meshDupVAO.Unbind();
		}
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

void MeshDrawer::setPerVertNormal(float* _normal_data, int _n)
{
	try
	{
		if (m_renderMode == PER_VERT)
		{
			m_meshVAO.Bind();
			m_normVBO.Bind(Buffer::Target::Array);

			m_normVBO.InvalidateData();
			Buffer::Data(Buffer::Target::Array, _n*3, _normal_data);
			m_simpleProg->Use();
			VertexAttribArray attr_nml(*m_simpleProg, "Normal");
			attr_nml.Setup<GLfloat>(3);
			attr_nml.Enable();

			m_normVBO.Unbind(Buffer::Target::Array);
			m_meshVAO.Unbind();
		}
		else
		{
			m_meshDupVAO.Bind();
			m_normDupVBO.Bind(Buffer::Target::Array);

			m_normDupVBO.InvalidateData();
			Buffer::Data(Buffer::Target::Array, _n*3, _normal_data);
			m_simpleProg->Use();
			VertexAttribArray attr_nml(*m_simpleProg, "Normal");
			attr_nml.Setup<GLfloat>(3);
			attr_nml.Enable();

			m_normDupVBO.Unbind(Buffer::Target::Array);
			m_meshDupVAO.Unbind();
		}
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

void MeshDrawer::setPerVertSaliency(float* _scalars, int _n)
{
	try
	{
		if (m_renderMode == PER_VERT)
		{
			m_meshVAO.Bind();
			m_saliencyVBO.Bind(Buffer::Target::Array);
			m_saliencyVBO.InvalidateData();
			Buffer::Data(Buffer::Target::Array, _n, _scalars);

			m_simpleProg->Use();
			VertexAttribArray attr_saliency(*m_simpleProg, "Saliency");
			attr_saliency.Setup<GLfloat>(1);
			attr_saliency.Enable();
			m_saliencyVBO.Unbind(Buffer::Target::Array);
			m_meshVAO.Unbind();
		}
		else
		{
			m_meshDupVAO.Bind();
			m_saliencyDupVBO.Bind(Buffer::Target::Array);
			m_saliencyDupVBO.InvalidateData();
			Buffer::Data(Buffer::Target::Array, _n, _scalars);

			m_simpleProg->Use();
			VertexAttribArray attr_saliency(*m_simpleProg, "Saliency");
			attr_saliency.Setup<GLfloat>(1);
			attr_saliency.Enable();
			m_saliencyDupVBO.Unbind(Buffer::Target::Array);
			m_meshDupVAO.Unbind();
		}
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

void MeshDrawer::setPerVertConstantSaliency(float _scalar, int _n)
{
	try
	{
		if (m_renderMode == PER_VERT)
		{
			m_meshVAO.Bind();
			m_saliencyVBO.Bind(Buffer::Target::Array);
			m_saliencyVBO.InvalidateData();
			float* scalars = new float[_n];
			for (int i = 0; i < _n; ++i)
				scalars[i] = _scalar;
			Buffer::Data(Buffer::Target::Array, _n, scalars);
			delete [] scalars;

			m_simpleProg->Use();
			VertexAttribArray attr_saliency(*m_simpleProg, "Saliency");
			attr_saliency.Setup<GLfloat>(1);
			attr_saliency.Enable();
			m_saliencyVBO.Unbind(Buffer::Target::Array);
			m_meshVAO.Unbind();
		}
		else
		{
			m_meshDupVAO.Bind();
			m_saliencyDupVBO.Bind(Buffer::Target::Array);
			m_saliencyDupVBO.InvalidateData();
			float* scalars = new float[_n];
			for (int i = 0; i < _n; ++i)
				scalars[i] = _scalar;
			Buffer::Data(Buffer::Target::Array, _n, scalars);
			delete [] scalars;

			m_simpleProg->Use();
			VertexAttribArray attr_saliency(*m_simpleProg, "Saliency");
			attr_saliency.Setup<GLfloat>(1);
			attr_saliency.Enable();
			m_saliencyDupVBO.Unbind(Buffer::Target::Array);
			m_meshDupVAO.Unbind();
		}
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

void MeshDrawer::setPointsPerFace(float* _dup_vts_coord, int _n)
{
	try{
		//cout << "setPerFaceVerts ..." << endl;

		m_meshDupVAO.Bind();
		m_vertsDupVBO.Bind(Buffer::Target::Array);
		m_vertsDupVBO.InvalidateData();

		Buffer::Data(Buffer::Target::Array, _n, _dup_vts_coord);

		m_simpleProg->Use();
		VertexAttribArray attr_pos(*m_simpleProg, "Position");
		attr_pos.Setup<GLfloat>(3);
		attr_pos.Enable();

		m_vertsDupVBO.Unbind(Buffer::Target::Array);
		m_meshDupVAO.Unbind();

		this->m_readyToDraw = true;
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

void MeshDrawer::setPointsPerFace(const vector<TriPoint>& _vts, const vector<TriFace>& _faces)
{
	try{
		m_nFaces2Draw = _faces.size();
		float* vts_data = new float[ _faces.size()*3*3 ];
		GLuint* indices = new GLuint[ _faces.size()*3 ];
		for (unsigned i = 0; i < _faces.size(); ++i)
		{
			auto& f = _faces[i];
			for (unsigned j = 0; j < 3; ++j)
			{
				auto& v = _vts[f[j]];
				vts_data[ (i*3 + j) * 3 + 0 ] = v[0];
				vts_data[ (i*3 + j) * 3 + 1 ] = v[1];
				vts_data[ (i*3 + j) * 3 + 2 ] = v[2];

				indices[ i*3 + j ] = i*3 + j;
			}
		}
		setPointsPerFace(vts_data, _faces.size() * 9);

		delete [] vts_data;
		vts_data = nullptr;

		m_meshDupVAO.Bind();
		m_indicesDup.Bind(Buffer::Target::ElementArray);
		Buffer::Data(Buffer::Target::ElementArray, m_nFaces2Draw * 3, indices);
		m_indicesDup.Unbind(Buffer::Target::ElementArray);
		m_meshDupVAO.Unbind();

		delete [] indices;
		indices = nullptr;
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

void MeshDrawer::setPerFaceColor(float* _color_data, int _n_face)
{
	//cout << "setPerFaceColor ..." << endl;

	unsigned n_faces = _n_face;
	float* color_data = new float[n_faces*3*3];
	for (unsigned fi = 0; fi < _n_face; ++fi)
	{
		for (unsigned vi = 0; vi < 3; ++vi)
		{
			color_data[(fi*3+vi)*3+0]=_color_data[fi*3+0];
			color_data[(fi*3+vi)*3+1]=_color_data[fi*3+1];
			color_data[(fi*3+vi)*3+2]=_color_data[fi*3+2];
		}
	}
	m_meshDupVAO.Bind();
	m_colorDupVBO.Bind(Buffer::Target::Array);
	m_colorDupVBO.InvalidateData();

	Buffer::Data(Buffer::Target::Array, 
		n_faces*3*3, color_data);
	delete [] color_data;

	m_simpleProg->Use();
	VertexAttribArray attr_color1(*m_simpleProg, "Color");
	attr_color1.Setup<GLfloat>(3);
	attr_color1.Enable();

	m_colorDupVBO.Unbind(Buffer::Target::Array);
	m_meshDupVAO.Unbind();

	//cout << "setPerFaceColor done!" << endl;
}

void MeshDrawer::setPerFaceNormal(float* _normal_data, int _n)
{
	//cout << "setPerFaceNormal ..." << endl;

	unsigned n_faces = _n;
	float* normal_data = new float[n_faces*3*3];
	for (unsigned fi = 0; fi < _n; ++fi)
	{
		for (unsigned vi = 0; vi < 3; ++vi)
		{
			normal_data[(fi*3+vi)*3+0]=_normal_data[fi*3+0];
			normal_data[(fi*3+vi)*3+1]=_normal_data[fi*3+1];
			normal_data[(fi*3+vi)*3+2]=_normal_data[fi*3+2];
		}
	}
	m_meshDupVAO.Bind();
	m_normDupVBO.Bind(Buffer::Target::Array);
	m_normDupVBO.InvalidateData();

	Buffer::Data(Buffer::Target::Array, 
		n_faces*3*3, normal_data);
	delete [] normal_data;

	m_simpleProg->Use();
	VertexAttribArray attr_nml(*m_simpleProg, "Normal");
	attr_nml.Setup<GLfloat>(3);
	attr_nml.Enable();

	m_diffuseEdgeProg->Use();
	VertexAttribArray attr_nml1(*m_diffuseEdgeProg, "Normal");
	attr_nml1.Setup<GLfloat>(3);
	attr_nml1.Enable();

	Program::UseNone();

	m_normDupVBO.Unbind(Buffer::Target::Array);
	m_meshDupVAO.Unbind();

	// cout << "setPerFaceNormal done!" << endl;
}

void MeshDrawer::setPerFaceSaliency(float* _scalars, int _n)
{
	//cout << "setting per face Saliency ..." << endl;

	unsigned n_faces = _n;
	float* saliency_data = new float[n_faces*3];
	for (unsigned fi = 0; fi < _n; ++fi)
	{
		for (unsigned vi = 0; vi < 3; ++vi)
		{
			saliency_data[fi*3+vi]=_scalars[fi];
		}
	}
	m_meshDupVAO.Bind();
	m_saliencyDupVBO.Bind(Buffer::Target::Array);
	m_saliencyDupVBO.InvalidateData();

	Buffer::Data(Buffer::Target::Array, 
		n_faces*3, saliency_data);
	delete [] saliency_data;

	m_simpleProg->Use();
	VertexAttribArray attr_saliency(*m_simpleProg, "Saliency");
	attr_saliency.Setup<GLfloat>(1);
	attr_saliency.Enable();

	m_saliencyDupVBO.Unbind(Buffer::Target::Array);
	m_meshDupVAO.Unbind();

	// cout << "setting per face Saliency done!" << endl;
}

///
/// screen vert picking related
///

void MeshDrawer::setPointsSVP(const vector<TriPoint>& _vts, const vector<unsigned>& _aux_data, const char* _name)
{
	unsigned m_nPoints = _vts.size();
	try {
		//cout << "MeshDrawer::setPointsSVP() setting up vertices... " << endl;
		// allocate data for vertices coordinates & aux data
		GLfloat * vert_data = new GLfloat[_vts.size() * 3];
		for (unsigned i = 0; i < _vts.size(); ++i)
		{
			TriPoint p = _vts[i];
			vert_data[3 * i + 0] = p[0];
			vert_data[3 * i + 1] = p[1];
			vert_data[3 * i + 2] = p[2];
		}
		GLuint * aux_data = new GLuint[_vts.size()];
		for (unsigned i = 0; i < _aux_data.size(); ++i)
		{
			unsigned datum = _aux_data[i];
			aux_data[i] = datum;	
		}

		// setup inputs for SVP program
		m_SVPProg->Use();
		m_meshVAO_SVP.Bind();

		m_vertsVBO_SVP.Bind(Buffer::Target::Array);
		Buffer::Data(Buffer::Target::Array, _vts.size()*3, vert_data);
		VertexAttribArray vert_attr(*m_SVPProg, "Position");
		vert_attr.Setup<GLfloat>(3);
		vert_attr.Enable();
		m_vertsVBO_SVP.Unbind(Buffer::Target::Array);

		delete [] vert_data;

		m_auxVBO_SVP.Bind(Buffer::Target::Array);
		Buffer::Data(Buffer::Target::Array, _aux_data.size(), aux_data);
		VertexAttribArray aux_attr(*m_SVPProg, _name);
		aux_attr.Setup<GLfloat>(1);
		aux_attr.Enable();
		m_auxVBO_SVP.Unbind(Buffer::Target::Array);

		delete [] aux_data;

		m_meshVAO_SVP.Unbind();

		m_readyToDraw = true;

		//cout << "MeshDrawer::setPointsSVP() done." << endl; // debug
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

std::shared_ptr<TriMesh> MeshDrawer::getMesh()
{
	return m;
}

void MeshDrawer::reshape(GLuint _w, GLuint _h)
{
	track_ball->setViewport(_w, _h);

	using namespace oglplus;

	gl.Viewport(_w, _h);
	// draw wire-frame?
	LazyUniform<Mat4f>* proj_mat;
	if (m_renderMode == SCREEN_PICKING)
	{
		m_SVPProg->Use();
		proj_mat = new LazyUniform<Mat4f>(*m_SVPProg, "ProjMat");
	}
	else
	{
		if (!m_drawEdge)
		{
			m_simpleProg->Use();
			proj_mat = new LazyUniform<Mat4f>(*m_simpleProg, "ProjMat");

			//cout << "m_simpleProg in Reshape." << endl;
		}
		else
		{
			// update m_diffuseEdgeProg projection & viewport
			m_diffuseEdgeProg->Use();
			proj_mat = new LazyUniform<Mat4f>(*m_diffuseEdgeProg, "ProjMat");
			LazyUniform<Vec2f> viewport_dim(*m_diffuseEdgeProg, "ViewportDim");
			viewport_dim.Set(Vec2f(_w, _h));

			//cout << "m_diffuseEdgeProg in Reshape." << endl;
		}
	}

	// adjust projection matrix
	if (m)
	{
		m->need_bbox();
		trimesh::point c = m->bbox.center();
		s = m->bbox.size().max();
	}
	proj_mat->Set(
		CameraMatrix<float>::PerspectiveX(
		Degrees(60), double(_w)/_h, 0.01f*s, 1000.0f*s
		)
		);

	delete proj_mat;
}

void MeshDrawer::render(double _time)
{
	if (!m_readyToDraw)
		return;

	//cout << "MeshDrawer::render()"<< endl;
	using namespace oglplus;

	try {
		LazyUniform<Mat4f>* cam_mat;
		LazyUniform<Mat4f>* model_mat;
		LazyUniform<Vec3f>* light_pos;
		LazyUniform<GLuint>* use_alpha;

		if (m_renderMode == SCREEN_PICKING)
		{
			/*cam_mat = new LazyUniform<Mat4f>(*m_simpleProg, "CamMat");
			model_mat = new LazyUniform<Mat4f>(*m_simpleProg, "ModelMat");
			model_mat->Set(track_ball->getModelMatrix());
			cam_mat->Set(track_ball->getViewMatrix()*track_ball->getCamMatrix());
			delete cam_mat;
			delete model_mat;

			m_meshVAO.Bind();
			m_indices.Bind(Buffer::Target::ElementArray);
			gl.DrawElements(
				PrimitiveType::Triangles, 
				m_nFaces2Draw * 3, 
				DataType::UnsignedInt
				);
			m_indices.Unbind(Buffer::Target::ElementArray);
			m_meshVAO.Unbind();*/
		}
		else
		{
			// draw wire-frame?
			if (!m_drawEdge)
			{
				m_simpleProg->Use();
				// update uniforms related to m_simpleProg
				cam_mat = new LazyUniform<Mat4f>(*m_simpleProg, "CamMat");
				model_mat = new LazyUniform<Mat4f>(*m_simpleProg, "ModelMat");
				light_pos = new LazyUniform<Vec3f>(*m_simpleProg, "LightPosition");
				use_alpha = new LazyUniform<GLuint>(*m_simpleProg, "UseAlpha");
				LazyUniform<GLuint> lighting(*m_simpleProg, "EnableLighting");
				lighting.Set(m_lightingEnabled);

				//cout << "m_simpleProg in Render." << endl;
			}
			else
			{
				// update uniforms related to m_diffuseEdgeProg
				m_diffuseEdgeProg->Use();
				cam_mat = new LazyUniform<Mat4f>(*m_diffuseEdgeProg, "CamMat");
				model_mat = new LazyUniform<Mat4f>(*m_diffuseEdgeProg, "ModelMat");
				light_pos = new LazyUniform<Vec3f>(*m_diffuseEdgeProg, "LightPosition");
				use_alpha = new LazyUniform<GLuint>(*m_diffuseEdgeProg, "UseAlpha");
				//light_pos->Set(Vec3f(10.0f, 10.0f, 7.0f));

				//cout << "m_diffuseEdgeProg in Render." << endl;
			}

			model_mat->Set(track_ball->getModelMatrix());
			cam_mat->Set(track_ball->getViewMatrix()*track_ball->getCamMatrix());
			//light_pos->Set(Vec3f(10.0f, 10.0f, 7.0f));
			light_pos->Set(Vec3f(0.0f, 0.0f, 50.0f*this->s));
			use_alpha->Set(m_use_alpha);

			delete cam_mat;
			delete model_mat;
			delete light_pos;

			if (m_renderMode == PER_FACE)
			{
				m_meshDupVAO.Bind();
				m_indicesDup.Bind(Buffer::Target::ElementArray);
				gl.DrawElements(
					PrimitiveType::Triangles, 
					m_nFaces2Draw * 3, 
					DataType::UnsignedInt
					);
				m_indicesDup.Unbind(Buffer::Target::ElementArray);
				m_meshDupVAO.Unbind();

				//cout << m_nFaces2Draw<<" faces rendered."<<endl;
			}
			else
			{
				m_meshVAO.Bind();
				m_indices.Bind(Buffer::Target::ElementArray);
				gl.DrawElements(
					PrimitiveType::Triangles, 
					m_nFaces2Draw * 3, 
					DataType::UnsignedInt
					);
				m_indices.Unbind(Buffer::Target::ElementArray);
				m_meshVAO.Unbind();
			}
		}
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

std::shared_ptr<TrackBall> MeshDrawer::getCamera()
{
	return track_ball;
}