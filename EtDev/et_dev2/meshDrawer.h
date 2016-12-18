#ifndef MESH_DRAWER_H
#define MESH_DRAWER_H

#include "all_gl.h"

#include <QMouseEvent>
#include <algorithm>

#include "allTypes.h"
#include "cameraTrackBall.h"
#include "drawable.h"

class MeshDrawer : public Drawable
{
public:
	enum RenderMode {PER_FACE, PER_VERT, SCREEN_PICKING};
public:
	MeshDrawer(
		GLuint _w, GLuint _h, 
		shared_ptr<oglplus::Program> _prog, 
		shared_ptr<oglplus::Program> _diffEdge_prog,
		shared_ptr<TrackBall> _track_ball = nullptr);
	MeshDrawer(
		GLuint _w, GLuint _h, 
		shared_ptr<oglplus::Program> _prog, 
		shared_ptr<oglplus::Program> _diffEdge_prog,
		shared_ptr<oglplus::Program> _svp_prog, /*the screen vert picking program*/
		shared_ptr<TrackBall> _track_ball = nullptr);
// 	MeshDrawer(
// 		shared_ptr<oglplus::Program> _prog, 
// 		shared_ptr<oglplus::Program> _diffEdge_prog, 
// 		shared_ptr<TrackBall> _track_ball);
	~MeshDrawer();

	void reset();

	// get/set lighting enabled or not
	void setLightingEnabled(bool _enable);
	bool getLightingEnabled();

	// regulate the transparency of the mesh 
	// by adding the given value to the current alpha used in rendering
	void setTransparencyEnabled(const bool _use_alpha);
	bool getTransparencyEnabled();

	void setScale(const float _s);
	/// position camera based on bbox of mesh.
	/// (use only when you have the mesh)
	void initCameraPosAndLens(void);
	// per-face rendering or per-vert rendering?
	RenderMode getRenderMode();
	void setRenderMode(RenderMode _mode);
	/// prepare for future rendering the given mesh object
	void setMeshToDraw(std::shared_ptr<TriMesh> _m);
	// only set the coodinates of the vertices
	void setPoints(const vector<TriPoint>& _pts);
	// set the faces to draw using indices into 
	// the face list of the currently bound mesh
	void setFaces(const vector<unsigned>& _face_indices);
	// set the faces to draw using indices into
	// the just uploaded faces. 
	// per face attributes are used.
	void setFaces(const vector<unsigned>& _face_indices, int _faces_must_be_already_uploaded);
	// set the given list of tri faces as faces to draw 
	// each face is a list of vertex indices 
	// into the already set-up vts list ( e.g. setPoints() ).
	// use per face attributes
	void setFaces(const vector<TriFace>& _faces, bool _use_per_face_attrib);
	// assign a value for each vert / face. vts/faces must be uploaded to buffer already.
	void setPerVertColor(float* _color_data, int _n);
	void setPerVertNormal(float* _normal_data, int _n);
	void setPerVertSaliency(float* _scalars, int _n);
	// set a constant saliency value for all vertices
	void setPerVertConstantSaliency(float _scalar, int _n);
	// set coordinates for vts in a per face manner, 
	// i.e. shared vts will be dup.ed
	void setPointsPerFace(float* _dup_vts_coord, int _n);
	void setPointsPerFace(const vector<TriPoint>& _vts, const vector<TriFace>& _faces);
	// upload per face attribute (data should be dup.ed properly). faces must be uploaded to buffer already.
	void setPerFaceColor(float* _color_data, int _n_face);
	void setPerFaceNormal(float* _normal_data, int _n_face);
	void setPerFaceSaliency(float* _scalars, int _n);

	/// screen vert picking SVP related 
	// set points with associated data to svp related buffer
	void setPointsSVP(const vector<TriPoint>& _pts, const vector<unsigned>& _aux_data, const char* _name);
	
	std::shared_ptr<TriMesh> getMesh();

	void reshape(GLuint _w, GLuint _h);
	void render(double _time);

	void MousePress(QMouseEvent* _evnt)
	{
		//cout << "mouse pressed: " << _evnt->x() << "," << track_ball->h()-_evnt->y() << endl;//debug

		if (_evnt->button() == Qt::LeftButton)
			track_ball->beginRotating(Vec2f(_evnt->x(), track_ball->h() - _evnt->y()));
		else if (_evnt->button() == Qt::RightButton)
			track_ball->beginPanning(Vec2f(_evnt->x(), track_ball->h() - _evnt->y()), 0.003f * s);
		else if (_evnt->button() == Qt::MiddleButton)
			track_ball->beginZooming(Vec2f(_evnt->x(), track_ball->h() - _evnt->y()), 0.003f * s);
	}

	void MouseMove(GLuint _x, GLuint _y, GLuint _w, GLuint _h)
	{
		//cout << "mouse moving: " << _x << "," << _y << endl;//debug

		track_ball->update(Vec2f(_x, _y));
	}

	void MouseRelease(QMouseEvent* _evnt)
	{
		// TODO
		track_ball->stopTracking();
	}

	// want to visualize edge on mesh?
	void setDrawEdge(bool draw_or_not)
	{
		this->m_drawEdge = draw_or_not;
	}

	std::shared_ptr<TrackBall> getCamera(void);

	/// helpers
private:

	/*void setupGLProgs();*/

private:
	// current gl context
	oglplus::Context gl;
	// shader related
	shared_ptr<oglplus::Program> m_simpleProg;
	oglplus::VertexShader m_simpleVS;
	oglplus::FragmentShader m_simpleFS;
	shared_ptr<oglplus::Program> m_diffuseEdgeProg; // diffuse with edge shader program
	
	// camera settings
	std::shared_ptr<TrackBall> track_ball;
	// gl-data
	oglplus::VertexArray m_meshVAO;
	oglplus::VertexArray m_meshDupVAO;
	oglplus::Buffer m_vertsVBO;
	oglplus::Buffer m_vertsDupVBO;
	oglplus::Buffer m_normVBO;
	oglplus::Buffer m_normDupVBO;
	oglplus::Buffer m_colorVBO;
	oglplus::Buffer m_colorDupVBO;
	oglplus::Buffer m_saliencyVBO;
	oglplus::Buffer m_saliencyDupVBO;
	oglplus::Buffer m_indices;
	oglplus::Buffer m_indicesDup;

	/// svp related
	shared_ptr<oglplus::Program> m_SVPProg; 
	oglplus::VertexArray m_meshVAO_SVP;
	oglplus::Buffer m_vertsVBO_SVP;
	oglplus::Buffer m_auxVBO_SVP;
	oglplus::Buffer m_indices_SVP;

	// model to draw ()
	std::shared_ptr<TriMesh> m;
	// flags & related parameters
	bool m_readyToDraw;
	bool m_drawEdge;
	float m_drawPointSize;
	unsigned m_nFaces2Draw;
	RenderMode m_renderMode;

	// geometry scale info
	float s;

	/* rendering related info*/
	// the transparency used in rendering the mesh
	GLuint m_use_alpha;
	GLuint m_lightingEnabled;
};

#endif