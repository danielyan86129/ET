#ifndef GL_AREA_H
#define GL_AREA_H

/*!
 * \file glArea.h
 *
 * \author yajieyan
 * \date May 2014
 * 
 * Declaration of class GLArea, which sits between QT interface and actual low-level geometric computation, 
 * e.g. SteinerGraph.

 * Rationale behind separating shader programs from drawers: 
 * Allow dynamic configuration of the shader program(s) that a drawer uses, 
 * & allow different drawer objects to share a same set of shader copies.
 * Maybe a bad idea since different drawers might want to maintain different shader states.
 * Solution? just reset the uniforms to whatever a shader likes?
 * 
 */

#include "all_gl.h"
#include "drawable.h"
#include "meshDrawer.h"
#include "lineDrawer.h"
#include "cameraTrackBall.h"
#include "geometryRepresentation.h"
#include "hybridSkeleton.h"
#include "surfaceFunc.h"
#include "isoSurfaceFrom2Manifold.h"
#include "IsoContourOnMA.h"
#include "ui_compact.h"
#include "graphDiffusion.h"

#include <QGLWidget>
#include <QMouseEvent>

#include <memory>

class GLArea : public QGLWidget
{
	Q_OBJECT

		/*static members*/
public:
	static QGLFormat getGLFormat(void);

	// assemble glsl program from shader files
	static shared_ptr<oglplus::Program> glslProgram(string _vert, string _frag);
	static shared_ptr<oglplus::Program> glslProgram(string _vert, string _geom, string _frag);

public:
	typedef SteinerSubdivision::SubdivideScheme SubdivScheme;
	// remember to also update m_lightingFlags if this is updated
#define N_OBJECT_TYPES 2
	enum ObjectType {OrigSurf=0, MA=1, FINE_MA_FOR_ISO_CONTOUR=2, FINE_MA_FOR_PRUNING=3};
	// possible scalar fields on verts of MA
	enum VertFieldType { BT2_v, BT3_v, DIFF_v, DIFF_REL_v };
	// possible scalar fields on faces of MA
	enum FaceFieldType {
		BT2_f, BT3_f, DIFF_f, DIFF_REL_f, ANGLE_f, LAMBDA_f, GEODESIC_f, FILE_MSURE
	};
	// visualize burned edges or dual edges
	enum BurntOrDualType { VIS_BURNT = 0x1, VIS_DUAL = 0x10 };
	// draw flags:
	enum DrawFlag {
		DRAW_MA=1, DRAW_ORIG=2, DRAW_MC=4, DRAW_HS=8, DRAW_MP=16, DRAW_ISOSURF=32, DRAW_ISOCONT=64, 
		DRAW_BURNT_EDGES=128, DRAW_POINTS = 256, DRAW_MA_FINE = 512, DRAW_QMAT = 1024, DRAW_MA_LINE = 2048
	};
	// possible scalar fields on vts of MC
	enum DistMC {BT3_MC, BT2_MC, BT2_BT3_MC, BT2_BT3_REL_MC, BT1_MC, BT1_BT2_MC, BT1_BT2_REL_MC};
public:
	GLArea(QWidget *parent = 0);
	~GLArea();

	/// overload
	QSize minimumSizeHint() const;
	QSize sizeHint() const;

	/// public interface
	bool event(QEvent *event);
	void reset();
	void saveView(std::shared_ptr<QSettings>& _qsetting);
	void loadView(std::shared_ptr<QSettings>& _qsetting);
	
	//void loadRadiiFromFile(std::string _r_file);
	void passToDrawable(
		std::string _medialMesh_file, 
		std::string _origMesh_file,
		std::string _r_file,
		bool _remove_dup_faces);
	void wireFrameOnMesh(bool _draw_or_not);
	void drawSteinerPoints(bool _draw_or_not, float _size);
	void drawMCProxy(bool _draw_or_not);
	// type:
	//	1 = burnt st edges
	//	2 = dual edges
	//	other = both
	// use_burnTree:
	//	1 = st edges from burn trees
	//	other = st edges from burn paths
	bool obtainDualLineStructure(
		bool _only_unburnt, 
		int _dual_method_idx, int _poly_dual_method_idx, 
		bool _stop_burn_at_junction);
	// burn the given curve network from boundary vertices
	// return distance and pre. vert of each vertex
	bool burnCurveNetwork(bool _stop_burn_at_junction, bool _only_unburnt, bool _protect_bt2, bool _underestimate_dist);
	// visualize the dual structure(medial curve) with the specified distance type
	bool setUpMCDistMetricForPruning(DistMC _metric, float& _min, float& _max);
	bool visDualLineStructure(DistMC _dist_type, int _type);
	// simply change distance field for visualization
	// update_opt: 1->only color, 2->only saliency, 3->both
	// Pre-condition: MC must have be uploaded!
	// (updage GL)
	bool colorMCEdgeBy(
		DistMC _dist_type, float& _min, float& _max, 
		bool _use_new_range, 
		float _min_alpha, float _alpha_exp, 
		int _update_opt);
	// compute the neighrborship of medial curves
	// usually called before pruning of the curves happens
	bool precomputeForMCPruning();
	// prune & vis. medial curve after pruning
	// can be used to show pruned MC that are exposed during iso-contour object's evolution
	// (for this feature, iso-contour object must have been setup properly)
	bool pruneAndVisMedialCurve(double _t, bool _preserve_topo, bool _constrained_by_iso_cont);
	/// prune helpers:
	// return prune value _v / bt2 range
	float toETcPruneRatio(float _v);
	// prune the medial curve based on the specified metric and threshold
	bool pruneMedialCurve( double _t, bool _preserve_topo );
	// upload vts and lines data to the specified line drawer
	bool uploadLinesToDrawer(const vector<TriPoint>& _vts, const vector<TriEdge>& _edges, shared_ptr<LineDrawer>& _line_drawer) const;
	// upload vts and lines of MC geometry to gpu for rendering
	bool uploadMCCompleteGeometry(const vector<TriPoint>& _vts, const vector<TriEdge>& _edges);
	// upload vts and line with adjacency info to gpu for the m_linesProxyDrawer
	bool uploadRemainedMCtoProxyDrawer();
	// show a subset of the medial curve structure already setup in line drawer
	bool uploadLineSubset(const vector<unsigned>& _edge_indices, const vector<TriEdge>& _edges);
	// print stats of remained MC to the screen
	// only call it after MC is generated
	void printRemainedMCStats();
	// set specified draw flag & re-paint
	void setDrawFlag(DrawFlag _obj_to_draw, bool _draw);
	// draw the specified faces of the current MA mesh
	void drawMeshSubset(
		std::shared_ptr<MeshDrawer> _mesh_drawer,
		const vector<unsigned>& _face_indices
		);
	// burn the loaded mesh and update gl to visualize the burn function
	bool burn( SteinerGraph::BurnScheme _burn_sch, SubdivScheme _scheme, double _steiner_param, int _edgeWeight_idx);
	// upload finner triangulation of MA before interactive session like iso-countour begins
	void uploadFinerMAStaticGeom();
	// upload colors for finner triangulation determined by a scalar field and min/max range
	void uploadMAFinerStaticColors(vector<vector<float>>& _scalar_per_sheet, float _min, float _max);
	void uploadMAFinerStaticColors(const TriColor& _c);
	// upload dynamic part of finer-tri of MA as a result of iso-contour interactive session
	void uploadFinerMADynamicGeom( const vector<TriPoint>& _vts_for_subFineFaces );	
	// upload color for the dynamic part of the finer-tri of MA (const color or scalar field color)
	void uploadFinerMADynamicColors(
		const vector<TriPoint>& _vts_for_subFineFaces, 
		const TriColor& _const_c
		);
	void uploadFinerMADynamicColors(
		const vector<TriPoint>& _vts_for_subFineFaces, 
		const vector<float>& vals_for_subFineFaces_vts, 
		const float cur_prune_min, const float cur_prune_max
		);
	// detect any topo. artifacts by burning the mesh, 
	// also visualize the burn function. 
	// a "cleaner" mesh will be output to a file.
	bool cleanTopo( std::string& _out_file, bool& _hasunburnt );
	// precompute stuff for coloring and pruning MA using different measures
	void precomputeForMeasures();
	// get the specified face distance metric into the list
	void getFaceDistMetricMA(
		FaceFieldType _face_field, 
		bool _do_face_diffusion,	
		vector<float>& _dist_per_face);
	void getFaceDistMetricMA(
		FaceFieldType _face_field, 
		bool _do_face_diffusion,
		vector<vector<float> >& _vert_dist_per_sheet);
	// get the specified per-sheet vertex distance metric into the list
	void getPerSheetDistMetricMA(VertFieldType _vert_field, vector<vector<float>>& _dist_per_sheet);
	// get the specified vertex distance metric into the list
	void getVertDistMetricMA(VertFieldType _vert_field, vector<float>& _dist_per_vert);
	// color each face according to the diff. dist of that face
	// if _use_new_range = true, [_min, _max] is used to rescale distance value
	// else [_min, _max] will be updated to current dist range
	// update opt: 1->only color, 2->only saliency, 3->both
	void colorMAFaceBy(
		FaceFieldType _type, float& _min, float& _max, 
		bool _use_new_range, 
		float _min_alpha, float _alpha_exp, 
		int _update_opt,
		bool _do_face_diffusion,
		bool _use_per_sheet);
	// color each vert according to the specified dist type
	// update opt: 1->only color, 2->only saliency, 3->both
	void colorMAVertBy(
		VertFieldType _type, float& _min, float& _max, 
		bool _use_new_range, 
		float _min_alpha, float _alpha_exp, 
		int _update_opt);
	// get the distance value of the specified dist metric type
	void getMCDistMetric(DistMC _type, vector<float>& _dist) const;
	// _dist_by_edge : {... ei<dist. at two ends resp. w.r.t. the edge> ... }
	void getMCDistMetric(DistMC _type, vector<trimesh::vec2>& _dist_by_edge) const;
	// get ready the distance metric for MA face pruning
	// _dist_i: among two dist metrics, which to be enabled?
	void setupMAFacePrune(	FaceFieldType _type, float& _min, float& _max, int _dist_i);
	// obtain absolute prune param for MA from a given relative prune value
	float getMAPruneValue(float _ratio);
	float toETmPruneRatio(float _v);
	// prune MA faces by the given diff. dist. threshold
	// using cur dist metric values
	void pruneMA(float _t1, float _t2);
	// prune MA-finner-triangulation by the given dist. threshold
	// using cur dist metric values
	void pruneFineMA(float _t1, float _t2);
	// finner triangulate faces & upload them to the drawer
	void initFinePruneMA();
	// tessellate finner faces crossed by iso-contour of the prune iso-value properly
	// & upload those faces to the drawer
	void updateFinePruneMA(float _iso);
	// use per-face attribute to render?
	void usePerFaceRender(bool _is_perFace);
	// visualize MP
	void visMP();
	// output MC and a specified measure on it
	void outputMCwMeasure( const std::string& _meas ) const;

	/*
	app: HS related
	*/

	// create HS
	void setHSFaceDegenerateThreshold(float _t);
	void setComponentFaceNumberThresh(float _t);
	// get hs param: given ratio -> absolute threshold value
	void getHSTheta2Abs(float _ratio, float& _abs_t);
	void getHSTheta1Abs(float _ratio, float& _abs_t);
	void getHSTheta2Ratio(float _abs_t, float& _ratio);
	void getHSTheta1Ratio(float _abs_t, float& _ratio);
	void createHS();
	void uploadHS();
	void pruneHS(
		float _f_diff_r, float _f_reldiff_r, float _l_diff_r, float _l_reldiff_r, 
		bool _use_inputs_directly,
		bool _remove_small_cmpnts
		);

	/*
	application surf func related
	eps_ratio specify the ratio of a smoothing factor
	*/
	enum SurfFuncCorrespScheme {SURF_FUNC_RANGE_CORRESP, SURF_FUNC_KNN_CORRESP};
	void computeSurfFuncCorrespondence(SurfFuncCorrespScheme _scheme, float _r_eps_ratio, int _k);
	void projectFunctionOnSurf(int _idx, float _eps_ratio, bool _do_diffuse = false);

	/*
	visualize iso-surface of 3d burning, iso-contour of 2d burning
	*/
	// setup isosurface
	void setupIsoSurface();
	// refine isosurface triangulation & corresp. onto MA/MC
	void refineIsoSurface(bool _is_recursive);
	// vis iso surface of 3d burning, i.e. orig 3d shape deformed towards MA
	void uploadIsoSurface(float _t, bool _hide_snapped_faces);
	void uploadIsoSurface(float _iso, bool _hide_snapped_faces, int);

	/// setup isocontour
	void setupIsoContour(float min_alpha, float exp_alpha, const float* const _new_min, const float* const _new_max);
	/// evolve & visualize isocontour as lines
	void uploadIsoContour(float _t, bool _hide_past_faces, bool _show_exposed_mc);
	void uploadIsoContour(float _t, bool _hide_past_faces, bool _show_exposed_mc, int);
	/// visualize MC exposed by current iso-contour
	void uploadFromRemainedMC(bool _exposed, bool _prune, float _iso, bool _fronts_on_MC);
	/// determine & upload the dynamic part of MC ( intersected by iso-contour at the moment )
	/// compare at each vertex the cur. chosen metric value against the given _iso value
	void uploadDynamicMC(float _iso_bt2, bool _prune_bt1, float __iso_bt1);

	/// iso-points related
	/// MC prune function extended by iso contour constrain:
	/// only display mc edges that are above _t and exposed by iso contour already 
	bool pruneMedialCurveConstrainedByIsoContour(double _t);

	/// debug related
	void printAndVisBurnTime();
	void visSphere(size_t _vid);
	void visSphere(vector<unsigned> _vid_list);
	void printMAOriginalVerts(vector<unsigned> _vid_list);
	void printMAOriginalVert(size_t _vid);

	// get stuff related to steiner graph
	float getRangeOfDiffDist();
	float getMinDiffDist();
	float getMaxDiffDist();
	const vector<TriPoint>& getOriginalVts() const;
	
	/// get stuff related to MA mesh
	int getNumFacesOfMA() const;

	// set & update
	bool changeOrigColor(float * _color);
	void setUseConstColorMA( bool _use_constColor_for_MA, const TriColor& _c );
	bool changeConstColorMA(float * _color, bool _color_vert);
	bool changeConstColorMC(float * _color);
	bool changeBGColor(float* _color);
	bool changeOrigTransparency(float _alpha);
	bool changeMAVertTransparency(VertFieldType _type, float _min_alpha, int _exp);
	bool changeMAFaceTransparency(FaceFieldType _type, float _min_alpha, int _exp);
	bool changeMCTransparency(DistMC _type, float _min_alpha, int _exp);
	bool changeMCLineWidth(float _width);
	void turnOffLighting(enum ObjectType _obj_t, bool _turn_off);
	void setTransparencyEnabled(enum ObjectType _obj_t, bool _enabled);
	void setTrueTransparencyRender(bool _is_enabled);
	void setDPMaxRenders(int _max_iter);
	void setLinesOffset(float _r);

	/// methods we compare to
	/// qmat related
	bool outputToWenping();
	bool readQMATFile(std::string _filename);
	void uploadQMAT();

	/* methods for exporting bt, skeleton, etc. */
	// export/import burn time per-vertex to/from file
	void exportPerVertexET();
	//bool importPerVertexBT(std::string _filename);
	// export/import burn time per-sector/sheet to/from file
	void exportPerSectorET();
	//bool importPerSectorBT(std::string _filename);
	/// export skeleton to file
	void exportSkeleton();

public:
	Ui::MainWindow* ui;

protected:
	void initializeGL();
	void paintGL();
	void resizeGL(int _w, int _h);
	void mousePressEvent(QMouseEvent *_evnt);
	void mouseMoveEvent(QMouseEvent *_evnt);
	void mouseReleaseEvent(QMouseEvent *_evnt);
	void mouseDoubleClickEvent(QMouseEvent *_evnt);

private:
	/// helpers
	// initialize shader programs that different drawers use
	void initShaders();

	/// rendering related general helpers
	void bindTexture2D(Texture& _name, int _unit,
		int _w, int _h,
		TextureMinFilter _min_filter = TextureMinFilter::Nearest,
		TextureMagFilter _mag_filter = TextureMagFilter::Nearest,
		TextureWrap _s_wrap = TextureWrap::ClampToEdge,
		TextureWrap _t_wrap = TextureWrap::ClampToEdge,
		PixelDataInternalFormat _internal = PixelDataInternalFormat::RGBA,
		PixelDataFormat _data_format = PixelDataFormat::RGBA,
		PixelDataType _data_type = PixelDataType::Float,
		const void* _data = nullptr
		);

	/// OIT related helpers
	// initialize OIT-related opengl states
	void initOIT();
	// resize texture etc. when window is resized
	void resizeOIT(int _w, int _h);
	// init required textures for OIT
	void initTexturesWithSizesOIT(int _w, int _h);
	// bind required textures for shaders to access in OIT
	void bindTexturesOIT();

	/// Depth-peeling (DP) related helpers
	// initialize DP opengl states
	void initDP();
	// init required textures for OIT
	void initTexturesWithSizesDP(int _w, int _h);
	// resize texture etc. when window is resized
	void resizeDP(int _w, int _h);

	/// debug tools: screen vertex picking (SVP)
	// initialize screen based vertex picking related texture, FBO, etc.
	void initVertPicking();
	// init required texture for SVP
	void initTexturesWithSizesSVP(int _w, int _h);
	// after render finished, the texture bound to fbo 
	// will hold the vert id for each texel
	// which can be queried by gl.ReadPixel() later
	void pickingPhaseSVP();

	// upload a simplicial complex (containing edges and triangle faces) to GPU for rendering
	void uploadSimplicialComplex( 
		const vector<TriPoint>& _vts, const vector<TriEdge>& _edges, const vector<TriFace>& _faces,
		const vector<TriColor>* _e_colors, const vector<TriColor>* _f_colors,
		shared_ptr<Drawable>& _edge_drawer, shared_ptr<Drawable>& _face_drawer );

private:
	std::shared_ptr<Drawable> m_origDrawer;
	std::shared_ptr<Drawable> m_MADrawer;
	std::shared_ptr<Drawable> m_MALineDrawer;
	std::shared_ptr<Drawable> m_FinerMAStaticDrawer;
	std::shared_ptr<Drawable> m_MAFinnerDynamicDrawer;
	std::shared_ptr<Drawable> m_pointDrawer;
	std::shared_ptr<Drawable> m_dualLinesDrawer;
	std::shared_ptr<Drawable> m_dynamic_dualLineDrawer;
	std::shared_ptr<Drawable> m_linesProxyDrawer;
	std::shared_ptr<Drawable> m_burntEdgesDrawer;
	std::shared_ptr<Drawable> m_hsFaceDrawer;
	std::shared_ptr<Drawable> m_hsLineDrawer;
	std::shared_ptr<Drawable> m_MPDrawer;
	std::shared_ptr<Drawable> m_isoSurfDrawer;
	vector<std::shared_ptr<Drawable> > m_drawers;

	// the drawers for drawing qmat MA
	std::shared_ptr<Drawable> m_qmatFaceDrawer;
	std::shared_ptr<Drawable> m_qmatEdgeDrawer;

	// the drawer for the screen vertex picking tool
	std::shared_ptr<Drawable> m_drawerSVP;

	// opengl context managed by oglplus
	oglplus::Context gl;

	// name of the MA file
	std::string m_medialAxisFile;

	// the mesh structure for MA and original surface
	std::shared_ptr<MyMesh> m_meshMA;
	std::shared_ptr<TriMesh> m_meshOrig;
	// the transformation used for orig & ma mesh
	trimesh::XForm<double> m_trans_mat;

	// a list of radius. empty means no radius available.
	vector<float> radii;
	std::shared_ptr<SteinerGraph> stg;
	std::shared_ptr<IsoSurfFrom2Manifold> m_iso_surf;
	std::shared_ptr<IsoContourOnMA> m_iso_cont;

	// cur distance over MA for visualization
	vector<float> cur_MA_visDist_f;
	// a diffuser for MA in case there are "scalar holes" to fill on MA
	GraphDiffuser m_diffuser_MA;
	// cur distance metric over MA used for pruning
	vector<float> cur_MA_pruneDist_f_1;
	vector<float> cur_MA_pruneDist_f_2;
	vector<vector<float>> cur_MA_pruneDist_persheet_1;
	vector<vector<float>> cur_MA_pruneDist_persheet_2;
	// cur distance metric over MC used for visualization
	vector<float> cur_MC_visDist;
	// cur distance over medial curve vertices
	vector<float> cur_medialCurve_distMetric;
	float minDist_MC;
	float maxDist_MC;
	// cache the neighbor-ship of medial curve structure 
	// usually computed (once) before pruning happens
	vector<vector<unsigned>> m_MC_vts_neighbors;
	vector<int> m_MC_vts_nbsCnt;
	// the remained MC structure ( a list of edge indices into stg->dual_edges)
	vector<unsigned> m_remained_MC;
	vector<int> m_remained_MC_nbsCnt;

	/// glsl programs for different drawers
	std::shared_ptr<oglplus::Program> m_simpProg;
	std::shared_ptr<oglplus::Program> m_edgeProg;
	std::shared_ptr<oglplus::Program> m_pointsProg;
	std::shared_ptr<oglplus::Program> m_linesProg;
	std::shared_ptr<oglplus::Program> m_linesProxyProg;
	// Depth-peeling (DP) programs
	std::shared_ptr<oglplus::Program> peel_prog, post_prog;
	// Bavoil's OIT programs
	std::shared_ptr<oglplus::Program> m_compositeProg;
	// debug tool screen vertex picking (SVP) program
	std::shared_ptr<oglplus::Program> m_svpProg;

	/// full-screen quad. Used for full-screen post-processing
	std::shared_ptr<oglplus::VertexArray> m_bigQuad_vao;
	std::shared_ptr<oglplus::Buffer> m_bigQuad_vts_vbo;

	/// OIT rendering related
	// Texture used for accumulating rgb (color) & a (alpha)
	// & counting # fragments for each texel in the accum_tex
	// & saving the background contributed by the opaque objects
	int accum_tex_unit_OIT, count_tex_unit_OIT, bg_tex_unit_OIT;
	std::shared_ptr<oglplus::Texture> accum_tex_OIT, count_tex_OIT, bg_tex_OIT, depth_tex_OIT;
	// FBO used for "rendering" to the two textures
	std::shared_ptr<oglplus::Framebuffer> m_accum_fbo_OIT;
	
	/// DP rendering related
	bool use_DP_transparency;
	int m_DPMaxRenders;
	// Textures (depths, foreground, background)
	int bg_tex_unit, fg_tex_unit, refDepth_tex_unit, fgDepth_tex_unit;
	enum {RefDepthTexID=0, FGDepthTexID=1, CompositeDepthTexID=2};
	std::shared_ptr<oglplus::Texture> bg_tex, fg_tex, depth_tex[3];
	// FBO attachment points used by above textures
	vector<oglplus::FramebufferAttachment> fbo_attachments;
	vector<FramebufferColorAttachment> fbo_drawBuffers;
	int bg_attachPoint, fg_attachPoint, compDepth_attachPoint;
	// the Depth-peeling FBO 
	std::shared_ptr<oglplus::Framebuffer> peel_blend_fbo;

	/// debug tool: screen vertex picking related
	bool SVP_activated; // use screen vertex picking
	// Texture (vertex id per pixel)
	int vertID_tex_unit;
	std::shared_ptr<oglplus::Texture> vertID_tex_SVP;
	std::shared_ptr<oglplus::Texture> depth_tex_SVP;
	// FBO attachment point used by texture
	oglplus::FramebufferAttachment fbo_attachment_SVP;
	FramebufferColorAttachment fbo_drawbuffer_SVP;
	// the SVP FBO
	std::shared_ptr<oglplus::Framebuffer> SVP_fbo;

	/// use track ball
	shared_ptr<TrackBall> track_ball;

	// extra info for going to full-screen and back
	QWidget* m_pParent;
	QSize m_pSize;
	Qt::WindowFlags m_enOrigWindowFlags;
	bool m_fullMode;

	/// render-call list of <drawer obj, is_face_draw? >
	vector<std::pair<std::shared_ptr<Drawable>,bool> > transparent_draws;
	vector<std::pair<std::shared_ptr<Drawable>,bool> > opaque_draws;
	/// flags
	bool m_drawOrig;
	bool m_drawMA, m_drawMALines;
	bool m_drawMAFinnerStatic;
	bool m_drawMC;
	bool m_drawBurntEdges;
	bool m_drawMP;
	bool m_drawPoints;
	bool m_drawLinesAsProxy;
	bool m_drawHS;
	bool m_drawIsoSurf;
	bool m_drawIsoCont;
	bool m_drawQMAT;
	bool m_MATransparent;
	bool m_OrigTransparent;
	bool m_lineTransparent;

	/// const color for MA
	bool m_use_constColor_for_MA;
	TriColor m_constColor_MA;

	bool m_radiiReady;

	TriColor bg_color;

	/// the surface function object (an app. of MC)
	std::shared_ptr<SurfaceFunc> m_surfF;
	/// hybrid skeleton
	std::shared_ptr<HybridSkeleton> m_hs;

    bool m_glInitialized = false; // workaround for force OpenGL initialize (xlzhang)
};

#endif