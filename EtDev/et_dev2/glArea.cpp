#include "glArea.h"
#include "meshDrawer.h"
#include "pointDrawer.h"
#include "lineDrawer.h"
#include "sphereDrawer.h"
#include "iso_contour.h"
#include "io.h"

#include <queue>
#include <unordered_set>
#include <float.h>
#include <QFileInfo>

QGLFormat GLArea::getGLFormat(void)
{
	QGLFormat qgl_format;
	qgl_format.setDepth(true); // depth buffer
	cout << "qgl depth size (before request): "<<qgl_format.depthBufferSize();
	qgl_format.setDepthBufferSize(32);
	cout << "qgl depth size (after request): "<<qgl_format.depthBufferSize();
	qgl_format.setDoubleBuffer(true); // double framebuffer
	return qgl_format;
}

shared_ptr<oglplus::Program> GLArea::glslProgram(string _vs, string _fs)
{
	try {
		cout << "Assembling shader program from: " << endl
			<< _vs << ",\n" << _fs << endl;
		using namespace oglplus;
		shared_ptr<Program> prog(new Program());
		prog->AttachShaders(
			Group<Shader>(
			VertexShader().Source(loadShaderSource(_vs.c_str())).Compile(),
			FragmentShader().Source(loadShaderSource(_fs.c_str())).Compile()
			)
			).Link();
		cout << "Done." << endl;
		return prog;
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

	exit(1);
}

shared_ptr<oglplus::Program> GLArea::glslProgram(string _vs, string _gs, string _fs)
{
	try
	{
		cout << "Assembling shader program from: " << endl
			<< _vs << ",\n" << _gs << ",\n" << _fs << endl;
		using namespace oglplus;
		shared_ptr<Program> prog(new Program());
		prog->AttachShaders(
			Group<Shader>(
			VertexShader().Source(loadShaderSource(_vs.c_str())).Compile(),
			GeometryShader().Source(loadShaderSource(_gs.c_str())).Compile(),
			FragmentShader().Source(loadShaderSource(_fs.c_str())).Compile()
			)
			).Link();
		cout << "Deon." << endl;
		return prog;
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

	exit(1);
}

GLArea::GLArea(QWidget *parent /* = 0 */)
	: QGLWidget(getGLFormat(), parent)
{
	cout << "creating GLArea." << endl; // debug

	// glarea should accept both mouse and tab focus
	setFocusPolicy(Qt::StrongFocus);

	bg_color = TriColor(0.85f, 0.85f, 0.85f);

	track_ball = std::shared_ptr<TrackBall>(new TrackBall());

	//m_lightingFlags.assign(N_OBJECT_TYPES, true);

	m_OrigTransparent = m_MATransparent = m_lineTransparent = false;
	m_use_constColor_for_MA = false;

	m_fullMode = false;
	m_pParent = nullptr;

	m_radiiReady = false;

	// create a surf func obj, an app of MC
	m_surfF = shared_ptr<SurfaceFunc>(new SurfaceFunc());

	// create a hybrid skeleton obj, an app of MC
	m_hs = shared_ptr<HybridSkeleton>(new HybridSkeleton());

	// an iso-surf obj for deforming the original 3D shape onto MA
	m_iso_surf = shared_ptr<IsoSurfFrom2Manifold>(new IsoSurfFrom2Manifold());

	// an iso-contour obj for displaying iso-contour on MA
	m_iso_cont = shared_ptr<IsoContourOnMA>(new IsoContourOnMA());

	/// setup rendering opengl related
	bg_tex_unit_OIT = 0;
	accum_tex_unit_OIT = 1;
	count_tex_unit_OIT = 2;

	/// setup DP related
	use_DP_transparency = false;
	m_DPMaxRenders = 5;
	fg_tex = nullptr;
	bg_tex = nullptr;
	depth_tex[0] = nullptr;
	depth_tex[1] = nullptr;
	depth_tex[2] = nullptr;
	peel_blend_fbo = nullptr;

	reset();
}

GLArea::~GLArea(void)
{

}

QSize GLArea::minimumSizeHint(void) const
{
	return QSize(500, 500);
}

QSize GLArea::sizeHint(void) const
{
	return QSize(800, 600);
}

bool GLArea::event(QEvent *_e)
{
	if (_e->type() == QEvent::KeyPress)
	{
		auto ke = static_cast<QKeyEvent*>(_e);
		if (ke->key() == Qt::Key_Space)
		{
			//cout << "key space pressed." << endl;

			this->track_ball->reset();
			updateGL();

			return true;
		}
	}

	return QGLWidget::event(_e);
}

void GLArea::reset()
{
#ifdef PROFILE_SPEED
	clock_t t_start, t_duration;
#endif

	//stg = NULL;
#ifdef PROFILE_SPEED
	t_start = clock();
#endif
	stg.reset();
	m_meshOrig.reset();
	m_meshMA.reset();
	m_surfF->reset();
	m_hs->reset();
	m_iso_surf->reset();
	m_iso_cont->reset();
#ifdef PROFILE_SPEED
	t_duration = clock() - t_start;
	cout << "PROFIELING SPEED: glarea.reset() took " << t_duration * 1000.0f / CLOCKS_PER_SEC << "ms."<<endl;
#endif

	this->cur_MA_visDist_f.clear();
	this->cur_MA_pruneDist_f_1.clear();
	this->cur_MA_pruneDist_f_2.clear();
	this->cur_MA_pruneDist_persheet_1.clear();
	this->cur_MA_pruneDist_persheet_2.clear();
	this->cur_medialCurve_distMetric.clear();
	this->m_MC_vts_neighbors.clear();

	// TODO: reset drawers
	m_drawOrig = false;
	m_drawPoints = false;
	m_drawMA = false;
	m_drawMALines = false;
	m_drawMAFinnerStatic = false;
	m_drawMC = false;
	m_drawBurntEdges = false;
	m_drawMP = false;
	m_drawHS = false;
	m_drawIsoSurf = false;
	m_drawIsoCont = false;
	m_drawQMAT = false;
	m_drawLinesAsProxy = false;
	m_OrigTransparent = m_MATransparent = m_lineTransparent = false;
	m_use_constColor_for_MA = false;

	for (auto it = m_drawers.begin(); it != m_drawers.end(); ++it)
	{
		if (*it)
		{
			(*it)->reset();
		}
	}

	use_DP_transparency = false;
	m_DPMaxRenders = 10;
}

void GLArea::saveView(std::shared_ptr<QSettings>& _qsetting)
{
	track_ball->saveToSetting(_qsetting);
}

void GLArea::loadView(std::shared_ptr<QSettings>& _qsetting)
{
	track_ball->loadFromSetting(_qsetting);
}

void GLArea::passToDrawable(
	std::string _medialMesh_file, 
	std::string _origMesh_file,
	std::string _r_file,
	bool _remove_dup_faces )
{
	/* pass MA file to its mesh drawer */
	/* and pass orig 3d surface file to its mesh drawer */
    while (!m_glInitialized) {
        QApplication::processEvents();
    }
    assert(m_MADrawer);

	std::shared_ptr<MyMesh> mesh_MA;
	float eps = (float)ui->epsilonSpin->value();
	float pert = (float)ui->pertSpin->value();
	std::shared_ptr<TriMesh> mesh_3d;

	bool loaded = et::loadMesh(
		_medialMesh_file.c_str(), mesh_MA, 
		_r_file.c_str(), &radii, 
		_origMesh_file.c_str(), mesh_3d, 
		m_trans_mat,
		eps, pert, _remove_dup_faces );

	if (!loaded)
	{
		cout << "draw nothing." << endl;
		m_radiiReady = false;
		return;
	}
	m_radiiReady = radii.size();
	this->m_meshMA = mesh_MA;
	this->m_medialAxisFile = _medialMesh_file;
	this->m_meshOrig = mesh_3d;

	printf( "Medial axis #V/#E/#F = %d/%d/%d\n", 
		m_meshMA->vertices.size(), m_meshMA->lines.size(), m_meshMA->faces.size() );

	std::dynamic_pointer_cast<MeshDrawer>(m_MADrawer)->setRenderMode(MeshDrawer::PER_VERT);
	std::dynamic_pointer_cast<MeshDrawer>(m_MADrawer)->setMeshToDraw(mesh_MA);
	std::dynamic_pointer_cast<MeshDrawer>(m_origDrawer)->setRenderMode(MeshDrawer::PER_VERT);
	std::dynamic_pointer_cast<MeshDrawer>(m_origDrawer)->setMeshToDraw(mesh_3d);
	std::dynamic_pointer_cast<MeshDrawer>(m_MADrawer)->initCameraPosAndLens();
	
	m_meshOrig->need_bbox();
	float s = m_meshOrig->bbox.size().max();
	
	/* setup line geometry of MA to render */
	auto ma_line_drawer_ptr = std::dynamic_pointer_cast<LineDrawer>( m_MALineDrawer );
	uploadLinesToDrawer( mesh_MA->vertices, mesh_MA->lines, ma_line_drawer_ptr );
	cout << "MA lines uploaded to GPU!" << endl;
	TriColor const_edge_color( 0.469, 0.469, 0.469 );
	auto color_data = new float[ mesh_MA->lines.size() * 2 * 3 ];
	for ( unsigned i = 0; i < mesh_MA->lines.size(); i++ )
	{
		const auto& edge_color = const_edge_color;
		color_data[ i * 2 * 3 + 0 ] = edge_color[ 0 ];
		color_data[ i * 2 * 3 + 1 ] = edge_color[ 1 ];
		color_data[ i * 2 * 3 + 2 ] = edge_color[ 2 ];
		color_data[ i * 2 * 3 + 3 + 0 ] = edge_color[ 0 ];
		color_data[ i * 2 * 3 + 3 + 1 ] = edge_color[ 1 ];
		color_data[ i * 2 * 3 + 3 + 2 ] = edge_color[ 2 ];
	}
	cout << "Setting color for MA lines..." << endl;
	ma_line_drawer_ptr->setPerVertColor( color_data, mesh_MA->lines.size() * 2 );
	delete[] color_data;
	
	// clear the history data in GPU managed by line/point drawers
	/*m_linesDrawer->clear();
	m_pointDrawer->clear();*/

	// need to re-compute steinergraph every time a new mesh loaded
	//this->stg = NULL; 

	// update the geometry scale of the other drawers
	std::dynamic_pointer_cast<PointDrawer>(m_pointDrawer)->setScale(s);
	std::dynamic_pointer_cast<LineDrawer>(m_dualLinesDrawer)->setScale(s);
	std::dynamic_pointer_cast<LineDrawer>(m_dynamic_dualLineDrawer)->setScale(s);
	std::dynamic_pointer_cast<LineDrawer>(m_burntEdgesDrawer)->setScale(s);
	std::dynamic_pointer_cast<MeshDrawer>(m_MADrawer)->setScale(s);
	std::dynamic_pointer_cast<LineDrawer>( m_MALineDrawer )->setScale( s );
	std::dynamic_pointer_cast<MeshDrawer>(m_FinerMAStaticDrawer)->setScale(s);
	std::dynamic_pointer_cast<MeshDrawer>(m_MAFinnerDynamicDrawer)->setScale(s);
	std::dynamic_pointer_cast<SphereDrawer>(m_MPDrawer)->setScale(s);
	std::dynamic_pointer_cast<MeshDrawer>(m_origDrawer)->setScale(s);
	std::dynamic_pointer_cast<MeshDrawer>(m_hsFaceDrawer)->setScale(s);
	std::dynamic_pointer_cast<LineDrawer>(m_hsLineDrawer)->setScale(s);
	std::dynamic_pointer_cast<MeshDrawer>(m_isoSurfDrawer)->setScale(s);
	m_iso_cont->setupGeometryScale(s);
	std::dynamic_pointer_cast<MeshDrawer>(m_qmatFaceDrawer)->setScale(s);
	std::dynamic_pointer_cast<LineDrawer>(m_qmatEdgeDrawer)->setScale(s);

	// set flags
	/*this->m_drawMA = true;
	this->m_drawOrig = true;*/

	this->resizeGL(width(), height());
	updateGL();
}

void GLArea::wireFrameOnMesh(bool _draw_or_not)
{
    if (!m_MADrawer) {
        return;
    }
	std::dynamic_pointer_cast<MeshDrawer>(m_MADrawer)->setDrawEdge(_draw_or_not);
	std::dynamic_pointer_cast<MeshDrawer>(m_origDrawer)->setDrawEdge(_draw_or_not);
	std::dynamic_pointer_cast<MeshDrawer>(m_isoSurfDrawer)->setDrawEdge(_draw_or_not);
	m_iso_cont->enableWireFrame(_draw_or_not);

	this->resizeGL(width(), height());
	updateGL();
}

void GLArea::drawSteinerPoints(bool _draw, float _size)
{
	if (_draw && stg)
	{
		this->m_drawPoints = true;
		std::dynamic_pointer_cast<PointDrawer>(m_pointDrawer)->setPointSize(_size);
		std::dynamic_pointer_cast<PointDrawer>(m_pointDrawer)->setPoints(stg->m_stSubdiv.getAllVts());
	}
	else
	{
		this->m_drawPoints = false;
	}

	this->resizeGL(width(), height());
	updateGL();
}

void GLArea::drawMCProxy(bool _draw_or_not)
{
	m_drawLinesAsProxy = _draw_or_not;

	resizeGL(width(), height());
	updateGL();
}

bool GLArea::obtainDualLineStructure( 
	bool _only_unburnt, int _edge_dual_method_idx, int _poly_dual_method_idx, bool _stop_burn_at_junction)
{
	if (stg)
	{
		try {
#ifdef PROFILE_SPEED
			auto t1 = clock();
#endif // PROFILE_SPEED

			shared_ptr<LineDrawer> line_drawer = std::dynamic_pointer_cast<LineDrawer>(m_dualLinesDrawer);
			/*line_drawer->setPoints(stg->st_vts);
			line_drawer->setLines(stg->st_edges);*/

			cout << "collecting burnt steiner edges... " << endl;
			vector<std::pair<TriEdge, int> > burntE_assocF_list;
			vector< std::pair<size_t,list<TriEdge>> > faces_w_crossed_paths;
			stg->getBurnedStEdgesFromBurnTree(burntE_assocF_list, nullptr/*&faces_w_crossed_paths*/);

			cout << "dualizing the planar subdivision... " << endl;
			m_meshMA->need_bbox();
			float eps = m_meshMA->bbox.size().max()*1e-3f;
			stg->dualizeFacesWithBurnPaths(
				burntE_assocF_list, 
				static_cast<SteinerGraph::DualOpt>(_edge_dual_method_idx),
				static_cast<SteinerGraph::DualOpt>(_poly_dual_method_idx),
				eps);
			cout << "Done! # dual vts/edges " 
				<<stg->dual_vts.size()<<"/"<<stg->dual_edges.size()<< endl;
			
#ifdef PROFILE_SPEED
			auto t2 = clock();
			auto t_ms = ( t2 - t1 ) * 1000.0 / CLOCKS_PER_SEC;
			cout << "Medial-curve creation took: " << t_ms << " ms." << endl;
#endif // PROFILE_SPEED

			// upload geometry & color of burnt st edges to drawer
			vector<TriEdge> burnt_edges;
			const auto& all_st_vts = stg->m_stSubdiv.getAllVts();
			burnt_edges.reserve(burntE_assocF_list.size());
			for (auto it = burntE_assocF_list.begin(); it != burntE_assocF_list.end(); ++it)
			{
				burnt_edges.push_back(it->first);
			}
			auto burntEdges_drawer = std::dynamic_pointer_cast<LineDrawer>(m_burntEdgesDrawer);
			//this->uploadLinesToDrawer(all_st_vts, burnt_edges, burntEdges_drawer);
			burntEdges_drawer->setPoints(all_st_vts);
			burntEdges_drawer->setLines(burnt_edges);
			TriColor burnt_edge_color(0.0f, 1.0f, 0.0f);
			float* color_data = new float[3 * all_st_vts.size()];
			for (unsigned i = 0; i < all_st_vts.size(); ++i)
			{
				color_data[3*i+0] = burnt_edge_color[0];
				color_data[3*i+1] = burnt_edge_color[1];
				color_data[3*i+2] = burnt_edge_color[2];
			}
			burntEdges_drawer->setPerVertColor(color_data, all_st_vts.size());
			delete [] color_data;

			// debug cycles in dual network
			char cycle_file[] = "dual_cycles.txt";
			std::ofstream out_file(cycle_file);
			if (!out_file.good())
				cout << "Error: couldn't open file "<< cycle_file << endl;
			else
			{
				stg->debug_dualline_cycles(
					stg->dual_vts.size(), stg->dual_edges, stg->m_dualEdge_from_which_face, true, /*&out_file*/nullptr);
				out_file.close();
				cout << "debug_dualline_cycles() finished!" << endl;
			}

			// debug conn. cmpnt.s in dual network
			char cc_file[] = "dual_cc.txt";
			out_file = std::ofstream(cc_file);
			if (!out_file.good())
				cout << "Error: couldn't open file "<< cc_file << endl;
			else
			{
				stg->debug_dualline_conn_comp(
					stg->dual_vts.size(), stg->dual_edges, stg->m_dualEdge_from_which_face, true, /*&out_file*/nullptr);
				out_file.close();
				cout << "debug_dualline_conn_comp() finished!" << endl;
			}

			// debug cycles in burn network
			//char file_name[] = "burn_cycles.txt";
			//out_file = std::ofstream(file_name);
			//if (!out_file.good())
			//	cout << "Error: couldn't open file "<< file_name << endl;
			//else
			//{
			//	stg->debug_burnline_cycles(
			//		stg->m_stSubdiv.getAllVts(), burntE_assocF_list, true, /*&out_file*/nullptr);
			//	out_file.close();
			//	cout << "debug_burnline_cycles() finished!" << endl;
			//}

			// clean up
			burntE_assocF_list.clear();burntE_assocF_list.shrink_to_fit();
			faces_w_crossed_paths.clear();faces_w_crossed_paths.shrink_to_fit();

			cout << "components & cycles analysis done!" << endl << endl;

			//this->m_drawMC = true;
			return true;
		}
		catch(const std::exception& se)
		{
			this->m_drawMC = false;
			std::cerr <<
				"General error: " <<
				se.what() << std::endl;
		}
	}
	return false;
}

bool GLArea::burnCurveNetwork(
	bool _stop_burn_at_junction, 
	bool _only_unburnt, 
	bool _protect_bt2, 
	bool _underestimate_dist)
{
#ifdef PROFILE_SPEED
	auto t1 = clock();
#endif // PROFILE_SPEED

	// carry out the burning process on the dual network
	stg->burnMedialCurveNetwork(true, _stop_burn_at_junction, _protect_bt2, _underestimate_dist);

	// upload line data, color(BT1), and get ready for MC pruning/viewing
	this->precomputeForMCPruning();
	float dummy = 0.0f;
	this->colorMCEdgeBy(DistMC::BT1_MC, dummy, dummy, false, 0.0f, 1.0f, 1);

#ifdef PROFILE_SPEED
	auto t2 = clock();
	auto t_ms = ( t2 - t1 ) * 1000.0 / CLOCKS_PER_SEC;
	cout << "Medial-curve burning took: " << t_ms << " ms." << endl;
#endif // PROFILE_SPEED

	return true;
}

bool GLArea::setUpMCDistMetricForPruning(DistMC _metric, float& _min, float& _max)
{
	if (!stg || stg->bt2_MC.empty() || stg->bt1_medialCurve.empty())
		return false;
	cur_medialCurve_distMetric.clear();
	getMCDistMetric( _metric, cur_medialCurve_distMetric );
	//switch (_metric)
	//{
	///*case DEFAULT_DIST_DUAL:
	//	break;*/
	//case BT2_MC:
	//	cur_medialCurve_distMetric = stg->bt2_MC;
	//	break;
	//case BT1_MC:
	//	cur_medialCurve_distMetric = stg->bt1_medialCurve;
	//	break;
	//case BT1_BT2_MC:
	//	cur_medialCurve_distMetric.clear();
	//	cur_medialCurve_distMetric.reserve(stg->dual_vts.size());
	//	for (unsigned i = 0; i < stg->bt1_medialCurve.size(); ++i)
	//	{
	//		cur_medialCurve_distMetric.push_back(
	//			stg->bt1_medialCurve[i] - stg->bt2_MC[i]
	//		);
	//	}
	//	break;
	//case BT1_BT2_REL_MC:
	//	cur_medialCurve_distMetric.clear();
	//	cur_medialCurve_distMetric.reserve(stg->dual_vts.size());
	//	for (unsigned i = 0; i < stg->bt1_medialCurve.size(); ++i)
	//	{
	//		cur_medialCurve_distMetric.push_back(
	//			1.0f - stg->bt2_MC[i] / stg->bt1_medialCurve[i]
	//		);
	//	}
	//	break;
	//}
	/*auto min_max = std::minmax_element(cur_medialCurve_distMetric.begin(),
		cur_medialCurve_distMetric.end());*/
	_min = numeric_limits<float>::max();
	_max = 0.0f;
	for (auto d_it = cur_medialCurve_distMetric.begin(); d_it != cur_medialCurve_distMetric.end(); ++d_it)
	{
		if (*d_it < numeric_limits<float>::max())
		{
			_min = std::min(_min, *d_it);
			_max = std::max(_max, *d_it);
		}
	}
	/*minDist_MC = _min = *min_max.first;
	maxDist_MC = _max = *min_max.second;*/
	minDist_MC = _min;
	maxDist_MC = _max;

	return true;
}

bool GLArea::visDualLineStructure(DistMC _dist_type, int _type)
{
	cout << "visDualLineStructure()..." << endl;
	if (stg)
	{
		try {
			const vector<float>* dist_to_vis = nullptr;
			bool mem_newed = false;
			switch (_dist_type)
			{
			case BT2_MC:
				dist_to_vis = &stg->bt2_MC;
				mem_newed = false;
				break;
			case BT1_MC:
				dist_to_vis = &stg->bt1_medialCurve;
				mem_newed = false;
				break;
			case BT1_BT2_MC:
				{
					auto* diff_dist = new vector<float>;
					mem_newed = true;
					diff_dist->reserve(stg->dual_vts.size());
					for (unsigned i = 0; i < stg->bt1_medialCurve.size(); ++i)
					{
						diff_dist->push_back(
							std::max(0.0f, stg->bt1_medialCurve[i] - stg->bt2_MC[i])
						);
					}
					dist_to_vis = diff_dist;
					break;
				}
			case BT1_BT2_REL_MC:
				{
					auto* rel_diff_dist = new vector<float>;
					mem_newed = true;
					rel_diff_dist->reserve(stg->dual_vts.size());
					for (unsigned i = 0; i < stg->bt1_medialCurve.size(); ++i)
					{
						rel_diff_dist->push_back(
							1.0f - stg->bt2_MC[i] / stg->bt1_medialCurve[i]
						);
					}
					dist_to_vis = rel_diff_dist;
					break;
				}
			}

			// add steiner points and dual points for visualizing burnt st. edges and dual edges
			vector<TriPoint> all_vts;
			all_vts.reserve(stg->m_stSubdiv.sizeOfVts() + stg->dual_vts.size());
			all_vts.insert(all_vts.end(), stg->m_stSubdiv.getAllVts().begin(), stg->m_stSubdiv.getAllVts().end());
			all_vts.insert(all_vts.end(), stg->dual_vts.begin(), stg->dual_vts.end());

			vector<TriEdge> all_edges;
			unsigned n_burntEdges = 0;
			if (_type > 0)
				all_edges.reserve(stg->m_stSubdiv.sizeOfStEdges() + stg->dual_edges.size());
			if (_type & VIS_BURNT)
			{
				// get all st edges that are burnt, add them to edges-to-visualize
				vector<std::pair<TriEdge, int> > burntE_assocF_list;
				stg->getBurnedStEdgesFromBurnTree(burntE_assocF_list);
				std::for_each(
					burntE_assocF_list.begin(), burntE_assocF_list.end(),
					[&all_edges] (const std::pair<TriEdge,int>& _e_f) {
						const TriEdge& e = _e_f.first;
						all_edges.push_back(e);
				} );
				n_burntEdges = burntE_assocF_list.size();
			}
			if (_type & VIS_DUAL)
			{
				// add all dual edges to edges-to-visualize
				all_edges.insert(all_edges.end(), stg->dual_edges.begin(), stg->dual_edges.end());
				// need to offset dual vert indices of dual edges 
				// because dual vts follow st. vts
				std::for_each(
					all_edges.begin() + n_burntEdges, all_edges.end(),
					[this] (TriEdge& _e) {
						_e[0] += this->stg->m_stSubdiv.sizeOfVts();
						_e[1] += this->stg->m_stSubdiv.sizeOfVts();
				} );
			}

			float* colors = new float[all_vts.size()*3];
			
			// give burnt edges' vertices green color
			for (unsigned i = 0; i < stg->m_stSubdiv.sizeOfVts(); ++i)
			{
				colors[i*3 + 0] = 0.0f;
				colors[i*3 + 1] = 1.0f;
				colors[i*3 + 2] = 0.0f;
			}

			// give dual vts blue color, or burn dist color
			if (dist_to_vis)
			{
				float min_dist = numeric_limits<float>::max();
				float max_dist = -1.0f;
				for (auto it = dist_to_vis->begin(); it != dist_to_vis->end(); ++it)
					if (*it < numeric_limits<float>::max())
					{
						min_dist = std::min(*it, min_dist);
						max_dist = std::max(*it, max_dist);
					}
				for (unsigned i = stg->m_stSubdiv.sizeOfVts(); i < all_vts.size(); ++i)
				{
					TriColor c = util::GetColour((*dist_to_vis)[i-stg->m_stSubdiv.sizeOfVts()], 
						min_dist, max_dist);
					colors[i*3 + 0] = c[0];
					colors[i*3 + 1] = c[1];
					colors[i*3 + 2] = c[2];
				}
				// debug
				cout << "medial curve dist range: ["<<min_dist<<","<<max_dist<<"]"<<endl;
			}
			else
			{
				for (unsigned i = stg->m_stSubdiv.sizeOfVts(); i < all_vts.size(); ++i)
				{
					colors[i*3 + 0] = 0.0f;
					colors[i*3 + 1] = 0.0f;
					colors[i*3 + 2] = 1.0f;
				}
			}
			if (mem_newed)
				delete dist_to_vis;

			// upload line data to GPU
			cout << "visDualLineStructure(): uploading to GPU for render... " << endl;
			shared_ptr<LineDrawer> line_drawer = std::dynamic_pointer_cast<LineDrawer>(m_dualLinesDrawer);
			line_drawer->setPoints(all_vts);
			line_drawer->setLines(all_edges);
			line_drawer->setPerVertColor(colors, all_vts.size());
			delete [] colors;
			cout << "Done!" << endl;

			//this->m_drawMC = true;
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
	else
	{
		this->m_drawMC = false;
	}

	this->resizeGL(width(), height());
	updateGL();

	cout << "visDualLineStructure() returns." << endl;

	return this->m_drawMC;
}

bool GLArea::colorMCEdgeBy(
	DistMC _dist_type, float& _min, float& _max, 
	bool _use_new_range, 
	float _min_alpha, float _alpha_exp, 
	int _update_opt)
{
	cout << "colorMCEdgeBy()... " << endl;
	if (stg)
	{
		try {
			//getMCDistMetric(_dist_type, cur_MC_visDist);
			vector<trimesh::vec2> cur_MC_visDist_byEdge;
			getMCDistMetric(_dist_type, cur_MC_visDist_byEdge);

			unsigned n_total_vts = stg->dual_edges.size()*2;
			if (_update_opt == 1 || _update_opt == 3)
			{
				// constant color within edge by duplicating color for shared vts
				float* colors = new float[n_total_vts*3];

				// color each edge by the smaller distance of two ends
				float min_dist = numeric_limits<float>::max();
				float max_dist = -1.0f;
				for (auto it = cur_MC_visDist_byEdge.begin(); it != cur_MC_visDist_byEdge.end(); ++it)
				{
					if ( it->at( 0 ) < stg->infiniteBurnDist() )
					{
						min_dist = std::min(it->at(0), min_dist);
						max_dist = std::max(it->at(0), max_dist);
					}
					if (it->at(1) < stg->infiniteBurnDist() )
					{
						min_dist = std::min(it->at(1), min_dist);
						max_dist = std::max(it->at(1), max_dist);
					}
				}

				if (_use_new_range)
				{
					min_dist = _min;
					max_dist = _max;
				}
				else
				{
					_min = min_dist;
					_max = max_dist;
				}
				TriColor c;
				for (auto it = cur_MC_visDist_byEdge.begin(); it != cur_MC_visDist_byEdge.end(); ++it)
				{
					auto ei = it - cur_MC_visDist_byEdge.begin();
					if ( it->at( 0 ) < stg->infiniteBurnDist() )
						c = util::GetColour( it->at( 0 ), min_dist, max_dist );
					else
						c = TriColor( 0, 0, 0 );
					colors[ei*2*3 + 0 + 0] = c[0];
					colors[ei*2*3 + 0 + 1] = c[1];
					colors[ei*2*3 + 0 + 2] = c[2];
					if ( it->at( 1 ) < stg->infiniteBurnDist() )
						c = util::GetColour( it->at( 1 ), min_dist, max_dist );
					else
						c = TriColor( 0, 0, 0 );
					colors[ei*2*3 + 3 + 0] = c[0];
					colors[ei*2*3 + 3 + 1] = c[1];
					colors[ei*2*3 + 3 + 2] = c[2];
				}
				/*for (unsigned ei = 0; ei < stg->dual_edges.size(); ei++)
				{
					const auto& e = stg->dual_edges[ei];
					float small_dist = std::min( cur_MC_visDist[e[0]], cur_MC_visDist[e[1]] );
					TriColor c = util::GetColour(small_dist, min_dist, max_dist);
					colors[ei*2*3 + 0 + 0] = c[0];
					colors[ei*2*3 + 0 + 1] = c[1];
					colors[ei*2*3 + 0 + 2] = c[2];

					colors[ei*2*3 + 3 + 0] = c[0];
					colors[ei*2*3 + 3 + 1] = c[1];
					colors[ei*2*3 + 3 + 2] = c[2];
				}*/
				// debug
				cout << "medial curve dist range: ["<<min_dist<<","<<max_dist<<"]"<<endl;

				// upload color data to GPU
				cout << "visDistOnMC(): uploading to GPU for render... " << endl;
				shared_ptr<LineDrawer> line_drawer = std::dynamic_pointer_cast<LineDrawer>(m_dualLinesDrawer);
				line_drawer->setPerVertColor(colors, n_total_vts);
				// also populate the color buffer for dynamic part of MC
				std::dynamic_pointer_cast<LineDrawer>(m_dualLinesDrawer)->setPerVertColor(colors, n_total_vts);
				delete [] colors;
			}
			if (_update_opt == 2 || _update_opt == 3)
			{
				// saliency of each face. can be used as alpha in shader.
				float* saliency_data = new float[n_total_vts];

				// color each edge by the smaller distance of two ends
				float min_dist = numeric_limits<float>::max();
				float max_dist = -1.0f;
				for (auto it = cur_MC_visDist_byEdge.begin(); it != cur_MC_visDist_byEdge.end(); ++it)
				{
					if ( it->at( 0 ) < stg->infiniteBurnDist() )
					{
						min_dist = std::min(it->at(0), min_dist);
						max_dist = std::max(it->at(0), max_dist);
					}
					if ( it->at( 1 ) < stg->infiniteBurnDist() )
					{
						min_dist = std::min(it->at(1), min_dist);
						max_dist = std::max(it->at(1), max_dist);
					}
				}
				for (auto it = cur_MC_visDist_byEdge.begin(); it != cur_MC_visDist_byEdge.end(); ++it)
				{
					auto ei = it - cur_MC_visDist_byEdge.begin();
					if ( it->at( 1 ) < stg->infiniteBurnDist() )
						saliency_data[ ei * 2 + 1 ] = util::rescale(
							it->at( 1 ), min_dist, max_dist, _alpha_exp, _min_alpha, 1.0f );
					else
						saliency_data[ ei * 2 + 1 ] = 1.0f;
					if ( it->at( 0 ) < stg->infiniteBurnDist() )
						saliency_data[ ei * 2 ] = util::rescale( 
							it->at( 0 ), min_dist, max_dist, _alpha_exp, _min_alpha, 1.0f );
					else
						saliency_data[ ei * 2 ] = 1.0f;
				}
				/*for (unsigned ei = 0; ei < stg->dual_edges.size(); ei++)
				{
					const auto& e = stg->dual_edges[ei];
					float small_dist = std::min( cur_MC_visDist[e[0]], cur_MC_visDist[e[1]] );

					saliency_data[ei*2 + 1] = saliency_data[ei*2] = 
						util::rescale(small_dist, min_dist, max_dist, _alpha_exp, _min_alpha, 1.0f);

				}*/
				// debug
				//cout << "medial curve dist range: ["<<min_dist<<","<<max_dist<<"]"<<endl;

				shared_ptr<LineDrawer> line_drawer = std::dynamic_pointer_cast<LineDrawer>(m_dualLinesDrawer);
				// upload saliency data
				line_drawer->setPerVertSaliency(saliency_data, n_total_vts);
				// also populate the saliency for dynamic part of MC
				std::dynamic_pointer_cast<LineDrawer>(m_dualLinesDrawer)->setPerVertSaliency(saliency_data, n_total_vts);
				delete [] saliency_data;
			}
			cout << "Done!" << endl;

			//this->m_drawMC = true;
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
	else
	{
		this->m_drawMC = false;
	}

	this->resizeGL(width(), height());
	updateGL();

	cout << "colorMCEdgeBy() returns" << endl;

	return this->m_drawMC;
}

bool GLArea::precomputeForMCPruning()
{
	if (!stg || stg->dual_vts.empty())
		return false;

	this->m_MC_vts_neighbors.clear();
	this->m_MC_vts_nbsCnt.clear();
	stg->build_neighborship(stg->dual_vts, stg->dual_edges, m_MC_vts_neighbors, m_MC_vts_nbsCnt);
	uploadMCCompleteGeometry(stg->dual_vts, stg->dual_edges);
	return true;
}

bool GLArea::pruneAndVisMedialCurve( double _t, bool _preserve_topo, bool _constrained_by_iso_cont )
{
	if (m_MC_vts_neighbors.empty())
		return false;

	bool success = (_constrained_by_iso_cont ? 
		pruneMedialCurveConstrainedByIsoContour(_t) :
		pruneMedialCurve(_t, _preserve_topo) );
	success = success && uploadLineSubset(m_remained_MC, stg->dual_edges);
	updateGL();

	return success;
}

//bool GLArea::

float GLArea::toETcPruneRatio(float _v)
{
	if (!m_meshOrig)
		return -999.0f;
	float bt2_range = stg->max_burnDistAllVts - stg->min_burnDistAllVts;
	return _v / bt2_range;
}

bool GLArea::pruneMedialCurve( double _t, bool _preserve_topo )
{
	if (!stg || stg->dual_vts.empty())
		return false;

	vector<bool> removed(stg->dual_vts.size(), false);
	set<TriEdge> pruned_edges;
	
	m_remained_MC.clear();
	m_remained_MC.reserve(stg->dual_edges.size() - pruned_edges.size());
	m_remained_MC_nbsCnt = m_MC_vts_nbsCnt;

	if (_preserve_topo)
	{
		// add those boundary vts whose distance smaller than threshold to queue
		int ends_not_in_q = 0;
		std::queue< std::pair<unsigned,unsigned> > q; // <vid, eid>
		for (unsigned ei = 0; ei < stg->dual_edges.size(); ++ei)
		{
			const auto& e = stg->dual_edges[ei];
			if (m_remained_MC_nbsCnt[e[0]] == 1 || m_remained_MC_nbsCnt[e[1]] == 1)
			{
				unsigned vi = m_remained_MC_nbsCnt[e[0]] == 1 ? e[0] : e[1];
				unsigned nbi = util::otherEnd(e, vi);
				if (std::min(
					cur_medialCurve_distMetric[e[0]],
					cur_medialCurve_distMetric[e[1]]) < _t 
					)
				{
					q.push(std::make_pair(vi,nbi));
				}
				else
				{
					ends_not_in_q ++;
				}
			}
		}
		/*cout << "# ends not in queue "<< ends_not_in_q<<endl;
		cout << "_t "<< _t<<endl;*/

		//	start pruning iterations:
		//
		//	while q not empty:
		//		cur_vi <- q.pop();
		//		removed[cur_vi] <- true;
		//
		//      // now cur_vi has been pruned. decrease its neighbors' ref count
		//		nbs <- cur_vi's neighbors
		//		for each nb in nbs:
		//			skip if nb is removed
		//			ref_cnt[nb] --;
		//			if ref_cnt[nb] now is 1 && dist[nb] < t:
		//				q.push(nb);
		//				save nb to output list
		//			end if
		//		end for
		//	end while
		while(!q.empty())
		{
			const auto& cur_vid_eid_pr = q.front();
			unsigned cur_vi = cur_vid_eid_pr.first;
			unsigned cur_nbi = cur_vid_eid_pr.second;
			q.pop();
			assert(m_remained_MC_nbsCnt[cur_vi] == 1 || m_remained_MC_nbsCnt[cur_vi] == 0);

			if (removed[cur_vi] )
				continue;

			removed[cur_vi] = true;

			// decrease cur_vi's neighbors' ref count
			m_remained_MC_nbsCnt[cur_nbi]--;
			// if cur nb's ref cnt is 1, that means it is the new boundary. 
			// find the corresponding adjacent edge of cur nb and add it to q
			if (m_remained_MC_nbsCnt[cur_nbi] == 1)
			{
				const auto& further_nbs = m_MC_vts_neighbors[cur_nbi];
				for (auto nb_it = further_nbs.begin(); nb_it != further_nbs.end(); ++nb_it)
				{
					if (!removed[*nb_it] && 
						std::min(
						cur_medialCurve_distMetric[cur_nbi],
						cur_medialCurve_distMetric[*nb_it]) < _t )
					{
						q.push(std::make_pair(cur_nbi, *nb_it));
						break;
					}
				}
			}
			// if cur nb's ref cnt is 0, that means it should be removed
			if (m_remained_MC_nbsCnt[cur_nbi] == 0)
			{
				removed[cur_nbi] = true;
			}
		}

		for (unsigned ei = 0; ei < stg->dual_edges.size(); ++ei)
		{
			const auto& e = stg->dual_edges[ei];
			if (removed[e[0]] || removed[e[1]])
				continue;
			m_remained_MC.push_back(ei);
		}
		//cout << "# remained MC edges: " << m_remained_MC.size() << endl;
	}
	else
	{
		for (unsigned ei = 0; ei < stg->dual_edges.size(); ++ei)
		{
			const auto& e = stg->dual_edges[ei];
			/*if (pruned_edges.count(util::makeEdge(e[0], e[1])) == 0)
			_remain_edges.push_back(e);*/
			if (cur_medialCurve_distMetric[e[0]] >= _t &&
				cur_medialCurve_distMetric[e[1]] >= _t)
			{
				m_remained_MC.push_back(ei);
			}
			else // e is discarded. modify the nbs cnt for remainied vts
			{
				m_remained_MC_nbsCnt[e[0]] -- ;
				m_remained_MC_nbsCnt[e[1]] -- ;
				assert(m_remained_MC_nbsCnt[e[0]] >= 0);
				assert(m_remained_MC_nbsCnt[e[1]] >= 0);
			}
		}
	}

	// debug: print # open ends after this prunning
	/*int cnt_openEnds = 0;
	for (auto it = m_remained_MC_nbsCnt.begin(); it != m_remained_MC_nbsCnt.end(); ++it)
	{
		if (*it == 1)
			cnt_openEnds++;
	}
	cout << "# open ends: "<<cnt_openEnds<<endl;*/

	return true;
}

bool GLArea::pruneMedialCurveConstrainedByIsoContour(double _t)
{
	if (m_MC_vts_neighbors.empty() || this->m_iso_cont == nullptr)
		return false;

	vector<bool> removed(stg->dual_vts.size(), false);
	set<TriEdge> pruned_edges;
	m_remained_MC.clear();
	m_remained_MC.reserve(stg->dual_edges.size() - pruned_edges.size());
	const auto& fine_tri_for_dual_edge = this->stg->m_fromFineTri_for_dualE;
	const auto& fine_tri_status = this->m_iso_cont->m_face_status_static;
	for (unsigned ei = 0; ei < stg->dual_edges.size(); ++ei)
	{
		const auto& e = stg->dual_edges[ei];
		/*if (pruned_edges.count(util::makeEdge(e[0], e[1])) == 0)
			_remain_edges.push_back(e);*/
		if ( 
			cur_medialCurve_distMetric[e[0]] >= _t &&
			cur_medialCurve_distMetric[e[1]] >= _t && 
			fine_tri_status[fine_tri_for_dual_edge[ei]] != IsoContourOnMA::AHEAD 
			)
			m_remained_MC.push_back(ei);
	}
	return true;
}

void GLArea::printAndVisBurnTime()
{
	vector<size_t> vts_indices;
	this->stg->printBurnTime(vts_indices);
	vector<TriPoint> vts;
	for (auto it = vts_indices.begin(); it != vts_indices.end(); ++it)
	{
		vts.push_back(stg->m_origG->vts[*it]);
	}
	auto sph_drawer = std::dynamic_pointer_cast<SphereDrawer>(this->m_MPDrawer);
	float r = 1.0f;
	m_meshMA->need_bbox();
	r = m_meshMA->bbox.radius() * 0.02f;
	sph_drawer->setCenter(vts, r);

	updateGL();
}

void GLArea::visSphere(size_t _vid)
{
	vector<TriPoint> vts;
	vts.push_back(stg->m_origG->vts[_vid]);

	auto sph_drawer = std::dynamic_pointer_cast<SphereDrawer>(this->m_MPDrawer);
	float r = 1.0f;
	m_meshMA->need_bbox();
	r = m_meshMA->bbox.radius() * 0.02f;
	sph_drawer->setCenter(vts, r);

	updateGL();
}

void GLArea::visSphere(vector<unsigned> _vid_list)
{
	vector<TriPoint> vts;
	for (auto it = _vid_list.begin(); it != _vid_list.end(); ++it)
		vts.push_back(stg->m_origG->vts[*it]);

	auto sph_drawer = std::dynamic_pointer_cast<SphereDrawer>(this->m_MPDrawer);
	float r = 1.0f;
	m_meshMA->need_bbox();
	r = m_meshMA->bbox.radius() * 0.02f;
	sph_drawer->setCenter(vts, r);

	updateGL();
}

void GLArea::printMAOriginalVert(size_t _vid)
{
	cout << "----vertex "<<_vid<<" info----"<<endl;
	cout << "topo type: "<<topo_type_str[stg->getTopoType(_vid)]<<endl;
	cout << "bt2: "<<stg->bt2MA_vert[_vid]<<endl;
}

void GLArea::printMAOriginalVerts(vector<unsigned> _vid_list)
{
	m_meshMA->need_bbox();
	cout << "MA mesh bounding volume r: "<<m_meshMA->bbox.radius()<<endl;
	float avg_len = 0.0f;
	for (auto it = stg->m_origG->edges.begin(); it != stg->m_origG->edges.end(); ++it)
	{
		const auto& u = stg->m_origG->vts[(*it)[0]];
		const auto& v = stg->m_origG->vts[(*it)[1]];
		avg_len += trimesh::len(u-v);
	}
	avg_len /= stg->m_origG->edges.size();
	cout << "MA mesh edges' average length: "<<avg_len<<endl;

	std::streamsize ss = std::cout.precision();
	auto flgs = cout.flags();
	cout.setf(ios::left, ios::adjustfield);
	cout.precision(numeric_limits<float>::max_digits10);
	for (size_t i = 0; i < _vid_list.size(); ++i)
	{
		auto vi = _vid_list[i];
		cout << "sample vert ";
		cout << setw(8) << vi;
		cout <<" ("<< topo_type_str[stg->getTopoType(vi)]<<") burnt time\t";
		cout <<stg->bt2MA_vert[vi]<<endl;
	}
	for (size_t i = 0; i < _vid_list.size(); ++i)
	{
		auto vi = _vid_list[i];
		cout << vi<<": "<<stg->m_origG->vts[vi]<<endl;
	}
	cout.precision(ss);
	cout.setf(flgs);
}

bool GLArea::uploadLinesToDrawer(const vector<TriPoint>& _vts, const vector<TriEdge>& _edges, 
	shared_ptr<LineDrawer>& _line_drawer) const
{
	try {
		/* give each edge a color by duplicating vts */
		vector<TriPoint> vts;
		vector<TriEdge> edges(_edges.size());
		vts.reserve(_edges.size()*2);
		for (unsigned i = 0; i < _edges.size(); ++i)
		{
			auto& e = _edges[i];
			vts.push_back(_vts[e[0]]);
			vts.push_back(_vts[e[1]]);
			edges[i] = TriEdge(vts.size()-2, vts.size()-1);
		}

		// upload line data to GPU
		/*cout << "showRemainedMedialCurve(): uploading to GPU for render... " << endl;*/
		shared_ptr<LineDrawer> line_drawer = std::dynamic_pointer_cast<LineDrawer>(_line_drawer);
		line_drawer->setPoints(vts);
		line_drawer->setLines(edges);

		return true;
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
	return false;
}

bool GLArea::uploadMCCompleteGeometry(const vector<TriPoint>& _vts, const vector<TriEdge>& _edges)
{
	if (!stg)
	{
		this->m_drawMC = false;
		return false;
	}

	uploadLinesToDrawer(_vts, _edges, std::dynamic_pointer_cast<LineDrawer>(m_dualLinesDrawer));

	vector<TriPoint> vts;
	vts.reserve(_edges.size()*2);
	for (unsigned i = 0; i < _edges.size(); ++i)
	{
		auto& e = _edges[i];
		vts.push_back(_vts[e[0]]);
		vts.push_back(_vts[e[1]]);
	}
	shared_ptr<LineDrawer> lineProxy_drawer = std::dynamic_pointer_cast<LineDrawer>(m_linesProxyDrawer);
	lineProxy_drawer->setPoints(vts);

	//this->m_drawMC = true;
	return true;
}

// pre: the vts data of lines are already residing in GPU, 
// with shared vts duplicated
bool GLArea::uploadRemainedMCtoProxyDrawer()
{
	try
	{
		/*vector<unsigned> adjacency;
		adjacency.reserve(4 * m_remained_MC.size());
		for (auto it = m_remained_MC.begin(); it != m_remained_MC.end(); ++it)
		{
			// in fact, use *fake* adjacency for edge <u, v>, {u, u, v, v}
			adjacency.push_back(2 * *it); 
			adjacency.push_back(2 * *it); 
			adjacency.push_back(2 * *it + 1); 
			adjacency.push_back(2 * *it + 1); 
		}

		dynamic_pointer_cast<LineDrawer>(m_linesProxyDrawer)->setLinesAdjacency(adjacency);*/

		//cout << "remained MC size: "<<m_remained_MC.size() << endl;
		std::dynamic_pointer_cast<LineDrawer>(m_linesProxyDrawer)->setLines(m_remained_MC);

		return true;
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

	return false;
}

bool GLArea::uploadLineSubset(const vector<unsigned>& _edge_indices, const vector<TriEdge>& _edges)
{
	if (!stg)
	{
		this->m_drawMC = false;
		return false;
	}

	std::dynamic_pointer_cast<LineDrawer>(m_dualLinesDrawer)->setLines(_edge_indices);

	this->resizeGL(width(), height());
	//updateGL();

	return true;
}

void GLArea::printRemainedMCStats()
{
	if (m_remained_MC.empty())
		return;

	set<unsigned> vts;
	for (unsigned i = 0; i < m_remained_MC.size(); ++i)
	{
		auto& e = stg->dual_edges[m_remained_MC[i]];
		vts.insert(e[0]);
		vts.insert(e[1]);
	}

	// compute average, std. dev, and range of cur distance value
	float min_val = numeric_limits<float>::max();
	float max_val = -1.0f;
	float avg = 0.0f;
	float std_dev = 0.0f;
	for (auto vit = vts.begin(); vit != vts.end(); ++vit)
	{
		float cur_val = this->cur_MC_visDist[*vit];
		min_val = std::min(min_val, cur_val);
		max_val = std::max(max_val, cur_val);
		avg += cur_val;
	}
	avg /= vts.size();
	for (auto vit = vts.begin(); vit != vts.end(); ++vit)
	{
		std_dev += std::pow(cur_MC_visDist[*vit] - avg, 2.0f);
	}
	std_dev = std::sqrt(std_dev / vts.size());

	cout << "stats of cur dist field on MC: "<<endl;
	cout << "range: ["<<min_val<<","<<max_val<<"]"<<endl;
	cout << "avg: "<<avg<<endl;
	cout << "std dev: "<<std_dev<<endl;
}

void GLArea::setDrawFlag(GLArea::DrawFlag _obj_to_draw, bool _draw)
{
	switch(_obj_to_draw)
	{
	case GLArea::DRAW_POINTS:
		m_drawPoints = _draw;
		break;
	case GLArea::DRAW_ORIG:
		m_drawOrig = _draw;
		break;
	case GLArea::DRAW_MA:
		m_drawMA = _draw;
		//m_drawMALines = _draw;
		break;
	case GLArea::DRAW_MA_LINE:
		m_drawMALines = _draw;
		break;
	case GLArea::DRAW_MA_FINE:
		m_drawMAFinnerStatic = _draw;		
		break;
	case GLArea::DRAW_MC:
		m_drawMC = _draw;
		break;
	case GLArea::DRAW_BURNT_EDGES:
		m_drawBurntEdges = _draw;
		break;
	case GLArea::DRAW_HS:
		m_drawHS = _draw;
		break;
	case GLArea::DRAW_MP:
		m_drawMP = _draw;
		break;
	case GLArea::DRAW_ISOSURF:
		m_drawIsoSurf = _draw;
		break;
	case GLArea::DRAW_ISOCONT:
		m_drawIsoCont = _draw;
		break;
	case GLArea::DRAW_QMAT:
		m_drawQMAT = _draw;
		break;
	default:
		break;
	}

	this->resizeGL(width(), height());
	updateGL();
}

void GLArea::drawMeshSubset(
	std::shared_ptr<MeshDrawer> _mesh_drawer,
	const vector<unsigned>& _face_indices)
{
	_mesh_drawer->setFaces(_face_indices);

	this->resizeGL(width(), height());
	updateGL();
}

bool GLArea::burn(SteinerGraph::BurnScheme _burn_sch, SubdivScheme _subd_sch, double _steiner_param, int _edgeWeight_idx)
{
	if (!this->m_meshMA)
		return false;

#ifdef PROFILE_SPEED
	clock_t t_start, t_duration;
#endif

	std::shared_ptr<MyGraph> graph_ptr(new MyGraph(m_meshMA->vertices, m_meshMA->lines, m_meshMA->faces));

	// convert steiner param to meaningful value
	if (_subd_sch == SteinerSubdivision::ADAPTIVE || 
		_subd_sch == SteinerSubdivision::ADAPTIVE_MIDPOINT)
	{
		m_meshOrig->need_bbox();
		auto boxdif = m_meshOrig->bbox.max - m_meshOrig->bbox.min;
		auto max_len = std::max(boxdif[0], std::max(boxdif[1],boxdif[2]));
		//_steiner_param = m_meshMA->bbox.radius() * _steiner_param;
		_steiner_param = 0.5f * max_len * _steiner_param;
	}

#ifdef PROFILE_SPEED
	t_start = clock();
#endif
	// clean other structs before cleaning stg.
	m_surfF->reset();
	m_hs->reset();
	stg.reset();
	//_burn_sch = SteinerGraph::STEINER_ONLY;
	if ( _burn_sch == SteinerGraph::STEINER_ONLY )
		cout << "burn-scheme: steiner-only." << endl;
	stg = shared_ptr<SteinerGraph>( new SteinerGraph(
		graph_ptr,
		radii/*vector<float>()*/,
		_subd_sch,
		_steiner_param, // # steiner points placed on tri edge according to _scheme
		_edgeWeight_idx,
		_burn_sch
		) );
#ifdef PROFILE_SPEED
	t_duration = clock() - t_start;
	cout << "PROFIELING SPEED: SteinerGraph() took " << t_duration * 1000.0f / CLOCKS_PER_SEC << " ms."<<endl;
#endif

	cout << "estimating radii for all faces..."<<endl;
	stg->computeDist2SurfForFaces(radii);
	cout << "done." << endl;
	cout << "computing difference measure for faces..."<<endl;
	stg->computeDiffDist();
	cout << "done." << endl;

	// stg is created, thus init the hybrid skeleton
	// - setup basic info to initialize hybrid skeleton 
	m_hs->setup(stg);
	// and surfFunc obj with it.
	m_surfF->setup(stg, m_meshOrig);

	if (stg->bt2MA_vert.empty())
		return false;

	// upload finner triangulation of MA before interactive session like iso-countour begins
	cout << "uploading finner MA triangulation data..."<<endl;
	uploadFinerMAStaticGeom();
	cout << "finner MA upload done!" <<endl;

	this->resizeGL(width(), height());
	updateGL();

	return true;
}

void GLArea::uploadFinerMAStaticGeom()
{
	auto drawer = std::dynamic_pointer_cast<MeshDrawer>(m_FinerMAStaticDrawer);
	drawer->setRenderMode(MeshDrawer::PER_VERT);
	std::dynamic_pointer_cast<MeshDrawer>(m_MAFinnerDynamicDrawer)->setRenderMode(MeshDrawer::PER_VERT);
	
	// upload vts & faces of the finner tessellation of MA
	vector<TriPoint> vts;
	vts.reserve(stg->m_MA_finer_tris.size() * 3);
	vector<float> scalar_field;
	scalar_field.reserve(stg->m_MA_finer_tris.size() * 3);
	vector<TriFace> faces(stg->m_MA_finer_tris.size());
	for (int fi = 0; fi < stg->m_MA_finer_tris.size(); ++fi)
	{
		const auto& f = stg->m_MA_finer_tris[fi].first;
		auto fi_MA = stg->m_MA_finer_tris[fi].second;
		// duplicate shared vertices
		for (int j = 0; j < 3; ++j)
		{
			int vi = f[j];
			TriPoint p = stg->m_stSubdiv.getAllVts()[vi];
			vts.push_back(p);
			auto tei = stg->mapTopo(vi, fi_MA);
			assert(tei >= 0);
			scalar_field.push_back( stg->bt2MA_vert_per_sheet[vi][tei] );
		}
		faces[fi] = util::makeFace(vts.size()-1, vts.size()-2, vts.size()-3);
	}
	drawer->setPoints(vts);
	drawer->setFaces(faces, false);

	// now set colors and normals info into drawer 
	float scalar;
	TriColor vert_color;
	trimesh::vec3 face_normal;
	float* color_data = new float[ vts.size()*3 ];
	float* normal_data = new float[ vts.size()*3 ];
	auto min_max_scalar = std::minmax_element(scalar_field.begin(), scalar_field.end());
	for (int fi = 0; fi < faces.size(); ++fi)
	{
		const auto& f = faces[fi];
		face_normal = trimesh::normalize(trimesh::trinorm(vts[f[0]], vts[f[1]], vts[f[2]]));
		// duplicate shared vertices
		for (int j = 0; j < 3; ++j)
		{
			int vi = f[j];
			
			// per vert color according to properly selected scalar value 
			scalar = scalar_field[vi];
			vert_color = util::GetColour( scalar, *min_max_scalar.first, *min_max_scalar.second );
			color_data[vi * 3 + 0] = vert_color[0];
			color_data[vi * 3 + 1] = vert_color[1];
			color_data[vi * 3 + 2] = vert_color[2];

			normal_data[3*vi + 0] = face_normal[0];
			normal_data[3*vi + 1] = face_normal[1];
			normal_data[3*vi + 2] = face_normal[2];
		}
	}
	drawer->setPerVertColor( color_data, vts.size() );
	drawer->setPerVertNormal( normal_data, vts.size() );

	delete [] normal_data;
	delete [] color_data;
}

void GLArea::uploadMAFinerStaticColors(
	vector<vector<float>>& _scalar_per_sheet, float _min, float _max
	)
{
	unsigned vts_size = stg->m_MA_finer_tris.size() * 3;
	float scalar;
	TriColor vert_color;
	float* color_data = new float[ vts_size * 3 ]; // to be uploaded to gpu
	for (int fi = 0; fi < stg->m_MA_finer_tris.size(); ++fi)
	{
		const auto& f = stg->m_MA_finer_tris[fi].first;
		auto fi_MA = stg->m_MA_finer_tris[fi].second;
		// assign a color for each vert for each face (shared vertices thus are duplicated)
		for (int j = 0; j < 3; ++j)
		{
			int vi = f[j];
			auto tei = stg->mapTopo( vi, fi_MA );
			assert( tei >= 0 );
			scalar = _scalar_per_sheet[ vi ][ tei ];
			if ( scalar == numeric_limits<float>::max() )
			{
				vert_color = TriColor( 0.0f, 0.0f, 0.0f ); 
			}
			else
			{
				vert_color = util::GetColour( scalar, _min, _max );
			}
			color_data[ vi * 3 + 0 ] = vert_color[ 0 ];
			color_data[ vi * 3 + 1 ] = vert_color[ 1 ];
			color_data[ vi * 3 + 2 ] = vert_color[ 2 ];
		}
	}

	auto drawer = std::dynamic_pointer_cast<MeshDrawer>(m_FinerMAStaticDrawer);
	drawer->setRenderMode(MeshDrawer::PER_VERT);
	drawer->setPerVertColor( color_data, vts_size );

	delete [] color_data;
}

void GLArea::uploadMAFinerStaticColors(const TriColor& _c)
{
	if (!stg)
		return;

	unsigned vts_size = stg->m_MA_finer_tris.size() * 3;
	TriColor vert_color;
	float* color_data = new float[ vts_size*3 ];
	for (int fi = 0; fi < stg->m_MA_finer_tris.size(); ++fi)
	{
		const auto& f = util::makeFace(3*fi+0, 3*fi+1, 3*fi+2);
		for (int j = 0; j < 3; ++j)
		{
			int vi = f[j];

			vert_color = _c;
			color_data[vi * 3 + 0] = vert_color[0];
			color_data[vi * 3 + 1] = vert_color[1];
			color_data[vi * 3 + 2] = vert_color[2];
		}
	}

	auto drawer = std::dynamic_pointer_cast<MeshDrawer>(m_FinerMAStaticDrawer);
	drawer->setRenderMode(MeshDrawer::PER_VERT);
	drawer->setPerVertColor( color_data, vts_size );
	delete [] color_data;
}

void GLArea::uploadFinerMADynamicGeom( const vector<TriPoint>& _vts_for_faces )
{
	auto dyn_drawer = std::dynamic_pointer_cast<MeshDrawer>(m_MAFinnerDynamicDrawer);
	dyn_drawer->setPoints(_vts_for_faces);
	vector<TriFace> subFine_faces;
	subFine_faces.reserve(_vts_for_faces.size() / 3);
	for ( unsigned vi = 0; vi < _vts_for_faces.size(); vi += 3 )
	{
		subFine_faces.push_back( util::makeFace(vi, vi+1, vi+2) );
	}
	dyn_drawer->setFaces(subFine_faces, false);

	// upload normal data
	float scalar;
	trimesh::vec3 face_normal;
	float* normal_data = new float[ _vts_for_faces.size()*3 ];
	for (int fi = 0; fi < _vts_for_faces.size()/3; ++fi)
	{
		auto f = util::makeFace(3*fi+0, 3*fi+1, 3*fi+2);
		face_normal = trimesh::normalize(
			trimesh::trinorm(_vts_for_faces[f[0]], _vts_for_faces[f[1]], _vts_for_faces[f[2]])
			);

		for ( int j = 0; j < 3; ++j )
		{
			normal_data[f[j] * 3 + 0] = face_normal[0];
			normal_data[f[j] * 3 + 1] = face_normal[1];
			normal_data[f[j] * 3 + 2] = face_normal[2];
		}
	}

	dyn_drawer->setPerVertNormal( normal_data, _vts_for_faces.size() );
	delete [] normal_data;
}

void GLArea::uploadFinerMADynamicColors(
	const vector<TriPoint>& _vts_for_subFineFaces, 
	const TriColor& _const_c
	)
{
	auto dyn_drawer = std::dynamic_pointer_cast<MeshDrawer>(m_MAFinnerDynamicDrawer);
	// upload color
	float scalar;
	TriColor vert_color;
	float* color_data = new float[ _vts_for_subFineFaces.size()*3 ];
	for (int fi = 0; fi < _vts_for_subFineFaces.size()/3; ++fi)
	{
		auto f = util::makeFace(3*fi+0, 3*fi+1, 3*fi+2);
		for ( int j = 0; j < 3; ++j )
		{
			vert_color = _const_c;
			color_data[f[j] * 3 + 0] = vert_color[0];
			color_data[f[j] * 3 + 1] = vert_color[1];
			color_data[f[j] * 3 + 2] = vert_color[2];
		}
	}

	dyn_drawer->setPerVertColor( color_data, _vts_for_subFineFaces.size() );
	delete color_data;
}
void GLArea::uploadFinerMADynamicColors(
	const vector<TriPoint>& _vts_for_faces, 
	const vector<float>& _vals_for_facesVts, 
	const float _scalar_min, const float _scalar_max
	)
{
	auto dyn_drawer = std::dynamic_pointer_cast<MeshDrawer>(m_MAFinnerDynamicDrawer);
	// upload color
	float scalar;
	TriColor vert_color;
	float* color_data = new float[ _vts_for_faces.size()*3 ];
	for (int fi = 0; fi < _vts_for_faces.size()/3; ++fi)
	{
		auto f = util::makeFace(3*fi+0, 3*fi+1, 3*fi+2);
		for ( int j = 0; j < 3; ++j )
		{
			vert_color = util::GetColour( _vals_for_facesVts[f[j]], _scalar_min, _scalar_max );
			color_data[f[j] * 3 + 0] = vert_color[0];
			color_data[f[j] * 3 + 1] = vert_color[1];
			color_data[f[j] * 3 + 2] = vert_color[2];
		}
	}

	dyn_drawer->setPerVertColor( color_data, _vts_for_faces.size() );
	delete color_data;
}

bool GLArea::cleanTopo(std::string& _out_file /* = "" */, bool& _hasunburnt)
{
	if ( !m_meshMA || !stg )
		return false; // check failed.

	vector<int> bad_faces = stg->checkUnburnt();	
	if ( bad_faces.empty() )
	{
		_hasunburnt = false;
		return true; // check succ.ed
	}

	_hasunburnt = true;
	vector<bool> face_flags(m_meshMA->faces.size(), false);
	for (auto itor = bad_faces.begin(); itor != bad_faces.end(); ++itor)
		face_flags[*itor] = true;
	trimesh::remove_faces(m_meshMA.get(), face_flags);

	// set again the content that meshDrawer will draw since faces have changed
	std::dynamic_pointer_cast<MeshDrawer>(m_MADrawer)->setMeshToDraw(m_meshMA);

	if (_out_file == "")
		_out_file = m_medialAxisFile;
	unsigned index = _out_file.find_last_of(".");
	_out_file.replace(index, 1, ".clean.");

	cout << "writing the topo-cleaned mesh to file: " << _out_file << endl;
	auto meshMA_cpy = *m_meshMA;
	auto trans_cpy = m_trans_mat;
	trimesh::invert(trans_cpy);
	trimesh::apply_xform(&meshMA_cpy, trans_cpy);
	//meshMA_cpy.write(_out_file);
	MyMesh::write( _out_file, &meshMA_cpy );
	cout << "writing done." << endl;

	return true;
}

void GLArea::precomputeForMeasures()
{
	if (!stg)
		return;
//	cout << "estimating radii for all MA points (steiner points)..."<<endl;
//	stg->computeDist2SurfForFaces(radii);
//	cout << "done." << endl;
//// 	vector<TriPoint> closest2pts;
//// 	find2Closest(m_meshMA, m_meshOrig, closest2pts);
//// 	cout << "Done." << endl;
//// 	stg->computeDist2SurfForFaces(closest2pts);
//	cout << "computing difference measure..."<<endl;
//	stg->computeDiffDist();
//	cout << "done." << endl;
	
	cout << "estimating foot points for all MA faces..."<<endl;
	stg->computeIntersectionsForFaces(m_meshOrig, radii);
	cout << "done." << endl;

	cout << "setting up diffuser for MA face scalar field..."<<endl;
	stg->setupDiffuserForMA(m_diffuser_MA);
	cout << "diffuser setup."<<endl;
}

void GLArea::getFaceDistMetricMA(
	FaceFieldType _face_field, 
	bool _do_face_diffusion,	
	vector<float>& _dist_per_face)
{
	_dist_per_face.assign( stg == nullptr ? m_meshMA->faces.size() : stg->m_origG->faces.size(), -1.0f );
	switch ( _face_field )
	{
	case BT2_f: case BT3_f: case DIFF_f: case DIFF_REL_f:
		//precomputeForMeasures();
		switch ( _face_field )
		{
		case BT2_f:
			_dist_per_face = stg->bt2MA_face;
			break;
		case BT3_f:
			_dist_per_face = stg->bt3MA_face;
			break;
		case DIFF_f:
			_dist_per_face = stg->bt2bt3diffMA_face;
			break;
		case DIFF_REL_f:
			_dist_per_face.resize( stg->bt2bt3diffMA_face.size() );
			for ( unsigned i = 0; i < _dist_per_face.size(); ++i )
				_dist_per_face[ i ] = !( stg->bt2MA_face[ i ] > 0.0f ) ?
				1.0f : stg->bt2bt3diffMA_face[ i ] / stg->bt2MA_face[ i ];
			break;
		}
		break;
	case ANGLE_f: case LAMBDA_f: case GEODESIC_f:
	{
		/*vector<TriPoint> isects;
				vector<bool> has_isects;
				stg->computeIntersectionsForFaces(this->m_meshOrig, radii);*/

		if ( _face_field == ANGLE_f )
		{
			cout << "computing angles..." << endl;
			stg->computeAngleMetricForFaces( stg->m_footPtsForMAFaces, _dist_per_face );
			cout << "done" << endl;
		}
		else if ( _face_field == LAMBDA_f )
		{
			cout << "computing lambda metric..." << endl;
			stg->computeLambdaMetricForFaces( stg->m_footPtsForMAFaces, _dist_per_face );
			cout << "done" << endl;
		}
		else
		{
			cout << "computing geodesic metric..." << endl;
			QString cache = QString( this->m_medialAxisFile.c_str() ).remove( ".clean" ).replace( ".off", "_geoDCache.txt" );
			stg->computeGeodMetricForFaces( this->m_meshOrig, stg->m_footPtsForMAFaces, _dist_per_face, cache.toStdString() );
			cout << "done" << endl;
		}

		// debug begin:
		cout << "getFaceDistMetricMA (per-face), before diffusion: " << endl;
		auto min_f_val = *_dist_per_face.begin();
		auto max_f_val = *_dist_per_face.begin();
		for ( unsigned fi = 0; fi < _dist_per_face.size(); ++fi )
		{
			if ( _dist_per_face[ fi ] < 0.0f )
				continue;
			min_f_val = std::min( min_f_val, _dist_per_face[ fi ] );
			max_f_val = std::max( max_f_val, _dist_per_face[ fi ] );
		}
		std::cout << "face field range: " << min_f_val << "," << max_f_val << endl;
		// debug end:

		if ( _do_face_diffusion )
		{
			m_diffuser_MA.solve( _dist_per_face );

			// debug begin:
			cout << "getFaceDistMetricMA (per-face), after diffusion: " << endl;
			min_f_val = *_dist_per_face.begin();
			max_f_val = *_dist_per_face.begin();
			for ( unsigned fi = 0; fi < _dist_per_face.size(); ++fi )
			{
				if ( _dist_per_face[ fi ] < 0.0f )
					continue;
				min_f_val = std::min( min_f_val, _dist_per_face[ fi ] );
				max_f_val = std::max( max_f_val, _dist_per_face[ fi ] );
			}
			std::cout << "face field range: " << min_f_val << "," << max_f_val << endl;
			// debug end:

			/*cout << "outputting diffusion system to Mathematica file..."<<endl;
			m_diffuser_MA.outputDiffuseSystem(_dist_per_face);
			cout << "Done: outputting diffusion system to Mathematica file"<<endl;*/
		}

		break;
	}
	case FILE_MSURE:
		_dist_per_face = m_meshMA->msure_face;
		break;
	}
}

void GLArea::getFaceDistMetricMA(
	FaceFieldType _face_field, 
	bool _do_face_diffusion,
	vector<vector<float> >& _vert_dist_per_sheet)
{
	cout << "face field: " << _face_field << endl;
	const auto& vts = m_meshMA->vertices;
	// by default BT2 per vert per sheet
	_vert_dist_per_sheet.assign(
		stg->bt2MA_vert_per_sheet.begin(), 
		//stg->bt2MA_vert_per_sheet.begin() + vts.size()
		stg->bt2MA_vert_per_sheet.end()
		);

	m_meshMA->need_bbox();
	auto convert_to_vert_field = [&] (
		const vector<float>& _per_face, vector<vector<float> >& _per_sheet
		) 
	{
		vector<int> all_assoc_faces; 
		vector<int> assoc_face_one_sheet;
		vector<TopoGraph::EdgeIdx> sheets;
		for (int i = 0; i < _per_sheet.size(); ++i)
		{
			auto cur_vert_per_sheet = _per_sheet[i];
			for (int tei = 0; tei < cur_vert_per_sheet.size(); ++tei)
			{
				all_assoc_faces.clear();
				sheets.clear();
				if ( stg->is_leading_topo_edge(i, (TopoGraph::EdgeIdx)tei) )
				{
					sheets.push_back( (TopoGraph::EdgeIdx)tei );
					
					// find all sheets (including tei) finalized by tei
					for ( int tej = 0; tej < cur_vert_per_sheet.size(); ++tej )
					{
						if ( stg->is_leading_topo_edge(i, (TopoGraph::EdgeIdx)tej) )
							continue;
						if ( stg->prev_tedge[i][tej] == (TopoGraph::EdgeIdx)tei )
						{
							sheets.push_back( (TopoGraph::EdgeIdx)tej );
						}
					}

					// find all faces assoc.ed with the sheets
					for ( auto it = sheets.begin(); it != sheets.end(); ++it)
					{
						stg->getAssocFaces(i, *it, assoc_face_one_sheet);
						all_assoc_faces.insert(
							all_assoc_faces.end(), 
							assoc_face_one_sheet.begin(), assoc_face_one_sheet.end()
							);
					}

					// average values of all assoc faces
					int valid_cnt = 0;
					float avg_val = 0.0f;
					for (int fi = 0; fi < all_assoc_faces.size(); ++fi)
					{
						if ( _per_face[all_assoc_faces[fi]] < 0.0f )
							continue;
						avg_val += _per_face[all_assoc_faces[fi]];
						valid_cnt ++;
					}
					avg_val /= std::max(1, valid_cnt);

					// all relevant sheets take on the same average value
					for ( auto it = sheets.begin(); it != sheets.end(); ++it )
					{
						cur_vert_per_sheet[*it] = avg_val;
					}
				}
			}
			_per_sheet[i] = cur_vert_per_sheet;
		}

		//// debug begin:
		//std::cout << "after conversion: "<<endl;
		//auto min_s_val = *( _per_sheet.begin()->begin() );
		//auto max_s_val = *( _per_sheet.begin()->begin() );
		//for (unsigned vi = 0; vi < _per_sheet.size(); ++vi)
		//{
		//	const auto& persheet_cur_vert = _per_sheet[vi];
		//	for (unsigned si = 0; si < persheet_cur_vert.size(); ++si)
		//	{
		//		if (persheet_cur_vert[si] < 0.0f)
		//			continue;
		//		min_s_val = std::min(min_s_val, persheet_cur_vert[si]);
		//		max_s_val = std::max(max_s_val, persheet_cur_vert[si]);
		//	}
		//}
		//std::cout << "per-sheet range: "<<min_s_val<<","<<max_s_val<<endl;
		//// debug end:
	};

	switch (_face_field)
	{
	case BT2_f: case BT3_f: case DIFF_f: case DIFF_REL_f:
		//precomputeForMeasures();
		switch (_face_field)
		{
		case BT2_f:
			break;
		case BT3_f:
			for (int i = 0; i < _vert_dist_per_sheet.size(); ++i)
			{
				_vert_dist_per_sheet[i].assign(_vert_dist_per_sheet[i].size(), stg->bt3MA_vert[i]);
			}
			break;
		case DIFF_f:
			for (int i = 0; i < _vert_dist_per_sheet.size(); ++i)
			{
				float bt3 = stg->bt3MA_vert[i];
				auto cur_vert_per_sheet = _vert_dist_per_sheet[i];
				for (int j = 0; j < cur_vert_per_sheet.size(); ++j)
				{
					float bt2_s = cur_vert_per_sheet[ j ];
					if ( bt2_s < stg->infiniteBurnDist() )
						cur_vert_per_sheet[ j ] = bt2_s - bt3;
					else
						cur_vert_per_sheet[ j ] = stg->infiniteBurnDist();
				}
				_vert_dist_per_sheet[i] = cur_vert_per_sheet;
			}
			break;
		case DIFF_REL_f:
			for (int i = 0; i < _vert_dist_per_sheet.size(); ++i)
			{
				float bt3 = stg->bt3MA_vert[i];
				auto cur_vert_per_sheet = _vert_dist_per_sheet[i];
				float bt2 = 1.0f;
				for (int j = 0; j < cur_vert_per_sheet.size(); ++j)
				{
					bt2 = cur_vert_per_sheet[j];
					cur_vert_per_sheet[j] = bt2 <= 0.0f ? 
						1.0f : (1.0f - bt3 / bt2);
				}
				_vert_dist_per_sheet[i] = cur_vert_per_sheet;
			}
			break;
		default:
			break;
		}
		break;
	case ANGLE_f: case LAMBDA_f: case GEODESIC_f:
		{
			/*vector<TriPoint> isects;
			stg->computeIntersectionsForFaces(this->m_meshOrig, radii, &isects, NULL);*/
			vector<float> dist_per_face;

			if (_face_field == ANGLE_f)
			{
				cout << "computing angles..." << endl;
				stg->computeAngleMetricForFaces(stg->m_footPtsForMAFaces, dist_per_face);
				cout << "done" << endl;
			}
			else if (_face_field == LAMBDA_f)
			{
				cout << "computing lambda metric..." << endl;
				stg->computeLambdaMetricForFaces(stg->m_footPtsForMAFaces, dist_per_face);
				cout << "done" << endl;
			}
			else
			{
				cout << "computing geodesic metric..." << endl;
				QString cache = QString(this->m_medialAxisFile.c_str()).remove(".clean").replace(".off", "_geoDCache.txt");
				stg->computeGeodMetricForFaces(this->m_meshOrig, stg->m_footPtsForMAFaces, dist_per_face, cache.toStdString());
				cout << "done" << endl;
			}

			// debug begin:
			cout << "getFaceDistMetricMA (per-face), before diffusion: "<<endl;
			auto min_f_val = *dist_per_face.begin();
			auto max_f_val = *dist_per_face.begin();
			for (unsigned fi = 0; fi < dist_per_face.size(); ++fi)
			{
				if (dist_per_face[fi] < 0.0f)
					continue;
				min_f_val = std::min(min_f_val, dist_per_face[fi]);
				max_f_val = std::max(max_f_val, dist_per_face[fi]);
			}
			std::cout << "face field range: "<<min_f_val<<","<<max_f_val<<endl;
			// debug end:

			m_diffuser_MA.solve(dist_per_face);

			// debug begin:
			cout << "getFaceDistMetricMA (per-face), after diffusion: "<<endl;
			min_f_val = *dist_per_face.begin();
			max_f_val = *dist_per_face.begin();
			for (unsigned fi = 0; fi < dist_per_face.size(); ++fi)
			{
				if (dist_per_face[fi] < 0.0f)
					continue;
				min_f_val = std::min(min_f_val, dist_per_face[fi]);
				max_f_val = std::max(max_f_val, dist_per_face[fi]);
			}
			std::cout << "face field range: "<<min_f_val<<","<<max_f_val<<endl;

			convert_to_vert_field(dist_per_face, _vert_dist_per_sheet);

			break;
		}
	default:
		break;
	}
}

void GLArea::getPerSheetDistMetricMA(VertFieldType _vert_field, vector<vector<float>>& _dist_per_sheet)
{
	switch (_vert_field)
	{
	case BT2_v:
		_dist_per_sheet = stg->bt2MA_vert_per_sheet;
		break;
	case BT3_v:
		_dist_per_sheet.resize(stg->bt2MA_vert_per_sheet.size());
		for (auto it = _dist_per_sheet.begin(); it != _dist_per_sheet.end(); ++it)
		{
			auto idx = it - _dist_per_sheet.begin();
			it->resize( stg->bt2MA_vert_per_sheet[idx].size(), stg->bt3MA_vert[idx] );
		}
		break;
	case DIFF_v:
		_dist_per_sheet = stg->bt2MA_vert_per_sheet;
		for (auto it = _dist_per_sheet.begin(); it != _dist_per_sheet.end(); ++it)
		{
			auto idx = it - _dist_per_sheet.begin();
			auto d_bt3 = stg->bt3MA_vert[idx];
			for (unsigned j = 0; j < it->size(); ++j)
			{
				float d_bt2 = ( *it )[ j ];
				if ( d_bt2 < stg->infiniteBurnDist() )
					( *it )[ j ] = d_bt2 - d_bt3;
				else
					( *it )[ j ] = stg->infiniteBurnDist();
			}
		}
		break;
	case DIFF_REL_v:
		_dist_per_sheet = stg->bt2MA_vert_per_sheet;
		for (auto it = _dist_per_sheet.begin(); it != _dist_per_sheet.end(); ++it)
		{
			auto idx = it - _dist_per_sheet.begin();
			auto bt3_val = stg->bt3MA_vert[idx];
			for (unsigned j = 0; j < it->size(); ++j)
			{
				(*it)[j] = 1.0f - bt3_val / (*it)[j];
			}
		}
		break;
	}
}

void GLArea::getVertDistMetricMA(VertFieldType _vert_field, vector<float>& _dist_per_vert)
{
	switch (_vert_field)
	{
	case BT2_v:
		_dist_per_vert = stg->bt2MA_vert;
		break;
	case BT3_v:
		_dist_per_vert = stg->bt3MA_vert; //radii;
		break;
	case DIFF_v:
		_dist_per_vert = stg->bt2bt3diffMA_vert;
		break;
	case DIFF_REL_v:
		_dist_per_vert.resize(stg->bt2bt3diffMA_vert.size());
		for (unsigned i = 0; i < _dist_per_vert.size(); ++i)
			_dist_per_vert[i] = !(stg->bt2MA_vert[i] > 0.0f) ?
			1.0f : stg->bt2bt3diffMA_vert[i] / stg->bt2MA_vert[i];
		break;
	}
}

void GLArea::colorMAFaceBy(
	FaceFieldType _type, float& _min, float& _max, 
	bool _use_new_range, 
	float _min_alpha, float _alpha_exp, 
	int _update_opt,
	bool _do_face_diffusion,
	bool _use_per_sheet)
{
	/*if (!stg)
		return;*/
	
	auto MA_drawer = std::dynamic_pointer_cast<MeshDrawer>(m_MADrawer);
	if (_use_per_sheet)
	{
		// we are really going to interpolate color across each face 
		// we'll use different color for different face for the same vert 
		// if it's a non-manifold vert.
		//MA_drawer->setRenderMode(MeshDrawer::PER_VERT);

		// update IBO of drawer to draw correct sequence of vertices
		//MA_drawer->setFaces(m_meshMA->faces, true);
		//MA_drawer->setRenderMode(MeshDrawer::PER_FACE);
		TriFace f;
		bool burnt;

		vector<vector<float> > vertDist_per_sheet;
		getFaceDistMetricMA(_type, _do_face_diffusion, vertDist_per_sheet);

		//// debug begin:
		//cout << "colorMAFaceBy() (per-sheet): "<<endl;
		//auto min_s_val = *( vertDist_per_sheet.begin()->begin() );
		//auto max_s_val = *( vertDist_per_sheet.begin()->begin() );
		//for (unsigned vi = 0; vi < vertDist_per_sheet.size(); ++vi)
		//{
		//	const auto& persheet_cur_vert = vertDist_per_sheet[vi];
		//	for (unsigned si = 0; si < persheet_cur_vert.size(); ++si)
		//	{
		//		if (persheet_cur_vert[si] < 0.0f)
		//			continue;
		//		min_s_val = std::min(min_s_val, persheet_cur_vert[si]);
		//		max_s_val = std::max(max_s_val, persheet_cur_vert[si]);
		//	}
		//}
		//std::cout << "per-sheet range: "<<min_s_val<<","<<max_s_val<<endl;
		//// debug end:

		float min_val = numeric_limits<float>::max();
		float max_val = numeric_limits<float>::min();

		for (unsigned vi = 0; vi < vertDist_per_sheet.size(); ++vi)
		{
			if ( ( _type == BT2_f || _type == DIFF_f || _type == DIFF_REL_f ) 
				&& !stg->burnt[ vi ] )
				continue;
			const auto & cur_vert_per_sheet = vertDist_per_sheet[vi];
			if ( !cur_vert_per_sheet.empty() )
			{
				auto min_max = minmax_element( cur_vert_per_sheet.begin(), cur_vert_per_sheet.end() );
				min_val = std::min( min_val, *min_max.first );
				max_val = std::max( max_val, *min_max.second );
			}
		}

		cout << "Visualize: using face dist measure "<<_type
			<<", min_val & max_val: " << min_val<<", "<<max_val<< endl; //debug
		if (_use_new_range)
		{
			min_val = _min;
			max_val = _max;
			cout << "rescale enabled: " << min_val <<", "<< max_val <<endl;
		}

		if (!_use_new_range) // then update dist range
		{
			_min = min_val; 
			_max = max_val;
		}

		if (_update_opt == 1 || _update_opt == 3)
		{
			cout << "updating color..." << endl;
			float* color_data = new float[m_meshMA->faces.size()*3*3];
			TriColor color;
			float s;

			for (unsigned fi = 0; fi < m_meshMA->faces.size(); ++fi)
			{
				f = m_meshMA->faces[fi];

				// visualize topo type: only if burnt && dist value valid
				for (unsigned j = 0; j < 3; ++j)
				{
					auto vi = f[ j ];
					int tei = stg->mapTopo( vi, fi);
					//assert(tei >= 0);
					if (tei < 0)
					{
						cout << "Fatal Error: cannot find topo sheet containing face "<<fi<<endl;
						exit(-1);
					}
					s = vertDist_per_sheet[ vi ][tei];
					bool invalid_color = ( _type == BT2_f || _type == DIFF_f || _type == DIFF_REL_f )
						&& !stg->burnt[ vi ];

					color = invalid_color ? TriColor( 0.0f, 0.0f, 0.0f ) 
						: util::GetColour( s, min_val, max_val );

					color_data[ 9*fi+3*j+0 ] = color[0];
					color_data[ 9*fi+3*j+1 ] = color[1];
					color_data[ 9*fi+3*j+2 ] = color[2];
				}
			}
			std::dynamic_pointer_cast<MeshDrawer>(m_MADrawer)->setPerVertColor(color_data, m_meshMA->faces.size()*3);

			delete [] color_data;
			cout << "updating color done." << endl;
		}
		if (_update_opt == 2 || _update_opt == 3)
		{
			// saliency of each face. can be used as alpha in shader.
			float scalar, saliency;

			float* saliency_data = new float[m_meshMA->faces.size()*3];
			for (unsigned fi = 0; fi < m_meshMA->faces.size(); ++fi)
			{
				f = m_meshMA->faces[fi];

				// visualize topo type: only if burnt && dist value valid
				for (unsigned j = 0; j < 3; ++j)
				{
					auto vi = f[ j ];
					if ( ( _type == BT2_f || _type == DIFF_f || _type == DIFF_REL_f )
						&& !stg->burnt[ vi ] )
					{
						saliency = 1.0f;
					}
					else
					{
						int tei = stg->mapTopo( f[ j ], fi );
						assert( tei >= 0 );
						scalar = vertDist_per_sheet[ f[ j ] ][ tei ];
						saliency = util::rescale( scalar, min_val, max_val, _alpha_exp, _min_alpha, 1.0f );
					}
					saliency_data[3*fi+j] = saliency;
				}
			}
			std::dynamic_pointer_cast<MeshDrawer>(m_MADrawer)->setPerVertSaliency(saliency_data, m_meshMA->faces.size()*3);
			delete [] saliency_data;
		}
		cout << "uploading MA finer static colors... " << endl;
		// don't forget to update colors of fine MA triangulation in case fine pruning is to perform
		this->uploadMAFinerStaticColors(vertDist_per_sheet, min_val, max_val);
		cout << "Done: uploading MA finer static colors." << endl;
	}
	else
	{
		// we are really going to interpolate color across each face 
		// we'll use different color for different face for the same vert 
		// if it's a non-manifold vert.
		MA_drawer->setRenderMode(MeshDrawer::PER_FACE);

		// update IBO of drawer to draw correct sequence of vertices
		//MA_drawer->setFaces(m_meshMA->faces, true);
		//MA_drawer->setRenderMode(MeshDrawer::PER_FACE);
		TriFace f;
		bool burnt;

		getFaceDistMetricMA(_type, _do_face_diffusion, cur_MA_visDist_f);

		//// debug begin:
		//cout << "colorMAFaceBy() (per-face): "<<endl;
		//auto min_f_val = *cur_MA_visDist_f.begin();
		//auto max_f_val = *cur_MA_visDist_f.begin();
		//for (unsigned fi = 0; fi < cur_MA_visDist_f.size(); ++fi)
		//{
		//	if (cur_MA_visDist_f[fi] < 0.0f)
		//		continue;
		//	min_f_val = std::min(min_f_val, cur_MA_visDist_f[fi]);
		//	max_f_val = std::max(max_f_val, cur_MA_visDist_f[fi]);
		//}
		//std::cout << "face field range: "<<min_f_val<<","<<max_f_val<<endl;
		//// debug end

		float min_val = numeric_limits<float>::max();
		float max_val = numeric_limits<float>::min();
		for (unsigned fi = 0; fi < cur_MA_visDist_f.size(); ++fi)
		{
			if (cur_MA_visDist_f[fi] < 0.0f)
				continue;
			min_val = std::min(min_val, cur_MA_visDist_f[fi]);
			max_val = std::max(max_val, cur_MA_visDist_f[fi]);
		}
		/*for ( auto s_iter = cur_MA_visDist_f.begin(); s_iter != cur_MA_visDist_f.end(); ++s_iter )
		{
			if (*s_iter < 0.0)
				continue;
			min_val = std::min(*s_iter, min_val);
			max_val = std::max(*s_iter, max_val);
		}*/
		
		cout << "Visualize: using face dist measure "<<_type
			<<", min_val & max_val: " << min_val<<", "<<max_val<< endl; //debug
		if (_use_new_range)
		{
			min_val = _min;
			max_val = _max;
			cout << "rescale enabled: " << min_val <<", "<< max_val <<endl;
		}

		if (!_use_new_range) // then update dist range
		{
			_min = min_val; 
			_max = max_val;
		}

		if (_update_opt == 1 || _update_opt == 3)
		{
			float* color_data = new float[m_meshMA->faces.size()*3];
			TriColor color;
			float s;			
			for (unsigned fi = 0; fi < m_meshMA->faces.size(); ++fi)
			{
			const auto& f = m_meshMA->faces[fi];
			s = cur_MA_visDist_f[fi];

			color = s >= 0.0f ? util::GetColour(
			s, min_val, max_val
			) : TriColor(0.0f, 0.0f, 0.0f);
			color_data[ 3*fi+0 ] = color[0];
			color_data[ 3*fi+1 ] = color[1];
			color_data[ 3*fi+2 ] = color[2];
			}
			std::dynamic_pointer_cast<MeshDrawer>(m_MADrawer)->setPerFaceColor(color_data, m_meshMA->faces.size());
			

			delete [] color_data;
		}
		if (_update_opt == 2 || _update_opt == 3)
		{
			// saliency of each face. can be used as alpha in shader.
			float* saliency_data = new float[m_meshMA->faces.size()];
			float scalar, saliency;
			
			for (unsigned fi = 0; fi < m_meshMA->faces.size(); ++fi)
			{
			// assign saliency to this face
			saliency_data[fi] = 
			util::rescale(cur_MA_visDist_f[fi], min_val, max_val, 
			_alpha_exp, _min_alpha, 1.0f);
			}
			std::dynamic_pointer_cast<MeshDrawer>(m_MADrawer)->setPerFaceSaliency(saliency_data, m_meshMA->faces.size());
			
			delete [] saliency_data;
		}

		// don't forget to update colors of fine MA triangulation in case fine pruning is to perform
		//this->uploadMAFinnerStaticColors(vertDist_per_sheet, min_val, max_val);
	}

	updateGL();
}

void GLArea::colorMAVertBy(
	VertFieldType _type, float& _min, float& _max, 
	bool _use_new_range, 
	float _min_alpha, float _alpha_exp, 
	int _update_opt)
{
	vector<float> field_vals;

	getVertDistMetricMA(_type, field_vals);
	float min_val = numeric_limits<float>::max();
	float max_val = numeric_limits<float>::min();
	for (unsigned i = 0; i < field_vals.size(); ++i)
	{
		if (!stg->burnt[i])
			continue;
		min_val = std::min(min_val, (field_vals)[i]);
		max_val = std::max(max_val, (field_vals)[i]);
	}

	//cout << "Visualize: using vert dist measure "<<_type
	//	<<", min_val & max_val: " << min_val<<", "<<max_val<< endl; //debug
	if (_use_new_range)
	{
		/*min_val = std::max(min_val, _min);
		max_val = std::min(max_val, _max);*/
		min_val = _min;
		max_val = _max;
		cout << "clamped enabled: " << min_val <<", "<< max_val <<endl;
	}

	if (!_use_new_range)
	{
		_min = min_val;
		_max = max_val;
	}

	std::dynamic_pointer_cast<MeshDrawer>(m_MADrawer)->setRenderMode(MeshDrawer::PER_VERT);
	if (_update_opt == 1 || _update_opt == 3) // update color
	{
		float* color_data = new float[m_meshMA->vertices.size() * 3];
		TriColor color;

		for (unsigned i = 0; i < m_meshMA->vertices.size(); ++i)
		{
			// visualize topo type
			color = stg->burnt[i] ? util::GetColour(
				(field_vals)[i], min_val, max_val
				) : TriColor(1.0f, 0.0f, 1.0f);
			//TriColor color = util::GetColour(stg->distF[i], stg->min_dist, stg->max_dist);
			color_data[3*i+0] = color[0];
			color_data[3*i+1] = color[1];
			color_data[3*i+2] = color[2];
		}
		m_MADrawer->setPerVertColor(color_data, m_meshMA->vertices.size());
		delete [] color_data;
	}
	if (_update_opt == 2 || _update_opt == 3) // update saliency
	{
		// saliency of each vert. can be used as alpha in shader.
		float* saliency_data = new float[m_meshMA->vertices.size()];

		for (unsigned i = 0; i < m_meshMA->vertices.size(); ++i)
		{
			// assign saliency to this vert
			saliency_data[i] = util::rescale(field_vals[i], min_val, max_val, _alpha_exp, _min_alpha, 1.0f);
		}
		std::dynamic_pointer_cast<MeshDrawer>(m_MADrawer)->setPerVertSaliency(saliency_data, m_meshMA->vertices.size());
		delete [] saliency_data;
	}

	updateGL();
}

void GLArea::getMCDistMetric(DistMC _type, vector<float>& _dist) const
{
	switch (_type)
	{
	case BT3_MC:
		_dist = stg->bt3_MC;
		break;
	case BT2_MC:
		_dist = stg->bt2_MC;
		break;
	case BT2_BT3_MC:
		_dist.clear();
		_dist.reserve(stg->dual_vts.size());
		for (unsigned i = 0; i < stg->bt2_MC.size(); ++i)
		{
			float d_bt2 = stg->bt2_MC[ i ], d_bt3 = stg->bt3_MC[ i ];
			if ( d_bt2 < stg->infiniteBurnDist() )
				_dist.push_back( d_bt2 - d_bt3 );
			else
				_dist.push_back( stg->infiniteBurnDist() );
		}
		break;
	case BT2_BT3_REL_MC:
		_dist.clear();
		_dist.reserve(stg->dual_vts.size());
		for (unsigned i = 0; i < stg->bt2_MC.size(); ++i)
		{
			_dist.push_back(
				1.0f - stg->bt3_MC[i] / stg->bt2_MC[i]
			);
		}
		break;
	case BT1_MC:
		_dist = stg->bt1_medialCurve;
		break;
	case BT1_BT2_MC:
		_dist.clear();
		_dist.reserve(stg->dual_vts.size());
		for (unsigned i = 0; i < stg->bt1_medialCurve.size(); ++i)
		{
			float d_bt2 = stg->bt2_MC[ i ], d_bt1 = stg->bt1_medialCurve[ i ];
			if ( d_bt2 < stg->infiniteBurnDist() )
				_dist.push_back( std::max( 0.0f, d_bt1 - d_bt2 ) );
			else
				_dist.push_back( stg->infiniteBurnDist() );
		}
		break;
	case BT1_BT2_REL_MC:
		_dist.clear();
		_dist.reserve(stg->dual_vts.size());
		for (unsigned i = 0; i < stg->bt1_medialCurve.size(); ++i)
		{
			_dist.push_back(
				1.0f - stg->bt2_MC[i] / stg->bt1_medialCurve[i]
			);
		}
		break;
	}
}

void GLArea::getMCDistMetric(DistMC _type, vector<trimesh::vec2>& _dist_by_edge) const
{
	_dist_by_edge.clear();
	trimesh::vec2 dist_on_e(0.0f, 0.0f), dist2_on_e(0.0f, 0.0f);
	trimesh::ivec2 query;// <dual v id, adjacent face id>
	
	auto get_bt3 = [ & ]( const TriEdge& _e, trimesh::vec2& _dist_on_e )
	{
		_dist_on_e[ 0 ] = stg->bt3_MC[ _e[ 0 ] ];
		_dist_on_e[ 1 ] = stg->bt3_MC[ _e[ 1 ] ];
	};
	auto get_bt2 = [&] (int _ei, const TriEdge& _e, trimesh::vec2& _dist_on_e)
	{
		query[1] = stg->m_dualEdge_from_which_face[_ei]; // the face where this dual edge resides
		if ( query[1] >= 0 )
		{
			query[ 0 ] = _e[ 0 ];
			_dist_on_e[ 0 ] = stg->m_bt2_MC_perFace.find( query )->second;
			query[ 0 ] = _e[ 1 ];
			_dist_on_e[ 1 ] = stg->m_bt2_MC_perFace.find( query )->second;
		}
		else
		{
			get_bt3( _e, _dist_on_e );
		}
	};
	auto get_bt1 = [this] (const TriEdge& _e, trimesh::vec2& _dist_on_e)
	{
		// needs to explicitly consider dist of _e b.c. our bt1_MC doesn't record
		// bt1 on a per-edge basis (i.e. the similar way per-sheet bt is recorded for MA)
		float d = trimesh::dist(this->stg->dual_vts[_e[0]], this->stg->dual_vts[_e[1]]);
		if (stg->burnNext_medialCurve[_e[0]] == _e[1])
		{
			_dist_on_e[0] = stg->bt1_medialCurve[_e[0]];
			_dist_on_e[1] = stg->bt1_medialCurve[_e[0]] + d;
		}
		else
		{
			_dist_on_e[1] = stg->bt1_medialCurve[_e[1]];
			_dist_on_e[0] = stg->bt1_medialCurve[_e[1]] + d;
		}
		if ( !( ::is_valid( _dist_on_e[ 0 ] ) && _dist_on_e[ 0 ] < stg->infiniteBurnDist() ) )
			_dist_on_e = vec2( stg->infiniteBurnDist() );
	};

	for (unsigned ei = 0; ei < stg->dual_edges.size(); ++ei)
	{
		const auto& e = stg->dual_edges[ei];
		switch (_type)
		{
		case BT3_MC:
			get_bt3(e, dist_on_e);
			break;
		case BT2_MC:
			get_bt2(ei, e, dist_on_e);
			break;
		case BT2_BT3_MC: 
			get_bt2(ei, e, dist_on_e);
			get_bt3(e, dist2_on_e);
			if ( dist_on_e[ 0 ] < stg->infiniteBurnDist() )
				dist_on_e = std::max( 0.0f, dist_on_e - dist2_on_e );
			else
				dist_on_e = vec2( stg->infiniteBurnDist() );
			break;
		case BT2_BT3_REL_MC:
			get_bt2(ei, e, dist_on_e);
			get_bt3(e, dist2_on_e);
			dist2_on_e = trimesh::vec2(1.0f) - dist_on_e / dist2_on_e;
			break;
		case BT1_MC:
			get_bt1(e, dist_on_e);
			break;
		case BT1_BT2_MC:
			get_bt1(e, dist_on_e);
			get_bt2(ei, e, dist2_on_e);
			if ( dist_on_e[ 0 ] < stg->infiniteBurnDist() )
				dist_on_e = std::max( 0.0f, dist_on_e - dist2_on_e );
			else
				dist_on_e = vec2( stg->infiniteBurnDist() );
			break;
		case BT1_BT2_REL_MC:
			get_bt1(e, dist_on_e);
			get_bt2(ei, e, dist2_on_e);
			dist2_on_e = trimesh::vec2(1.0f) - dist_on_e / dist2_on_e;
			break;
		}
		if ( !( ::is_valid( dist_on_e[ 0 ] ) && dist_on_e[ 0 ] < stg->infiniteBurnDist() ) )
			dist_on_e = vec2( stg->infiniteBurnDist() );
		_dist_by_edge.push_back(dist_on_e);
	}
}

void GLArea::setupMAFacePrune( FaceFieldType _type, float& _min, float& _max, int _dist_i )
{
	/*auto& dist_field = _dist_i == 1 
		? cur_MA_pruneDist_f_1 : cur_MA_pruneDist_f_2;
	getFaceDistMetricMA(_type, dist_field);*/

	// begin: used for face-pruning by inspecting vertex values
	auto& dist_field = _dist_i == 1 
		? cur_MA_pruneDist_persheet_1 : cur_MA_pruneDist_persheet_2;
	getFaceDistMetricMA(_type, true, dist_field);
	// end: used for face-pruning by inspecting vertex values

	float min_val = numeric_limits<float>::max();
	float max_val = numeric_limits<float>::min();
	/*for (unsigned fi = 0; fi < dist_field.size(); ++fi)
	{
		TriFace f = m_meshMA->faces[fi];
		if (!stg->burnt[f[0]] || !stg->burnt[f[1]] || !stg->burnt[f[2]] || 
			dist_field[fi] < 0.0f)
			continue;
		min_val = std::min(min_val, dist_field[fi]);
		max_val = std::max(max_val, dist_field[fi]);
	}*/
	// begin: used for face-pruning by inspecting vertex values
	for (unsigned vi = 0; vi < dist_field.size(); ++vi)
	{
		for (int tei = 0; tei < dist_field[vi].size(); ++tei)
		{
			if (!stg->burnt[vi] || dist_field[vi][tei] < 0.0f)
				continue;
			min_val = std::min(min_val, dist_field[vi][tei]);
			max_val = std::max(max_val, dist_field[vi][tei]);
		}
	}
	// end: used for face-pruning by inspecting vertex values
	_min = min_val; 
	_max = max_val;

	cout << "Prune: using face dist measure "<<_type
		<<"[" << _min<<", "<<_max<<"]"<< endl; //debug
}

float GLArea::getMAPruneValue(float _ratio)
{
	if (!m_meshOrig)
		return -999.0f;

	m_meshOrig->need_bbox();
	auto d = m_meshOrig->bbox.size().max();
	float prune_val = d * _ratio;

	return prune_val;
}

float GLArea::toETmPruneRatio(float _v)
{
	if (!m_meshOrig)
		return -999.0f;

	m_meshOrig->need_bbox();
	return _v / m_meshOrig->bbox.size().max();
}

void GLArea::pruneMA(float _t1, float _t2)
{
	if (!stg)
		return;
	// do nothing if no threshold is valid
	if ( _t1 < 0 && _t2 < 0 )
		return;

	vector<unsigned> _faces;
	TriFace f;
	// begin: used for face-pruning by inspecting vertex values
	TopoGraph::EdgeIdx tei;
	// end: used for face-pruning by inspecting vertex values
	for (unsigned fi = 0; fi < stg->m_origG->faces.size(); ++fi)
	{
		stg->setFacePruned(fi, false);
		f = stg->m_origG->faces[fi];

		if ( !(_t1 < 0) ) // t1 valid
		{
			if ( !(_t2 < 0) ) // t2 also valid
			{
				/*if (cur_MA_pruneDist_f_1[fi] >= _t1 && cur_MA_pruneDist_f_2[fi] >= _t2)
					_faces.push_back(fi);
				else
					stg->setFacePruned(fi, true);*/

				// begin: used for face-pruning by inspecting vertex values
				int v_remain_cnt = 0;
				for (int i = 0; i < 3; ++i)
				{
					int ret = stg->mapTopo(f[i], fi);
					assert(ret != -1);
					tei = (TopoGraph::EdgeIdx)ret;
					if (cur_MA_pruneDist_persheet_1[f[i]][tei] >= _t1 && 
						cur_MA_pruneDist_persheet_2[f[i]][tei] >= _t2)
						v_remain_cnt ++;
				}
				if (v_remain_cnt > 0)
					_faces.push_back(fi);
				else
					stg->setFacePruned(fi, true);
				// end: used for face-pruning by inspecting vertex values
			}
			else // only t1 valid i.e. only compare to t1
			{
				/*if (cur_MA_pruneDist_f_1[fi] >= _t1)
					_faces.push_back(fi);
				else
					stg->setFacePruned(fi, true);*/

				// begin: used for face-pruning by inspecting vertex values
				int v_remain_cnt = 0;
				for (int i = 0; i < 3; ++i)
				{
					int ret = stg->mapTopo(f[i], fi);
					assert(ret != -1);
					tei = (TopoGraph::EdgeIdx)ret;
					if ( cur_MA_pruneDist_persheet_1[f[i]][tei] >= _t1 )
						v_remain_cnt ++;
				}
				if (v_remain_cnt > 0)
					_faces.push_back(fi);
				else
					stg->setFacePruned(fi, true);
				// end: used for face-pruning by inspecting vertex values
			}
		}
		else // only t2 valid
		{
			/*if (cur_MA_pruneDist_f_2[fi] >= _t2)
				_faces.push_back(fi);
			else
				stg->setFacePruned(fi, true);*/

			// begin: used for face-pruning by inspecting vertex values
			int v_remain_cnt = 0;
			for (int i = 0; i < 3; ++i)
			{
				int ret = stg->mapTopo(f[i], fi);
				assert(ret != -1);
				tei = (TopoGraph::EdgeIdx)ret;
				if ( cur_MA_pruneDist_persheet_2[f[i]][tei] >= _t1 )
					v_remain_cnt ++;
			}
			if (v_remain_cnt >= 2)
				_faces.push_back(fi);
			else
				stg->setFacePruned(fi, true);
			// end: used for face-pruning by inspecting vertex values
		}
	}

	cout << "MA faces left: " << _faces.size() << endl;
	this->drawMeshSubset(std::dynamic_pointer_cast<MeshDrawer>(m_MADrawer), _faces);
}

void GLArea::pruneFineMA( float _t1, float _t2 )
{
	if (!stg)
		return;

	vector<unsigned> faces_ahead;
	vector<unsigned> faces_crossed;
	TriFace f;
	TopoGraph::EdgeIdx tei;
	// stats about current prune measure
	float cur_prune_min = std::numeric_limits<float>::max();
	float cur_prune_max = 0.0f;

	for (auto fine_it = stg->m_MA_finer_tris.begin(); fine_it != stg->m_MA_finer_tris.end(); ++fine_it)
	{
		const auto& fine_f = fine_it->first; 
		auto fi_MA = fine_it->second;
		int v_above_cnt = 0;
		for ( int i = 0; i < 3; ++i )
		{
			auto vi = fine_f[i];
			int ret = stg->mapTopo(vi, fi_MA);
			assert(ret >= 0);
			tei = (TopoGraph::EdgeIdx)(ret);
			float cur_val = cur_MA_pruneDist_persheet_1[vi][tei];
			cur_prune_max = std::max(cur_val, cur_prune_max);
			cur_prune_min = std::min(cur_val, cur_prune_min);

			v_above_cnt += IsoContourExtractor::above( 
				cur_val, _t1 );
		}
		if ( v_above_cnt == 3 )
		{
			faces_ahead.push_back( fine_it - stg->m_MA_finer_tris.begin() );
		}
		else if ( v_above_cnt > 0 )
		{
			faces_crossed.push_back( fine_it - stg->m_MA_finer_tris.begin() );
		}
	}

	// generate dynamic part of the MA within those fine triangles crossed by iso-contour
	vector<TriPoint> vts_of_crossed_faces;
	vector<float> vals_at_vts;
	for ( auto cross_it = faces_crossed.begin(); cross_it != faces_crossed.end(); ++cross_it )
	{
		const auto& f = stg->m_MA_finer_tris[*cross_it].first;
		auto fi_MA = stg->m_MA_finer_tris[*cross_it].second;
		for ( int i = 0; i < 3; ++i )
		{
			vts_of_crossed_faces.push_back( stg->m_stSubdiv.getAllVts()[f[i]] );
			int ret = stg->mapTopo(f[i], fi_MA);
			assert (ret >= 0);
			tei = (TopoGraph::EdgeIdx)(ret);
			vals_at_vts.push_back(cur_MA_pruneDist_persheet_1[f[i]][tei]);
		}
	}
	vector<TriPoint> iso_cont_vts;
	vector<TriPoint> vts_for_subFineFaces;
	vector<float> vals_for_subFineFaces_vts;
	IsoContourExtractor iso_extractor;
	iso_extractor.marchTriangle(
		vts_of_crossed_faces, vals_at_vts, _t1, 
		iso_cont_vts, vts_for_subFineFaces, vals_for_subFineFaces_vts
		);

	// visualize the remaining part of fine-MA: static part & the dynamic part
	// this is static part
	std::dynamic_pointer_cast<MeshDrawer>(m_FinerMAStaticDrawer)->setFaces(faces_ahead, false);
	faces_ahead.clear();
	// this is dynamic part
	uploadFinerMADynamicGeom(vts_for_subFineFaces);
	// color for the dynamic part
	if ( !m_use_constColor_for_MA )
	{
		uploadFinerMADynamicColors(vts_for_subFineFaces, vals_for_subFineFaces_vts, cur_prune_min, cur_prune_max);
	}
	else
	{
		uploadFinerMADynamicColors(vts_for_subFineFaces, m_constColor_MA);
	}
	
	this->resizeGL(width(), height());
	updateGL();
}

void GLArea::updateFinePruneMA(float _iso)
{
	// st edge -> the iso vert idx on it
	map<TriEdge, unsigned> stE_isoV_map;
	// order vts on a face (in ccw or cw)
	vector<unsigned> vts_2_traverse;
	// iso verts and edges
	vector<TriPoint> iso_vts;
	vector<TriEdge> iso_edges;
	// fine faces resulted from tessellation
	vector<TriFace> fine_faces;

	// for each face, sort original & steiner vts on it in a certain order, 
	// identify iso vert for each st edge following that order, 
	// & connect iso verts for the face s.t. vts with "-" sign are connected with iso edges
	// and vts with "+" sign are isolated by the edges
	const auto& faces = stg->m_origG->faces;
	TriEdge e;
	vector<unsigned> stVts_e;
	for (size_t fi = 0; fi < faces.size(); ++fi)
	{
		vts_2_traverse.clear();
		const auto& f = faces[fi];

		// order vts on f consistently
		for (size_t i = 0; i < 3; ++i)
		{
			e = util::makeEdge(f[i], f[(i+1)%3]);
			stVts_e.clear();
			stg->m_stSubdiv.getStVertIndicesOnTriEdge(e, stVts_e);
			if ( e[0] != f[i] )
			{// need to reverse st-vts-on-edge
				vts_2_traverse.insert(vts_2_traverse.end(), stVts_e.rbegin(), stVts_e.rend());
			}
			else
			{
				vts_2_traverse.insert(vts_2_traverse.end(), stVts_e.begin(), stVts_e.end());
			}
		}

		//  do a quick filtering to see if this face is worth looking into
		int crossed_number = 0;
		for (auto it = vts_2_traverse.begin(); it != vts_2_traverse.end(); ++it)
		{
			
		}
	}
}

void GLArea::visMP()
{
	if (!stg)
		return;
	// 0. preprocess: find the biggest conn. component
	vector<vector<TriEdge>> cc_elist;
	vector<vector<unsigned>> cc_vlist;
	stg->find_conn_cmpnts(stg->dual_vts.size(), stg->dual_edges, cc_elist, cc_vlist);
	int max_cc_i = 0;
	int max_cc_vtsCnt = cc_vlist[0].size();
	for (int i = 0; i < cc_vlist.size(); ++i)
	{
		if (cc_vlist[i].size() > max_cc_vtsCnt)
		{
			max_cc_vtsCnt = cc_vlist[i].size();
			max_cc_i = i;
		}
	}
	const auto& max_cc = cc_vlist[max_cc_i];

	// 1. find the MP, i.e. the one with next burn vertex -1 in the maximum conn. cmpnt.
	auto next_burn = stg->burnNext_medialCurve;
	unsigned i;
	for (i = 0; i < max_cc.size(); ++i)
	{
		if (next_burn[max_cc[i]] == -1)
		{
			const auto& mp = stg->dual_vts[max_cc[i]];
			float r = 1.0f;
			m_meshMA->need_bbox();
			r = m_meshMA->bbox.radius() * 0.03f;

			auto mp_drawer = std::dynamic_pointer_cast<SphereDrawer>(m_MPDrawer);
			//mp_drawer->setCenter(mp);
			vector<TriPoint> pts;
			pts.push_back(mp);
			mp_drawer->setCenter(pts, r);

			cout << "MP drawn: center "<< mp << "radius "<<r<<endl;

			break;
		}
	}

	// MP must be found!
	assert(i < max_cc.size()); 
	
	updateGL();
}

void GLArea::outputMCwMeasure( const std::string & _meas_name ) const
{
	auto meas_name = _meas_name;
	std::unordered_set<string> valid_meas( { "bt2", "bt3" } );
	if ( valid_meas.count( meas_name ) == 0 )
		return; // meas_name = "BT3";
	std::string mc_file_name = 
		m_medialAxisFile.substr( 0, m_medialAxisFile.find_last_of( '.' ) ) + ".mc";
	std::string meas_file_name =
		m_medialAxisFile.substr( 0, m_medialAxisFile.find_last_of( '.' ) ) + "." + meas_name + ".msure";
	std::string mc_v_order_filename =
		m_medialAxisFile.substr( 0, m_medialAxisFile.find_last_of( '.' ) ) + ".mcorder";

	// output mc
	//stg->outputMC( mc_file_name );
	auto trans_cpy = m_trans_mat;
	trimesh::invert( trans_cpy );
	const auto& mc_vts = stg->getDualVts();
	const auto& mc_edges = stg->dual_edges;
	std::ofstream mc_os( mc_file_name );
	mc_os << mc_vts.size() << std::endl;
	for ( auto v : mc_vts )
	{
		v = trans_cpy * v;
		mc_os << v[ 0 ] << " " << v[ 1 ] << " " << v[ 2 ] << std::endl;
	}
	mc_os.close();
	cout << "Done: MC geometry output to file " << mc_file_name << endl;

	// output measure
	auto meas_type = DistMC::BT2_MC;
	vector<float> meas_per_mc_vert;
	if ( meas_name == "bt2" )
		meas_type = DistMC::BT2_MC;
	else if ( meas_name == "bt3" )
		meas_type = DistMC::BT3_MC;
	else if ( meas_name == "bt1" )
		meas_type = DistMC::BT1_MC;
	getMCDistMetric( meas_type, meas_per_mc_vert );
	std::ofstream meas_os( meas_file_name );
	meas_os << meas_per_mc_vert.size() << std::endl;
	for ( auto s : meas_per_mc_vert )
	{
		meas_os << s << std::endl;
	}
	meas_os.close();
	cout << "Done: MC measure output to file " << meas_file_name << endl;

	// output mc vert order defined by burning
	const auto& mc_v_order = stg->burnNext_medialCurve;
	std::ofstream mc_v_order_os( mc_v_order_filename );
	//mc_v_order_os << "# each pair of two verts a, b defines an order (from a to b)" << std::endl;
	mc_v_order_os << mc_v_order.size() << std::endl;
	for ( auto i = 0; i < mc_v_order.size(); ++i )
	{
		mc_v_order_os << i << " " << ( mc_v_order[ i ] < 0 ? i : mc_v_order[ i ] ) << std::endl;
	}
	mc_v_order_os.close();
	cout << "Done: MC verts order output to file " << mc_v_order_filename << endl;
}

void GLArea::usePerFaceRender(bool _is_perFace)
{
	//_is_perFace ? 
	//	cout << "use per face render" << endl : 
	//	cout << "use per vert render" << endl; // debug
	std::dynamic_pointer_cast<MeshDrawer>(m_MADrawer)->setRenderMode(
		_is_perFace? MeshDrawer::PER_FACE : MeshDrawer::PER_VERT);
}

void GLArea::setHSFaceDegenerateThreshold(float _t)
{
	this->m_hs->setPolyFaceDegenerateThresh(_t);
}

void GLArea::setComponentFaceNumberThresh(float _t)
{
	this->m_hs->setComponentFaceNumberThresh((int)_t);
}

void GLArea::getHSTheta2Abs(float _ratio, float& _abs_t)
{
	/*m_meshOrig->need_bbox();
	_abs_t = m_meshOrig->bbox.size().max() * _ratio;*/
	m_hs->ratioToAbs(_ratio, &_abs_t);
}
void GLArea::getHSTheta1Abs(float _ratio, float& _abs_t)
{
	/*m_meshOrig->need_bbox();
	_abs_t = m_meshOrig->bbox.size().max() * _ratio;*/
	m_hs->ratioToAbs(0.0f, nullptr, 0.0f, nullptr, _ratio, &_abs_t);
}
void GLArea::getHSTheta2Ratio(float _abs_t, float& _ratio)
{
	/*m_meshOrig->need_bbox();
	_ratio = _abs_t / m_meshOrig->bbox.size().max();*/
	m_hs->absToRatio(_abs_t, &_ratio);
}	 
void GLArea::getHSTheta1Ratio(float _abs_t, float& _ratio)
{
	/*m_meshOrig->need_bbox();
	_ratio = _abs_t / m_meshOrig->bbox.size().max();*/
	m_hs->absToRatio(0.0f, nullptr, 0.0f, nullptr, _abs_t, &_ratio);
}

void GLArea::createHS()
{
	if (!stg || stg->dual_vts.empty() || stg->dual_edges.empty())
	{
		return;
	}

#ifdef PROFILE_SPEED
	auto t1 = clock();
#endif // PROFILE_SPEED

	auto& hs = m_hs;
	cout << "Preprocessing data for hybrid skeleton..." << endl;
	hs->preprocess();
	cout << "Done: hybrid skeleton preprocess." << endl;
	/*hs->prepareFaceExtraction(stg->dual_vts, stg->dual_edges);
	cout << "Done: preparing for face extraction" << endl;
	hs->extractPolyFaces();
	cout << "Done: extracting poly faces" << endl;*/

	// now assign values to HS elements 
	vector<vector<float>> orig_persheet_bt2bt3diff;
	vector<vector<float>> orig_persheet_bt2bt3reldiff;
	vector<float> orig_vts_bt2bt3diff;
	vector<float> orig_vts_bt2bt3reldiff;
	vector<float> dual_vts_bt2bt3diff;
	vector<float> dual_vts_bt2bt3reldiff;
	vector<float> dual_vts_bt1bt2diff;
	vector<float> dual_vts_bt1bt2reldiff;
	getPerSheetDistMetricMA(DIFF_v, orig_persheet_bt2bt3diff);
	//cout << "orig_persheet_bt2bt3diff obtained" << endl;
	getPerSheetDistMetricMA(DIFF_REL_v, orig_persheet_bt2bt3reldiff);
	//cout << "orig_persheet_bt2bt3reldiff obtained" << endl;
	getVertDistMetricMA(DIFF_v, orig_vts_bt2bt3diff);
	//cout << "orig_vts_bt2bt3diff obtained" << endl;
	getVertDistMetricMA(DIFF_REL_v, orig_vts_bt2bt3reldiff);
	//cout << "orig_vts_bt2bt3reldiff obtained" << endl;
	getMCDistMetric(BT1_BT2_MC, dual_vts_bt1bt2diff);
	//cout << "dual_vts_bt1bt2diff obtained" << endl;
	getMCDistMetric(BT1_BT2_REL_MC, dual_vts_bt1bt2reldiff);
	//cout << "dual_vts_bt1bt2reldiff obtained" << endl;
	getMCDistMetric(BT2_BT3_MC, dual_vts_bt2bt3diff);
	getMCDistMetric(BT2_BT3_REL_MC, dual_vts_bt2bt3reldiff);

	hs->assignElementValues(
		orig_persheet_bt2bt3diff, orig_persheet_bt2bt3reldiff,
		orig_vts_bt2bt3diff, orig_vts_bt2bt3reldiff,
		dual_vts_bt2bt3diff, dual_vts_bt2bt3reldiff,
		dual_vts_bt1bt2diff, dual_vts_bt1bt2reldiff);
	cout << "Done: assigning element values" << endl;

#ifdef PROFILE_SPEED
	auto t2 = clock();
	auto t_ms = ( t2 - t1 ) * 1000.0 / CLOCKS_PER_SEC;
	cout << "Skeleton preprocessing took: " << t_ms << " ms." << endl;
#endif // PROFILE_SPEED
}

void GLArea::uploadHS()
{
	// obtain edges and faces of HS
	auto& vts_hs = m_hs->getVts();
	vector<TriEdge> edges_hs;
	m_hs->getRemainedEdges(edges_hs);
	//m_hs->getDualEdges( edges_hs );
	vector<TriFace> tri_faces_hs;
	m_hs->getRemainedFaces(tri_faces_hs);
	/*vector<TriColor> edge_colors;
	m_hs->getRemainedEdgesColor( edge_colors );*/

	//this->m_drawHS = true;
	cout << "hs faces left: " << tri_faces_hs.size() << endl;
	vector<TriColor> edge_color( 1, TriColor( 0, 0, 0 ) );
	vector<TriColor> face_color( 1, m_constColor_MA );
	uploadSimplicialComplex(vts_hs, edges_hs, tri_faces_hs, &edge_color/*nullptr*/, &face_color/*nullptr*/, m_hsLineDrawer, m_hsFaceDrawer);
}

//void GLArea::uploadHS()
//{
//	// obtain edges and faces of HS
//	auto& vts_hs = m_hs->getVts();
//	vector<TriEdge> edges_hs;
//	m_hs->getRemainedEdges(edges_hs);
//	vector<TriFace> tri_faces_hs;
//	m_hs->getRemainedTriFaces(tri_faces_hs);
//
//	auto hs_face_drawer_ptr = std::dynamic_pointer_cast<MeshDrawer>(this->m_hsFaceDrawer);
//	hs_face_drawer_ptr->setPointsPerFace(vts_hs, tri_faces_hs);
//	cout << "setPointsPerFace done!" << endl;
//
//	// TODO: color the HS face and line according to the cur dist metrics
//	// now just use default uniform color
//	
//	const TriColor face_color(1.0f, 98.0f/255.0f, 0.0f);
//	const TriColor edge_color(0.0f, 0.0f, 0.0f);
//	// setup HS faces
//	float* color_data = new float[tri_faces_hs.size()*3];
//	for (unsigned i = 0; i < tri_faces_hs.size()*3; i += 3)
//	{
//		color_data[i+0] = face_color[0];
//		color_data[i+1] = face_color[1];
//		color_data[i+2] = face_color[2];
//	}
//	hs_face_drawer_ptr->setPerFaceColor(color_data, tri_faces_hs.size());
//	delete [] color_data;
//	color_data = nullptr;
//	
//	float* normal_data = new float[tri_faces_hs.size()*3];
//	for (unsigned i = 0; i < tri_faces_hs.size(); ++i)
//	{
//		auto& f = tri_faces_hs[i];
//		auto& nml = (vts_hs[f[0]] - vts_hs[f[1]]).cross(vts_hs[f[0]] - vts_hs[f[2]]);
//		trimesh::normalize(nml);
//		normal_data[i*3 + 0] = nml[0];
//		normal_data[i*3 + 1] = nml[1];
//		normal_data[i*3 + 2] = nml[2];
//	}
//	hs_face_drawer_ptr->setPerFaceNormal(normal_data, tri_faces_hs.size());
//	delete [] normal_data;
//	normal_data = nullptr;
//	
//	float* saliency = new float[tri_faces_hs.size()];
//	for (unsigned i = 0; i < tri_faces_hs.size(); ++i)
//	{
//		saliency[i] = 0.2f;
//	}
//	hs_face_drawer_ptr->setPerFaceSaliency(saliency, tri_faces_hs.size());
//	delete [] saliency;
//	saliency = nullptr;
//
//	this->m_drawHS = true;
//	hs_face_drawer_ptr->setRenderMode(MeshDrawer::PER_FACE);
//
//	// setup HS lines
//	auto hs_line_drawer_ptr = std::dynamic_pointer_cast<LineDrawer>(this->m_hsLineDrawer);
//	/*hs_line_drawer_ptr->setPoints(vts_hs);
//	hs_line_drawer_ptr->setLines(edges_hs);*/
//	uploadLinesToDrawer(vts_hs, edges_hs, hs_line_drawer_ptr);
//	cout << "setPoints & setLines done!"<<endl;
//
//	color_data = new float[edges_hs.size() * 2 * 3];
//	for (unsigned i = 0; i < edges_hs.size() * 2 * 3; i += 3)
//	{
//		color_data[i + 0] = edge_color[0];
//		color_data[i + 1] = edge_color[1];
//		color_data[i + 2] = edge_color[2];
//	}
//	//stg->m_hs.getEdgeDiffValColor(color_data);
//	cout << "GLArea::setting line color ..."<<endl;
//	hs_line_drawer_ptr->setPerVertColor(color_data, edges_hs.size() * 2);
//	delete [] color_data;
//	color_data = nullptr;
//	cout << "LineDrawer::setPerVertColor done!"<<endl;
//
//	updateGL();
//}


void GLArea::pruneHS(
	float _f_diff_r, float _f_reldiff_r, float _l_diff_r, float _l_reldiff_r, 
	bool _use_inputs_directly,
	bool _remove_small_cmpnts)
{
#ifdef PROFILE_SPEED
	auto t1 = clock();
#endif // PROFILE_SPEED

	m_hs->prune(_f_diff_r, _f_reldiff_r, _l_diff_r, _l_reldiff_r, 
		_use_inputs_directly, _remove_small_cmpnts
		);

#ifdef PROFILE_SPEED
	auto t2 = clock();
	auto t_ms = ( t2 - t1 ) * 1000.0 / CLOCKS_PER_SEC;
	cout << "Skeleton creation took: " << t_ms << " ms." << endl;
#endif // PROFILE_SPEED
}

void GLArea::computeSurfFuncCorrespondence(SurfFuncCorrespScheme _scheme, float _r_eps_ratio, int _k)
{
	if (_scheme == SURF_FUNC_RANGE_CORRESP)
	{
		this->m_surfF->computeCorrespondenceWithRangeSearch(_r_eps_ratio, _k);
	}
	else if (_scheme == SURF_FUNC_KNN_CORRESP)
	{
		this->m_surfF->computeCorrespondenceWithKNNSearch(_k);
	}

	/*cout << "setting up diffusion system..."<<endl;
	m_surfF->setupDiffusionSystem();
	cout << "diffusion system setup."<<endl;*/
}

void GLArea::projectFunctionOnSurf(int _idx, float _smooth_ratio, bool _do_diffuse/* = false*/)
{
	SurfaceFunc::SurfFuncType surf_func_type = SurfaceFunc::SHAPE_DIAM;
	switch(_idx)
	{
	case 0: surf_func_type = SurfaceFunc::SHAPE_DIAM; break;
	case 1: surf_func_type = SurfaceFunc::SHAPE_WIDTH; break;
	case 2: surf_func_type = SurfaceFunc::SHAPE_WIDTH_EXTREMITY; break;
	case 3: surf_func_type = SurfaceFunc::SHAPE_LENGTH_EXTREMITY; break;
	default: surf_func_type = SurfaceFunc::SHAPE_WIDTH; break;
	}

	vector<float> scalars_per_vert;
	m_surfF->getSurfaceFunction(surf_func_type, scalars_per_vert, _smooth_ratio, _do_diffuse);
	
	float max_val = numeric_limits<float>::min();
	float min_val = numeric_limits<float>::max();
	for (unsigned i = 0; i < scalars_per_vert.size(); ++i)
	{
		float val = scalars_per_vert[i];
		if (val == SurfaceFunc::INVALID_VAL || val < 0.0f)
		{
			continue;
		}
		min_val = std::min(min_val, scalars_per_vert[i]);
		max_val = std::max(max_val, scalars_per_vert[i]);
	}
	cout << "project value on surf. range of values: "<<min_val<<", "<<max_val<<endl;

	// color each vert of the orig 3d mesh by the scalars obtained above
	float* color = new float[scalars_per_vert.size() * 3];
	TriColor c;
	for (unsigned i = 0; i < scalars_per_vert.size(); ++i)
	{
		float val = scalars_per_vert[i];
		if (val < 0.0f || val == SurfaceFunc::INVALID_VAL)
		{
			c = TriColor(0.0f, 0.0f, 0.0f);
		}
		else
		{
			c = util::GetColour(val, min_val, max_val);
		}
		color[3*i + 0] = c[0];
		color[3*i + 1] = c[1];
		color[3*i + 2] = c[2];
	}
	std::dynamic_pointer_cast<MeshDrawer>(this->m_origDrawer)->setPerVertColor(color, scalars_per_vert.size());
	delete [] color;

	cout << "Surface function projection done."<<endl;

	updateGL();
}

void GLArea::setupIsoSurface()
{
	//m_iso_surf = shared_ptr<IsoSurfFrom2Manifold>(new IsoSurfFrom2Manifold());
	m_iso_surf->setup(stg, m_meshOrig, m_isoSurfDrawer);
	m_iso_surf->precompute();

	cout << "iso-surface setup." << endl;
}

void GLArea::refineIsoSurface(bool _is_recursive)
{
	m_iso_surf->refineTriangulation(_is_recursive); // only one iteration of refinement.
}

void GLArea::uploadIsoSurface(float _t, bool _hide_snapped_faces)
{
	m_iso_surf->updateIsoSurf(_t, _hide_snapped_faces);
}

void GLArea::uploadIsoSurface(float _iso, bool _hide_snapped_faces, int)
{
	m_iso_surf->updateIsoSurf(_iso, _hide_snapped_faces, 1);
}

void GLArea::setupIsoContour(float min_alpha, float exp_alpha, const float* const _new_min, const float* const _new_max)
{
	//this->precomputeForMeasures();
	bool is_transparent = min_alpha < 0.98f;
	m_iso_cont->setup(stg, is_transparent, min_alpha, exp_alpha, _new_min, _new_max);
	//m_iso_cont->precompute();

	/*
	// set current dist metric to be BT2. 
	// dynamic part of MC SWEPT BY ISO CONTOUR will use such info.
	float min_val, max_val;
	this->setUpMCDistMetricForPruning(GLArea::BT2_MC, min_val, max_val);
	*/

	cout << "iso-contour object setup." << endl;
}

void GLArea::uploadIsoContour(float _t, bool _hide_past_faces, bool _show_exposed_mc)
{
	float iso = m_iso_cont->getIsoValue(_t);
	m_iso_cont->genContour(iso, 1);
	if (_show_exposed_mc)
	{
		uploadFromRemainedMC(true, false, iso, false);
	}
}

void GLArea::uploadIsoContour(float _iso, bool _hide_past_faces, bool _show_exposed_mc, int)
{
	m_iso_cont->genContour(_iso, 1);
	if (_show_exposed_mc)
	{
		uploadFromRemainedMC(true, false, _iso, false);
	}
}

void GLArea::uploadFromRemainedMC(bool _only_show_exposed, bool _prune, float _iso, bool _fronts_on_MC)
{
	if (m_MC_vts_neighbors.empty() || this->m_iso_cont == nullptr)
		return;

	vector<unsigned> MC_edges_to_upload;
	vector<int> MC_nbsCnt_cpy(m_remained_MC_nbsCnt);
	vector<TriPoint> open_ends; 
	const auto& fine_tri_for_dual_edge = this->stg->m_fromFineTri_for_dualE;
	const auto& fine_tri_status = this->m_iso_cont->m_face_status_static;
	bool pruned, unexposed;

	for (unsigned i = 0; i < m_remained_MC.size(); ++i)
	{
		unsigned ei = m_remained_MC[i];
		const auto& e = stg->dual_edges[ei];
		/*if (pruned_edges.count(util::makeEdge(e[0], e[1])) == 0)
			_remain_edges.push_back(e);*/
		pruned = _prune && (
			cur_medialCurve_distMetric[e[0]] < _iso ||
			cur_medialCurve_distMetric[e[1]] < _iso 
			);
		unexposed = _only_show_exposed && (
			fine_tri_status[fine_tri_for_dual_edge[ei]] != IsoContourOnMA::BEHIND 
			);
		if ( pruned || unexposed )
		{
			if (pruned)
			{
				MC_nbsCnt_cpy[e[0]] --;
				MC_nbsCnt_cpy[e[1]] --;
				assert(MC_nbsCnt_cpy[e[1]] >= 0);
				assert(MC_nbsCnt_cpy[e[0]] >= 0);
			}
		}
		else
		{
			MC_edges_to_upload.push_back(ei);
		}
	}

	// collect open ends: 
	// only those open ends of the edges that are to upload
	for (unsigned i = 0; i < MC_edges_to_upload.size(); ++i)
	{
		unsigned ei = MC_edges_to_upload[i];
		const auto& e = stg->dual_edges[ei];

		if (MC_nbsCnt_cpy[e[0]] == 1)
			open_ends.push_back(stg->dual_vts[e[0]]);
		if (MC_nbsCnt_cpy[e[1]] == 1)
			open_ends.push_back(stg->dual_vts[e[1]]);
	}

	// upload remaining MC subset
	uploadLineSubset(MC_edges_to_upload, stg->dual_edges);
	if(_fronts_on_MC)
	{
		// upload burning fronts (as red points)
		std::dynamic_pointer_cast<PointDrawer>(m_pointDrawer)->setPoints(open_ends);
		float* color_data = new float[open_ends.size()*3];
		for (unsigned i = 0; i < open_ends.size(); ++i)
		{
			color_data[3*i+0] = 1.0f;
			color_data[3*i+1] = 0.0f;
			color_data[3*i+2] = 0.0f;
		}
		std::dynamic_pointer_cast<PointDrawer>(m_pointDrawer)->setPerVertColor(color_data, open_ends.size());
		delete [] color_data;
	}

	// upload the dynamic part of MC. 
	// use the curr iso as the prune threshold which will be used if _prune == true
	uploadDynamicMC(_iso, _prune, _iso);
}

void GLArea::uploadDynamicMC(float _iso_bt2, bool _prune_bt1, float _pruneThresh)
{
	vector<TriPoint> dynamic_vts;
	vector<TriEdge> dynamic_edges;
	float u_val, v_val;
	for (unsigned i = 0; i < m_remained_MC.size(); ++i)
	{
		unsigned ei = m_remained_MC[i];
		auto e = stg->dual_edges[ei];
		
		if (_prune_bt1 && 
			(stg->bt1_medialCurve[e[0]] < _pruneThresh || 
			stg->bt1_medialCurve[e[1]] < _pruneThresh) 
			)
		{
			continue;
		}

		// make sure u_val <= v_val (values at two ends of cur edge e)
		u_val = stg->bt2_MC[e[0]]; /*cur_medialCurve_distMetric[e[0]];*/
		v_val = stg->bt2_MC[e[1]]; /*cur_medialCurve_distMetric[e[1]];*/
		if (u_val > v_val)
		{
			std::swap(u_val, v_val);
			std::swap(e[0], e[1]);
		}
		// store dynamic part of the line behind the iso value ( i.e. <e[0], interpolated vert> )
		if ( (u_val - _iso_bt2) * (v_val - _iso_bt2) < 0.0f )
		{
			dynamic_vts.push_back(stg->dual_vts[e[0]]);
			dynamic_vts.push_back(
				trimesh::mix( stg->dual_vts[e[0]], stg->dual_vts[e[1]], 
				std::abs(u_val - _iso_bt2) / (v_val - u_val) ) 
				);
			dynamic_edges.push_back(util::makeEdge(dynamic_vts.size()-2, dynamic_vts.size()-1));
		}
	}

	auto dynamic_line_drawer = std::dynamic_pointer_cast<LineDrawer>(m_dynamic_dualLineDrawer);
	dynamic_line_drawer->setPoints(dynamic_vts);
	dynamic_line_drawer->setLines(dynamic_edges);
}

float GLArea::getRangeOfDiffDist()
{
	if (!stg || stg->min_diffDistOrigFaces == numeric_limits<float>::max())
		return 1.0f;
	return stg->max_diffDistOrigFaces - stg->min_diffDistOrigFaces;
}

float GLArea::getMinDiffDist()
{
	return stg->min_diffDistOrigFaces;
}

float GLArea::getMaxDiffDist()
{
	return stg->max_diffDistOrigFaces;
}

const vector<TriPoint>& GLArea::getOriginalVts() const
{
	return stg->m_origG->vts;
}

int GLArea::getNumFacesOfMA() const
{
	if (m_meshMA)
		return m_meshMA->faces.size();
	return -1;
}

bool GLArea::changeOrigColor(float * _color)
{
	if ( !m_meshOrig )
		return false;

	auto& orig_drawer = std::dynamic_pointer_cast<MeshDrawer>( m_origDrawer );
	//float* color_data = new float[m_meshOrig->vertices.size() * 3];
	auto color_data = new float[ m_meshOrig->faces.size() * 3 ];
	for ( unsigned i = 0; i < m_meshOrig->faces.size(); ++i )
	{
		color_data[ 3 * i + 0 ] = _color[ 0 ];
		color_data[ 3 * i + 1 ] = _color[ 1 ];
		color_data[ 3 * i + 2 ] = _color[ 2 ];
	}
	orig_drawer->setRenderMode( MeshDrawer::PER_FACE );
	orig_drawer->setPerFaceColor( color_data, m_meshOrig->faces.size() );
	/*for (unsigned i = 0; i < m_meshOrig->vertices.size(); ++i)
	{
		color_data[3*i+0] = _color[0];
		color_data[3*i+1] = _color[1];
		color_data[3*i+2] = _color[2];
	}
	m_origDrawer->setPerVertColor(color_data, m_meshOrig->vertices.size());*/
	delete [] color_data;

	updateGL();
	return true;
}

void GLArea::setUseConstColorMA( bool _use_constColor_for_MA, const TriColor& _c )
{
	m_use_constColor_for_MA = _use_constColor_for_MA;
	m_constColor_MA = _c;
}

bool GLArea::changeConstColorMA(float * _color, bool _color_vert)
{
	if ( !m_meshMA )
		return false;

	auto& MA_drawer = std::dynamic_pointer_cast<MeshDrawer>(m_MADrawer);
	float* color_data;
	if (_color_vert)
	{
		color_data = new float[m_meshMA->vertices.size() * 3];
		for (unsigned i = 0; i < m_meshMA->vertices.size(); ++i)
		{
			color_data[3*i+0] = _color[0];
			color_data[3*i+1] = _color[1];
			color_data[3*i+2] = _color[2];
		}
		MA_drawer->setRenderMode(MeshDrawer::PER_VERT);
		MA_drawer->setPerVertColor(color_data, m_meshMA->vertices.size());
	}
	else
	{
		color_data = new float[m_meshMA->faces.size() * 3];
		for (unsigned i = 0; i < m_meshMA->faces.size(); ++i)
		{
			color_data[3*i+0] = _color[0];
			color_data[3*i+1] = _color[1];
			color_data[3*i+2] = _color[2];
		}
		MA_drawer->setRenderMode(MeshDrawer::PER_FACE);
		MA_drawer->setPerFaceColor(color_data, m_meshMA->faces.size());
	}
	delete [] color_data;

	// change color of the 1-cell in input MA
	color_data = new float[ m_meshMA->lines.size() * 2 * 3 ];
	for ( auto ei = 0; ei < m_meshMA->lines.size(); ++ei )
	{
		const auto& e = m_meshMA->lines[ ei ];
		color_data[ ei * 2 * 3 + 0 + 0 ] = _color[ 0 ];
		color_data[ ei * 2 * 3 + 0 + 1 ] = _color[ 1 ];
		color_data[ ei * 2 * 3 + 0 + 2 ] = _color[ 2 ];

		color_data[ ei * 2 * 3 + 3 + 0 ] = _color[ 0 ];
		color_data[ ei * 2 * 3 + 3 + 1 ] = _color[ 1 ];
		color_data[ ei * 2 * 3 + 3 + 2 ] = _color[ 2 ];
	}
	std::dynamic_pointer_cast<LineDrawer>( m_MALineDrawer )->setPerVertColor( color_data, m_meshMA->lines.size() * 2 );
	delete[] color_data;

	// also update color of finer MA 
	uploadMAFinerStaticColors( TriColor(_color) );

	updateGL();
	return true;
}

bool GLArea::changeConstColorMC(float * _color)
{
	if ( !stg || stg->dual_vts.empty() )
		return false;
	float *color_data = new float[ stg->dual_edges.size() * 2 * 3 ];
	for (unsigned ei = 0; ei < stg->dual_edges.size(); ei++)
	{
		const auto& e = stg->dual_edges[ei];
		color_data[ei*2*3 + 0 + 0] = _color[0];
		color_data[ei*2*3 + 0 + 1] = _color[1];
		color_data[ei*2*3 + 0 + 2] = _color[2];

		color_data[ei*2*3 + 3 + 0] = _color[0];
		color_data[ei*2*3 + 3 + 1] = _color[1];
		color_data[ei*2*3 + 3 + 2] = _color[2];
	}

	// change color of MC
	std::dynamic_pointer_cast<LineDrawer>(m_dualLinesDrawer)->setPerVertColor(color_data, stg->dual_edges.size() * 2);
	// also change color of the dynamic part of MC
	std::dynamic_pointer_cast<LineDrawer>(m_dynamic_dualLineDrawer)->setPerVertColor(color_data, stg->dual_edges.size()*2 );
	delete[] color_data; 
	color_data = nullptr;

	return true;
}

bool GLArea::changeBGColor(float* _color)
{
	bg_color[0] = _color[0];
	bg_color[1] = _color[1];
	bg_color[2] = _color[2];
	gl.ClearColor(_color[0], _color[1], _color[2], 1.0f);

	updateGL();
	return true;
}

bool GLArea::changeMAVertTransparency(VertFieldType _type, float _min_alpha, int _exp)
{
	if (!m_MADrawer || ui->visMADistCombo->currentIndex() < 0)
		return false;

	if (_min_alpha > 0.9f)
	{
		m_MATransparent = false;
	}
	else
	{
		m_MATransparent = true;
		float dummy;
		colorMAVertBy(_type, dummy, dummy, false, _min_alpha, _exp, 2);
	}
	m_MADrawer->setTransparencyEnabled(m_MATransparent);
	
	updateGL();
	return true;
}

bool GLArea::changeMAFaceTransparency(FaceFieldType _type, float _min_alpha, int _exp)
{
	if (!m_MADrawer || ui->visMADistCombo->currentIndex() < 0)
		return false;

	if (_min_alpha > 0.9f)
	{
		m_MATransparent = false;
	}
	else
	{
		m_MATransparent = true;
		float dummy;
		colorMAFaceBy( 
			_type, 
			dummy, dummy, false, 
			_min_alpha, _exp, 2, 
			ui->doMAFaceScalarDiffusion->isChecked(),
			ui->usePerSheetBox->isChecked() );
	}
	m_MADrawer->setTransparencyEnabled(m_MATransparent);

	updateGL();
	return true;
}

bool GLArea::changeOrigTransparency(float _alpha)
{
	if (!m_origDrawer)
		return false;

	auto drawer = std::dynamic_pointer_cast<MeshDrawer>(m_origDrawer);
	if (_alpha > 0.9f)
	{
		m_OrigTransparent = false;
	}
	else
	{
		m_OrigTransparent = true;
		/*drawer->setPerVertConstantSaliency(_alpha, this->m_meshOrig->vertices.size());*/
		drawer->setPerVertConstantSaliency( _alpha, this->m_meshOrig->faces.size()*3 );
	}
	drawer->setTransparencyEnabled(m_OrigTransparent);

	updateGL();
	return true;
}

bool GLArea::changeMCTransparency(DistMC _type, float _min_alpha, int _exp)
{
	if (!m_dualLinesDrawer)
		return false;

	if (_min_alpha > 0.9f)
	{
		m_lineTransparent = false;
	}
	else
	{
		m_lineTransparent = true;
		float dummy = 0.0f;
		this->colorMCEdgeBy(_type, dummy, dummy, false, _min_alpha, _exp, 2);
	}
	m_dualLinesDrawer->setTransparencyEnabled(m_lineTransparent);
	m_dynamic_dualLineDrawer->setTransparencyEnabled(m_lineTransparent);
	
	updateGL();
	return true;
}

bool GLArea::changeMCLineWidth(float _width)
{
	if (!m_dualLinesDrawer)
		return false;

	gl.LineWidth(std::max(_width, 1.0f));

	updateGL();
	return true;
}

void GLArea::turnOffLighting(enum ObjectType _obj_t, bool _turn_off) 
{
	if (_obj_t == FINE_MA_FOR_ISO_CONTOUR)
	{
		if (this->m_iso_cont != nullptr)
		{
			this->m_iso_cont->setLightingEnabled(!_turn_off);
		}
	}
	else if (_obj_t == FINE_MA_FOR_PRUNING)
	{
		if (this->m_FinerMAStaticDrawer != nullptr && 
			this->m_MAFinnerDynamicDrawer != nullptr)
		{
			this->m_FinerMAStaticDrawer->setLightingEnabled(!_turn_off);
			this->m_MAFinnerDynamicDrawer->setLightingEnabled(!_turn_off);
		}
	}
	else if (_obj_t == MA)
	{
		if (this->m_MADrawer)
			m_MADrawer->setLightingEnabled(!_turn_off);
	}
	
	updateGL();
};

void GLArea::setTransparencyEnabled(enum ObjectType _obj_t, bool _enabled)
{
	if (_obj_t == FINE_MA_FOR_ISO_CONTOUR)
	{
		this->m_iso_cont->setTransparencyEnabled(_enabled);
	}
	else
	{
		m_drawers[_obj_t]->setTransparencyEnabled(_enabled);
	}

	updateGL();
};

void GLArea::setTrueTransparencyRender(bool _is_enabled)
{
	this->use_DP_transparency = _is_enabled;
}

void GLArea::setDPMaxRenders(int _max_iter)
{
	this->m_DPMaxRenders = _max_iter;
}

void GLArea::setLinesOffset(float _r)
{
	std::dynamic_pointer_cast<LineDrawer>(m_dualLinesDrawer)->setZFightOffset(_r);
	std::dynamic_pointer_cast<LineDrawer>(m_hsLineDrawer)->setZFightOffset(_r);
	std::dynamic_pointer_cast<LineDrawer>( m_MALineDrawer )->setZFightOffset( _r );
	std::dynamic_pointer_cast<LineDrawer>(m_dynamic_dualLineDrawer)->setZFightOffset(_r);

	if (m_iso_cont)
		m_iso_cont->setZFightOffset(_r);
	if (m_dualLinesDrawer)
		std::dynamic_pointer_cast<LineDrawer>(m_dualLinesDrawer)->setZFightOffset(_r);

	updateGL();
	this->resizeGL(width(), height());
}

bool GLArea::outputToWenping()
{
	bool ret_code;
	auto temp = QString(this->m_medialAxisFile.c_str()).remove(".clean").replace(".off", ".ma");
	auto _file_name = temp.toStdString();
	if (radii.empty())
	{
		cout << "all radii are 0s" << endl;
		ret_code = stg->outputToWenping(_file_name.c_str(), vector<float>(m_meshMA->vertices.size(), 0.0f));
	}
	else
	{
		ret_code = stg->outputToWenping(_file_name.c_str(), radii);
	}
	return ret_code;
}

bool GLArea::readQMATFile(std::string _filename)
{
	cout << "reading medial structure from qmat file "<<_filename<<endl;
	vector<float> r;
	auto mesh = MyMesh::readQMATFile( _filename, r );
	if ( !mesh )
	{
		cout << "Error: cannot read qmat file " << _filename << endl;
		return false;
	}
	cout << "Done: reading qmat file."<<endl;

	// draw this simplicial complex
	cout << "uploading qmat data to GPU for rendering..."<<endl;
	this->m_drawQMAT = true;
	uploadSimplicialComplex(mesh->vertices, mesh->lines, mesh->faces, nullptr, nullptr, m_qmatEdgeDrawer, m_qmatFaceDrawer);
	cout << "uploading qmat data to GPU done!"<<endl;
	
	delete mesh;

	return true;
}

void GLArea::exportPerVertexET()
{
	QString bt_name = QString(m_medialAxisFile.c_str());
	bt_name.remove( ".clean" );
	bt_name = QFileInfo( bt_name ).baseName() + ".et";
	stg->exportPerVertexET(bt_name.toStdString());
}

//bool GLArea::importPerVertexBT(std::string _filename)
//{
//	if (stg->setPerVertexBT(_filename) == false)
//	{
//		cout << "failed to set per-vertex BT from file "<< _filename<<endl;
//		return false;
//	}
//
//	return true;
//}

void GLArea::exportPerSectorET()
{
	QString bt_name = QString(m_medialAxisFile.c_str());
	QString ext = ".off";
	bt_name.remove(".clean").replace(".off", ".etps");
	stg->exportPerSectorET(bt_name.toStdString());
}

void GLArea::exportSkeleton()
{
	QString skel_name = QString(m_medialAxisFile.c_str());
	skel_name.remove( ".clean" );
	auto skel_w_msure_name = skel_name;

	if ( skel_w_msure_name.indexOf( ".off" ) >= 0 )
		skel_w_msure_name.replace( ".off", "_skel.skMsure" );
	else if ( skel_w_msure_name.indexOf( ".ply" ) >= 0 )
		skel_w_msure_name.replace( ".ply", "_skel.skMsure" );
	if ( skel_name.indexOf( ".off" ) >= 0 )
		skel_name.replace( ".off", "_skel.ply" );
	else if ( skel_name.indexOf( ".ply" ) >= 0 )
		skel_name.replace( ".ply", "_skel.ply" );

	auto trans_cpy = m_trans_mat;
	trimesh::invert(trans_cpy);
	m_hs->exportSkeleton( { skel_name.toStdString(), skel_w_msure_name.toStdString() }, trans_cpy );
}

void GLArea::initializeGL(void)
{
	cout << "Initializing GL... " << endl; // debug

	// init gl wrangler, e.g. glew
	oglplus::GLAPIInitializer gl_init_obj;

	// set up gl states and enable various capabilities
	using oglplus::Capability;
	gl.ClearColor(1.0f, 1.0f, 1.0f, 1.0f);
	gl.ClearDepth(1.0f);
	//gl.Enable(Capability::DepthTest);
	//gl.Enable(Capability::Blend);
	//gl.BlendFunc(BlendFn::SrcAlpha, BlendFn::One);
	gl.BlendFunc(BlendFn::SrcAlpha, BlendFn::OneMinusSrcAlpha);
	//gl.BlendFunc(BlendFn::One, BlendFn::One);
	// line smooth and width
	gl.Enable(Capability::LineSmooth);
	//gl.LineWidth(1.5f);
	
	// specify polygon offset params
	//gl.PolygonOffset(1.0f, 1.0f);

	initShaders();
	m_MADrawer = std::shared_ptr<MeshDrawer>(
		new MeshDrawer(width(), height(), m_simpProg, m_edgeProg)
		);
	track_ball = dynamic_pointer_cast<MeshDrawer>(m_MADrawer)->getCamera();
	m_MALineDrawer = std::shared_ptr<LineDrawer>(
		new LineDrawer( m_linesProg, track_ball )
		);
	m_origDrawer = std::shared_ptr<MeshDrawer>(
		new MeshDrawer(width(), height(), m_simpProg, m_edgeProg, track_ball)
		);
	m_FinerMAStaticDrawer = std::shared_ptr<MeshDrawer>(
		new MeshDrawer( width(), height(), m_simpProg, m_edgeProg, track_ball )
		);
	m_MAFinnerDynamicDrawer = std::shared_ptr<MeshDrawer>(
		new MeshDrawer( width(), height(), m_simpProg, m_edgeProg, track_ball )
		);
	m_pointDrawer = std::shared_ptr<PointDrawer>(
		new PointDrawer(m_pointsProg, track_ball)
		);
	m_dualLinesDrawer = std::shared_ptr<LineDrawer>(
		new LineDrawer(m_linesProg, track_ball)
		);
	m_dynamic_dualLineDrawer = std::shared_ptr<LineDrawer>(
		new LineDrawer(m_linesProg, track_ball)
		);
	m_burntEdgesDrawer = std::shared_ptr<LineDrawer>(
		new LineDrawer(m_linesProg, track_ball)
		);
	m_linesProxyDrawer = std::shared_ptr<LineDrawer>(
		new LineDrawer(m_linesProxyProg, track_ball)
		);
	m_hsFaceDrawer = std::shared_ptr<MeshDrawer>(
		new MeshDrawer(width(), height(), m_simpProg, m_edgeProg, track_ball)
		);
	m_hsLineDrawer = std::shared_ptr<LineDrawer>(
		new LineDrawer(m_linesProg, track_ball)
		);
	m_MPDrawer = std::shared_ptr<SphereDrawer>(
		new SphereDrawer(track_ball)
		);
	m_isoSurfDrawer = std::shared_ptr<MeshDrawer>(
		new MeshDrawer(width(), height(), m_simpProg, m_edgeProg, track_ball)
		);
	m_iso_cont->setupDrawers(width(), height(), m_simpProg, m_edgeProg, m_linesProg, track_ball);

	m_drawers.clear();
	m_drawers.push_back(m_origDrawer);
	m_drawers.push_back(m_MADrawer);
	m_drawers.push_back( m_MALineDrawer );
	m_drawers.push_back(m_pointDrawer);
	m_drawers.push_back(m_dualLinesDrawer);
	m_drawers.push_back(m_dynamic_dualLineDrawer);
	m_drawers.push_back(m_linesProxyDrawer);
	m_drawers.push_back(m_burntEdgesDrawer);
	m_drawers.push_back(m_hsFaceDrawer);
	m_drawers.push_back(m_hsLineDrawer);
	m_drawers.push_back(m_MPDrawer);
	m_drawers.push_back(m_isoSurfDrawer);
	m_drawers.push_back(m_FinerMAStaticDrawer);
	m_drawers.push_back(m_MAFinnerDynamicDrawer);
	m_drawers.push_back(m_qmatEdgeDrawer);
	m_drawers.push_back(m_qmatFaceDrawer);

	// SVP related
	/*m_drawerSVP = std::shared_ptr<MeshDrawer>(
		new MeshDrawer(width(), height(), m_simpProg, m_edgeProg, m_svpProg, track_ball)
		);*/

	//initOIT();
	initDP();

	// init qmat's drawers
	m_qmatFaceDrawer = std::shared_ptr<MeshDrawer>(
		new MeshDrawer(width(), height(), m_simpProg, m_edgeProg, track_ball)
		);
	m_qmatEdgeDrawer = std::shared_ptr<LineDrawer>(
		new LineDrawer(m_linesProg, track_ball)
		);
    m_glInitialized = true;
	cout << "Done!" << endl; // debug
    
}

void GLArea::resizeGL(int _w, int _h)
{
	_w = max(_w, 0);
	_h = max(_h, 0);

	//assert(m_MADrawer);
	//cout << "GLArea::resizeGL()"<< endl;
	if (m_drawOrig)
		m_origDrawer->reshape(_w, _h);
	if (m_drawMA)
	{
		m_MADrawer->reshape( _w, _h );
		m_MALineDrawer->reshape( _w, _h );
	}
	if (m_drawMAFinnerStatic)
	{
		m_FinerMAStaticDrawer->reshape(_w, _h);
		m_MAFinnerDynamicDrawer->reshape(_w, _h);
	}
	if (m_drawPoints)
		m_pointDrawer->reshape(_w, _h);
	if (m_drawMC)
		m_dualLinesDrawer->reshape(_w, _h);
	if (m_drawBurntEdges)
		m_burntEdgesDrawer->reshape(_w, _h);
	if (m_drawLinesAsProxy)
	{	
		m_linesProxyDrawer->reshape(_w, _h);
	}
	if (m_drawHS)
	{
		m_hsLineDrawer->reshape(_w, _h);
		m_hsFaceDrawer->reshape(_w, _h);
	}
	if (m_drawMP)
	{
		m_MPDrawer->reshape(_w, _h);
	}
	if (m_drawIsoSurf)
	{
		m_isoSurfDrawer->reshape(_w, _h);
		m_dynamic_dualLineDrawer->reshape(_w, _h);
	}
	if (m_drawIsoCont)
	{
		m_iso_cont->reshape(_w, _h);
	}
	if (m_drawQMAT)
	{
	m_qmatEdgeDrawer->reshape(_w, _h);
	m_qmatFaceDrawer->reshape(_w, _h);
	}

	//resizeOIT(_w, _h);
	resizeDP(_w, _h);
	//cout << "GLArea::resizeGL() done." << endl;
}

void GLArea::paintGL(void)
{
	try
	{
		assert(m_MADrawer);
		//cout << "GLArea::paintGL()"<< endl;
		transparent_draws.clear();
		opaque_draws.clear();

		if (m_drawMA)
		{
			if (m_MATransparent)
			{
				/*gl.Enable(Capability::Blend);
				gl.Disable(Capability::DepthTest);*/
				transparent_draws.push_back( std::make_pair(m_MADrawer, true) );
			}
			else
			{
				opaque_draws.push_back( std::make_pair(m_MADrawer, true) );
				//m_MADrawer->render(0.0);
			}
			if ( m_drawMALines )
			{
				if ( m_lineTransparent )
				{
					transparent_draws.push_back( std::make_pair( m_MALineDrawer, false ) );
				}
				else
				{
					opaque_draws.push_back( std::make_pair( m_MALineDrawer, false ) );
				}
			}
			//m_MADrawer->render(0.0);
		}

		if (m_drawOrig)
		{
			if (m_OrigTransparent)
			{
				/*gl.Enable(Capability::Blend);
				gl.Disable(Capability::DepthTest);*/
				transparent_draws.push_back(std::make_pair(m_origDrawer, true));
			}
			else
			{
				/*gl.Enable(Capability::DepthTest);
				gl.Disable(Capability::Blend);*/
				opaque_draws.push_back(std::make_pair(m_origDrawer, true));
				//m_origDrawer->render(0.0);
			}
			//m_origDrawer->render(0.0);
		}

		if (m_drawMAFinnerStatic)
		{
			opaque_draws.push_back(std::make_pair(m_FinerMAStaticDrawer, true));
			opaque_draws.push_back(std::make_pair(m_MAFinnerDynamicDrawer, true));
		}

		if (m_drawPoints || m_drawBurntEdges || m_drawMC || m_drawLinesAsProxy)
		{
			if (m_drawMC)
			{
				if (m_lineTransparent)
				{
					/*gl.Enable(Capability::DepthTest);
					gl.Disable(Capability::Blend);*/
					transparent_draws.push_back(std::make_pair(m_dualLinesDrawer, false));
				}
				else
				{
					/*gl.Enable(Capability::Blend);
					gl.Disable(Capability::DepthTest);*/
					opaque_draws.push_back(std::make_pair(m_dualLinesDrawer, false));
					//m_linesDrawer->render(0.0);
				}
			}

			if (m_drawBurntEdges)
			{
				opaque_draws.push_back(std::make_pair(m_burntEdgesDrawer, false));
				// m_burntEdgesDrawer->render(0.0);
			}

			if (m_drawPoints)
			{
				std::dynamic_pointer_cast<PointDrawer>(m_pointDrawer)->setPointSize(5);
				opaque_draws.push_back(std::make_pair(m_pointDrawer, false));
				//m_pointDrawer->render(0.0);
			}
			if (m_drawLinesAsProxy)
			{
				//cout << "GLArea::paintGL(), drawing lines as proxy..."<< endl;

				m_linesProxyProg->Use();

				// radius of ball-stick
				LazyUniform<float> radius(*m_linesProxyProg, "Radius");
				this->m_meshMA->need_bbox();
				radius.Set(this->m_meshMA->bbox.radius() * 0.01f);

				// screen dimension
				LazyUniform<oglplus::Vec2f> screen(*m_linesProxyProg, "Screen");
				screen.Set(oglplus::Vec2f(width(), height()));
				// near & far 
				const auto& near_far = dynamic_pointer_cast<LineDrawer>(m_linesProxyDrawer)->getNearFar();
				LazyUniform<float> z_near(*m_linesProxyProg, "Znear");
				LazyUniform<float> z_far(*m_linesProxyProg, "Zfar");
				z_near.Set( near_far[0] );
				z_far.Set( near_far[1] );

				// eye position & camera basis axes
				const auto& track_ball = dynamic_pointer_cast<LineDrawer>(m_linesProxyDrawer)->getCamera();
				auto model_view = track_ball->getViewMatrix() * track_ball->getCamMatrix();
				const auto& col4 = model_view.Col(3);
				const auto& X = model_view.Col(0);
				const auto& Y = model_view.Col(1);
				const auto& Z = model_view.Col(2);

				LazyUniform<oglplus::Vec3f> eye(*m_linesProxyProg, "EyePos");
				using oglplus::Dot;
				eye.Set(oglplus::Vec3f(-Dot(X, col4), -Dot(Y, col4), -Dot(Z, col4)));

				float ratio = (float)height() / width();
				float tan30 = oglplus::Tan(Degrees(30));
				LazyUniform<oglplus::Vec3f> U(*m_linesProxyProg, "U");
				U.Set(X.xyz() * z_near * tan30);
				LazyUniform<oglplus::Vec3f> V(*m_linesProxyProg, "V");
				V.Set(Y.xyz() * z_near * tan30 * ratio );
				LazyUniform<oglplus::Vec3f> W(*m_linesProxyProg, "W");
				W.Set(-Z.xyz()*z_near);

				//std::dynamic_pointer_cast<LineDrawer>(m_linesProxyDrawer)->renderAsProxy(0.0);
				opaque_draws.push_back(std::make_pair(m_linesProxyDrawer, false));
				//m_linesProxyDrawer->render(0.0);

				Program::UseNone();

				//cout << "GLArea::paintGL(), drawing lines as proxy done."<< endl;
			}
		}

		if (m_drawHS)
		{
			/*gl.Enable(Capability::DepthTest);
			gl.Disable(Capability::Blend);*/
			/*m_hsFaceDrawer->render(0.0);
			m_hsLineDrawer->render(0.0);*/
			opaque_draws.push_back(std::make_pair(m_hsFaceDrawer, true));
			opaque_draws.push_back(std::make_pair(m_hsLineDrawer, false));

		}

		if (m_drawMP)
		{
			/*gl.Enable(Capability::DepthTest);
			gl.Disable(Capability::Blend);*/
			//m_MPDrawer->render(0.0);
			opaque_draws.push_back(std::make_pair(m_MPDrawer, true));

			//cout << "m_drawMP = true"<<endl;
		}

		if (m_drawIsoSurf)
		{
			//m_isoSurfDrawer->render(0.0);
			opaque_draws.push_back(std::make_pair(m_isoSurfDrawer, true));
			//cout << "GLArea::paintGL() - iso surf rendered!" << endl; //debug
		}

		if (m_drawIsoCont)
		{
			//m_iso_cont->render(true, transparent_draws, opaque_draws);
			m_iso_cont->submitDrawCalls(transparent_draws, opaque_draws);
			opaque_draws.push_back( std::make_pair(m_dynamic_dualLineDrawer, false) );
		}

		if (m_drawQMAT)
		{
			//cout << "qmat rendered." << endl;
			opaque_draws.push_back( std::make_pair(m_qmatEdgeDrawer, false));
			opaque_draws.push_back( std::make_pair(m_qmatFaceDrawer, true));
		}

		/// From here we implement depth-peeling based transparency
		int fg_depth_idx = FGDepthTexID, ref_depth_idx = RefDepthTexID, comp_depth_idx = CompositeDepthTexID;
		if (use_DP_transparency) // use depth peeling for transparency
		{
			/* render opaque object to blend FBO */
			/// to init bg texture & init depth buffer.
			/// use normal depth compare function.
			peel_blend_fbo->Bind(FramebufferTarget::Draw);
			peel_blend_fbo->AttachTexture(
				Framebuffer::Target::Draw,
				FramebufferAttachment::Depth,
				*depth_tex[RefDepthTexID],
				0
				);
			gl.DrawBuffer(fbo_drawBuffers[bg_attachPoint]);
			gl.Enable(Capability::DepthTest);
			gl.DepthMask(true);
			gl.DepthFunc(CompareFunction::Less);
			gl.ClearDepth(1.0f);
			gl.Disable(Capability::Blend);
			gl.ClearColor(bg_color[0], bg_color[1], bg_color[2], 0.0f);
			gl.Clear().ColorBuffer().DepthBuffer();
			ProgramUniform<GLuint>(*m_simpProg, "Peel").Set(0u);
			ProgramUniform<GLuint>(*m_linesProg, "Peel").Set(0u);
			ProgramUniform<GLuint>(*m_edgeProg, "Peel").Set(0u);
			for (auto iter = opaque_draws.begin(); iter != opaque_draws.end(); ++iter)
			{
// 				if ( (iter)->second ) // this render-call draws faces
// 				{
// 					gl.Enable(Capability::PolygonOffsetFill);
// 				}
// 				else
// 				{
// 					gl.Disable(Capability::PolygonOffsetFill);
// 				}
				(iter->first)->render(0.0);
			}
			
			/* second, render transparent object to peel_blend fbo. */
			/// by repeatedly peeling off the farthest geometry layer
			//
			int i = 0;
			while (i != m_DPMaxRenders)
			{
				/* Stage 1: render to peel fbo to get new depth info */
				peel_blend_fbo->Bind(FramebufferTarget::Draw);
				// peeling prog needs to access ref depth texture
				Texture::Active(refDepth_tex_unit);
				depth_tex[ref_depth_idx]->Bind(Texture::Target::_2D);
				// attach fg depth tex to Depth point
				peel_blend_fbo->AttachTexture(
					Framebuffer::Target::Draw,
					FramebufferAttachment::Depth,
					*depth_tex[fg_depth_idx],
					0);
				gl.Enable(Capability::DepthTest);
				gl.DepthMask(true);
				gl.DepthFunc(CompareFunction::Greater);
				gl.ClearDepth(0.0f);
				gl.ClearColor(0.0f, 0.0f, 0.0f, 0.0f);
				gl.DrawBuffer( fbo_drawBuffers[fg_attachPoint] );
				gl.Clear().ColorBuffer().DepthBuffer();
				ProgramUniform<GLuint>(*m_simpProg, "Peel").Set(1u);
				ProgramUniform<GLuint>(*m_linesProg, "Peel").Set(1u);
				ProgramUniform<GLuint>(*m_edgeProg, "Peel").Set(1u);
				for (auto iter = transparent_draws.begin(); iter != transparent_draws.end(); ++iter)
				{
// 					if ( (iter)->second ) // this render-call draws faces
// 					{
// 						gl.Enable(Capability::PolygonOffsetFill);
// 					}
// 					else
// 					{
// 						gl.Disable(Capability::PolygonOffsetFill);
// 					}
					(iter->first)->render(0.0);
				}

				/* Stage 2: Composite stages */
				m_bigQuad_vao->Bind();
				{
					/* Stage 2.1: composite ref. and fg. depth. */
					// post processing prog needs to access ref & fg depth textures 
					Texture::Active(refDepth_tex_unit);
					depth_tex[ref_depth_idx]->Bind(Texture::Target::_2D);
					Texture::Active(fgDepth_tex_unit);
					depth_tex[fg_depth_idx]->Bind(Texture::Target::_2D);
					// post prog will write to ref depth tex simply
					gl.Enable(Capability::DepthTest);
					gl.DepthFunc(CompareFunction::Always);
					gl.DepthMask(true);
					peel_blend_fbo->AttachTexture(
						Framebuffer::Target::Draw,
						FramebufferAttachment::Depth,
						*depth_tex[ref_depth_idx],
						0);
					gl.ColorMask(false, false, false, false);
					post_prog->Use();
					ProgramUniform<GLuint>(*post_prog, "CompositeStage").Set(1u);
					gl.DrawArrays(PrimitiveType::Triangles, 0, 6);

					/* Stage 2.2: composite fg and bg color to default framebuffer */
					// post prog needs fg color tex
					Texture::Active(fg_tex_unit);
					fg_tex->Bind(Texture::Target::_2D);
					Texture::Active(bg_tex_unit);
					bg_tex->Bind(Texture::Target::_2D);

					// run the shader program
					ProgramUniform<GLuint>(*post_prog, "CompositeStage").Set(2u);
					gl.Disable(Capability::DepthTest);
					gl.DepthMask(false);
					gl.ColorMask(true, true, true, true);
					Framebuffer::BindDefault(Framebuffer::Target::Draw);
					gl.DrawBuffer(ColorBuffer(GL_BACK));
					gl.DrawArrays(PrimitiveType::Triangles, 0, 6);

					// save the color buffer content to the bg_attachment of peel_fbo
					Framebuffer::BindDefault(Framebuffer::Target::Read);
					gl.ReadBuffer(ColorBuffer::Back);
					peel_blend_fbo->Bind(Framebuffer::Target::Draw);
					gl.DrawBuffer(fbo_drawBuffers[bg_attachPoint]);
					gl.BlitFramebuffer(
						0, 0, width(), height(), 0, 0, width(), height(),
						Bitfield<BufferSelectBit>(GL_COLOR_BUFFER_BIT), BlitFilter::Nearest
						);
				}
				m_bigQuad_vao->Unbind();

				i++;
			}
		}
		else // opengl built-in transparency
		{
			ProgramUniform<GLuint>(*m_simpProg, "Peel").Set(0u);
			ProgramUniform<GLuint>(*m_linesProg, "Peel").Set(0u);
			ProgramUniform<GLuint>(*m_edgeProg, "Peel").Set(0u);

			Framebuffer::BindDefault(Framebuffer::Target::Draw);
			gl.DrawBuffer(ColorBuffer::Back);
			gl.Enable(Capability::DepthTest);
			gl.Disable(Capability::Blend);
			gl.DepthMask(true);
			gl.DepthFunc(CompareFunction::Less);
			gl.ClearDepth(1.0f);
			gl.ClearColor(bg_color[0], bg_color[1], bg_color[2], 1.0);
			gl.Clear().ColorBuffer().DepthBuffer();

			for (auto iter = opaque_draws.begin(); iter != opaque_draws.end(); ++iter)
			{
				//if ( (iter)->second ) // this render-call draws faces
				//{
				//	gl.Enable(Capability::PolygonOffsetFill);
				//}
				//else
				//{
				//	gl.Disable(Capability::PolygonOffsetFill);
				//}
				(iter->first)->render(0.0);
			}
			/* Then draw transparent stuff to update Accum and Count texture
			disable depth-writing, keep accum fbo bound */
			gl.Enable(Capability::DepthTest);
			gl.DepthMask(false);
			gl.Enable(Capability::Blend);
			for (auto iter = transparent_draws.begin(); iter != transparent_draws.end(); ++iter)
			{
				//if ( (iter)->second ) // this render-call draws faces
				//{
				//	gl.Enable(Capability::PolygonOffsetFill);
				//}
				//else
				//{
				//	gl.Disable(Capability::PolygonOffsetFill);
				//}
				(iter->first)->render(0.0);
			}
		} // END: opengl alpha blending rendering 

		///* Lastly, composite the bg_tex with the accum_tex 
		//by drawing a full-screen quad
		//*/
		//gl.Disable(Capability::Blend);
		//gl.Disable(Capability::DepthTest);
		//gl.DepthMask(true);
		//m_bigQuad_vao->Bind();
		//{
		//	Framebuffer::BindDefault(Framebuffer::Target::Draw);
		//	gl.ClearColor(0.0f, 0.0f, 0.0f, 0.0f);
		//	gl.ClearDepth(1.0f);
		//	gl.Clear().ColorBuffer().DepthBuffer();
		//	m_compositeProg->Use();
		//	gl.DrawArrays(PrimitiveType::Triangles, 0, 6);

		//	/*opaque_vao.Bind();
		//	gl.DrawArrays(PrimitiveType::Triangles, 0, 3);
		//	opaque_vao.Unbind();*/
		//}
		//m_bigQuad_vao->Unbind();
		//cout << "GLArea::paintGL() done." << endl;
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

void GLArea::mouseMoveEvent(QMouseEvent *_evnt)
{
	assert(m_MADrawer && m_origDrawer);
	
	m_MADrawer->MouseMove(
		_evnt->x(), 
		height() - _evnt->y(),
		width(),
		height()
		);

	updateGL();
	/*cout << "qmouse: " << height() - _evnt->y() << endl;*/
}

void GLArea::mousePressEvent(QMouseEvent *_evnt)
{
	m_MADrawer->MousePress(_evnt);
	updateGL();
}

void GLArea::mouseReleaseEvent(QMouseEvent *_evnt)
{
	m_MADrawer->MouseRelease(_evnt);
	updateGL();
}

void GLArea::mouseDoubleClickEvent(QMouseEvent *_evnt)
{
	//if (m_fullMode== false)
	//{
	//	cout << "setting to full..." << endl;
	//	m_enOrigWindowFlags = this->windowFlags();
	//	m_pSize = this->size();
	//	m_pParent = (QWidget*)this->parent();
	//	cout << "changing parent to 0..." << endl;
	//	this->setParent(0);
	//	/*cout << "inserting window flags..." << endl;
	//	this->setWindowFlags( Qt::FramelessWindowHint|Qt::WindowStaysOnTopHint);*/
	//	cout << "going full..." << endl;
	//	this->showMaximized();
	//	m_fullMode = true;
	//}
	//else
	//{
	//	cout << "setting back to normal size..." << endl;
	//	cout << "changing parent back..." << endl;
	//	this->setParent(m_pParent);
	//	cout << "changing sizes back..." << endl;
	//	this ->resize(m_pSize);
	//	/*cout << "changing window flags back..." << endl;
	//	this->overrideWindowFlags(m_enOrigWindowFlags);*/
	//	cout << "showing..." << endl;
	//	this->show();
	//	m_fullMode =  false;
	//}
	//this->showFullScreen();
}


void GLArea::initShaders()
{
	try {
		m_simpProg = glslProgram("./shaders/DiffuseOnly.vert.glsl", 
			"./shaders/DiffuseOnly.frag.glsl");
		m_edgeProg = glslProgram("./shaders/DiffuseWithEdges.vert.glsl", 
			"./shaders/DiffuseWithEdges.geom.glsl", 
			"./shaders/DiffuseWithEdges.frag.glsl");
		m_pointsProg = glslProgram("./shaders/DrawPoints.vert.glsl", 
			"./shaders/DrawPoints.frag.glsl");
		m_linesProg = glslProgram("./shaders/DrawLines.vert.glsl", 
			"./shaders/DrawLines.frag.glsl");
		m_linesProxyProg = glslProgram("./shaders/realCylinder.vert.glsl",
			"./shaders/realCylinder.geom.glsl",
			"./shaders/realCylinder.frag.glsl");
		peel_prog = glslProgram("./shaders/simpleDP.vert.glsl",
			"./shaders/peelDP.frag.glsl");
		post_prog = glslProgram("./shaders/simpleDP.vert.glsl",
			"./shaders/postDP.frag.glsl");
		/*m_compositeProg = glslProgram("./shaders/simple.vert.glsl",
			"./shaders/composite.frag.glsl");*/

		/*m_svpProg = glslProgram("./shaders/SVP.vert.glsl",
		"./shaders/SVP.geom.glsl", "./shaders/SVP.frag.glsl");*/
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

void GLArea::bindTexture2D(Texture& _name, int _unit,
	int _w, int _h,
	TextureMinFilter _min_filter,
	TextureMagFilter _mag_filter,
	TextureWrap _s_wrap,
	TextureWrap _t_wrap,
	PixelDataInternalFormat _internal,
	PixelDataFormat _data_format,
	PixelDataType _data_type,
	const void* _data
	)
{
	Texture::Active(_unit);
	_name.Bind(Texture::Target::_2D);
	_name.MinFilter(Texture::Target::_2D, _min_filter);
	_name.MagFilter(Texture::Target::_2D, _mag_filter);
	_name.WrapS(Texture::Target::_2D, _s_wrap);
	_name.WrapT(Texture::Target::_2D, _t_wrap);
	_name.Image2D(
		Texture::Target::_2D,
		0,
		_internal,
		_w, _h,
		0,
		_data_format,
		_data_type,
		_data
		);
	Texture::Unbind(Texture::Target::_2D);
}

void GLArea::initOIT()
{
	using namespace oglplus;
	try 
	{
		depth_tex_OIT = std::shared_ptr<Texture>(new Texture());
		bg_tex_OIT = std::shared_ptr<Texture>(new Texture());
		accum_tex_OIT = std::shared_ptr<Texture>(new Texture());
		count_tex_OIT = std::shared_ptr<Texture>(new Texture());

		m_accum_fbo_OIT = std::shared_ptr<Framebuffer>(new Framebuffer());

		m_bigQuad_vao = std::shared_ptr<VertexArray>(new VertexArray());
		m_bigQuad_vts_vbo = std::shared_ptr<Buffer>(new Buffer());

		/// Setup textures
		initTexturesWithSizesOIT(width(), height());

		/// Setup fbo's
		m_accum_fbo_OIT->Bind(Framebuffer::Target::Draw);
		m_accum_fbo_OIT->AttachTexture(
			Framebuffer::Target::Draw,
			FramebufferAttachment::Depth,
			*depth_tex_OIT,
			0
			);
		m_accum_fbo_OIT->AttachTexture(
			Framebuffer::Target::Draw,
			FramebufferAttachment::Color,
			*bg_tex_OIT,
			0
			);
		m_accum_fbo_OIT->AttachTexture(
			Framebuffer::Target::Draw,
			FramebufferAttachment::Color1,
			*accum_tex_OIT,
			0
			);
		m_accum_fbo_OIT->AttachTexture(
			Framebuffer::Target::Draw,
			FramebufferAttachment::Color2,
			*count_tex_OIT,
			0
			);
		assert(m_accum_fbo_OIT->IsComplete(Framebuffer::Target::Draw));
		Framebuffer::Unbind(FramebufferTarget::Draw);

		// geometry for full screen quad
		GLfloat bigQuad_vts[24] = {
			-1.0f,	1.0f,	0.1f,	1.0f, 
			-1.0f,	-1.0f,	0.1f,	1.0f,
			1.0f,	-1.0f,	0.1f,	1.0f,
			-1.0f,	1.0f,	0.1f,	1.0f,
			1.0f,	-1.0f,	0.1f,	1.0f,
			1.0f,	1.0f,	0.1f,	1.0f
		};
		m_bigQuad_vao->Bind();
		m_bigQuad_vts_vbo->Bind(Buffer::Target::Array);
		Buffer::Data(
			Buffer::Target::Array,
			24,
			bigQuad_vts
			);
		auto vert_attr = VertexAttribArray(*m_compositeProg, "Vertex");
		vert_attr.Setup<GLfloat>(4);
		vert_attr.Enable();
		m_bigQuad_vts_vbo->Unbind(Buffer::Target::Array);
		m_bigQuad_vao->Unbind();

		// Setup the uniforms & inputs
		ProgramUniformSampler(
			*m_compositeProg, "BGTex"
			).Set(bg_tex_unit_OIT);
		ProgramUniformSampler(
			*m_compositeProg, "AccumTex"
			).Set(accum_tex_unit_OIT);
		ProgramUniformSampler(
			*m_compositeProg, "CountTex"
			).Set(count_tex_unit_OIT);
		ProgramUniform<Vec2f>(
			*m_compositeProg, "ViewportSizes"
			).Set(width(), height());
	}
	catch(oglplus::Error& err)
	{
		std::cerr << 
			"GLArea::initOIT()" << endl <<
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
			"GLArea::initOIT()" << endl <<
			"General error: " <<
			se.what() << std::endl;
	}
}

void GLArea::resizeOIT(int _w, int _h)
{
	if (m_compositeProg == nullptr)
		return;
	ProgramUniform<Vec2f>(
		*m_compositeProg, "ViewportSizes"
		).Set(_w, _h);

	if (depth_tex_OIT == nullptr)
		return;
	initTexturesWithSizesOIT(_w, _h);
}

void GLArea::initTexturesWithSizesOIT(int _w, int _h)
{
	depth_tex_OIT->Bind(Texture::Target::_2D);
	depth_tex_OIT->MinFilter(Texture::Target::_2D, TextureMinFilter::Linear);
	depth_tex_OIT->MagFilter(Texture::Target::_2D, TextureMagFilter::Linear);
	depth_tex_OIT->WrapS(Texture::Target::_2D, TextureWrap::ClampToEdge);
	depth_tex_OIT->WrapT(Texture::Target::_2D, TextureWrap::ClampToEdge);
	depth_tex_OIT->Image2D(
		Texture::Target::_2D,
		0,
		PixelDataInternalFormat::DepthComponent32,
		_w, _h,
		0,
		PixelDataFormat::DepthComponent,
		PixelDataType::Float,
		nullptr
		);
	Texture::Unbind(Texture::Target::_2D);
	Texture::Active(bg_tex_unit_OIT);
	bg_tex_OIT->Bind(Texture::Target::_2D);
	bg_tex_OIT->MinFilter(Texture::Target::_2D, TextureMinFilter::Linear);
	bg_tex_OIT->MagFilter(Texture::Target::_2D, TextureMagFilter::Linear);
	bg_tex_OIT->WrapS(Texture::Target::_2D, TextureWrap::ClampToEdge);
	bg_tex_OIT->WrapT(Texture::Target::_2D, TextureWrap::ClampToEdge);
	bg_tex_OIT->Image2D(
		Texture::Target::_2D,
		0,
		PixelDataInternalFormat::RGBA,
		_w, _h,
		0,
		PixelDataFormat::RGBA,
		PixelDataType::Float,
		nullptr
		);
	Texture::Unbind(Texture::Target::_2D);
	Texture::Active(accum_tex_unit_OIT);
	accum_tex_OIT->Bind(Texture::Target::_2D);
	accum_tex_OIT->MinFilter(Texture::Target::_2D, TextureMinFilter::Linear);
	accum_tex_OIT->MagFilter(Texture::Target::_2D, TextureMagFilter::Linear);
	accum_tex_OIT->WrapS(Texture::Target::_2D, TextureWrap::ClampToEdge);
	accum_tex_OIT->WrapT(Texture::Target::_2D, TextureWrap::ClampToEdge);
	accum_tex_OIT->Image2D(
		Texture::Target::_2D,
		0,
		PixelDataInternalFormat::RGBA,
		_w, _h,
		0,
		PixelDataFormat::RGBA,
		PixelDataType::Float,
		nullptr
		);
	Texture::Unbind(Texture::Target::_2D);
	Texture::Active(count_tex_unit_OIT);
	count_tex_OIT->Bind(Texture::Target::_2D);
	count_tex_OIT->MinFilter(Texture::Target::_2D, TextureMinFilter::Linear);
	count_tex_OIT->MagFilter(Texture::Target::_2D, TextureMagFilter::Linear);
	count_tex_OIT->WrapS(Texture::Target::_2D, TextureWrap::ClampToEdge);
	count_tex_OIT->WrapT(Texture::Target::_2D, TextureWrap::ClampToEdge);
	count_tex_OIT->Image2D(
		Texture::Target::_2D,
		0,
		PixelDataInternalFormat::RGBA16F, // this makes += 1 right
		_w, _h,
		0,
		PixelDataFormat::RGBA,
		PixelDataType::Float,
		nullptr
		);
	Texture::Unbind(Texture::Target::_2D);
}

void GLArea::bindTexturesOIT()
{
	Texture::Active(bg_tex_unit_OIT);
	bg_tex_OIT->Bind(Texture::Target::_2D);
	Texture::Active(accum_tex_unit_OIT);
	accum_tex_OIT->Bind(Texture::Target::_2D);
	Texture::Active(count_tex_unit_OIT);
	count_tex_OIT->Bind(Texture::Target::_2D);
}

void GLArea::initDP()
{
	try
	{
		m_bigQuad_vao = std::shared_ptr<VertexArray>(new VertexArray());
		m_bigQuad_vts_vbo = std::shared_ptr<Buffer>(new Buffer());
		bg_tex = shared_ptr<Texture>(new Texture(ObjectDesc("DP background texture")));
		fg_tex = shared_ptr<Texture>(new Texture(ObjectDesc("DP foreground texture")));
		depth_tex[RefDepthTexID] = shared_ptr<Texture>(new Texture(ObjectDesc("DP ref.(B.G.) depth texture")));
		depth_tex[FGDepthTexID] = shared_ptr<Texture>(new Texture(ObjectDesc("DP F.G. depth texture")));
		depth_tex[CompositeDepthTexID] = shared_ptr<Texture>(new Texture(ObjectDesc("DP composite depth texture")));
		peel_blend_fbo = shared_ptr<Framebuffer>(new Framebuffer(ObjectDesc("DP FBO")));

		// prepare FBO attachment points
		bg_attachPoint = 0; 
		fbo_attachments.push_back(FramebufferAttachment::Color);
		fbo_drawBuffers.push_back(FramebufferColorAttachment::_0);
		fg_attachPoint = 1;
		fbo_attachments.push_back(FramebufferAttachment::Color1);
		fbo_drawBuffers.push_back(FramebufferColorAttachment::_1);
		compDepth_attachPoint = 2;
		fbo_attachments.push_back(FramebufferAttachment::Color2);
		fbo_drawBuffers.push_back(FramebufferColorAttachment::_2);
		bg_tex_unit = 0;
		fg_tex_unit = 1;
		refDepth_tex_unit = 2;
		fgDepth_tex_unit = 3;

		initTexturesWithSizesDP(width(), height());

		/// Now textures are ready. Finish setting up FBO
		peel_blend_fbo->Bind(Framebuffer::Target::Draw);
		// need one effective depth buffer
		peel_blend_fbo->AttachTexture(
			Framebuffer::Target::Draw,
			FramebufferAttachment::Depth,
			*depth_tex[RefDepthTexID],
			0
			);
		// foreground texture as one color attach
		peel_blend_fbo->AttachTexture(
			Framebuffer::Target::Draw,
			fbo_attachments[fg_attachPoint],
			*fg_tex,
			0
			);
		// background texture as one color attach
		peel_blend_fbo->AttachTexture(
			Framebuffer::Target::Draw,
			fbo_attachments[bg_attachPoint],
			*bg_tex,
			0
			);
		assert(peel_blend_fbo->IsComplete(Framebuffer::Target::Draw));
		Framebuffer::Unbind(FramebufferTarget::Draw);

		/* Setup inputs and uniforms */
		GLfloat bigQuad_vts[24] = {
			-1.0f, 1.0f, 0.1f, 1.0f,
			-1.0f, -1.0f, 0.1f, 1.0f,
			1.0f, -1.0f, 0.1f, 1.0f,
			-1.0f, 1.0f, 0.1f, 1.0f,
			1.0f, -1.0f, 0.1f, 1.0f,
			1.0f, 1.0f, 0.1f, 1.0f
		};
		m_bigQuad_vao->Bind();
		m_bigQuad_vts_vbo->Bind(Buffer::Target::Array);
		Buffer::Data(
			Buffer::Target::Array,
			24,
			bigQuad_vts
			);
		auto vert_attr = VertexAttribArray(*peel_prog, "Vertex");
		vert_attr.Setup<GLfloat>(4);
		vert_attr.Enable();
		m_bigQuad_vts_vbo->Unbind(Buffer::Target::Array);

		// Setup the uniforms & inputs
		// depth samplers
		ProgramUniformSampler(
			*m_simpProg, "DepthRefTex"
			).Set(refDepth_tex_unit);
		ProgramUniformSampler(
			*m_linesProg, "DepthRefTex"
			).Set(refDepth_tex_unit);
		ProgramUniformSampler(
			*m_edgeProg, "DepthRefTex"
			).Set(refDepth_tex_unit);
		ProgramUniformSampler(
			*post_prog, "DepthRefTex"
			).Set(refDepth_tex_unit);
		ProgramUniformSampler(
			*post_prog, "DepthFGTex"
			).Set(fgDepth_tex_unit);
		// color samplers
		ProgramUniformSampler(
			*post_prog, "FGColorTex"
			).Set(fg_tex_unit);
		ProgramUniformSampler(
			*post_prog, "BGColorTex"
			).Set(bg_tex_unit);

		cout << "DP uniform samplers setup." << endl;
	}
	catch(oglplus::Error& err)
	{
		std::cerr << 
			"GLArea::initOIT()" << endl <<
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
			"GLArea::initOIT()" << endl <<
			"General error: " <<
			se.what() << std::endl;
	}
}

void GLArea::initTexturesWithSizesDP(int _w, int _h)
{
	/// Setup textures
	// first 3 depth textures: ref depth, fg depth, and composite depth
	bindTexture2D(*depth_tex[RefDepthTexID], refDepth_tex_unit, _w, _h,
		TextureMinFilter::Linear, TextureMagFilter::Linear,
		TextureWrap::ClampToEdge, TextureWrap::ClampToEdge,
		PixelDataInternalFormat::DepthComponent32, PixelDataFormat::DepthComponent, PixelDataType::Float, nullptr);
	bindTexture2D(*depth_tex[FGDepthTexID], fgDepth_tex_unit, _w, _h,
		TextureMinFilter::Linear, TextureMagFilter::Linear,
		TextureWrap::ClampToEdge, TextureWrap::ClampToEdge,
		PixelDataInternalFormat::DepthComponent32, PixelDataFormat::DepthComponent, PixelDataType::Float, nullptr);
	int dummy_unit = 0;
	bindTexture2D(*depth_tex[CompositeDepthTexID], dummy_unit, _w, _h,
		TextureMinFilter::Linear, TextureMagFilter::Linear,
		TextureWrap::ClampToEdge, TextureWrap::ClampToEdge,
		PixelDataInternalFormat::DepthComponent32, PixelDataFormat::DepthComponent, PixelDataType::Float, nullptr);
	// then fg color tex, and bg color tex
	bindTexture2D(*fg_tex, fg_tex_unit, _w, _h,
		TextureMinFilter::Nearest, TextureMagFilter::Nearest,
		TextureWrap::ClampToEdge, TextureWrap::ClampToEdge,
		PixelDataInternalFormat::RGBA, PixelDataFormat::RGBA, PixelDataType::Float, nullptr);
	bindTexture2D(*bg_tex, bg_tex_unit, _w, _h,
		TextureMinFilter::Nearest, TextureMagFilter::Nearest,
		TextureWrap::ClampToEdge, TextureWrap::ClampToEdge,
		PixelDataInternalFormat::RGBA, PixelDataFormat::RGBA, PixelDataType::Float, nullptr);
}

void GLArea::resizeDP(int _w, int _h)
{
	if (fg_tex == nullptr)
		return;
	initTexturesWithSizesDP(_w, _h);
}

void GLArea::initVertPicking()
{
	vertID_tex_SVP = shared_ptr<Texture>(new Texture(ObjectDesc("Vertex ID texture (screen picking)")));
	vertID_tex_unit = 0;

	fbo_drawbuffer_SVP = FramebufferColorAttachment::_0;

	initTexturesWithSizesSVP(width(), height());

	SVP_fbo->Bind(Framebuffer::Target::Draw);
	SVP_fbo->AttachTexture(
		Framebuffer::Target::Draw,
		FramebufferAttachment::Color,
		*vertID_tex_SVP,
		0
		);
	SVP_fbo->AttachTexture(
		Framebuffer::Target::Draw,
		FramebufferAttachment::Depth,
		*depth_tex_SVP,
		0
		);
	assert(SVP_fbo->IsComplete(Framebuffer::Target::Draw));
	Framebuffer::Unbind(FramebufferTarget::Draw);

	auto drawer = std::dynamic_pointer_cast<MeshDrawer>(m_drawerSVP);
	auto rendermode = drawer->getRenderMode();
	drawer->setRenderMode(MeshDrawer::PER_VERT);
	vector<unsigned int> vertIDs;
	for (size_t i = 0; i < stg->m_origG->vts.size(); ++i)
		vertIDs.push_back(i);
	drawer->setPointsSVP(stg->m_origG->vts, vertIDs, "VertexID");
}

void GLArea::initTexturesWithSizesSVP(int _w, int _h)
{

}

void GLArea::pickingPhaseSVP()
{
	// bind SVP fbo
	SVP_fbo->Bind(FramebufferTarget::Draw);
	
	// set up correct gl render params
	gl.DrawBuffer(fbo_drawbuffer_SVP);
	gl.Enable(Capability::DepthTest);
	gl.Disable(Capability::Blend);
	gl.DepthMask(true);
	gl.DepthFunc(CompareFunction::Less);
	gl.ClearDepth(1.0f);
	gl.ClearColor(0.0f, 0.0f, 0.0f, 0.0f);
	gl.Clear().ColorBuffer().DepthBuffer();

	// render all visible geometry
	if (m_drawMA)
	{
		m_svpProg->Use();
		auto drawer = std::dynamic_pointer_cast<MeshDrawer>(m_drawerSVP);
		drawer->render(0.0);
		oglplus::Program::UseNone();
	}

	// do normal gl draw call...
	updateGL();
}

void GLArea::uploadSimplicialComplex( 
	const vector<TriPoint>& _vts, const vector<TriEdge>& _edges, const vector<TriFace>& _faces,
	const vector<TriColor>* _e_colors, const vector<TriColor>* _f_colors,
	shared_ptr<Drawable>& _edge_drawer, shared_ptr<Drawable>& _face_drawer )
{
	cout << "uploading simplicial complex to GPU..." << endl;
	// obtain edges and faces of HS
	const auto& vts_hs = _vts;
	const auto& edges_hs = _edges;
	const auto& tri_faces_hs = _faces;

	auto hs_face_drawer_ptr = std::dynamic_pointer_cast<MeshDrawer>(_face_drawer);
	hs_face_drawer_ptr->setPointsPerFace(vts_hs, tri_faces_hs);
	cout << "setPointsPerFace done!" << endl;

	// TODO: color the HS face and line according to the cur dist metrics
	// now just use default uniform color
	
	const TriColor const_face_color(1.0f, 98.0f/255.0f, 0.0f);
	const TriColor const_edge_color(0.0f, 0.0f, 0.0f);
	// setup HS faces
	float* color_data = new float[tri_faces_hs.size()*3];
	for (unsigned i = 0; i < tri_faces_hs.size()*3; i += 3)
	{
		const auto& face_color = _f_colors ? ( *_f_colors )[ _f_colors->size() == 1 ? 0 : i ] : const_face_color;
		color_data[i+0] = face_color[0];
		color_data[i+1] = face_color[1];
		color_data[i+2] = face_color[2];
	}
	hs_face_drawer_ptr->setPerFaceColor(color_data, tri_faces_hs.size());
	delete [] color_data;
	color_data = nullptr;
	
	float* normal_data = new float[tri_faces_hs.size()*3];
	for (unsigned i = 0; i < tri_faces_hs.size(); ++i)
	{
		auto& f = tri_faces_hs[i];
		auto& nml = (vts_hs[f[0]] - vts_hs[f[1]]).cross(vts_hs[f[0]] - vts_hs[f[2]]);
		trimesh::normalize(nml);
		normal_data[i*3 + 0] = nml[0];
		normal_data[i*3 + 1] = nml[1];
		normal_data[i*3 + 2] = nml[2];
	}
	hs_face_drawer_ptr->setPerFaceNormal(normal_data, tri_faces_hs.size());
	delete [] normal_data;
	normal_data = nullptr;
	
	float* saliency = new float[tri_faces_hs.size()];
	for (unsigned i = 0; i < tri_faces_hs.size(); ++i)
	{
		saliency[i] = 0.2f;
	}
	hs_face_drawer_ptr->setPerFaceSaliency(saliency, tri_faces_hs.size());
	delete [] saliency;
	saliency = nullptr;

	hs_face_drawer_ptr->setRenderMode(MeshDrawer::PER_FACE);

	// setup HS lines
	auto hs_line_drawer_ptr = std::dynamic_pointer_cast<LineDrawer>(_edge_drawer);
	/*hs_line_drawer_ptr->setPoints(vts_hs);
	hs_line_drawer_ptr->setLines(edges_hs);*/
	uploadLinesToDrawer(vts_hs, edges_hs, hs_line_drawer_ptr);
	cout << "setPoints & setLines done!"<<endl;

	color_data = new float[ edges_hs.size() * 2 * 3 ];
	for ( unsigned i = 0; i < edges_hs.size(); i++ )
	{
		const auto& edge_color = _e_colors ? ( *_e_colors )[ _e_colors->size() == 1 ? 0 : i ] : const_edge_color;
		color_data[ i * 2 * 3 + 0 ] = edge_color[ 0 ];
		color_data[ i * 2 * 3 + 1 ] = edge_color[ 1 ];
		color_data[ i * 2 * 3 + 2 ] = edge_color[ 2 ];
		color_data[ i * 2 * 3 + 3 + 0 ] = edge_color[ 0 ];
		color_data[ i * 2 * 3 + 3 + 1 ] = edge_color[ 1 ];
		color_data[ i * 2 * 3 + 3 + 2 ] = edge_color[ 2 ];
	}
	//stg->m_hs.getEdgeDiffValColor(color_data);
	cout << "GLArea::setting line color ..."<<endl;
	hs_line_drawer_ptr->setPerVertColor(color_data, edges_hs.size() * 2);
	delete [] color_data;
	color_data = nullptr;
	cout << "LineDrawer::setPerVertColor done!"<<endl;

	cout << "uploading simplicial complex done!" << endl;

	updateGL();
}