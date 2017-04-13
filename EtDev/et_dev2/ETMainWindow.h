#ifndef QTPUREGL_H
#define QTPUREGL_H

#include "glArea.h"

#include <QMainWindow>
#include <QtOpenGL/QGLWidget>
#include <QMessageBox>
#include <QFileDialog>
#include <string>
#include "ui_compact.h"
#include "geometryRepresentation.h"

class ETMainWindow : public QMainWindow
{
	Q_OBJECT

public:
	ETMainWindow(QWidget *parent = 0, Qt::WindowFlags flags = 0);
	~ETMainWindow();
	
	//
	// set the program into debug mode
	//
	void setDebugMode(bool _enabled);
	//
	// set the program into console mode
	//
	void setConsoleMode( bool _enabled );
	//
	// remember the inputs, including file names, parameter settings, etc.
	//
	void setInputs( 
		const std::string& _ma_file, const std::string& _shape_file, const std::string& _r_file,
		float _omega,
		int _burn_sch_id = 0 /*default: steinergraph::orig-and-steiner*/ );

public slots:
	void onSaveView();
	void onLoadView();
	void onLoadFilesClicked();
	void onBrowseMAClicked();
	void onBrowseOrigClicked();
	void onBrowseRadiusClicked();

	void onTrueTransparencyChanged(int _s);
	void onDPMaxRendersSpun(int _v);

	void onBurnBtnClicked();
	void onCleanTopoBtnClicked();
	void onDrawEdgeChecked(int _state);
	void onDrawStPointsChecked(int _state);
	void onPointSizeSpun(double _value);
	void onHideMAChanged(int _state);
	void onHideOrigChanged(int _state);
	void onHideMCChanged(int _state);
	void onVisBurntEdgesChanged(int _state);
	void onHideHSChanged(int _state);
	void onHideIsoSurfChanged(int _state);
	void onHideIsoContourChanged(int _state);

	void onTurnOffLightingForMAChanged(int _state);

	void onComputeAllDistClicked();
	void onEnableFinePrune(int _state);
	void onPruneMASlided_1(int _value);
	void onPruneMASpun_1(double _value);
	void onPruneMASlided_2(int _value);
	void onPruneMASpun_2(double _value);
	void onOrigTransparentSlided(int _value);
	// listen to either min_alpha or exp_alpha change
	void onMATransparentSlided(int _value);
	void onMCTransparentSlided(int _value);
	void onIsoContTransparentSlided(int _value);
	void onMCWidthSpinChanged(double _value);
	
	//void onColorPerFaceChecked(int _state);

	void onOrigColorClicked();
	void onMAColorClicked();
	void onMCColorClicked();
	void onChangeBGColorClicked();
	void onColorPerVertBtnClicked();
	void onColorPerFaceBtnClicked();
	void onVisDistMAComboChanged(int _idx);
	void onClampMADistBoxChanged(int _state);
	void onMaxVisDistMASpun(double _value);
	void onMinVisDistMASpun(double _value);
	void onPruneDistMAComboChanged_1(int _idx);
	void onPruneDistMAComboChanged_2(int _idx);

	/* medial curve creating, burning, pruning... */
	void onDualizeClicked();
	void onBurnMCClicked();
	void onVisDistMCComboChanged(int _idx);
	void onPruneCurveSlided(int _value);
	void onPruneCurveSpinChanged(double _val);
	bool onPruneDistMCComboChanged(int _idx);
	void onClampMCDistBoxChanged(int _state);
	void onMaxVisDistMCSpun(double _value);
	void onMinVisDistMCSpun(double _value);
	void onVisBallStickChecked(int _state);
	void onPrintStatsClicked();

	/* hybrid skeleton exploration */
	void onCreateHSClicked();
	void onVisHSClicked();
	void onHSRemoveSmallCmpntsChecked(int _state);
	void onDegenerateFaceThreshChanged(int _val);
	void onHSTheta2Changed(int _val);
	void onHSTheta2Spun(double _val);
	void onHSTheta1Changed(int _val);
	void onHSTheta1Spun(double _val);

	/* surf func exploration */
	void onSurfFuncCorrespSetupClicked();
	void onSurfFuncCorrespSchemeComboChanged(int _idx);
	void onSurfFuncComboChanged(int _idx);
	void onSurfFuncSlided(int _value);
	void onDiffuseSurfFuncChecked(int _state);
	void onSurfFuncProjectClicked();

	/* medial point exploration */
	void onShowMPChecked(int _state);

	/* iso surface evolution */
	void onIsoSurfPrecomputeClicked();
	void onIsoSurfRefineClicked();
	void onIsoSurfSlided(int _value);
	void onHideSnappedFacesChanged(int _state);

	/* iso contour evolution */
	void onIsoContPrecomputeClicked();
	void onIsoContSlided(int _value);
	void onHidePastSurfChanged(int _state);
	void onIsoContShowMCChanged(int _state);

	/* choreograph all iso-objects */
	void onEvolveAllPrecomputeClicked();
	void onEvolveAllSlided(int _value);

	/* MISC */
	void onPrintBTClicked();
	void onShowSpheresChecked(int _state);
	void onVisSphereSlided(int _val);
	void onVisSphereSpun(int _val);
	// print info about the selected original vertex in MA mesh
	void onPrintSelectVertInfoClicked();
	void onVisSpheresFromFileClicked();
	bool readVertsFromFile(QString _file, vector<unsigned>& _vid_list);
	// set lines-rendering offset
	void onLinesRenderOffsetSlided(int _val);

	// output current MA to wenping's qmat input format
	void onOutputWenpingBtnClicked();
	// read in qmat's medial axis file and display the structure
	void onReadQMATBtnClicked();
	// 
	void onHideQMATChanged(int _state);

	/* EXPORT */
	void onExportETBtnClicked();
	void onExportHSBtnClicked();
	void onExportBtnClicked();

private:
	//
	// steps that are usually referred to internally
	//
	enum FineStep {
		OPEN_FILE, POST_OPEN_FILE, BURN_MA, PRUNE_MA, 
		DUALIZE, BURN_MC, POST_BURN_MC, COMPUTE_HS, PRUNE_HS, WRITE_BT_SKEL
	};
	//
	// steps that make more sense to users
	//
	enum HighlevelStep {IMPORT, BURN, SKELETON, EXPORT};
	//
	// return the set of widgets/action crucial for the specified step
	//
	void getRelevantWidgets(FineStep _stp, QList<QObject*>& _list);
	//
	// prepare program for the specified step
	//
	void prepareForStep(FineStep _stp);
	//
	// disable widgets relevant to the steps at & after the *specified step*
	//
	void disableWidgetsForStep(FineStep _stp);
	//
	// initialize UI to a proper state (e.g. hiding extra details)
	//
	void initUI();

private:
	GLArea::FaceFieldType get_MA_face_dist_type(int _idx);
	GLArea::VertFieldType get_MA_vert_dist_type(int _idx);
	GLArea::DistMC get_MC_dist_type(int _idx);

	QColor getIdealTextColor(const QColor& rBackgroundColor) const;
	void resetWidgetsStates();
	void resetParams(FineStep _stp);

	void resetPruneMC();
	void enablePruneMC();

	/// surf func related
	void resetSurfFuncWidgets();
	void enableSurfFuncWidgets();
	void enableSmoothSurfFuncWidgets();

	void getPruneMAThreshold(float& _v1, float& _v2);

	float get_cur_isoCont_t();

private:
	Ui::MainWindow ui_compact;
	GLArea *glarea;

	/* maintain some widgets' values here */
	// values for MA visualization
	float visDistOnMA_min, visDistOnMA_max;
	// values for MA pruning
	float pruneDistOnMA_min_1, pruneDistOnMA_max_1;
	float pruneDistOnMA_min_2, pruneDistOnMA_max_2;

	// values for medial curve pruning
	// current range of medial curve's distance metric
	float distMin_MC, distMax_MC;

	// values for hs pruning
	float ET_ma_min, ET_ma_max;
	float ET_mc_min, ET_mc_max;

	// values for iso-objects evolution
	float iso_surf_max;
	float iso_cont_max;
	float iso_point_max;

	/* program states */
	bool m_debugMode;
	bool m_consoleMode;

	/*related to console mode.
	**TODO: process console mode in another class!!!*/
	// input file names
	std::string m_ma_file;
	std::string m_shape_file;
	std::string m_r_file;
	// sampling rate
	float m_omega;
	// burning scheme
	SteinerGraph::BurnScheme m_burn_sch;
};

#endif // QTPUREGL_H
