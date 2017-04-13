#include "ETMainWindow.h"
#include <QColorDialog>
#include <QSettings>
#ifdef PRINT_MEM_USAGE
#include "windows.h"
#include "psapi.h"
#endif // PRINT_MEM_USAGE

ETMainWindow::ETMainWindow(QWidget *parent, Qt::WindowFlags flags)
	: QMainWindow(parent, flags)
{
	ui_compact.setupUi(this);
	
	QObject::connect( ui_compact.saveView, SIGNAL(triggered()), this, SLOT(onSaveView()) );
	QObject::connect( ui_compact.loadView, SIGNAL(triggered()), this, SLOT(onLoadView()) );
	//QObject::connect( ui_compact.readRadii, SIGNAL(triggered()), this, SLOT(onReadRadii()) );
	QObject::connect( ui_compact.loadFilesBtn, SIGNAL(clicked()), this, SLOT(onLoadFilesClicked()) );
	QObject::connect( ui_compact.saveViewBtn, SIGNAL(clicked()), this, SLOT(onSaveView()) );
	QObject::connect( ui_compact.loadViewBtn, SIGNAL(clicked()), this, SLOT(onLoadView()) );
	QObject::connect( ui_compact.browseMAFileBtn, SIGNAL(clicked()), this, SLOT(onBrowseMAClicked()) );
	QObject::connect( ui_compact.browse3DShapeFileBtn, SIGNAL(clicked()), this, SLOT(onBrowseOrigClicked()) );
	QObject::connect( ui_compact.browseRadiusFileBtn, SIGNAL(clicked()), this, SLOT(onBrowseRadiusClicked()) );
					  
	QObject::connect( ui_compact.useTrueTransparencyBox, SIGNAL(stateChanged(int)), this, SLOT(onTrueTransparencyChanged(int)) );
	QObject::connect( ui_compact.DPMaxRendersSpin, SIGNAL(valueChanged(int)), this, SLOT(onDPMaxRendersSpun(int)) );
					  
	QObject::connect( ui_compact.burnBtn, SIGNAL(clicked()), this, SLOT(onBurnBtnClicked()) );
	QObject::connect( ui_compact.cleanTopoBtn, SIGNAL(clicked()), this, SLOT(onCleanTopoBtnClicked()) );
	QObject::connect( ui_compact.drawEdge, SIGNAL(stateChanged(int)), this, SLOT(onDrawEdgeChecked(int)) );
	QObject::connect( ui_compact.drawStPoints, SIGNAL(stateChanged(int)), this, SLOT(onDrawStPointsChecked(int)) );
	QObject::connect( ui_compact.pointSizeSpin, SIGNAL(valueChanged(double)), this, SLOT(onPointSizeSpun(double)) );
	QObject::connect( ui_compact.dualizeBtn, SIGNAL(clicked()), this, SLOT(onDualizeClicked()) );
	QObject::connect( ui_compact.hideOrig, SIGNAL(stateChanged(int)), this, SLOT(onHideOrigChanged(int)) );
	QObject::connect( ui_compact.hideMA, SIGNAL(stateChanged(int)), this, SLOT(onHideMAChanged(int)) );
	QObject::connect( ui_compact.hideMC, SIGNAL(stateChanged(int)), this, SLOT(onHideMCChanged(int)) );
	QObject::connect( ui_compact.visBurnt, SIGNAL(stateChanged(int)), this, SLOT(onVisBurntEdgesChanged(int)) );
	QObject::connect( ui_compact.hideHS, SIGNAL(stateChanged(int)), this, SLOT(onHideHSChanged(int)) );
	QObject::connect( ui_compact.hideIsoSurf, SIGNAL(stateChanged(int)), this, SLOT(onHideIsoSurfChanged(int)) );
	QObject::connect( ui_compact.hideIsoCont, SIGNAL(stateChanged(int)), this, SLOT(onHideIsoContourChanged(int)) );
					  
	QObject::connect( ui_compact.turnOffLightingForMABox, SIGNAL(stateChanged(int)), this, SLOT(onTurnOffLightingForMAChanged(int)));

	QObject::connect( ui_compact.computeAllDistBtn, SIGNAL(clicked()), this, SLOT(onComputeAllDistClicked()) );
	QObject::connect( ui_compact.enableFinePruneMA, SIGNAL(stateChanged(int)), this, SLOT(onEnableFinePrune(int)) );
	//QObject::connect(ui.colorPerFace, SIGNAL(stateChanged(int)), this, SLOT(onColorPerFaceChecked(int)));
	QObject::connect( ui_compact.pruneMADistCombo1, SIGNAL(currentIndexChanged(int)), this, SLOT(onPruneDistMAComboChanged_1(int)) );
	QObject::connect( ui_compact.pruneMASlider1, SIGNAL(valueChanged(int)), this, SLOT(onPruneMASlided_1(int)) );
	QObject::connect( ui_compact.pruneMASpin1, SIGNAL(valueChanged(double)), this, SLOT(onPruneMASpun_1(double)) );
	QObject::connect( ui_compact.pruneMADistCombo2, SIGNAL(currentIndexChanged(int)), this, SLOT(onPruneDistMAComboChanged_2(int)) );
	QObject::connect( ui_compact.pruneMASlider2, SIGNAL(valueChanged(int)), this, SLOT(onPruneMASlided_2(int)) );
	QObject::connect( ui_compact.pruneMASpin2, SIGNAL(valueChanged(double)), this, SLOT(onPruneMASpun_2(double)) );
					  
	QObject::connect( ui_compact.clampMADistBox, SIGNAL(stateChanged(int)), this, SLOT(onClampMADistBoxChanged(int)) );
	QObject::connect( ui_compact.minVisDistMASpin, SIGNAL(valueChanged(double)), this, SLOT(onMinVisDistMASpun(double)) );
	QObject::connect( ui_compact.maxVisDistMASpin, SIGNAL(valueChanged(double)), this, SLOT(onMaxVisDistMASpun(double)) );
					  
	QObject::connect( ui_compact.colorOrigBtn, SIGNAL(clicked()), this, SLOT(onOrigColorClicked()) );
	QObject::connect( ui_compact.colorBGBtn, SIGNAL(clicked()), this, SLOT(onChangeBGColorClicked()) );
	QObject::connect( ui_compact.colorMABtn, SIGNAL(clicked()), this, SLOT(onMAColorClicked()) );
	QObject::connect( ui_compact.colorMCBtn, SIGNAL(clicked()), this, SLOT(onMCColorClicked()) );
	QObject::connect( ui_compact.colorPerVert, SIGNAL(toggled(bool)), this, SLOT(onColorPerVertBtnClicked()) );
	//QObject::connect( ui_compact.colorPerVert, SIGNAL(clicked()), this, SLOT(onColorPerVertBtnClicked()) );
	QObject::connect( ui_compact.colorPerFace, SIGNAL(toggled(bool)), this, SLOT(onColorPerFaceBtnClicked()) );
	//QObject::connect( ui_compact.colorPerFace, SIGNAL(clicked()), this, SLOT(onColorPerFaceBtnClicked()) );
	QObject::connect( ui_compact.visMADistCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(onVisDistMAComboChanged(int)) );
	//QObject::connect( ui_compact.visMADistCombo, SIGNAL(activated(int)), this, SLOT(onVisDistMAComboChanged(int)) );
	QObject::connect( ui_compact.origTransparentSlider, SIGNAL(valueChanged(int)), this, SLOT(onOrigTransparentSlided(int)) );
	QObject::connect( ui_compact.MATransparentSlider, SIGNAL(valueChanged(int)), this, SLOT(onMATransparentSlided(int)) );
	QObject::connect( ui_compact.MATransExpSlider, SIGNAL(valueChanged(int)), this, SLOT(onMATransparentSlided(int)) );
	QObject::connect( ui_compact.MCAlphaSlider, SIGNAL(valueChanged(int)), this, SLOT(onMCTransparentSlided(int)) );
	QObject::connect( ui_compact.MCTransExpSlider, SIGNAL(valueChanged(int)), this, SLOT(onMCTransparentSlided(int)) );
	QObject::connect( ui_compact.MCWidthSpin, SIGNAL(valueChanged(double)), this, SLOT(onMCWidthSpinChanged(double)) );
					  
	QObject::connect( ui_compact.burnMedialCurveBtn, SIGNAL(clicked()), this, SLOT(onBurnMCClicked()));
	QObject::connect( ui_compact.visMCDistCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(onVisDistMCComboChanged(int)));
	//QObject::connect( ui_compact.visMCDistCombo, SIGNAL(activated(int)), this, SLOT(onVisDistMCComboChanged(int)));
	QObject::connect( ui_compact.clampMCDistBox, SIGNAL(stateChanged(int)), this, SLOT(onClampMCDistBoxChanged(int)));
	QObject::connect( ui_compact.minVisDistMCSpin, SIGNAL(valueChanged(double)), this, SLOT(onMinVisDistMCSpun(double)));
	QObject::connect( ui_compact.maxVisDistMCSpin, SIGNAL(valueChanged(double)), this, SLOT(onMaxVisDistMCSpun(double)));
	QObject::connect( ui_compact.pruneMCSlider1, SIGNAL(valueChanged(int)), this, SLOT(onPruneCurveSlided(int)));
	QObject::connect( ui_compact.pruneMCDistSpin1, SIGNAL(valueChanged(double)), this, SLOT(onPruneCurveSpinChanged(double)));
	QObject::connect( ui_compact.pruneMCDistCombo1, SIGNAL(currentIndexChanged(int)), this, SLOT(onPruneDistMCComboChanged(int)));
	QObject::connect( ui_compact.visBallStick, SIGNAL(stateChanged(int)), this, SLOT(onVisBallStickChecked(int)));
	QObject::connect( ui_compact.printMCStatsBtn, SIGNAL(clicked()), this, SLOT(onPrintStatsClicked()));
					  
	QObject::connect( ui_compact.createHSBtn, SIGNAL(clicked()), this, SLOT(onCreateHSClicked()));
	QObject::connect( ui_compact.visHSBtn, SIGNAL(clicked()), this, SLOT(onVisHSClicked()));
	QObject::connect( ui_compact.hsBt2Bt3Slider, SIGNAL(valueChanged(int)), this, SLOT(onHSTheta2Changed(int)) );
	QObject::connect( ui_compact.hsBt1Bt2Slider, SIGNAL(valueChanged(int)), this, SLOT(onHSTheta1Changed(int)) );
	QObject::connect( ui_compact.hsFaceDiffSpin, SIGNAL(valueChanged(double)), this, SLOT(onHSTheta2Spun(double)) );
	QObject::connect( ui_compact.hsCurveDiffSpin, SIGNAL(valueChanged(double)), this, SLOT(onHSTheta1Spun(double)) );
	QObject::connect( ui_compact.removeSmallCmpnts, SIGNAL(stateChanged(int)), this, SLOT(onHSRemoveSmallCmpntsChecked(int)));
					  
	QObject::connect( ui_compact.surfFuncCorrespSchemeCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(onSurfFuncCorrespSchemeComboChanged(int)));
	QObject::connect( ui_compact.surfFuncCorrespSchemeCombo, SIGNAL(activated(int)), this, SLOT(onSurfFuncCorrespSchemeComboChanged(int)));
	QObject::connect( ui_compact.surfFuncSetupCorrespBtn, SIGNAL(clicked()), this, SLOT(onSurfFuncCorrespSetupClicked()));
	QObject::connect( ui_compact.pruneSMetricCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(onSurfFuncComboChanged(int)));
	QObject::connect( ui_compact.pruneSMetricCombo, SIGNAL(activated(int)), this, SLOT(onSurfFuncComboChanged(int)));
	QObject::connect( ui_compact.smoothSMetricSlider, SIGNAL(valueChanged(int)), this, SLOT(onSurfFuncSlided(int)));
	QObject::connect( ui_compact.diffuseSurfFunc, SIGNAL(stateChanged(int)), this, SLOT(onDiffuseSurfFuncChecked(int)));
	QObject::connect( ui_compact.surfFuncProjectBtn, SIGNAL(clicked()), this, SLOT(onSurfFuncProjectClicked()));
					  
	QObject::connect( ui_compact.showMPBox, SIGNAL(stateChanged(int)), this, SLOT(onShowMPChecked(int)));
					  
	QObject::connect( ui_compact.isoSurfPrecomputeBtn, SIGNAL(clicked()), this, SLOT(onIsoSurfPrecomputeClicked()));
	QObject::connect( ui_compact.isoSurfRefineBtn, SIGNAL(clicked()), this, SLOT(onIsoSurfRefineClicked()));
	QObject::connect( ui_compact.isoSurfSlider, SIGNAL(valueChanged(int)), this, SLOT(onIsoSurfSlided(int)));
	QObject::connect( ui_compact.hideSnappedFaceBox, SIGNAL(stateChanged(int)), this, SLOT(onHideSnappedFacesChanged(int)));
					  
	QObject::connect( ui_compact.isoContPrecomputeBtn, SIGNAL(clicked()), this, SLOT(onIsoContPrecomputeClicked()));
	QObject::connect( ui_compact.isoContSlider, SIGNAL(valueChanged(int)), this, SLOT(onIsoContSlided(int)));
	QObject::connect( ui_compact.hidePastSurfBox, SIGNAL(stateChanged(int)), this, SLOT(onHidePastSurfChanged(int)));
	QObject::connect( ui_compact.isoContShowMC, SIGNAL(stateChanged(int)), this, SLOT(onIsoContShowMCChanged(int)));
					  
	QObject::connect( ui_compact.evolveAllPrecomputeBtn, SIGNAL(clicked()), this, SLOT(onEvolveAllPrecomputeClicked()));
	QObject::connect( ui_compact.evolveAllSlider, SIGNAL(valueChanged(int)), this, SLOT(onEvolveAllSlided(int)));

	QObject::connect( ui_compact.printBTBtn, SIGNAL(clicked()), this, SLOT(onPrintBTClicked()) );
	QObject::connect( ui_compact.showSphereBox, SIGNAL(stateChanged(int)), this, SLOT(onShowSpheresChecked(int)) );
	QObject::connect( ui_compact.visSphereSlider, SIGNAL(valueChanged(int)), this, SLOT(onVisSphereSlided(int)) );
	QObject::connect( ui_compact.visSphereSpin, SIGNAL(valueChanged(int)), this, SLOT(onVisSphereSpun(int)) );
	QObject::connect( ui_compact.printSelectVertInfoBtn, SIGNAL(clicked()), this, SLOT(onPrintSelectVertInfoClicked()) );
	QObject::connect( ui_compact.visSphereFromFileBtn, SIGNAL(clicked()), this, SLOT(onVisSpheresFromFileClicked()) );

	QObject::connect( ui_compact.linesRenderOffsetSlider, SIGNAL(valueChanged(int)), this, SLOT(onLinesRenderOffsetSlided(int)) );
	
	QObject::connect( ui_compact.outputWenpingBtn, SIGNAL(clicked()), this, SLOT(onOutputWenpingBtnClicked()) );
	QObject::connect( ui_compact.readQMatMA, SIGNAL(clicked()), this, SLOT(onReadQMATBtnClicked()) );
	QObject::connect( ui_compact.hideQMatBox, SIGNAL(stateChanged(int)), this, SLOT(onHideQMATChanged(int)) );

	QObject::connect( ui_compact.exportBTBtn, SIGNAL(clicked()), this, SLOT(onExportETBtnClicked()) );
	QObject::connect( ui_compact.exportHSBtn, SIGNAL(clicked()), this, SLOT(onExportHSBtnClicked()) );
	QObject::connect( ui_compact.exportBtn, SIGNAL(clicked()), this, SLOT(onExportBtnClicked()) );

	/*resetParams();*/
	
	cout << "ui setup!" << endl;

	//ui.splitter->setMinimumSize(QSize(800, 800));
	glarea = new GLArea(this);
	//glarea->reset();
	glarea->ui = &ui_compact; // In case glarea want to access ui. Maybe a bad idea.
	ui_compact.glLayout->addWidget(glarea, 1);
	glarea->updateGeometry();

	initUI();
	prepareForStep(OPEN_FILE);

	glarea->updateGL();
}

ETMainWindow::~ETMainWindow()
{
	if (glarea)
		delete glarea;
}

void ETMainWindow::initUI()
{
	// hide some pages of tab (containing details user don't need to know)
	vector<QWidget*> pages_to_remove;
	ui_compact.toolTabs->setCurrentIndex(0);
	pages_to_remove.push_back(ui_compact.SMetricTab);
	//pages_to_remove.push_back(ui_compact.VideoTab);
	pages_to_remove.push_back(ui_compact.MiscTab);
	pages_to_remove.push_back(ui_compact.hiddenTab);
	pages_to_remove.push_back(ui_compact.ExportTab);
	pages_to_remove.push_back(ui_compact.MCTab);
	for (auto it = pages_to_remove.begin(); it != pages_to_remove.end(); ++it)
	{
		auto index = ui_compact.toolTabs->indexOf(*it);
		assert(index >= 0);
		ui_compact.toolTabs->removeTab(
			index
			);
	}
}

void ETMainWindow::setDebugMode(bool _enabled)
{
	m_debugMode = _enabled;
}

void ETMainWindow::setConsoleMode( bool _enabled )
{
	m_consoleMode = _enabled;
}

void ETMainWindow::setInputs(
	const std::string& _ma_file, const std::string& _shape_file, const std::string& _r_file,
	float _omega,
	const std::string& _mc_meas,
	double _theta_2, double _theta_1 )
{
	m_ma_file = _ma_file;
	m_shape_file = _shape_file;
	m_r_file = _r_file;
	m_omega = _omega;
	m_mc_meas_to_output = _mc_meas;
	m_skel_theta2 = _theta_2;
	m_skel_theta1 = _theta_1;
}

void ETMainWindow::onTrueTransparencyChanged(int _s)
{
	glarea->setTrueTransparencyRender(_s == Qt::Checked);
	ui_compact.DPMaxRendersSpin->setEnabled(_s == Qt::Checked);
	glarea->setDPMaxRenders(ui_compact.DPMaxRendersSpin->value());
	glarea->updateGL();
}

void ETMainWindow::onDPMaxRendersSpun(int _v)
{
	int v = std::max(_v, 1);
	ui_compact.DPMaxRendersSpin->blockSignals(true);
	ui_compact.DPMaxRendersSpin->setValue(v);
	ui_compact.DPMaxRendersSpin->blockSignals(false);

	glarea->setDPMaxRenders(v);
	glarea->updateGL();
}

void ETMainWindow::onBurnBtnClicked()
{
	ui_compact.statusbar->showMessage("burning medial axis...");
	bool success = glarea->burn( 
		(GLArea::SubdivScheme)ui_compact.stSubdivCombo->currentIndex(), 
		m_consoleMode ? m_omega : ui_compact.nFixedSteinerSpin->value(),
		ui_compact.edgeWeightCombo->currentIndex() );
	if (success)
	{
		prepareForStep(PRUNE_MA);
		prepareForStep(DUALIZE);
		ui_compact.statusbar->showMessage("Done: burning medial axis.");
		/*ui.maMeasureGroup->setEnabled(true);
		ui.mcDualizeGroup->setEnabled(true);
		ui.computeAllDistBtn->setEnabled(true);*/
	}
}

void ETMainWindow::onCleanTopoBtnClicked()
{
	bool has_unburnt = false;
	glarea->cleanTopo( std::string( "" ), has_unburnt );
	if ( has_unburnt ) // need to burn again
		this->onBurnBtnClicked();
}


void ETMainWindow::onSaveView()
{
	/*auto sizes = this->size();
	this->resize(sizes.width() * 0.8f, sizes.height() * 0.8f);*/
	QString ini_file_name = QFileDialog::getSaveFileName(
		this, "Save current view to a .ini file", "", "", NULL, QFileDialog::DontUseNativeDialog
		);
	if ( ini_file_name.isNull() )
	{
		QMessageBox::warning(this, "Warning: no file specified.", "Current view was not saved.");
	}
	else // use the specified file as the .ini file for saving settings
	{
		typedef std::shared_ptr<QSettings> setting_ptr;
		setting_ptr view_settings = setting_ptr(
			new QSettings(ini_file_name, QSettings::IniFormat)
			);
		view_settings->setValue("QtPureGL/sizes", this->size());

		// propagate down the saving command
		glarea->saveView(view_settings);

		QMessageBox::information(this, "success", "Current view saved.");
	}
}

void ETMainWindow::onLoadView()
{
	/*auto sizes = this->size();
	this->resize(sizes.width() * 0.8f, sizes.height() * 0.8f);*/
	QString ini_file_name = QFileDialog::getOpenFileName(
		this, "Open a .ini file to load view", "", "", NULL, QFileDialog::DontUseNativeDialog
		);

	if ( ini_file_name.isNull() )
	{
		QMessageBox::warning(this, "Warning: no file specified.", "View not loaded.");
	}
	else // use the specified file as the .ini file for loading settings
	{
		typedef std::shared_ptr<QSettings> setting_ptr;
		setting_ptr view_settings( new QSettings(ini_file_name, QSettings::IniFormat) );

		// propagate down the saving command
		glarea->loadView(view_settings);

		auto window_sizes = view_settings->value("QtPureGL/sizes").toSize();
		this->resize(window_sizes.width(), window_sizes.height());

		//QMessageBox::information(this, "success", "Current view saved.");
	}
}

void ETMainWindow::onLoadFilesClicked()
{
	assert(glarea);
	// load MA mesh, orig 3d mesh, and radius files
	glarea->reset();
	//glarea->passToDrawable(mesh_file.toStdString(), origMesh_file.toStdString());

	std::string ma_file, origshape_file, r_file;
	if ( m_consoleMode )
	{
		ma_file = m_ma_file;
		origshape_file = m_shape_file;
		r_file = m_r_file;
	}
	else
	{
		ma_file = ui_compact.maFileEdit->text().toStdString();
		origshape_file = ui_compact.shapeFileEdit->text().toStdString();
		r_file = ui_compact.radiusFileEdit->text().toStdString();
	}

	if (ma_file == "")
	{
		cout << "onOpenFile(): no valid medial axis file specified. "<< endl;
		ui_compact.statusbar->showMessage("no valid medial axis file specified.");
		return;
	}

	ui_compact.statusbar->showMessage("Loading original shape and its medial axis...");
	prepareForStep(OPEN_FILE);
	// load MA as the 3d shape if the user somehow doesn't specify the 3d shape
	if (origshape_file == "")
		origshape_file = ma_file;
	glarea->reset();
	glarea->passToDrawable(ma_file, origshape_file, r_file);

	prepareForStep(POST_OPEN_FILE);
	prepareForStep(BURN_MA);
	
	// save the states in case other MA related geometry needs this color (i.e. finer MA)
	TriColor color(1.0f, 0.2f, 0.2f);
	glarea->setUseConstColorMA( true, color );
	glarea->changeConstColorMA(color, ui_compact.colorPerVert->isChecked());

	ui_compact.statusbar->showMessage("Files loaded.");
}

void ETMainWindow::onBrowseMAClicked()
{
	QString ma_file = QFileDialog::getOpenFileName( this, "Select medial axis file", "", "*.off *.ply",
		NULL/*, QFileDialog::DontUseNativeDialog*/ );
	ui_compact.maFileEdit->setText( ma_file );
}
void ETMainWindow::onBrowseOrigClicked()
{
	QString shape_file = QFileDialog::getOpenFileName( this, "Select original 3d shape file", "", "*.off",
		NULL/*, QFileDialog::DontUseNativeDialog*/ );
	ui_compact.shapeFileEdit->setText( shape_file );
}
void ETMainWindow::onBrowseRadiusClicked()
{
	QString radius_file = QFileDialog::getOpenFileName( this, "Select radius file", "", "*.r",
		NULL/*, QFileDialog::DontUseNativeDialog*/ );
	ui_compact.radiusFileEdit->setText( radius_file );
}

void ETMainWindow::onDrawEdgeChecked(int _state)
{
	glarea->wireFrameOnMesh(_state == Qt::Checked);
}

void ETMainWindow::onDrawStPointsChecked(int _state)
{
	ui_compact.pointSizeSpin->setEnabled(_state == Qt::Checked);
	glarea->drawSteinerPoints(_state == Qt::Checked, (float)ui_compact.pointSizeSpin->value());
}

void ETMainWindow::onPointSizeSpun(double _value)
{
	//cout << "slider: " << _value << endl;
	glarea->drawSteinerPoints(true, _value);
}

void ETMainWindow::onDualizeClicked()
{
	int type = 0;
	/*type |= ui.visBurnt->isChecked() ? GLArea::VIS_BURNT : 0;
	type |= ui.visDual->isChecked() ? GLArea::VIS_DUAL : 0;*/
	bool only_unburnt = ui_compact.onlyUnburntBox->isChecked();
	int edge_dual_method_idx = ui_compact.edgeDualOptCombo->currentIndex();
	int poly_dual_method_idx = ui_compact.polyDualOptCombo->currentIndex();
	bool stop_burn_at_junc = ui_compact.stopAtJuncBox->isChecked();

	ui_compact.statusbar->showMessage("computing medial curve...");

#ifdef PRINT_MEM_USAGE
	PROCESS_MEMORY_COUNTERS pmc;
	GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc));
	SIZE_T memPreUsed, memCurUsed;
	memPreUsed = memCurUsed = pmc.WorkingSetSize;
#endif // PRINT_MEM_USAGE

	bool success = glarea->obtainDualLineStructure( 
		only_unburnt, 
		edge_dual_method_idx, poly_dual_method_idx, stop_burn_at_junc);
	cout << "obtainDualLineStructure() finished." << endl << endl;

#ifdef PRINT_MEM_USAGE
	GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc));
	memCurUsed = pmc.WorkingSetSize;
	cout << "mem. usage: " << memCurUsed / (1024*1024) << endl;
	cout << "mem. usage (by dualization() ) in MB: " << (memCurUsed - memPreUsed) / (1024*1024) << endl;
	memPreUsed = memCurUsed;
#endif // PRINT_MEM_USAGE

	ui_compact.statusbar->showMessage("Done: computing medial curve.");

	if (success)
	{
		prepareForStep(BURN_MC);
		if (!m_debugMode)
		{
			ui_compact.burnMedialCurveBtn->click();
		}
		/*enablePruneMC();
		enableSurfFuncWidgets();*/
	}
}

void ETMainWindow::onBurnMCClicked()
{
	bool only_unburnt = ui_compact.onlyUnburntBox->isChecked();
	bool stop_burn_at_junc = ui_compact.stopAtJuncBox->isChecked();
	bool protect_bt2 = ui_compact.protectBT2Box->isChecked();
	bool underestimate_dist = /*false;*/ ui_compact.underEstimateBox->isChecked();

	ui_compact.statusbar->showMessage("Burning the medial curve...");
	bool success = 
		glarea->burnCurveNetwork(stop_burn_at_junc, only_unburnt, protect_bt2, underestimate_dist);

	if (success)
	{
		cout << "Done: burning the medial curve." << endl << endl;
		prepareForStep(POST_BURN_MC);
		prepareForStep(COMPUTE_HS);
		ui_compact.statusbar->showMessage("Done: burning the medial curve.");
		//enablePruneMC();
		handle_extra_params();
	}
}

void ETMainWindow::onVisDistMCComboChanged(int _idx)
{
	bool clamp = ui_compact.clampMCDistBox->isChecked();
	float visDistOnMC_min, visDistOnMC_max;
	if (clamp)
	{
		visDistOnMC_min = ui_compact.minVisDistMCSpin->value();
		visDistOnMC_max = ui_compact.maxVisDistMCSpin->value();
	}

	GLArea::DistMC field_type = get_MC_dist_type(_idx);
	glarea->colorMCEdgeBy(field_type, visDistOnMC_min, visDistOnMC_max, clamp, 0.0f, 1.0f, 1);

	// update range of the cur distance metric
	if (!clamp)
	{
		ui_compact.minVisDistMCSpin->setValue(visDistOnMC_min);
		ui_compact.maxVisDistMCSpin->setValue(visDistOnMC_max);
	}
}

void ETMainWindow::onPruneCurveSlided(int _value)
{
	// convert to abs value
	double range = distMax_MC - distMin_MC;
	double t = distMin_MC + range * ui_compact.pruneMCSlider1->value() / 
		(ui_compact.pruneMCSlider1->maximum() - ui_compact.pruneMCSlider1->minimum());

	ui_compact.mcPruneRatio->setValue( glarea->toETcPruneRatio(t) );

	bool success = glarea->pruneAndVisMedialCurve(t, ui_compact.preserveMCTopo->isChecked(), false);
	if ( !success )
	{
		cout << "ERROR: failed to prune medial curves." << endl;
	}

	// update digit spin
	ui_compact.pruneMCDistSpin1->blockSignals(true);
	ui_compact.pruneMCDistSpin1->setValue(t);
	ui_compact.pruneMCDistSpin1->blockSignals(false);
}

void ETMainWindow::onPruneCurveSpinChanged(double _val)
{
	_val = std::min(distMax_MC,
		std::max(distMin_MC, (float)_val));
	bool success = glarea->pruneAndVisMedialCurve( _val, ui_compact.preserveMCTopo->isChecked(), false );
	if ( !success )
	{
		cout << "ERROR: failed to prune medial curves." << endl;
	}

	// update the slider
	ui_compact.pruneMCSlider1->blockSignals( true );
	float range = distMax_MC - distMin_MC;
	ui_compact.pruneMCSlider1->setValue( ui_compact.pruneMCSlider1->minimum() +
		(ui_compact.pruneMCSlider1->maximum() - ui_compact.pruneMCSlider1->minimum()) * 
		(_val - distMin_MC) / range );
	ui_compact.pruneMCSlider1->blockSignals( false );

}

bool ETMainWindow::onPruneDistMCComboChanged(int _idx)
{
	GLArea::DistMC metric;
	switch (_idx)
	{
	case 0:
		metric = GLArea::BT2_MC;
		break;
	case 1:
		metric = GLArea::BT1_MC;
		break;
	case 2:
		metric = GLArea::BT1_BT2_MC;
		break;
	case 3:
		metric = GLArea::BT1_BT2_REL_MC;
		break;
	}
	if ( !glarea->setUpMCDistMetricForPruning(metric, distMin_MC, distMax_MC) )
	{
		cout << "ERROR: failed to get current distance metric to be used in pruning." << endl;
		return false;
	}

	/*int line_type = 0;
	line_type |= ui.visBurnt->isChecked() ? GLArea::VIS_BURNT : 0;
	line_type |= ui.visDual->isChecked() ? GLArea::VIS_DUAL : 0;*/

	float dummy = 0.0f;
	if ( glarea->colorMCEdgeBy(metric, dummy, dummy, false, 0.0f, 1.0f, 1) )
	{
		ui_compact.pruneMCSlider1->setEnabled(true);
		return true;
	}
	else 
		return false;
}

void ETMainWindow::onVisBallStickChecked(int _state)
{
	if ( _state == Qt::Checked )
	{
		if ( glarea->uploadRemainedMCtoProxyDrawer() )
		{
			ui_compact.pruneMCSlider1->setEnabled(false);
			ui_compact.pruneMCDistSpin1->setEnabled(false);
			glarea->drawMCProxy(true);
		}
		else
		{
			glarea->drawMCProxy(false);
		}
	}
	else
	{
		ui_compact.pruneMCSlider1->setEnabled(true);
		ui_compact.pruneMCDistSpin1->setEnabled(true);

		glarea->drawMCProxy(false);
	}
}

void ETMainWindow::onPrintStatsClicked()
{
	glarea->printRemainedMCStats();
}

void ETMainWindow::onHideMAChanged(int _state)
{
	glarea->setDrawFlag( GLArea::DRAW_MA, _state != Qt::Checked );
}

void ETMainWindow::onHideOrigChanged(int _state)
{
	glarea->setDrawFlag( GLArea::DRAW_ORIG, _state != Qt::Checked );
}

void ETMainWindow::onHideMCChanged(int _state)
{
	glarea->setDrawFlag( GLArea::DRAW_MC, _state != Qt::Checked );
}

void ETMainWindow::onVisBurntEdgesChanged(int _state)
{
	glarea->setDrawFlag(GLArea::DRAW_BURNT_EDGES, _state == Qt::Checked);
}

void ETMainWindow::onHideHSChanged(int _state)
{
	glarea->setDrawFlag( GLArea::DRAW_HS, _state != Qt::Checked );
}

void ETMainWindow::onHideIsoSurfChanged(int _state)
{
	glarea->setDrawFlag( GLArea::DRAW_ISOSURF, _state != Qt::Checked);
}

void ETMainWindow::onHideIsoContourChanged(int _state)
{
	glarea->setDrawFlag( GLArea::DRAW_ISOCONT, _state != Qt::Checked );
}

void ETMainWindow::onTurnOffLightingForMAChanged(int _state)
{
	glarea->turnOffLighting(GLArea::MA, _state == Qt::Checked);
	glarea->turnOffLighting(GLArea::FINE_MA_FOR_ISO_CONTOUR, _state == Qt::Checked);
	glarea->turnOffLighting(GLArea::FINE_MA_FOR_PRUNING, _state == Qt::Checked);
}

void ETMainWindow::onMAColorClicked()
{
	const QString COLOR_STYLE("QPushButton { background-color : %1; color : %2; }");

	QColor chosen_color = QColorDialog::getColor();
	QColor text_color = getIdealTextColor(chosen_color);
	ui_compact.colorMABtn->setStyleSheet(COLOR_STYLE.arg(chosen_color.name()).arg(text_color.name()));

	float color[3] = { 
		(float)chosen_color.redF(), 
		(float)chosen_color.greenF(), 
		(float)chosen_color.blueF() 
	};
	// save the states in case other MA related geometry needs this color (i.e. finer MA)
	glarea->setUseConstColorMA( true, TriColor(color) );
	
	glarea->changeConstColorMA(color, ui_compact.colorPerVert->isChecked());
}

void ETMainWindow::onMCColorClicked()
{
	QColor chosen_color = QColorDialog::getColor();
	QColor text_color = getIdealTextColor(chosen_color);

	const QString COLOR_STYLE("QPushButton { background-color : %1; color : %2; }");
	ui_compact.colorMCBtn->setStyleSheet(COLOR_STYLE.arg(chosen_color.name()).arg(text_color.name()));

	float color[3] = { 
		(float)chosen_color.redF(), 
		(float)chosen_color.greenF(), 
		(float)chosen_color.blueF() };
	glarea->changeConstColorMC(color);
}

void ETMainWindow::onOrigColorClicked()
{
	const QString COLOR_STYLE("QPushButton { background-color : %1; color : %2; }");

	QColor chosen_color = QColorDialog::getColor();
	QColor text_color = getIdealTextColor(chosen_color);
	ui_compact.colorOrigBtn->setStyleSheet(COLOR_STYLE.arg(chosen_color.name()).arg(text_color.name()));

	float color[3] = { 
		(float)chosen_color.redF(), 
		(float)chosen_color.greenF(), 
		(float)chosen_color.blueF() };
		glarea->changeOrigColor(color);
}

void ETMainWindow::onChangeBGColorClicked()
{
	const QString COLOR_STYLE("QPushButton { background-color : %1; color : %2; }");

	QColor chosen_color = QColorDialog::getColor();
	QColor text_color = getIdealTextColor(chosen_color);
	ui_compact.colorBGBtn->setStyleSheet(COLOR_STYLE.arg(chosen_color.name()).arg(text_color.name()));

	float color[3] = {
		(float)chosen_color.redF(),
		(float)chosen_color.greenF(),
		(float)chosen_color.blueF()
	};
	glarea->changeBGColor(color);
}

void ETMainWindow::onComputeAllDistClicked()
{
	cout << "precomputing for various measures..."<<endl;
	glarea->precomputeForMeasures();
	cout << "precomputation done."<<endl;
	ui_compact.pruneMASlider2->setEnabled(true);
	ui_compact.pruneMASpin2->setEnabled(true);
}

void ETMainWindow::onEnableFinePrune(int _state)
{
	glarea->setDrawFlag(GLArea::DRAW_MA_FINE, _state == Qt::Checked);
	if (_state == Qt::Checked)
	{
		// hide coarse ma
		ui_compact.hideMA->setChecked(true);
		// force prune ma slider 1 to "move"
		auto cur_value = ui_compact.pruneMASlider1->value();
		ui_compact.pruneMASlider1->blockSignals(true);
		ui_compact.pruneMASlider1->setValue(cur_value == ui_compact.pruneMASlider1->maximum() ? cur_value - 1 : cur_value + 1);
		ui_compact.pruneMASlider1->blockSignals(false);
		ui_compact.pruneMASlider1->setValue(cur_value);
	}
}

void ETMainWindow::onPruneMASlided_1(int _value)
{
	//cout << "prune MA slider changed!" << endl;

	float val1, val2;
	getPruneMAThreshold(val1, val2);

	if ( ui_compact.enableFinePruneMA->isChecked() )
		glarea->pruneFineMA(val1, val2);
	else
		glarea->pruneMA(val1, val2);
}

void ETMainWindow::onPruneMASpun_1(double _value)
{
	float range = pruneDistOnMA_max_1 - pruneDistOnMA_min_1;

	_value = std::min(pruneDistOnMA_max_1,
		std::max(pruneDistOnMA_min_1, (float)_value));
	ui_compact.pruneMASlider1->setValue( ui_compact.pruneMASlider1->minimum() + 
		(ui_compact.pruneMASlider1->maximum() - ui_compact.pruneMASlider1->minimum()) * 
		(_value - pruneDistOnMA_min_1) / range );
}

void ETMainWindow::onPruneMASlided_2(int _value)
{
	float value_1, value_2;

	getPruneMAThreshold(value_1, value_2);

	glarea->pruneMA(value_1, value_2);

}

void ETMainWindow::onPruneMASpun_2(double _value)
{
	float range = pruneDistOnMA_max_2 - pruneDistOnMA_min_2;

	_value = std::min(pruneDistOnMA_max_2,
		std::max(pruneDistOnMA_min_2, (float)_value));
	ui_compact.pruneMASlider2->setValue( ui_compact.pruneMASlider2->minimum() + 
		(ui_compact.pruneMASlider2->maximum() - ui_compact.pruneMASlider2->minimum()) * 
		(_value - pruneDistOnMA_min_2) / range );
}

void ETMainWindow::onOrigTransparentSlided(int _value)
{
	glarea->changeOrigTransparency( 
		((float)_value - ui_compact.origTransparentSlider->minimum()) / 
		(ui_compact.origTransparentSlider->maximum() - ui_compact.origTransparentSlider->minimum()) 
		);
}

void ETMainWindow::onMATransparentSlided(int _value)
{
	if (
		(ui_compact.colorPerVert->isChecked() || ui_compact.colorPerFace->isChecked()) && 
		ui_compact.visMADistCombo->currentIndex() >= 0
		)
	{
		float min_alpha = ((float)ui_compact.MATransparentSlider->value() - ui_compact.MATransparentSlider->minimum()) / 
			(ui_compact.MATransparentSlider->maximum() - ui_compact.MATransparentSlider->minimum());
		int exp = ui_compact.MATransExpSlider->value();
		if (ui_compact.colorPerVert->isChecked())
		{
			glarea->changeMAVertTransparency( 
				this->get_MA_vert_dist_type(ui_compact.visMADistCombo->currentIndex()), 
				min_alpha, exp
				);
		}
		else if (ui_compact.colorPerFace->isChecked())
		{
			glarea->changeMAFaceTransparency( 
				this->get_MA_face_dist_type(ui_compact.visMADistCombo->currentIndex()), 
				min_alpha, exp
				);
		}
	}
}

void ETMainWindow::onMCTransparentSlided(int _value)
{
	float min_alpha = ((float)ui_compact.MCAlphaSlider->value() - ui_compact.MCAlphaSlider->minimum()) / 
		(ui_compact.MCAlphaSlider->maximum() - ui_compact.MCAlphaSlider->minimum()) ;
	int exp = ui_compact.MCTransExpSlider->value();
	glarea->changeMCTransparency( 
		get_MC_dist_type(ui_compact.visMCDistCombo->currentIndex()), 
		min_alpha, exp
		);
}

void ETMainWindow::onIsoContTransparentSlided(int _value)
{
	float min_alpha = ((float)_value - ui_compact.MCAlphaSlider->minimum()) / 
		(ui_compact.MCAlphaSlider->maximum() - ui_compact.MCAlphaSlider->minimum()) ;
}

void ETMainWindow::onMCWidthSpinChanged(double _value)
{
	glarea->changeMCLineWidth(_value);
}

void ETMainWindow::onClampMADistBoxChanged(int _state)
{
	if (_state == Qt::Checked)
	{
		ui_compact.maxVisDistMASpin->setEnabled(true);
		ui_compact.maxVisDistMASpin->blockSignals(false);
		ui_compact.minVisDistMASpin->setEnabled(true);
		ui_compact.minVisDistMASpin->blockSignals(false);
	}
	else
	{
		ui_compact.maxVisDistMASpin->setEnabled(false);
		ui_compact.maxVisDistMASpin->blockSignals(true);
		ui_compact.minVisDistMASpin->setEnabled(false);
		ui_compact.minVisDistMASpin->blockSignals(true);

		ui_compact.visMADistCombo->setCurrentIndex(ui_compact.visMADistCombo->currentIndex());
	}
}

void ETMainWindow::onMaxVisDistMASpun(double _value)
{
	if (!ui_compact.maxVisDistMASpin->isEnabled())
		return;

	float min_val = (float)ui_compact.minVisDistMASpin->value();
	float max_val = (float)ui_compact.maxVisDistMASpin->value();

	GLArea::FaceFieldType face_field;
	GLArea::VertFieldType vert_field;
	ui_compact.colorPerFace->isChecked() 
		? (
		face_field = get_MA_face_dist_type(ui_compact.visMADistCombo->currentIndex()),
		glarea->colorMAFaceBy(
		face_field, 
		min_val, max_val, 
		ui_compact.clampMADistBox->isChecked(), 0.0f, 0.0f, 1, 
		ui_compact.doMAFaceScalarDiffusion->isChecked(),
		ui_compact.usePerSheetBox->isChecked() ) 
		)
		: (
		vert_field = get_MA_vert_dist_type(ui_compact.visMADistCombo->currentIndex()),
		glarea->colorMAVertBy(
		vert_field, 
		min_val, max_val, 
		ui_compact.clampMADistBox->isChecked(), 0.0f, 0.0f, 1)
		);
}

void ETMainWindow::onMinVisDistMASpun(double _value)
{
	if (!ui_compact.minVisDistMASpin->isEnabled())
		return;

	float min_val = (float)ui_compact.minVisDistMASpin->value();
	float max_val = (float)ui_compact.maxVisDistMASpin->value();

	GLArea::FaceFieldType face_field;
	GLArea::VertFieldType vert_field;
	ui_compact.colorPerFace->isChecked() 
		? (
		face_field = get_MA_face_dist_type(ui_compact.visMADistCombo->currentIndex()),
		glarea->colorMAFaceBy(
		face_field, 
		min_val, max_val, 
		ui_compact.clampMADistBox->isChecked(), 
		0.0f, 0.0f, 1, 
		ui_compact.doMAFaceScalarDiffusion->isChecked(),
		ui_compact.usePerSheetBox->isChecked() ) 
		)
		: (
		vert_field = get_MA_vert_dist_type(ui_compact.visMADistCombo->currentIndex()),
		glarea->colorMAVertBy(
		vert_field, 
		min_val, max_val, 
		ui_compact.clampMADistBox->isChecked(),
		0.0f, 0.0f, 1)
		);
}

void ETMainWindow::onClampMCDistBoxChanged(int _state)
{
	if (_state == Qt::Checked)
	{
		ui_compact.maxVisDistMCSpin->setEnabled(true);
		ui_compact.maxVisDistMCSpin->blockSignals(false);
		ui_compact.minVisDistMCSpin->setEnabled(true);
		ui_compact.minVisDistMCSpin->blockSignals(false);
	}
	else
	{
		ui_compact.maxVisDistMCSpin->setEnabled(false);
		ui_compact.maxVisDistMCSpin->blockSignals(true);
		ui_compact.minVisDistMCSpin->setEnabled(false);
		ui_compact.minVisDistMCSpin->blockSignals(true);

		ui_compact.visMCDistCombo->setCurrentIndex(ui_compact.visMCDistCombo->currentIndex());
	}
}

void ETMainWindow::onMaxVisDistMCSpun(double _value)
{
	if (!ui_compact.maxVisDistMCSpin->isEnabled())
		return;

	float min_val = (float)ui_compact.minVisDistMCSpin->value();
	float max_val = (float)ui_compact.maxVisDistMCSpin->value();

	GLArea::DistMC MC_field;
	MC_field = get_MC_dist_type(ui_compact.visMCDistCombo->currentIndex());
	glarea->colorMCEdgeBy(
		MC_field, 
		min_val, max_val, 
		ui_compact.clampMCDistBox->isChecked(), 0.0f, 0.0f, 1);
}

void ETMainWindow::onMinVisDistMCSpun(double _value)
{
	if (!ui_compact.minVisDistMCSpin->isEnabled())
		return;

	float min_val = (float)ui_compact.minVisDistMCSpin->value();
	float max_val = (float)ui_compact.maxVisDistMCSpin->value();

	GLArea::DistMC MC_field;
	MC_field = get_MC_dist_type(ui_compact.visMCDistCombo->currentIndex());
	glarea->colorMCEdgeBy(
		MC_field, 
		min_val, max_val, 
		ui_compact.clampMCDistBox->isChecked(), 0.0f, 0.0f, 1);
}

void ETMainWindow::onColorPerVertBtnClicked()
{
	// load the items into combo 
	ui_compact.visMADistCombo->blockSignals(true);

	ui_compact.visMADistCombo->clear();
	ui_compact.visMADistCombo->addItem("BT_M");
	ui_compact.visMADistCombo->addItem("R");
	ui_compact.visMADistCombo->addItem("ET");
	//ui_compact.visMADistCombo->addItem("ET(rel)");
	ui_compact.visMADistCombo->setCurrentIndex(-1);

	ui_compact.visMADistCombo->blockSignals(false);
}

void ETMainWindow::onColorPerFaceBtnClicked()
{
	// load the items into combo 
	ui_compact.visMADistCombo->blockSignals(true);

	ui_compact.visMADistCombo->clear();
	ui_compact.visMADistCombo->addItem("BT_M");
	ui_compact.visMADistCombo->addItem("R");
	ui_compact.visMADistCombo->addItem("ET");
	ui_compact.visMADistCombo->addItem("ET(rel)");
	ui_compact.visMADistCombo->addItem("Angle");
	ui_compact.visMADistCombo->addItem("Lambda");
	ui_compact.visMADistCombo->addItem("MGF");
	ui_compact.visMADistCombo->setCurrentIndex(-1);

	ui_compact.visMADistCombo->setItemData(
		0, 
		"burn time on MA", 
		Qt::ToolTipRole);
	ui_compact.visMADistCombo->setItemData(
		1, 
		"radius field", 
		Qt::ToolTipRole);
	ui_compact.visMADistCombo->setItemData(
		2, 
		"erosion thickness", 
		Qt::ToolTipRole);
	ui_compact.visMADistCombo->setItemData(
		3, 
		"relative erosion thickness", 
		Qt::ToolTipRole);
	ui_compact.visMADistCombo->setItemData(
		4, 
		"angle btw. vectors formed by center of face & 3-sphere-intersections", 
		Qt::ToolTipRole);
	ui_compact.visMADistCombo->setItemData(
		5, 
		"circumradius", 
		Qt::ToolTipRole);
	ui_compact.visMADistCombo->setItemData(
		6, 
		"medial geodesic function", 
		Qt::ToolTipRole);

	ui_compact.visMADistCombo->blockSignals(false);
}

void ETMainWindow::onVisDistMAComboChanged(int _idx)
{
	//cout << "current index: " << _idx << endl; //debug
	if (_idx < 0)
		return;

	bool clamp = ui_compact.clampMADistBox->isChecked();
	float visDistOnMA_min, visDistOnMA_max;
	if (clamp)
	{
		visDistOnMA_min = ui_compact.minVisDistMASpin->value();
		visDistOnMA_max = ui_compact.maxVisDistMASpin->value();
	}

	// grab transparency params
	float min_alpha = ((float)ui_compact.MATransparentSlider->value() - ui_compact.MATransparentSlider->minimum()) / 
		(ui_compact.MATransparentSlider->maximum() - ui_compact.MATransparentSlider->minimum());
	int exp = ui_compact.MATransExpSlider->value();
	min_alpha = exp == 0 ? 1.0f : min_alpha;
	int update_option = 3;

	glarea->setUseConstColorMA( false, TriColor()/*dummy value*/ );
	if (ui_compact.colorPerFace->isChecked())
	{
		// process per-face distance
		glarea->usePerFaceRender(true);
		GLArea::FaceFieldType face_field = get_MA_face_dist_type(_idx);
		glarea->colorMAFaceBy(
			face_field, 
			visDistOnMA_min, visDistOnMA_max, clamp, 
			min_alpha, exp, 
			update_option, 
			ui_compact.doMAFaceScalarDiffusion->isChecked(),
			ui_compact.usePerSheetBox->isChecked());
	}
	else
	{
		// process per-vert distance
		glarea->usePerFaceRender(false);
		GLArea::VertFieldType vert_field = get_MA_vert_dist_type(_idx);
		glarea->colorMAVertBy(vert_field, visDistOnMA_min, visDistOnMA_max, clamp, min_alpha, exp, update_option);
	}

	// update range of the cur distance metric
	if (!clamp)
	{
		ui_compact.minVisDistMASpin->setValue(visDistOnMA_min);
		ui_compact.maxVisDistMASpin->setValue(visDistOnMA_max);
	}
}

void ETMainWindow::onPruneDistMAComboChanged_1(int _idx)
{
	// only prune based on distance per face
	// prepare the distance metric used for pruning
	GLArea::FaceFieldType face_field_1 = get_MA_face_dist_type(_idx);
	glarea->setupMAFacePrune(face_field_1, pruneDistOnMA_min_1, pruneDistOnMA_max_1, 1);

	// update range of the cur distance metric
	// but don't propagate the event
	ui_compact.minPruneDistMASpin1->blockSignals(true);
	ui_compact.minPruneDistMASpin1->setValue(pruneDistOnMA_min_1);
	ui_compact.maxPruneDistMASpin1->setValue(pruneDistOnMA_max_1);
	ui_compact.minPruneDistMASpin1->blockSignals(false);

	ui_compact.pruneMASlider1->setEnabled(true);
	ui_compact.pruneMASpin1->setEnabled(true);
	ui_compact.enableFinePruneMA->setEnabled(true);

	// update current pruning result
	// by forcing slider to *move*
	auto cur_value = ui_compact.pruneMASlider1->value();
	ui_compact.pruneMASlider1->blockSignals(true);
	ui_compact.pruneMASlider1->setValue(cur_value == ui_compact.pruneMASlider1->maximum() ? cur_value - 1 : cur_value + 1);
	ui_compact.pruneMASlider1->blockSignals(false);
	ui_compact.pruneMASlider1->setValue(cur_value);

}

void ETMainWindow::onPruneDistMAComboChanged_2(int _idx)
{
	// only prune based on distance per face
	// prepare the distance metric used for pruning
	GLArea::FaceFieldType face_field = get_MA_face_dist_type(_idx);
	glarea->setupMAFacePrune(face_field, pruneDistOnMA_min_2, pruneDistOnMA_max_2, 2);

	// update range of the cur distance metric
	// but don't propagate the event
	ui_compact.minPruneDistMASpin2->blockSignals(true);
	ui_compact.minPruneDistMASpin2->setValue(pruneDistOnMA_min_2);
	ui_compact.maxPruneDistMASpin2->setValue(pruneDistOnMA_max_2);
	ui_compact.minPruneDistMASpin2->blockSignals(false);

	ui_compact.pruneMASlider2->setEnabled(true);
	ui_compact.pruneMASpin2->setEnabled(true);

	// update current pruning result
	// by forcing slider to *move*
	auto cur_value = ui_compact.pruneMASlider2->value();
	ui_compact.pruneMASlider2->blockSignals(true);
	ui_compact.pruneMASlider2->setValue(cur_value == ui_compact.pruneMASlider2->maximum() ? cur_value - 1 : cur_value + 1);
	ui_compact.pruneMASlider2->blockSignals(false);
	ui_compact.pruneMASlider2->setValue(cur_value);

}

/////////////////////////////////////////////////////
///////////////// private helpers ///////////////////
/////////////////////////////////////////////////////

void ETMainWindow::prepareForStep(FineStep _stp)
{
	switch (_stp)
	{
	case OPEN_FILE:
		/// disable all widgets
		disableWidgetsForStep(OPEN_FILE);

		/// enable following widgets/actions:
		ui_compact.openFile->setEnabled(true);
		resetParams(OPEN_FILE);
		break;
	case POST_OPEN_FILE:
		disableWidgetsForStep(POST_OPEN_FILE);
		ui_compact.readRadii->setEnabled(true);
		ui_compact.origSurfViewGroup->setEnabled(true);
		ui_compact.saveView->setEnabled(true);
		ui_compact.loadView->setEnabled(true);
		resetParams(POST_OPEN_FILE);
		ui_compact.colorPerFace->setChecked(true);
		break;
	case BURN_MA:
		disableWidgetsForStep(BURN_MA);
		ui_compact.maViewGroup->setEnabled(true);
		ui_compact.maBurnGroup->setEnabled(true);
		resetParams(BURN_MA);
		if (!m_debugMode)
		{
			// adaptive subdivision
			ui_compact.stSubdivCombo->setCurrentIndex(1);
		}
		// show both orig and MA
		ui_compact.hideOrig->blockSignals(true);
		ui_compact.hideOrig->setChecked(true);
		ui_compact.hideOrig->blockSignals(false);
		ui_compact.hideOrig->setChecked(false);
		ui_compact.hideMA->blockSignals(true);
		ui_compact.hideMA->setChecked(true);
		ui_compact.hideMA->blockSignals(false);
		ui_compact.hideMA->setChecked(false);
		break;
	case PRUNE_MA:
		disableWidgetsForStep(PRUNE_MA);
		ui_compact.maMeasureGroup->setEnabled(true);
		ui_compact.distVisBtnGroup->setEnabled(true);
		resetParams(PRUNE_MA);
		if (!m_debugMode)
		{
			ui_compact.computeAllDistBtn->click();
			ui_compact.doMAFaceScalarDiffusion->setChecked(true);
			ui_compact.usePerSheetBox->setChecked(true);
		}
		ui_compact.colorPerFace->setChecked(true);
		// forcefully to emit a signal to set the measure to visualize on ma to be "2"
		ui_compact.visMADistCombo->blockSignals(true);
		ui_compact.visMADistCombo->setCurrentIndex(1);
		ui_compact.visMADistCombo->blockSignals(false);
		ui_compact.visMADistCombo->setCurrentIndex(2);
		// forcefully to emit a signal to set the measure to prune on ma to be "2"
		ui_compact.pruneMADistCombo1->blockSignals(true);
		ui_compact.pruneMADistCombo1->setCurrentIndex(1);
		ui_compact.pruneMADistCombo1->blockSignals(false);
		ui_compact.pruneMADistCombo1->setCurrentIndex(2);
		// allow exporting
		ui_compact.exportBTGroup->setEnabled(true);
		break;
	case DUALIZE:
		disableWidgetsForStep(DUALIZE);
		ui_compact.mcDualizeGroup->setEnabled(true);
		resetParams(DUALIZE);
		// only show MA
		ui_compact.hideMA->blockSignals(true);
		ui_compact.hideMA->setChecked(true);
		ui_compact.hideMA->blockSignals(false);
		ui_compact.hideMA->setChecked(false);
		if (!m_debugMode)
		{
			ui_compact.hsCreateGroup->setEnabled(true);
		}
		break;
	case BURN_MC:
		disableWidgetsForStep(BURN_MC);
		ui_compact.mcBurnGroup->setEnabled(true);
		resetParams(BURN_MC);
		if (!m_debugMode)
		{
			ui_compact.protectBT2Box->setChecked(true);
			ui_compact.underEstimateBox->setChecked(true);
		}
		break;
	case POST_BURN_MC:
		disableWidgetsForStep(POST_BURN_MC);
		ui_compact.mcViewGroup->setEnabled(true);
		ui_compact.mcPruneGroup->setEnabled(true);
		ui_compact.mpGroup->setEnabled(true);
		resetParams(POST_BURN_MC);
		// visualize measure # 5 on MC
		ui_compact.visMCDistCombo->blockSignals(true);
		ui_compact.visMCDistCombo->setCurrentIndex(0);
		ui_compact.visMCDistCombo->blockSignals(false);
		ui_compact.visMCDistCombo->setCurrentIndex(2);
		ui_compact.pruneMCDistCombo1->setCurrentIndex( 2 );
		//onVisDistMCComboChanged(ui_compact.visMCDistCombo->currentIndex());
		break;
	case COMPUTE_HS:
		disableWidgetsForStep(COMPUTE_HS);
		ui_compact.hsCreateGroup->setEnabled(true);
		resetParams(COMPUTE_HS);
		// only show MC
		ui_compact.hideMC->blockSignals(true);
		ui_compact.hideMC->setChecked(true);
		ui_compact.hideMC->blockSignals(false);
		ui_compact.hideMC->setChecked(false);
		break;
	case PRUNE_HS:
		disableWidgetsForStep(PRUNE_HS);
		ui_compact.hsViewGroup->setEnabled(true);
		ui_compact.hsPruneGroup->setEnabled(true);
		ui_compact.exportHSGroup->setEnabled(true);
		ui_compact.ExportSkelBox->setEnabled(true);
		resetParams(PRUNE_HS);
		break;
	case WRITE_BT_SKEL:
		disableWidgetsForStep(WRITE_BT_SKEL);
		ui_compact.exportGroup->setEnabled(true);
		ui_compact.ExportSkelBox->setEnabled(true);
		resetParams(WRITE_BT_SKEL);
		// only show HS
		ui_compact.hideHS->blockSignals(true);
		ui_compact.hideHS->setChecked(true);
		ui_compact.hideHS->blockSignals(false);
		ui_compact.hideHS->setChecked(false);
		break;
	default:
		break;
	}
	if (_stp > BURN_MA)
	{
		ui_compact.MAAlphaGroup->setEnabled(true);
	}
	if (_stp >= PRUNE_MA)
	{
		ui_compact.exportGroup->setEnabled(true);
		ui_compact.exportBTGroup->setEnabled(true);
		ui_compact.exportBtn->setEnabled(true);
	}
}

void ETMainWindow::disableWidgetsForStep(FineStep _stp)
{
	switch(_stp)
	{
	case OPEN_FILE:
		ui_compact.saveView->setDisabled(true);
		ui_compact.loadView->setDisabled(true);
		ui_compact.readRadii->setDisabled(true);
	case POST_OPEN_FILE:
		ui_compact.origSurfViewGroup->setDisabled(true);
		ui_compact.MAAlphaGroup->setDisabled(true);
	case BURN_MA:
		ui_compact.maViewGroup->setDisabled(true);
		ui_compact.maBurnGroup->setDisabled(true);
		ui_compact.distVisBtnGroup->setDisabled(true);
		ui_compact.exportBTGroup->setDisabled(true);
	case PRUNE_MA:
		ui_compact.maMeasureGroup->setDisabled(true);
	case DUALIZE:
		ui_compact.mcDualizeGroup->setDisabled(true);
	case BURN_MC:
		ui_compact.mcBurnGroup->setDisabled(true);
	case POST_BURN_MC:
		ui_compact.mcViewGroup->setDisabled(true);
		ui_compact.mcPruneGroup->setDisabled(true);
		ui_compact.mpGroup->setDisabled(true);
	case COMPUTE_HS:
		ui_compact.hsCreateGroup->setDisabled(true);
		ui_compact.exportHSGroup->setDisabled(true);
	case PRUNE_HS:
		ui_compact.hsViewGroup->setDisabled(true);
		ui_compact.hsPruneGroup->setDisabled(true);
	case WRITE_BT_SKEL:
		ui_compact.exportBTGroup->setDisabled(true);
		ui_compact.ExportSkelBox->setDisabled(true);
		ui_compact.exportGroup->setDisabled(true);
	default:
		// surf func section
		ui_compact.surfFuncGroup->setDisabled(true);
		// misc section
		ui_compact.preprocessGroup->setDisabled(true);
		ui_compact.miscBurntimeGroup->setDisabled(true);
		ui_compact.renderParamGroup->setDisabled(true);
		ui_compact.qmatGroup->setDisabled(true);
	}
}

GLArea::FaceFieldType ETMainWindow::get_MA_face_dist_type(int _idx)
{
	GLArea::FaceFieldType face_field;
	switch (_idx) // burn dist 
	{
	case 0:
		face_field = GLArea::BT2_f;
		break;
	case 1: 
		face_field = GLArea::BT3_f;
		break;
	case 2:
		face_field = GLArea::DIFF_f;
		break;
	case 3:
		face_field = GLArea::DIFF_REL_f;
		break;
	case 4:
		face_field = GLArea::ANGLE_f;
		break;
	case 5:
		face_field = GLArea::LAMBDA_f;
		break;
	case 6:
		face_field = GLArea::GEODESIC_f;
		break;
	default:
		face_field = GLArea::BT2_f;
		break;
	}

	return face_field;
}

GLArea::VertFieldType ETMainWindow::get_MA_vert_dist_type(int _idx)
{
	GLArea::VertFieldType vert_field;
	switch (_idx)
	{
	case 0: // burn dist
		vert_field = GLArea::BT2_v;
		break;
	case 1: // dist2surf
		vert_field = GLArea::BT3_v;
		break;
	case 2: // diff
		vert_field = GLArea::DIFF_v;
		break;
	//case 3: // relative diff
	//	vert_field = GLArea::DIFF_REL_v;
	//	break;
	default:
		vert_field = GLArea::BT2_v;
		break;
	}

	return vert_field;
}

GLArea::DistMC ETMainWindow::get_MC_dist_type(int _idx)
{
	GLArea::DistMC dist_type;
	switch (ui_compact.visMCDistCombo->currentIndex())
	{
	/*case 0:
		dist_type = GLArea::BT3_MC;
		break;
	case 1:
		dist_type = GLArea::BT2_MC;
		break;
	case 2:
		dist_type = GLArea::BT2_BT3_MC;
		break;
	case 3:
		dist_type = GLArea::BT2_BT3_REL_MC;
		break;
	case 4:
		dist_type = GLArea::BT1_MC;
		break;
	case 5:
		dist_type = GLArea::BT1_BT2_MC;
		break;
	case 6:
		dist_type = GLArea::BT1_BT2_REL_MC;
		break;*/
	case 0:
		dist_type = GLArea::BT3_MC;
		break;
	case 1:
		dist_type = GLArea::BT2_MC;
		break;
	case 2:
		dist_type = GLArea::BT2_BT3_MC;
		break; 
	case 3:
		dist_type = GLArea::BT1_MC;
		break;
	case 4:
		dist_type = GLArea::BT1_BT2_MC;
		break;
	default:
		dist_type = GLArea::BT1_BT2_MC;
		break;
	}

	return dist_type;
}

QColor ETMainWindow::getIdealTextColor(const QColor& rBackgroundColor) const
{
	const int THRESHOLD = 105;
	int BackgroundDelta = 
		(rBackgroundColor.red() * 0.299) 
		+ (rBackgroundColor.green() * 0.587) 
		+ (rBackgroundColor.blue() * 0.114);
	return QColor((255- BackgroundDelta < THRESHOLD) ? Qt::black : Qt::white);
}

//
// H.S. exploration
//

void ETMainWindow::onCreateHSClicked()
{
	if (!m_debugMode)
	{
		// click the button for users
		ui_compact.dualizeBtn->click();
	}

#ifdef PRINT_MEM_USAGE
	PROCESS_MEMORY_COUNTERS pmc;
	GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc));
	SIZE_T memPreUsed, memCurUsed;
	memPreUsed = memCurUsed = pmc.WorkingSetSize;
#endif // PRINT_MEM_USAGE

	ui_compact.statusbar->showMessage("Creating the skeleton...");
	glarea->createHS();
	ui_compact.statusbar->showMessage("Done: creating skeleton.");

#ifdef PRINT_MEM_USAGE
	GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc));
	memCurUsed = pmc.WorkingSetSize;
	cout << "mem. usage: " << memCurUsed / (1024*1024) << endl;
	cout << "mem. usage (by skeleton creation() ) in MB: " << (memCurUsed - memPreUsed) / (1024*1024) << endl;
	memPreUsed = memCurUsed;
#endif // PRINT_MEM_USAGE
	
	prepareForStep(PRUNE_HS);
	prepareForStep(WRITE_BT_SKEL);
	if (!m_debugMode)
	{
		// click "update" button for user
		ui_compact.visHSBtn->click();
	}
}

void ETMainWindow::onVisHSClicked()
{
	ui_compact.statusbar->showMessage("Pruning the skeleton...");

	// update the threshold for degenerate poly faces
	glarea->setHSFaceDegenerateThreshold( ui_compact.hsFaceDegenThreshold->value() );
	glarea->setComponentFaceNumberThresh( ui_compact.hsCmpntFaceNumTSpin->value() );
	// setup pruning params
	float face_diff_ratio;
	float face_reldiff_ratio;
	float line_diff_ratio;
	float line_reldiff_ratio;

	if ( ui_compact.hsUseTypedInput->isChecked() )
	{
		face_diff_ratio = (float)ui_compact.hsFaceDiffSpin->value();
		face_reldiff_ratio = (float)ui_compact.hsFaceRelDiffSpin->value();
		line_diff_ratio = (float)ui_compact.hsCurveDiffSpin->value();
		line_reldiff_ratio = (float)ui_compact.hsCurveRelDiffSpin->value();
	}
	else
	{
		face_diff_ratio = (float)ui_compact.hsBt2Bt3Slider->value() / ui_compact.hsBt2Bt3Slider->maximum();
		face_reldiff_ratio = (float)ui_compact.hsBt2Bt3RelSlider->value() / ui_compact.hsBt2Bt3RelSlider->maximum();
		line_diff_ratio = (float)ui_compact.hsBt1Bt2Slider->value() / ui_compact.hsBt1Bt2Slider->maximum();
		line_reldiff_ratio = (float)ui_compact.hsBt1Bt2RelSlider->value() / ui_compact.hsBt1Bt2RelSlider->maximum();
	}

	if ( m_consoleMode ) 
	{
		// override some parameters when in console mode
		face_diff_ratio = m_skel_theta2;
		line_diff_ratio = m_skel_theta1;
		glarea->pruneHS(
			face_diff_ratio, face_reldiff_ratio,
			line_diff_ratio, line_reldiff_ratio,
			true, /*treat values as absolute*/
			ui_compact.removeSmallCmpnts->isChecked()
		);
	}
	else
	{
		// perform pruning of HS
		glarea->pruneHS(
			face_diff_ratio, face_reldiff_ratio,
			line_diff_ratio, line_reldiff_ratio,
			ui_compact.hsUseTypedInput->isChecked(),
			ui_compact.removeSmallCmpnts->isChecked()
		);
		// upload hs geometry and associated attributes
		glarea->uploadHS();
	}

	ui_compact.statusbar->showMessage("Done: pruning the skeleton. Skeleton is updated.");
}

void ETMainWindow::onHSRemoveSmallCmpntsChecked(int _state)
{
	ui_compact.hsFaceDegenThreshold->setEnabled(_state == Qt::Checked);
}

void ETMainWindow::onDegenerateFaceThreshChanged(int _val)
{
	
}

void ETMainWindow::onHSTheta2Changed(int _val)
{
	float ratio = (float)(_val) / 
		(ui_compact.hsBt2Bt3Slider->maximum() - ui_compact.hsBt2Bt3Slider->minimum());
	float abs_t;
	glarea->getHSTheta2Abs(ratio, abs_t);
	ui_compact.hsFaceDiffSpin->blockSignals(true);
	ui_compact.hsFaceDiffSpin->setValue(abs_t);
	ui_compact.hsFaceDiffSpin->blockSignals( false );
}

void ETMainWindow::onHSTheta2Spun(double _val)
{
	float abs_t = (float)_val;
	float ratio;
	glarea->getHSTheta2Ratio(abs_t, ratio);
	ui_compact.hsBt2Bt3Slider->blockSignals(true);
	ui_compact.hsBt2Bt3Slider->setValue(
		ratio * (ui_compact.hsBt2Bt3Slider->maximum() - ui_compact.hsBt2Bt3Slider->minimum())
	);
	ui_compact.hsBt2Bt3Slider->blockSignals( false );
}

void ETMainWindow::onHSTheta1Changed(int _val)
{
	float ratio = (float)(_val) / 
		(ui_compact.hsBt2Bt3Slider->maximum() - ui_compact.hsBt2Bt3Slider->minimum());
	float abs_t;
	glarea->getHSTheta1Abs(ratio, abs_t);
	ui_compact.hsCurveDiffSpin->blockSignals(true);
	ui_compact.hsCurveDiffSpin->setValue(abs_t);
	ui_compact.hsCurveDiffSpin->blockSignals(false);
}

void ETMainWindow::onHSTheta1Spun(double _val)
{
	float abs_t = (float)_val;
	float ratio;
	glarea->getHSTheta1Ratio(abs_t, ratio);
	ui_compact.hsBt1Bt2Slider->blockSignals(true);
	ui_compact.hsBt1Bt2Slider->setValue(
		ratio * (ui_compact.hsBt1Bt2Slider->maximum() - ui_compact.hsBt1Bt2Slider->minimum())
		);
	ui_compact.hsBt1Bt2Slider->blockSignals( false );
}

void ETMainWindow::onSurfFuncCorrespSchemeComboChanged(int _idx)
{
	if ( _idx == 0 ) // enable range search params
	{
		ui_compact.surfFuncRadiusEps->setEnabled(true);
		ui_compact.surfFuncKNN->setEnabled(true);
	}
	else // enable KNN search params 
	{
		ui_compact.surfFuncKNN->setEnabled(true);
		ui_compact.surfFuncRadiusEps->setEnabled(false);
	}
}

void ETMainWindow::onSurfFuncCorrespSetupClicked()
{
	auto corresp_idx = ui_compact.surfFuncCorrespSchemeCombo->currentIndex();
	if ( corresp_idx == 0 )
	{
		glarea->computeSurfFuncCorrespondence(
			GLArea::SURF_FUNC_RANGE_CORRESP, 
			ui_compact.surfFuncRadiusEps->value(),
			ui_compact.surfFuncKNN->value()
			);
		enableSmoothSurfFuncWidgets();
	}
	else if ( corresp_idx == 1 )
	{
		glarea->computeSurfFuncCorrespondence(
			GLArea::SURF_FUNC_KNN_CORRESP,
			0.0f, // dummy value
			ui_compact.surfFuncKNN->value()
			);
		enableSmoothSurfFuncWidgets();
	}
	else
	{
		QMessageBox::warning(this, "Error", "Please choose a valid correspondence scheme!");
	}
}

void ETMainWindow::onSurfFuncProjectClicked()
{
	ui_compact.hideMA->setChecked(true);
	ui_compact.hideMC->setChecked(true);
	ui_compact.hideHS->setChecked(true);
	ui_compact.hideOrig->setChecked(false);
	float ratio = 
		(float)(ui_compact.smoothSMetricSlider->value()) / ui_compact.smoothSMetricSlider->maximum();
	int smetric_idx = ui_compact.pruneSMetricCombo->currentIndex();
	bool do_diffuse = ui_compact.diffuseSurfFunc->isChecked();
	glarea->projectFunctionOnSurf(smetric_idx, ratio, do_diffuse);
}

void ETMainWindow::onSurfFuncComboChanged(int _idx)
{
}

void ETMainWindow::onSurfFuncSlided(int _value)
{
}

void ETMainWindow::onDiffuseSurfFuncChecked(int _state)
{
}

void ETMainWindow::onShowMPChecked(int _state)
{
	glarea->setDrawFlag(GLArea::DRAW_MP, _state == Qt::Checked);
	if (_state == Qt::Checked)
	{
		glarea->visMP();
		//cout << "showMP turned on."<<endl;
	}
	else
	{
		//cout << "showMP turned off."<<endl;
	}
}

void ETMainWindow::onIsoSurfPrecomputeClicked()
{
	glarea->setupIsoSurface();
}

void ETMainWindow::onIsoSurfRefineClicked()
{
	glarea->refineIsoSurface(ui_compact.recursiveRefineBox->isChecked());

	float t = (float)ui_compact.isoSurfSlider->minimum() + 
		((float)ui_compact.isoSurfSlider->value() - ui_compact.isoSurfSlider->minimum()) / 
		((float)ui_compact.isoSurfSlider->maximum() - ui_compact.isoSurfSlider->minimum());
	glarea->uploadIsoSurface(t, ui_compact.hideSnappedFaceBox->isChecked());

	glarea->updateGL();
}

void ETMainWindow::onIsoSurfSlided(int _value)
{
	float t = (float)ui_compact.isoSurfSlider->minimum() + 
		((float)_value - ui_compact.isoSurfSlider->minimum()) / 
		((float)ui_compact.isoSurfSlider->maximum() - ui_compact.isoSurfSlider->minimum());
	glarea->uploadIsoSurface(t, ui_compact.hideSnappedFaceBox->isChecked());

	glarea->updateGL();
}

void ETMainWindow::onHideSnappedFacesChanged(int _state)
{
	onIsoSurfSlided(ui_compact.isoSurfSlider->value());
}

void ETMainWindow::onIsoContPrecomputeClicked()
{
	float MA_alpha = (float)(ui_compact.MATransparentSlider->value() - ui_compact.MATransparentSlider->minimum()) / 
		(ui_compact.MATransparentSlider->maximum() - ui_compact.MATransparentSlider->minimum());
	float MA_exp = ui_compact.MATransExpSlider->value();
	float new_min = ui_compact.minVisDistMASpin->value();
	float new_max = ui_compact.maxVisDistMASpin->value();
	bool use_rescale = ui_compact.clampMADistBox->isChecked();
	glarea->setupIsoContour(
		MA_alpha, MA_exp, 
		use_rescale ? &new_min : nullptr, use_rescale ? &new_max : nullptr
		);
	float t = get_cur_isoCont_t();
	glarea->uploadIsoContour(t, ui_compact.hidePastSurfBox->isChecked(), ui_compact.isoContShowMC->isChecked());
	glarea->setDrawFlag(GLArea::DRAW_ISOCONT, !ui_compact.hideIsoCont->isChecked());
}

void ETMainWindow::onIsoContSlided(int _value)
{
	float t = get_cur_isoCont_t();
	glarea->uploadIsoContour(t, ui_compact.hidePastSurfBox->isChecked(), ui_compact.isoContShowMC->isChecked());

	glarea->updateGL();
}

void ETMainWindow::onHidePastSurfChanged(int _state)
{
	onIsoContSlided(ui_compact.isoContSlider->value());
}

void ETMainWindow::onIsoContShowMCChanged(int _state)
{
	//if (_state != Qt::Checked) // show exposed mc
	//{
		glarea->uploadFromRemainedMC(_state != Qt::Checked, false, 0.0f, false);
	//}

	this->onHideMCChanged(_state != Qt::Checked);
}

void ETMainWindow::onEvolveAllPrecomputeClicked()
{
	vector<float> mc_bt1;
	glarea->getMCDistMetric(GLArea::BT1_MC, mc_bt1);
	iso_point_max = *( std::max_element(mc_bt1.begin(), mc_bt1.end()) );

	// don't forget to pre-compute for evolving MC
	if ( !glarea->setUpMCDistMetricForPruning(GLArea::BT1_MC, distMin_MC, distMax_MC) )
	{
		cout << "ERROR: failed to get current distance metric to be used in pruning." << endl;
	}
	else
	{
		float dummy = 0.0f;
		glarea->colorMCEdgeBy(GLArea::BT1_MC, dummy, dummy, false, 0.0f, 1.0f, 1);
	}
}

void ETMainWindow::onEvolveAllSlided(int _value)
{
	float t = (float)ui_compact.evolveAllSlider->minimum() + 
		((float)_value - ui_compact.evolveAllSlider->minimum()) / 
		((float)ui_compact.evolveAllSlider->maximum() - ui_compact.evolveAllSlider->minimum());
	float iso = (iso_point_max - 0.0f) * t + 0.0f;

	// evolve iso-surf
	glarea->uploadIsoSurface(iso, ui_compact.hideSnappedFaceBox->isChecked(), false);
	// evolve iso-cont
	glarea->uploadIsoContour(iso, ui_compact.hidePastSurfBox->isChecked(), false, 1);
	// evolve MC
	glarea->uploadFromRemainedMC(!ui_compact.evolveAllShowEntireMC->isChecked(), true, iso, true);

	glarea->setDrawFlag(GLArea::DRAW_POINTS, true);
	//glarea->updateGL();
}

void ETMainWindow::onPrintBTClicked()
{
	glarea->printAndVisBurnTime();
}

void ETMainWindow::onShowSpheresChecked(int _state)
{
	// right now all spheres are controled by draw_mp flag
	glarea->setDrawFlag(GLArea::DRAW_MP, _state == Qt::Checked);
}

void ETMainWindow::onVisSphereSlided(int _val)
{
	int n_vts = glarea->getOriginalVts().size();
	ui_compact.visSphereSlider->setMaximum(n_vts - 1);
	ui_compact.visSphereSpin->setMaximum(n_vts - 1);
	size_t vert_id = (size_t)_val;
	glarea->visSphere(vert_id);

	ui_compact.visSphereSpin->blockSignals(true);
	ui_compact.visSphereSpin->setValue(vert_id);
	ui_compact.visSphereSpin->blockSignals(false);
}

void ETMainWindow::onVisSphereSpun(int _val)
{
	int n_vts = glarea->getOriginalVts().size();
	ui_compact.visSphereSpin->setMaximum(n_vts - 1);
	ui_compact.visSphereSlider->setMaximum(n_vts - 1);
	size_t vert_id = (size_t)_val;
	glarea->visSphere(vert_id);

	ui_compact.visSphereSlider->blockSignals(true);
	ui_compact.visSphereSlider->setValue(vert_id);
	ui_compact.visSphereSlider->blockSignals(false);
}

void ETMainWindow::onPrintSelectVertInfoClicked()
{
	glarea->printMAOriginalVert(ui_compact.visSphereSpin->value());
}

void ETMainWindow::onVisSpheresFromFileClicked()
{
	QString vts_file_name = QFileDialog::getOpenFileName(this, "Open File", "", "", NULL, QFileDialog::DontUseNativeDialog);
	vector<unsigned> vts_list;
	if ( readVertsFromFile(vts_file_name, vts_list) )
	{
		glarea->visSphere(vts_list);
		glarea->printMAOriginalVerts(vts_list);
	}
}

bool ETMainWindow::readVertsFromFile(QString _file, vector<unsigned>& _vid_list)
{
	std::ifstream ifile(_file.toStdString());
	if ( !ifile.good() )
	{
		cout << "Error: cannot open file "<<_file.toStdString()<<endl;
		return false;
	}

	unsigned vid;
	while ( !ifile.eof() )
	{
		ifile >> vid;
		_vid_list.push_back(vid);
	}

	return true;
}

void ETMainWindow::onLinesRenderOffsetSlided(int _val)
{
	float r = 0.0001f + (-0.0001f - 0.0001f) * _val / ui_compact.linesRenderOffsetSlider->maximum();
	glarea->setLinesOffset(r);
}

void ETMainWindow::onOutputWenpingBtnClicked()
{
	glarea->outputToWenping();
	cout << "Current MA output to qmat .ma file." << endl;
}

void ETMainWindow::onReadQMATBtnClicked()
{
	QString filename = QFileDialog::getOpenFileName(this, 
		tr("Select a qmat (.ma) file to open"), 
		NULL, 
		tr("qmat file (*.ma)")
		);	
	if (!filename.isEmpty())
	{
		glarea->readQMATFile(filename.toStdString());
	}
}

void ETMainWindow::onHideQMATChanged(int _state)
{
	glarea->setDrawFlag( GLArea::DRAW_QMAT, _state != Qt::Checked );
}

void ETMainWindow::onExportETBtnClicked()
{
	if (ui_compact.exportPerSectorBTBox->isChecked())
	{
		cout << endl;
		cout << "Exporting per-sector ET to file ..."<<endl;
		ui_compact.statusbar->showMessage("Exporting per-sector ET to file ...");
		glarea->exportPerSectorET();
		cout << "Done."<<endl;
		ui_compact.statusbar->showMessage("Done: exporting per-sector ET.");
		cout << endl;
	}
	else
	{
		cout << endl;
		cout << "Exporting per-vertex ET to file ..."<<endl;
		ui_compact.statusbar->showMessage("Exporting per-vertex ET to file ...");
		glarea->exportPerVertexET();
		cout << "Done."<<endl;
		ui_compact.statusbar->showMessage("Done: exporting per-vertex ET.");
		cout << endl;
	}
}

void ETMainWindow::onExportHSBtnClicked()
{
	cout << endl;
	cout << "Exporting skeleton to file ..."<<endl;
	ui_compact.statusbar->showMessage("Exporting skeleton to file ...");
	glarea->exportSkeleton();
	cout << "Done."<<endl;
	ui_compact.statusbar->showMessage("Done: exporting skeleton.");
	cout << endl;
}

void ETMainWindow::onExportBtnClicked()
{
	/*if (ui_compact.exportPerVertexBTBox->isChecked())
	{
		cout << endl;
		cout << "Exporting per-vertex BT to file ..."<<endl;
		glarea->exportPerVertexET();
		cout << "Done."<<endl;
		cout << endl;
	}
	if (ui_compact.exportPerSectorBTBox->isChecked())
	{
		cout << endl;
		cout << "Exporting per-sector BT to file ..."<<endl;
		glarea->exportPerSectorET();
		cout << "Done."<<endl;
		cout << endl;
	}
	if (ui_compact.ExportSkelBox->isChecked())
	{
		cout << endl;
		cout << "Exporting skeleton to file ..."<<endl;
		glarea->exportSkeleton();
		cout << "Done."<<endl;
		cout << endl;
	}*/
}

void ETMainWindow::handle_extra_params()
{
	if ( m_mc_meas_to_output != "null" )
		glarea->outputMCwMeasure( m_mc_meas_to_output );
}

void ETMainWindow::resetWidgetsStates()
{
	ui_compact.toolTabs->setCurrentIndex(0);

	/*ui.hideOrig->blockSignals(true);
	ui.hideMA->blockSignals(true);
	ui.hideMC->blockSignals(true);
	ui.hideHS->blockSignals(true);
	ui.drawEdge->blockSignals(true);
	ui.drawStPoints->blockSignals(true);
	ui.showMPBox->blockSignals(true);*/
	ui_compact.hideMA->setCheckState(Qt::Unchecked);
	ui_compact.drawEdge->setCheckState(Qt::Unchecked);
	ui_compact.drawStPoints->setCheckState(Qt::Unchecked);
	ui_compact.hideOrig->setCheckState(Qt::Unchecked);
	ui_compact.hideMC->setCheckState(Qt::Unchecked);
	ui_compact.showMPBox->setCheckState(Qt::Unchecked);
	/*ui.hideOrig->blockSignals(false);
	ui.hideMA->blockSignals(false);
	ui.hideMC->blockSignals(false);
	ui.drawEdge->blockSignals(false);
	ui.drawStPoints->blockSignals(false);
	ui.hideHS->blockSignals(false);
	ui.showMPBox->blockSignals(false);*/

	//ui.pointSizeSpin->setEnabled(false);

	// reset steiner subdivision widgets
	ui_compact.stSubdivCombo->blockSignals(true);
	ui_compact.stSubdivCombo->clear();
	ui_compact.stSubdivCombo->addItems(QStringList() <<"Fixed"<<"Adaptive"<<"Midpoint");
	ui_compact.stSubdivCombo->blockSignals(false);

	// reset burning MA related widgets
	ui_compact.edgeWeightCombo->blockSignals(true);
	ui_compact.edgeWeightCombo->clear();
	ui_compact.edgeWeightCombo->addItems(QStringList() << "euclidean"<<"diff");
	ui_compact.edgeWeightCombo->blockSignals(false);

	// reset prune MA related widgets
	//ui.computeAllDistBtn->setEnabled(false);
	//ui.enableFinePruneMA->setEnabled(false);
	//ui.pruneMASlider2->setEnabled(false);
	//ui.pruneMASpin2->setEnabled(false);
	//ui.pruneMASlider1->setEnabled(false);
	//ui.pruneMASpin1->setEnabled(false);
	// reset and disable distance measure group
	ui_compact.colorPerVert->setChecked(true);
	ui_compact.clampMADistBox->setChecked(false);
	//ui.minVisDistMASpin->setEnabled(false);
	//ui.maxVisDistMASpin->setEnabled(false);
	ui_compact.clampMCDistBox->setChecked(false);
	//ui.minVisDistMCSpin->setEnabled(false);
	//ui.maxVisDistMCSpin->setEnabled(false);

	//ui.enableFinePruneMA->blockSignals(true);
	ui_compact.enableFinePruneMA->setChecked(false);
	//ui.enableFinePruneMA->blockSignals(false);

	ui_compact.pruneMADistCombo1->blockSignals(true);
	ui_compact.pruneMADistCombo1->clear();
	ui_compact.pruneMADistCombo1->addItems(
		QStringList() <<"BT2"<<"BT3"<<"diff"<<"rel diff"<<"angle"<<"lambda"<<"geodesic"
		);
	ui_compact.pruneMADistCombo1->setCurrentIndex(-1);
	ui_compact.pruneMADistCombo1->blockSignals(false);
	//ui.maMeasureGroup->setEnabled(false);
	ui_compact.pruneMADistCombo2->blockSignals(true);
	ui_compact.pruneMADistCombo2->clear();
	ui_compact.pruneMADistCombo2->addItems(
		QStringList() <<"BT2"<<"BT3"<<"diff"<<"rel diff"<<"angle"<<"lambda"<<"geodesic"
		);
	ui_compact.pruneMADistCombo2->setCurrentIndex(-1);
	ui_compact.pruneMADistCombo2->blockSignals(false);
	//ui.maMeasureGroup->setEnabled(false);

	int filterByDiffSlider_max = 10000;
	ui_compact.pruneMASlider2->setMaximum(filterByDiffSlider_max);

	int transparent_slider_max = 100;
	ui_compact.MATransparentSlider->setMaximum(transparent_slider_max);

	/* reset stEdgesBox */
	ui_compact.hideMC->setChecked(false);
	ui_compact.visBurnt->setChecked(false);

	//ui.mcDualizeGroup->setEnabled(false);

	ui_compact.edgeDualOptCombo->blockSignals(true);
	ui_compact.edgeDualOptCombo->clear();
	ui_compact.edgeDualOptCombo->addItems(
		QStringList() 
		<<"simple"
		<<"pseudo inverse"
		<<"weight closeness to center"
		<<"new weighting" );
	ui_compact.edgeDualOptCombo->setCurrentIndex(0);
	ui_compact.edgeDualOptCombo->blockSignals(false);

	ui_compact.polyDualOptCombo->blockSignals(true);
	ui_compact.polyDualOptCombo->clear();
	ui_compact.polyDualOptCombo->addItems(
		QStringList() 
		<<"simple"
		<<"pseudo inverse"
		<<"weight closeness to center"
		<<"new weighting" );
	ui_compact.polyDualOptCombo->setCurrentIndex(0);
	ui_compact.polyDualOptCombo->blockSignals(false);

	ui_compact.visBurnt->setChecked(false);
	ui_compact.hideMC->setChecked(false);

	// reset the dist metric combo for medial curve vis.
	ui_compact.visMCDistCombo->blockSignals(true);
	ui_compact.visMCDistCombo->clear();
	ui_compact.visMCDistCombo->addItems(QStringList() <<"BT3"<<"BT2"<<"BT2-BT3"<<"1-BT3/BT2"<<"BT1"<<"BT1-BT2"<<"1-BT2/BT1");
	ui_compact.visMCDistCombo->setCurrentIndex(0);
	ui_compact.visMCDistCombo->blockSignals(false);

	resetPruneMC();

	// reset widgets for surf func app.
	ui_compact.pruneSMetricCombo->blockSignals(true);
	ui_compact.pruneSMetricCombo->clear();
	ui_compact.pruneSMetricCombo->addItems(
		QStringList() <<"diameter"<<"width"<<"width extremity"<<"length extremity"
		);
	ui_compact.pruneSMetricCombo->setCurrentIndex(-1);
	ui_compact.pruneSMetricCombo->blockSignals(false);
	resetSurfFuncWidgets();

	// reset iso-surface evolution slider
	ui_compact.isoSurfSlider->blockSignals(true);
	ui_compact.isoSurfSlider->setValue(0);
	ui_compact.isoSurfSlider->blockSignals(false);

	// reset DP transparency 
	if (glarea)
	{
		ui_compact.DPMaxRendersSpin->blockSignals(true);
		ui_compact.DPMaxRendersSpin->setValue(3);
		ui_compact.DPMaxRendersSpin->blockSignals(false);
		ui_compact.useTrueTransparencyBox->blockSignals(true);
		ui_compact.useTrueTransparencyBox->setChecked(false);
		ui_compact.useTrueTransparencyBox->blockSignals(false);
	}

	//resetParams();
}

void ETMainWindow::resetParams(FineStep _stp)
{
	switch (_stp)
	{
	case OPEN_FILE:
		ui_compact.toolTabs->setCurrentIndex(0);
		// perturbation params
		ui_compact.epsilonSpin->setValue(0.0);
		ui_compact.pertSpin->setValue(0.0);
	case POST_OPEN_FILE:
		
	case BURN_MA:
		// reset steiner subdivision widgets
		ui_compact.nFixedSteinerSpin->setValue(0.004);
		ui_compact.stSubdivCombo->blockSignals(true);
		ui_compact.stSubdivCombo->clear();
		ui_compact.stSubdivCombo->addItems(QStringList() <<"Fixed"<<"Adaptive"<<"Midpoint");
		ui_compact.stSubdivCombo->blockSignals(false);
		// reset burning MA related widgets
		ui_compact.edgeWeightCombo->blockSignals(true);
		ui_compact.edgeWeightCombo->clear();
		ui_compact.edgeWeightCombo->addItems(QStringList() << "euclidean"<<"diff");
		ui_compact.edgeWeightCombo->blockSignals(false);
	case PRUNE_MA:
		/*ui_compact.clampMADistBox->setChecked(true);*/
		ui_compact.clampMADistBox->setChecked(false);
		onClampMADistBoxChanged(Qt::Unchecked);
		ui_compact.enableFinePruneMA->setChecked(false);

		ui_compact.pruneMADistCombo1->blockSignals(true);
		ui_compact.pruneMADistCombo1->clear();
		ui_compact.pruneMADistCombo1->addItems(
			QStringList() <<"BT_M"<<"R"<<"ET"<<"ET(rel)"<<"Angle"<<"Lambda"<<"MGF"
			);
		ui_compact.pruneMADistCombo1->setCurrentIndex(-1);
		ui_compact.pruneMADistCombo1->blockSignals(false);
		ui_compact.pruneMADistCombo2->blockSignals(true);
		ui_compact.pruneMADistCombo2->clear();
		ui_compact.pruneMADistCombo2->addItems(
			QStringList() <<"BT_M"<<"R"<<"ET"/*<<"ET(rel)"<<"Angle"<<"Lambda"<<"MGF"*/
			);
		ui_compact.pruneMADistCombo2->setCurrentIndex(-1);
		ui_compact.pruneMADistCombo2->blockSignals(false);
	case DUALIZE:
		ui_compact.edgeDualOptCombo->blockSignals(true);
		ui_compact.edgeDualOptCombo->clear();
		ui_compact.edgeDualOptCombo->addItems(
			QStringList() 
			<<"simple"
			<<"pseudo inverse"
			<<"weight closeness to center"
			<<"new weighting" );
		ui_compact.edgeDualOptCombo->setCurrentIndex(0);
		ui_compact.edgeDualOptCombo->blockSignals(false);

		ui_compact.polyDualOptCombo->blockSignals(true);
		ui_compact.polyDualOptCombo->clear();
		ui_compact.polyDualOptCombo->addItems(
			QStringList() 
			<<"simple"
			<<"pseudo inverse"
			<<"weight closeness to center"
			<<"new weighting" );
		ui_compact.polyDualOptCombo->setCurrentIndex(0);
		ui_compact.polyDualOptCombo->blockSignals(false);
		ui_compact.clampMCDistBox->setChecked(true);
		ui_compact.clampMCDistBox->setChecked(false);
	case BURN_MC:
	case POST_BURN_MC:
		ui_compact.visMCDistCombo->blockSignals( true );
		ui_compact.visMCDistCombo->clear();
		//ui_compact.visMCDistCombo->addItems(QStringList() <<"r"<<"BT_M"<<"ET_M"<<"ET_M(rel)"<<"BT_C"<<"ET_C"<<"ET_C(rel)");
		//ui_compact.visMCDistCombo->addItems(QStringList() <<"BT_M"<<"BT_C"<<"ET_C");
		ui_compact.visMCDistCombo->addItems( QStringList() << "r" << "BT_M" <<"ET_M" << "BT_C" << "ET_C" );
		ui_compact.visMCDistCombo->setCurrentIndex( 0 );
		ui_compact.visMCDistCombo->blockSignals(false);
		ui_compact.pruneMCSlider1->setMaximum(10000);
		ui_compact.pruneMCDistCombo1->blockSignals(true);
		ui_compact.pruneMCDistCombo1->clear();
		//ui_compact.pruneMCDistCombo1->addItems(QStringList() <<"BT_M"<<"BT_C"<<"ET_C"<<"ET_C(rel)");
		ui_compact.pruneMCDistCombo1->addItems(QStringList() <<"BT_M"<<"BT_C"<<"ET_C");
		ui_compact.pruneMCDistCombo1->setCurrentIndex( 0 );
		ui_compact.pruneMCDistCombo1->blockSignals(false);
		//ui_compact.pruneMCDistSpin1->setMaximum( 10000 );
	case PRUNE_HS:
		{
			// make sure nothing is pruned in the beginning. 
			ui_compact.hsBt2Bt3Slider->setValue( -1 );
			// setting threshold for degenerate face cluster.
			ui_compact.hsFaceDegenThreshold->setValue(0.05);
			int facenum_t = 2000;
			if (glarea)
			{
				int numfaces = glarea->getNumFacesOfMA();
				facenum_t = numfaces > 0 ? numfaces * 0.1f : facenum_t;
			}
			ui_compact.hsCmpntFaceNumTSpin->setValue(facenum_t);
		}
		ui_compact.removeSmallCmpnts->setChecked(false);
		ui_compact.removeSmallCmpnts->setChecked(true);
	default:

		break;
	}

	// hide all geometry by default
	ui_compact.hideMA->blockSignals(true);
	ui_compact.hideMA->setChecked(false);
	ui_compact.hideMA->blockSignals(false);
	ui_compact.hideMA->setChecked(true);
	ui_compact.hideOrig->blockSignals(true);
	ui_compact.hideOrig->setChecked(false);
	ui_compact.hideOrig->blockSignals(false);
	ui_compact.hideOrig->setChecked(true);
	ui_compact.hideMC->blockSignals(true);
	ui_compact.hideMC->setChecked(false);
	ui_compact.hideMC->blockSignals(false);
	ui_compact.hideMC->setChecked(true);
	ui_compact.hideHS->blockSignals(true);
	ui_compact.hideHS->setChecked(false);
	ui_compact.hideHS->blockSignals(false);
	ui_compact.hideHS->setChecked(true);
	ui_compact.drawEdge->blockSignals(true);
	ui_compact.drawEdge->setChecked(true);
	ui_compact.drawEdge->blockSignals(false);
	ui_compact.drawEdge->setChecked(false);
	ui_compact.drawStPoints->blockSignals(true);
	ui_compact.drawStPoints->setChecked(true);
	ui_compact.drawStPoints->blockSignals(false);
	ui_compact.drawStPoints->setChecked(false);
	ui_compact.showMPBox->blockSignals(true);
	ui_compact.showMPBox->setChecked(true);
	ui_compact.showMPBox->blockSignals(false);
	ui_compact.showMPBox->setChecked(false);
	ui_compact.enableFinePruneMA->blockSignals(true);
	ui_compact.enableFinePruneMA->setChecked(true);
	ui_compact.enableFinePruneMA->blockSignals(false);
	ui_compact.enableFinePruneMA->setChecked(false);
}

void ETMainWindow::resetPruneMC()
{
	// medial curve related widget reset/setup
	ui_compact.mcPruneGroup->setEnabled(false);
	ui_compact.pruneMCSlider1->setMaximum(10000);
	ui_compact.pruneMCDistCombo1->blockSignals(true);
	ui_compact.pruneMCDistCombo1->clear();
	ui_compact.pruneMCDistCombo1->addItems(QStringList() <<"BT2"<<"BT1"<<"diff"<<"rel diff");
	ui_compact.pruneMCDistCombo1->setCurrentIndex(0);
	ui_compact.pruneMCDistCombo1->blockSignals(false);

	//ui_compact.pruneMCDistSpin1->setMaximum(10000);
}

void ETMainWindow::enablePruneMC()
{
	ui_compact.mcPruneGroup->setEnabled(true);
	ui_compact.pruneMCSlider1->setEnabled(true);
	ui_compact.visBallStick->blockSignals(true);
	ui_compact.visBallStick->setChecked(false);
	ui_compact.visBallStick->blockSignals(false);
	ui_compact.visBallStick->setEnabled(true);
}

void ETMainWindow::resetSurfFuncWidgets()
{
	ui_compact.surfFuncCorrespSchemeCombo->blockSignals(true);
	ui_compact.surfFuncRadiusEps->blockSignals(true);
	ui_compact.surfFuncKNN->blockSignals(true);
	ui_compact.pruneSMetricCombo->blockSignals(true);
	ui_compact.smoothSMetricSlider->blockSignals(true);
	ui_compact.diffuseSurfFunc->blockSignals(true);

	ui_compact.surfFuncCorrespSchemeCombo->clear();
	ui_compact.surfFuncCorrespSchemeCombo->addItems(
		QStringList() <<"range-based"<<"knn-based"
		);
	ui_compact.surfFuncCorrespSchemeCombo->setCurrentIndex(-1);
	ui_compact.surfFuncRadiusEps->setValue(0.0);
	ui_compact.surfFuncKNN->setValue(1);
	ui_compact.pruneSMetricCombo->setCurrentIndex(-1);
	ui_compact.smoothSMetricSlider->setValue(0);
	ui_compact.diffuseSurfFunc->setChecked(false);

	ui_compact.surfFuncCorrespSchemeCombo->blockSignals(false);
	ui_compact.surfFuncRadiusEps->blockSignals(false);
	ui_compact.surfFuncKNN->blockSignals(false);
	ui_compact.pruneSMetricCombo->blockSignals(false);
	ui_compact.smoothSMetricSlider->blockSignals(false);
	ui_compact.diffuseSurfFunc->blockSignals(false);

	ui_compact.surfFuncGroup->setEnabled(false);
}

void ETMainWindow::enableSurfFuncWidgets()
{
	ui_compact.surfFuncGroup->setEnabled(true);
	ui_compact.surfFuncCorrespSchemeCombo->setEnabled(true);
	ui_compact.surfFuncSetupCorrespBtn->setEnabled(true);
	ui_compact.surfFuncKNN->setEnabled(true);
	ui_compact.surfFuncRadiusEps->setEnabled(true);
}

void ETMainWindow::enableSmoothSurfFuncWidgets()
{
	ui_compact.pruneSMetricCombo->setEnabled(true);
	ui_compact.smoothSMetricSlider->setEnabled(true);
	ui_compact.diffuseSurfFunc->setEnabled(true);
	ui_compact.surfFuncProjectBtn->setEnabled(true);
}

void ETMainWindow::getPruneMAThreshold(float& _v1, float& _v2)
{
	if (ui_compact.pruneMADistCombo1->currentIndex() >= 0)
	{
		float range_1 = pruneDistOnMA_max_1 - pruneDistOnMA_min_1;
		_v1 = pruneDistOnMA_min_1 + 
			range_1 * ui_compact.pruneMASlider1->value() / ui_compact.pruneMASlider1->maximum();
		ui_compact.pruneMASpin1->blockSignals(true);
		ui_compact.pruneMASpin1->setValue(_v1);
		ui_compact.pruneMASpin1->blockSignals(false);
	}
	else
		_v1 = -999.0f;// effectively will be ignored

	// reflect current ratio: prune value / bbox size
	ui_compact.maPruneRatio->setValue( glarea->toETmPruneRatio(_v1) );

	if (ui_compact.pruneMADistCombo2->currentIndex() >= 0)
	{
		float range_2 = pruneDistOnMA_max_2 - pruneDistOnMA_min_2;
		_v2 = pruneDistOnMA_min_2 + 
			range_2 * ui_compact.pruneMASlider2->value() / ui_compact.pruneMASlider2->maximum();

		ui_compact.pruneMASpin2->blockSignals(true);
		ui_compact.pruneMASpin2->setValue(_v2);
		ui_compact.pruneMASpin2->blockSignals(false);
	}
	else 
		_v2 = -999.0f; // will be ignored
}

float ETMainWindow::get_cur_isoCont_t()
{
	return 1.1f * (float)ui_compact.isoContSlider->minimum() + 
		((float)ui_compact.isoContSlider->value() - ui_compact.isoContSlider->minimum()) / 
		((float)ui_compact.isoContSlider->maximum() - ui_compact.isoContSlider->minimum());
}