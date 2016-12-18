/********************************************************************************
** Form generated from reading UI file 'compact.ui'
**
** Created by: Qt User Interface Compiler version 5.7.0
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_COMPACT_H
#define UI_COMPACT_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QCheckBox>
#include <QtWidgets/QComboBox>
#include <QtWidgets/QDoubleSpinBox>
#include <QtWidgets/QFormLayout>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QGroupBox>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QRadioButton>
#include <QtWidgets/QSlider>
#include <QtWidgets/QSpacerItem>
#include <QtWidgets/QSpinBox>
#include <QtWidgets/QSplitter>
#include <QtWidgets/QStatusBar>
#include <QtWidgets/QTabWidget>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_MainWindow
{
public:
    QAction *openFile;
    QAction *saveView;
    QAction *loadView;
    QAction *readRadii;
    QAction *actionHelp;
    QWidget *centralwidget;
    QGridLayout *gridLayout_11;
    QSplitter *splitter;
    QTabWidget *toolTabs;
    QWidget *fileTab;
    QGridLayout *gridLayout_17;
    QGroupBox *loadFileGroup;
    QGridLayout *gridLayout_40;
    QLabel *label_11;
    QLineEdit *maFileEdit;
    QPushButton *browseMAFileBtn;
    QLabel *label_12;
    QLineEdit *shapeFileEdit;
    QPushButton *browse3DShapeFileBtn;
    QLabel *label_13;
    QLineEdit *radiusFileEdit;
    QPushButton *browseRadiusFileBtn;
    QPushButton *loadFilesBtn;
    QSpacerItem *verticalSpacer_3;
    QWidget *origSurfTab;
    QGridLayout *gridLayout_2;
    QGroupBox *groupBox;
    QGridLayout *gridLayout_30;
    QSpinBox *DPMaxRendersSpin;
    QPushButton *colorBGBtn;
    QCheckBox *drawEdge;
    QPushButton *saveViewBtn;
    QPushButton *loadViewBtn;
    QCheckBox *useTrueTransparencyBox;
    QGroupBox *origSurfViewGroup;
    QGridLayout *gridLayout_16;
    QPushButton *colorOrigBtn;
    QCheckBox *hideOrig;
    QLabel *label_4;
    QSlider *origTransparentSlider;
    QGroupBox *maViewGroup;
    QGridLayout *gridLayout_6;
    QCheckBox *hideMA;
    QCheckBox *turnOffLightingForMABox;
    QCheckBox *drawStPoints;
    QPushButton *colorMABtn;
    QGroupBox *distVisBtnGroup;
    QGridLayout *gridLayout_37;
    QComboBox *visMADistCombo;
    QGridLayout *gridLayout_38;
    QCheckBox *clampMADistBox;
    QDoubleSpinBox *minVisDistMASpin;
    QDoubleSpinBox *maxVisDistMASpin;
    QGroupBox *MAAlphaGroup;
    QGridLayout *gridLayout_10;
    QSlider *MATransparentSlider;
    QSlider *MATransExpSlider;
    QGroupBox *mcViewGroup;
    QGridLayout *gridLayout_4;
    QCheckBox *visBurnt;
    QPushButton *colorMCBtn;
    QCheckBox *hideMC;
    QLabel *label_10;
    QDoubleSpinBox *MCWidthSpin;
    QGroupBox *groupBox_2;
    QGridLayout *gridLayout_8;
    QComboBox *visMCDistCombo;
    QGridLayout *gridLayout_39;
    QCheckBox *clampMCDistBox;
    QDoubleSpinBox *minVisDistMCSpin;
    QDoubleSpinBox *maxVisDistMCSpin;
    QGroupBox *MCAlphaGroup;
    QGridLayout *gridLayout_22;
    QSlider *MCAlphaSlider;
    QSlider *MCTransExpSlider;
    QGroupBox *hsViewGroup;
    QGridLayout *gridLayout_9;
    QCheckBox *hideHS;
    QWidget *MATab;
    QGridLayout *gridLayout_12;
    QGroupBox *maBurnGroup;
    QGridLayout *gridLayout_5;
    QDoubleSpinBox *nFixedSteinerSpin;
    QPushButton *cleanTopoBtn;
    QPushButton *burnBtn;
    QGroupBox *maMeasureGroup;
    QGridLayout *gridLayout_3;
    QSlider *pruneMASlider1;
    QComboBox *pruneMADistCombo1;
    QDoubleSpinBox *pruneMASpin1;
    QCheckBox *enableFinePruneMA;
    QGroupBox *exportBTGroup;
    QGridLayout *gridLayout_32;
    QPushButton *exportBTBtn;
    QCheckBox *exportPerSectorBTBox;
    QSpacerItem *verticalSpacer;
    QWidget *HSTab;
    QGridLayout *gridLayout_34;
    QGroupBox *mcPruneGroup;
    QGridLayout *gridLayout_15;
    QSlider *pruneMCSlider1;
    QDoubleSpinBox *pruneMCDistSpin1;
    QComboBox *pruneMCDistCombo1;
    QGroupBox *exportHSGroup;
    QGridLayout *gridLayout_35;
    QPushButton *exportHSBtn;
    QGroupBox *hsPruneGroup;
    QGridLayout *gridLayout_18;
    QLabel *label_5;
    QCheckBox *removeSmallCmpnts;
    QSlider *hsBt2Bt3Slider;
    QLabel *label_7;
    QDoubleSpinBox *hsFaceDiffSpin;
    QSlider *hsBt1Bt2Slider;
    QDoubleSpinBox *hsCurveDiffSpin;
    QDoubleSpinBox *hsFaceDegenThreshold;
    QPushButton *visHSBtn;
    QGroupBox *hsCreateGroup;
    QGridLayout *gridLayout_33;
    QPushButton *createHSBtn;
    QSpacerItem *verticalSpacer_8;
    QWidget *MCTab;
    QGridLayout *gridLayout_36;
    QSpacerItem *verticalSpacer_4;
    QGroupBox *mcDualizeGroup;
    QGridLayout *gridLayout_7;
    QPushButton *dualizeBtn;
    QGroupBox *mpGroup;
    QGridLayout *gridLayout_21;
    QCheckBox *showMPBox;
    QWidget *ExportTab;
    QGridLayout *gridLayout_29;
    QGroupBox *exportGroup;
    QGridLayout *gridLayout_13;
    QCheckBox *ExportSkelBox;
    QPushButton *exportBtn;
    QSpacerItem *verticalSpacer_7;
    QWidget *SMetricTab;
    QGridLayout *gridLayout_20;
    QGroupBox *surfFuncGroup;
    QGridLayout *gridLayout_19;
    QComboBox *pruneSMetricCombo;
    QPushButton *surfFuncProjectBtn;
    QComboBox *surfFuncCorrespSchemeCombo;
    QSlider *smoothSMetricSlider;
    QCheckBox *diffuseSurfFunc;
    QDoubleSpinBox *surfFuncRadiusEps;
    QPushButton *surfFuncSetupCorrespBtn;
    QSpinBox *surfFuncKNN;
    QSpacerItem *verticalSpacer_5;
    QWidget *VideoTab;
    QGridLayout *gridLayout_23;
    QGroupBox *isoSurfGroup;
    QGridLayout *gridLayout_24;
    QPushButton *isoSurfPrecomputeBtn;
    QCheckBox *hideIsoSurf;
    QCheckBox *recursiveRefineBox;
    QSlider *isoSurfSlider;
    QCheckBox *hideSnappedFaceBox;
    QPushButton *isoSurfRefineBtn;
    QGroupBox *groupBox_7;
    QGridLayout *gridLayout_25;
    QPushButton *isoContPrecomputeBtn;
    QCheckBox *hideIsoCont;
    QCheckBox *isoContShowMC;
    QCheckBox *hidePastSurfBox;
    QSlider *isoContSlider;
    QGroupBox *groupBox_9;
    QGridLayout *gridLayout_27;
    QSlider *isoContMinAlphaSlider;
    QSlider *isoContExpAlphaSlider;
    QSpacerItem *verticalSpacer_6;
    QGroupBox *groupBox_8;
    QGridLayout *gridLayout_26;
    QSlider *evolveAllSlider;
    QPushButton *evolveAllPrecomputeBtn;
    QCheckBox *evolveAllShowEntireMC;
    QWidget *MiscTab;
    QGridLayout *gridLayout_31;
    QSpacerItem *verticalSpacer_2;
    QGroupBox *preprocessGroup;
    QGridLayout *gridLayout;
    QLabel *label_2;
    QDoubleSpinBox *pertSpin;
    QLabel *label;
    QDoubleSpinBox *epsilonSpin;
    QGroupBox *miscBurntimeGroup;
    QGridLayout *gridLayout_28;
    QPushButton *visSphereFromFileBtn;
    QSpinBox *visSphereSpin;
    QPushButton *randBTBtn;
    QPushButton *printBTBtn;
    QSlider *visSphereSlider;
    QCheckBox *showSphereBox;
    QPushButton *printSelectVertInfoBtn;
    QGroupBox *renderParamGroup;
    QFormLayout *formLayout_3;
    QLabel *label_3;
    QSlider *linesRenderOffsetSlider;
    QGroupBox *qmatGroup;
    QFormLayout *formLayout_4;
    QLabel *label_9;
    QPushButton *outputWenpingBtn;
    QPushButton *readQMatMA;
    QCheckBox *hideQMatBox;
    QWidget *hiddenTab;
    QComboBox *stSubdivCombo;
    QComboBox *edgeWeightCombo;
    QDoubleSpinBox *pointSizeSpin;
    QCheckBox *doMAFaceScalarDiffusion;
    QPushButton *computeAllDistBtn;
    QCheckBox *usePerSheetBox;
    QDoubleSpinBox *maxPruneDistMASpin2;
    QDoubleSpinBox *maPruneRatio;
    QDoubleSpinBox *minPruneDistMASpin1;
    QDoubleSpinBox *pruneMASpin2;
    QDoubleSpinBox *maxPruneDistMASpin1;
    QSlider *pruneMASlider2;
    QComboBox *pruneMADistCombo2;
    QDoubleSpinBox *minPruneDistMASpin2;
    QComboBox *edgeDualOptCombo;
    QComboBox *polyDualOptCombo;
    QCheckBox *underEstimateBox;
    QCheckBox *onlyUnburntBox;
    QCheckBox *stopAtJuncBox;
    QCheckBox *protectBT2Box;
    QSlider *hsBt1Bt2RelSlider;
    QDoubleSpinBox *hsCurveRelDiffSpin;
    QDoubleSpinBox *hsFaceRelDiffSpin;
    QSlider *hsBt2Bt3RelSlider;
    QLabel *label_8;
    QLabel *label_6;
    QDoubleSpinBox *hsCmpntFaceNumTSpin;
    QGroupBox *mcBurnGroup;
    QGridLayout *gridLayout_14;
    QPushButton *burnMedialCurveBtn;
    QPushButton *printMCStatsBtn;
    QSlider *pruneMCSlider2;
    QDoubleSpinBox *mcPruneRatio;
    QComboBox *pruneMCDistCombo2;
    QDoubleSpinBox *pruneMCDistSpin2;
    QRadioButton *colorPerVert;
    QRadioButton *colorPerFace;
    QCheckBox *hsUseTypedInput;
    QCheckBox *visBallStick;
    QCheckBox *preserveMCTopo;
    QWidget *layoutWidget;
    QHBoxLayout *glLayout;
    QSpacerItem *horizontalSpacer;
    QStatusBar *statusbar;

    void setupUi(QMainWindow *MainWindow)
    {
        if (MainWindow->objectName().isEmpty())
            MainWindow->setObjectName(QStringLiteral("MainWindow"));
        MainWindow->resize(1002, 846);
        QSizePolicy sizePolicy(QSizePolicy::Minimum, QSizePolicy::Minimum);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(MainWindow->sizePolicy().hasHeightForWidth());
        MainWindow->setSizePolicy(sizePolicy);
        MainWindow->setMinimumSize(QSize(800, 600));
        openFile = new QAction(MainWindow);
        openFile->setObjectName(QStringLiteral("openFile"));
        saveView = new QAction(MainWindow);
        saveView->setObjectName(QStringLiteral("saveView"));
        loadView = new QAction(MainWindow);
        loadView->setObjectName(QStringLiteral("loadView"));
        readRadii = new QAction(MainWindow);
        readRadii->setObjectName(QStringLiteral("readRadii"));
        actionHelp = new QAction(MainWindow);
        actionHelp->setObjectName(QStringLiteral("actionHelp"));
        centralwidget = new QWidget(MainWindow);
        centralwidget->setObjectName(QStringLiteral("centralwidget"));
        gridLayout_11 = new QGridLayout(centralwidget);
        gridLayout_11->setObjectName(QStringLiteral("gridLayout_11"));
        splitter = new QSplitter(centralwidget);
        splitter->setObjectName(QStringLiteral("splitter"));
        splitter->setOrientation(Qt::Horizontal);
        toolTabs = new QTabWidget(splitter);
        toolTabs->setObjectName(QStringLiteral("toolTabs"));
        toolTabs->setEnabled(true);
        sizePolicy.setHeightForWidth(toolTabs->sizePolicy().hasHeightForWidth());
        toolTabs->setSizePolicy(sizePolicy);
        toolTabs->setMaximumSize(QSize(280, 16777215));
        fileTab = new QWidget();
        fileTab->setObjectName(QStringLiteral("fileTab"));
        gridLayout_17 = new QGridLayout(fileTab);
        gridLayout_17->setObjectName(QStringLiteral("gridLayout_17"));
        loadFileGroup = new QGroupBox(fileTab);
        loadFileGroup->setObjectName(QStringLiteral("loadFileGroup"));
        sizePolicy.setHeightForWidth(loadFileGroup->sizePolicy().hasHeightForWidth());
        loadFileGroup->setSizePolicy(sizePolicy);
        gridLayout_40 = new QGridLayout(loadFileGroup);
        gridLayout_40->setObjectName(QStringLiteral("gridLayout_40"));
        label_11 = new QLabel(loadFileGroup);
        label_11->setObjectName(QStringLiteral("label_11"));

        gridLayout_40->addWidget(label_11, 0, 0, 1, 1);

        maFileEdit = new QLineEdit(loadFileGroup);
        maFileEdit->setObjectName(QStringLiteral("maFileEdit"));

        gridLayout_40->addWidget(maFileEdit, 0, 1, 1, 1);

        browseMAFileBtn = new QPushButton(loadFileGroup);
        browseMAFileBtn->setObjectName(QStringLiteral("browseMAFileBtn"));

        gridLayout_40->addWidget(browseMAFileBtn, 0, 2, 1, 1);

        label_12 = new QLabel(loadFileGroup);
        label_12->setObjectName(QStringLiteral("label_12"));

        gridLayout_40->addWidget(label_12, 1, 0, 1, 1);

        shapeFileEdit = new QLineEdit(loadFileGroup);
        shapeFileEdit->setObjectName(QStringLiteral("shapeFileEdit"));

        gridLayout_40->addWidget(shapeFileEdit, 1, 1, 1, 1);

        browse3DShapeFileBtn = new QPushButton(loadFileGroup);
        browse3DShapeFileBtn->setObjectName(QStringLiteral("browse3DShapeFileBtn"));

        gridLayout_40->addWidget(browse3DShapeFileBtn, 1, 2, 1, 1);

        label_13 = new QLabel(loadFileGroup);
        label_13->setObjectName(QStringLiteral("label_13"));

        gridLayout_40->addWidget(label_13, 2, 0, 1, 1);

        radiusFileEdit = new QLineEdit(loadFileGroup);
        radiusFileEdit->setObjectName(QStringLiteral("radiusFileEdit"));

        gridLayout_40->addWidget(radiusFileEdit, 2, 1, 1, 1);

        browseRadiusFileBtn = new QPushButton(loadFileGroup);
        browseRadiusFileBtn->setObjectName(QStringLiteral("browseRadiusFileBtn"));

        gridLayout_40->addWidget(browseRadiusFileBtn, 2, 2, 1, 1);

        loadFilesBtn = new QPushButton(loadFileGroup);
        loadFilesBtn->setObjectName(QStringLiteral("loadFilesBtn"));

        gridLayout_40->addWidget(loadFilesBtn, 3, 2, 1, 1);


        gridLayout_17->addWidget(loadFileGroup, 0, 0, 1, 1);

        verticalSpacer_3 = new QSpacerItem(20, 40, QSizePolicy::Minimum, QSizePolicy::Expanding);

        gridLayout_17->addItem(verticalSpacer_3, 1, 0, 1, 1);

        toolTabs->addTab(fileTab, QString());
        origSurfTab = new QWidget();
        origSurfTab->setObjectName(QStringLiteral("origSurfTab"));
        gridLayout_2 = new QGridLayout(origSurfTab);
        gridLayout_2->setObjectName(QStringLiteral("gridLayout_2"));
        groupBox = new QGroupBox(origSurfTab);
        groupBox->setObjectName(QStringLiteral("groupBox"));
        gridLayout_30 = new QGridLayout(groupBox);
        gridLayout_30->setObjectName(QStringLiteral("gridLayout_30"));
        DPMaxRendersSpin = new QSpinBox(groupBox);
        DPMaxRendersSpin->setObjectName(QStringLiteral("DPMaxRendersSpin"));
        DPMaxRendersSpin->setMaximum(99999);

        gridLayout_30->addWidget(DPMaxRendersSpin, 3, 1, 1, 1);

        colorBGBtn = new QPushButton(groupBox);
        colorBGBtn->setObjectName(QStringLiteral("colorBGBtn"));

        gridLayout_30->addWidget(colorBGBtn, 0, 0, 1, 1);

        drawEdge = new QCheckBox(groupBox);
        drawEdge->setObjectName(QStringLiteral("drawEdge"));
        QSizePolicy sizePolicy1(QSizePolicy::Preferred, QSizePolicy::Fixed);
        sizePolicy1.setHorizontalStretch(0);
        sizePolicy1.setVerticalStretch(0);
        sizePolicy1.setHeightForWidth(drawEdge->sizePolicy().hasHeightForWidth());
        drawEdge->setSizePolicy(sizePolicy1);

        gridLayout_30->addWidget(drawEdge, 0, 1, 1, 1);

        saveViewBtn = new QPushButton(groupBox);
        saveViewBtn->setObjectName(QStringLiteral("saveViewBtn"));

        gridLayout_30->addWidget(saveViewBtn, 1, 0, 1, 1);

        loadViewBtn = new QPushButton(groupBox);
        loadViewBtn->setObjectName(QStringLiteral("loadViewBtn"));

        gridLayout_30->addWidget(loadViewBtn, 1, 1, 1, 1);

        useTrueTransparencyBox = new QCheckBox(groupBox);
        useTrueTransparencyBox->setObjectName(QStringLiteral("useTrueTransparencyBox"));

        gridLayout_30->addWidget(useTrueTransparencyBox, 3, 0, 1, 1);


        gridLayout_2->addWidget(groupBox, 0, 0, 1, 1);

        origSurfViewGroup = new QGroupBox(origSurfTab);
        origSurfViewGroup->setObjectName(QStringLiteral("origSurfViewGroup"));
        QSizePolicy sizePolicy2(QSizePolicy::Preferred, QSizePolicy::Preferred);
        sizePolicy2.setHorizontalStretch(0);
        sizePolicy2.setVerticalStretch(0);
        sizePolicy2.setHeightForWidth(origSurfViewGroup->sizePolicy().hasHeightForWidth());
        origSurfViewGroup->setSizePolicy(sizePolicy2);
        origSurfViewGroup->setMinimumSize(QSize(250, 0));
        origSurfViewGroup->setMaximumSize(QSize(200, 16777215));
        gridLayout_16 = new QGridLayout(origSurfViewGroup);
        gridLayout_16->setObjectName(QStringLiteral("gridLayout_16"));
        colorOrigBtn = new QPushButton(origSurfViewGroup);
        colorOrigBtn->setObjectName(QStringLiteral("colorOrigBtn"));
        QSizePolicy sizePolicy3(QSizePolicy::Expanding, QSizePolicy::Fixed);
        sizePolicy3.setHorizontalStretch(0);
        sizePolicy3.setVerticalStretch(0);
        sizePolicy3.setHeightForWidth(colorOrigBtn->sizePolicy().hasHeightForWidth());
        colorOrigBtn->setSizePolicy(sizePolicy3);
        colorOrigBtn->setMaximumSize(QSize(500, 23));

        gridLayout_16->addWidget(colorOrigBtn, 0, 1, 1, 1);

        hideOrig = new QCheckBox(origSurfViewGroup);
        hideOrig->setObjectName(QStringLiteral("hideOrig"));
        sizePolicy3.setHeightForWidth(hideOrig->sizePolicy().hasHeightForWidth());
        hideOrig->setSizePolicy(sizePolicy3);

        gridLayout_16->addWidget(hideOrig, 0, 0, 1, 1);

        label_4 = new QLabel(origSurfViewGroup);
        label_4->setObjectName(QStringLiteral("label_4"));
        label_4->setAlignment(Qt::AlignCenter);

        gridLayout_16->addWidget(label_4, 1, 0, 1, 1);

        origTransparentSlider = new QSlider(origSurfViewGroup);
        origTransparentSlider->setObjectName(QStringLiteral("origTransparentSlider"));
        origTransparentSlider->setOrientation(Qt::Horizontal);

        gridLayout_16->addWidget(origTransparentSlider, 1, 1, 1, 1);


        gridLayout_2->addWidget(origSurfViewGroup, 1, 0, 1, 1);

        maViewGroup = new QGroupBox(origSurfTab);
        maViewGroup->setObjectName(QStringLiteral("maViewGroup"));
        maViewGroup->setMinimumSize(QSize(250, 0));
        maViewGroup->setMaximumSize(QSize(200, 16777215));
        gridLayout_6 = new QGridLayout(maViewGroup);
        gridLayout_6->setObjectName(QStringLiteral("gridLayout_6"));
        hideMA = new QCheckBox(maViewGroup);
        hideMA->setObjectName(QStringLiteral("hideMA"));
        QSizePolicy sizePolicy4(QSizePolicy::Minimum, QSizePolicy::Fixed);
        sizePolicy4.setHorizontalStretch(0);
        sizePolicy4.setVerticalStretch(0);
        sizePolicy4.setHeightForWidth(hideMA->sizePolicy().hasHeightForWidth());
        hideMA->setSizePolicy(sizePolicy4);

        gridLayout_6->addWidget(hideMA, 0, 0, 1, 1);

        turnOffLightingForMABox = new QCheckBox(maViewGroup);
        turnOffLightingForMABox->setObjectName(QStringLiteral("turnOffLightingForMABox"));
        QSizePolicy sizePolicy5(QSizePolicy::MinimumExpanding, QSizePolicy::Fixed);
        sizePolicy5.setHorizontalStretch(0);
        sizePolicy5.setVerticalStretch(0);
        sizePolicy5.setHeightForWidth(turnOffLightingForMABox->sizePolicy().hasHeightForWidth());
        turnOffLightingForMABox->setSizePolicy(sizePolicy5);

        gridLayout_6->addWidget(turnOffLightingForMABox, 0, 1, 1, 1);

        drawStPoints = new QCheckBox(maViewGroup);
        drawStPoints->setObjectName(QStringLiteral("drawStPoints"));
        sizePolicy1.setHeightForWidth(drawStPoints->sizePolicy().hasHeightForWidth());
        drawStPoints->setSizePolicy(sizePolicy1);

        gridLayout_6->addWidget(drawStPoints, 1, 0, 1, 1);

        colorMABtn = new QPushButton(maViewGroup);
        colorMABtn->setObjectName(QStringLiteral("colorMABtn"));
        sizePolicy1.setHeightForWidth(colorMABtn->sizePolicy().hasHeightForWidth());
        colorMABtn->setSizePolicy(sizePolicy1);
        colorMABtn->setMaximumSize(QSize(500, 23));

        gridLayout_6->addWidget(colorMABtn, 1, 1, 1, 1);

        distVisBtnGroup = new QGroupBox(maViewGroup);
        distVisBtnGroup->setObjectName(QStringLiteral("distVisBtnGroup"));
        sizePolicy.setHeightForWidth(distVisBtnGroup->sizePolicy().hasHeightForWidth());
        distVisBtnGroup->setSizePolicy(sizePolicy);
        gridLayout_37 = new QGridLayout(distVisBtnGroup);
        gridLayout_37->setObjectName(QStringLiteral("gridLayout_37"));
        visMADistCombo = new QComboBox(distVisBtnGroup);
        visMADistCombo->setObjectName(QStringLiteral("visMADistCombo"));
        sizePolicy4.setHeightForWidth(visMADistCombo->sizePolicy().hasHeightForWidth());
        visMADistCombo->setSizePolicy(sizePolicy4);

        gridLayout_37->addWidget(visMADistCombo, 0, 0, 1, 1);

        gridLayout_38 = new QGridLayout();
        gridLayout_38->setObjectName(QStringLiteral("gridLayout_38"));
        clampMADistBox = new QCheckBox(distVisBtnGroup);
        clampMADistBox->setObjectName(QStringLiteral("clampMADistBox"));
        sizePolicy3.setHeightForWidth(clampMADistBox->sizePolicy().hasHeightForWidth());
        clampMADistBox->setSizePolicy(sizePolicy3);
        clampMADistBox->setMaximumSize(QSize(100, 16777215));

        gridLayout_38->addWidget(clampMADistBox, 0, 0, 1, 1);

        minVisDistMASpin = new QDoubleSpinBox(distVisBtnGroup);
        minVisDistMASpin->setObjectName(QStringLiteral("minVisDistMASpin"));
        sizePolicy3.setHeightForWidth(minVisDistMASpin->sizePolicy().hasHeightForWidth());
        minVisDistMASpin->setSizePolicy(sizePolicy3);
        minVisDistMASpin->setDecimals(8);

        gridLayout_38->addWidget(minVisDistMASpin, 0, 1, 1, 1);

        maxVisDistMASpin = new QDoubleSpinBox(distVisBtnGroup);
        maxVisDistMASpin->setObjectName(QStringLiteral("maxVisDistMASpin"));
        sizePolicy3.setHeightForWidth(maxVisDistMASpin->sizePolicy().hasHeightForWidth());
        maxVisDistMASpin->setSizePolicy(sizePolicy3);
        maxVisDistMASpin->setMaximumSize(QSize(500, 16777215));
        maxVisDistMASpin->setDecimals(8);
        maxVisDistMASpin->setMaximum(1e+7);

        gridLayout_38->addWidget(maxVisDistMASpin, 0, 2, 1, 1);


        gridLayout_37->addLayout(gridLayout_38, 1, 0, 1, 1);


        gridLayout_6->addWidget(distVisBtnGroup, 2, 0, 1, 2);

        MAAlphaGroup = new QGroupBox(maViewGroup);
        MAAlphaGroup->setObjectName(QStringLiteral("MAAlphaGroup"));
        gridLayout_10 = new QGridLayout(MAAlphaGroup);
        gridLayout_10->setObjectName(QStringLiteral("gridLayout_10"));
        MATransparentSlider = new QSlider(MAAlphaGroup);
        MATransparentSlider->setObjectName(QStringLiteral("MATransparentSlider"));
        MATransparentSlider->setMaximum(1000);
        MATransparentSlider->setOrientation(Qt::Horizontal);

        gridLayout_10->addWidget(MATransparentSlider, 0, 0, 1, 1);

        MATransExpSlider = new QSlider(MAAlphaGroup);
        MATransExpSlider->setObjectName(QStringLiteral("MATransExpSlider"));
        MATransExpSlider->setOrientation(Qt::Horizontal);

        gridLayout_10->addWidget(MATransExpSlider, 0, 1, 1, 1);


        gridLayout_6->addWidget(MAAlphaGroup, 3, 0, 1, 2);


        gridLayout_2->addWidget(maViewGroup, 2, 0, 1, 1);

        mcViewGroup = new QGroupBox(origSurfTab);
        mcViewGroup->setObjectName(QStringLiteral("mcViewGroup"));
        gridLayout_4 = new QGridLayout(mcViewGroup);
        gridLayout_4->setObjectName(QStringLiteral("gridLayout_4"));
        visBurnt = new QCheckBox(mcViewGroup);
        visBurnt->setObjectName(QStringLiteral("visBurnt"));
        sizePolicy4.setHeightForWidth(visBurnt->sizePolicy().hasHeightForWidth());
        visBurnt->setSizePolicy(sizePolicy4);

        gridLayout_4->addWidget(visBurnt, 0, 0, 1, 2);

        colorMCBtn = new QPushButton(mcViewGroup);
        colorMCBtn->setObjectName(QStringLiteral("colorMCBtn"));
        sizePolicy3.setHeightForWidth(colorMCBtn->sizePolicy().hasHeightForWidth());
        colorMCBtn->setSizePolicy(sizePolicy3);

        gridLayout_4->addWidget(colorMCBtn, 0, 2, 1, 1);

        hideMC = new QCheckBox(mcViewGroup);
        hideMC->setObjectName(QStringLiteral("hideMC"));
        sizePolicy4.setHeightForWidth(hideMC->sizePolicy().hasHeightForWidth());
        hideMC->setSizePolicy(sizePolicy4);

        gridLayout_4->addWidget(hideMC, 1, 0, 1, 1);

        label_10 = new QLabel(mcViewGroup);
        label_10->setObjectName(QStringLiteral("label_10"));

        gridLayout_4->addWidget(label_10, 2, 1, 1, 1);

        MCWidthSpin = new QDoubleSpinBox(mcViewGroup);
        MCWidthSpin->setObjectName(QStringLiteral("MCWidthSpin"));
        sizePolicy3.setHeightForWidth(MCWidthSpin->sizePolicy().hasHeightForWidth());
        MCWidthSpin->setSizePolicy(sizePolicy3);

        gridLayout_4->addWidget(MCWidthSpin, 2, 2, 1, 1);

        groupBox_2 = new QGroupBox(mcViewGroup);
        groupBox_2->setObjectName(QStringLiteral("groupBox_2"));
        gridLayout_8 = new QGridLayout(groupBox_2);
        gridLayout_8->setObjectName(QStringLiteral("gridLayout_8"));
        visMCDistCombo = new QComboBox(groupBox_2);
        visMCDistCombo->setObjectName(QStringLiteral("visMCDistCombo"));
        sizePolicy3.setHeightForWidth(visMCDistCombo->sizePolicy().hasHeightForWidth());
        visMCDistCombo->setSizePolicy(sizePolicy3);

        gridLayout_8->addWidget(visMCDistCombo, 0, 0, 1, 1);

        gridLayout_39 = new QGridLayout();
        gridLayout_39->setObjectName(QStringLiteral("gridLayout_39"));
        clampMCDistBox = new QCheckBox(groupBox_2);
        clampMCDistBox->setObjectName(QStringLiteral("clampMCDistBox"));

        gridLayout_39->addWidget(clampMCDistBox, 0, 0, 1, 1);

        minVisDistMCSpin = new QDoubleSpinBox(groupBox_2);
        minVisDistMCSpin->setObjectName(QStringLiteral("minVisDistMCSpin"));
        sizePolicy3.setHeightForWidth(minVisDistMCSpin->sizePolicy().hasHeightForWidth());
        minVisDistMCSpin->setSizePolicy(sizePolicy3);
        minVisDistMCSpin->setDecimals(8);

        gridLayout_39->addWidget(minVisDistMCSpin, 0, 1, 1, 1);

        maxVisDistMCSpin = new QDoubleSpinBox(groupBox_2);
        maxVisDistMCSpin->setObjectName(QStringLiteral("maxVisDistMCSpin"));
        sizePolicy3.setHeightForWidth(maxVisDistMCSpin->sizePolicy().hasHeightForWidth());
        maxVisDistMCSpin->setSizePolicy(sizePolicy3);
        maxVisDistMCSpin->setDecimals(8);

        gridLayout_39->addWidget(maxVisDistMCSpin, 0, 2, 1, 1);


        gridLayout_8->addLayout(gridLayout_39, 1, 0, 1, 1);


        gridLayout_4->addWidget(groupBox_2, 3, 0, 1, 3);

        MCAlphaGroup = new QGroupBox(mcViewGroup);
        MCAlphaGroup->setObjectName(QStringLiteral("MCAlphaGroup"));
        gridLayout_22 = new QGridLayout(MCAlphaGroup);
        gridLayout_22->setObjectName(QStringLiteral("gridLayout_22"));
        MCAlphaSlider = new QSlider(MCAlphaGroup);
        MCAlphaSlider->setObjectName(QStringLiteral("MCAlphaSlider"));
        MCAlphaSlider->setMaximum(1000);
        MCAlphaSlider->setOrientation(Qt::Horizontal);

        gridLayout_22->addWidget(MCAlphaSlider, 1, 0, 1, 1);

        MCTransExpSlider = new QSlider(MCAlphaGroup);
        MCTransExpSlider->setObjectName(QStringLiteral("MCTransExpSlider"));
        MCTransExpSlider->setOrientation(Qt::Horizontal);

        gridLayout_22->addWidget(MCTransExpSlider, 1, 1, 1, 1);


        gridLayout_4->addWidget(MCAlphaGroup, 4, 0, 1, 3);


        gridLayout_2->addWidget(mcViewGroup, 3, 0, 1, 1);

        hsViewGroup = new QGroupBox(origSurfTab);
        hsViewGroup->setObjectName(QStringLiteral("hsViewGroup"));
        sizePolicy.setHeightForWidth(hsViewGroup->sizePolicy().hasHeightForWidth());
        hsViewGroup->setSizePolicy(sizePolicy);
        gridLayout_9 = new QGridLayout(hsViewGroup);
        gridLayout_9->setObjectName(QStringLiteral("gridLayout_9"));
        hideHS = new QCheckBox(hsViewGroup);
        hideHS->setObjectName(QStringLiteral("hideHS"));

        gridLayout_9->addWidget(hideHS, 0, 0, 1, 1);


        gridLayout_2->addWidget(hsViewGroup, 4, 0, 1, 1);

        toolTabs->addTab(origSurfTab, QString());
        MATab = new QWidget();
        MATab->setObjectName(QStringLiteral("MATab"));
        gridLayout_12 = new QGridLayout(MATab);
        gridLayout_12->setObjectName(QStringLiteral("gridLayout_12"));
        maBurnGroup = new QGroupBox(MATab);
        maBurnGroup->setObjectName(QStringLiteral("maBurnGroup"));
        QSizePolicy sizePolicy6(QSizePolicy::Expanding, QSizePolicy::Minimum);
        sizePolicy6.setHorizontalStretch(0);
        sizePolicy6.setVerticalStretch(0);
        sizePolicy6.setHeightForWidth(maBurnGroup->sizePolicy().hasHeightForWidth());
        maBurnGroup->setSizePolicy(sizePolicy6);
        gridLayout_5 = new QGridLayout(maBurnGroup);
        gridLayout_5->setObjectName(QStringLiteral("gridLayout_5"));
        nFixedSteinerSpin = new QDoubleSpinBox(maBurnGroup);
        nFixedSteinerSpin->setObjectName(QStringLiteral("nFixedSteinerSpin"));
        sizePolicy4.setHeightForWidth(nFixedSteinerSpin->sizePolicy().hasHeightForWidth());
        nFixedSteinerSpin->setSizePolicy(sizePolicy4);
        nFixedSteinerSpin->setDecimals(8);

        gridLayout_5->addWidget(nFixedSteinerSpin, 0, 0, 1, 1);

        cleanTopoBtn = new QPushButton(maBurnGroup);
        cleanTopoBtn->setObjectName(QStringLiteral("cleanTopoBtn"));
        sizePolicy4.setHeightForWidth(cleanTopoBtn->sizePolicy().hasHeightForWidth());
        cleanTopoBtn->setSizePolicy(sizePolicy4);
        cleanTopoBtn->setMinimumSize(QSize(40, 23));
        cleanTopoBtn->setMaximumSize(QSize(100, 23));

        gridLayout_5->addWidget(cleanTopoBtn, 0, 1, 1, 1);

        burnBtn = new QPushButton(maBurnGroup);
        burnBtn->setObjectName(QStringLiteral("burnBtn"));
        burnBtn->setEnabled(true);
        sizePolicy4.setHeightForWidth(burnBtn->sizePolicy().hasHeightForWidth());
        burnBtn->setSizePolicy(sizePolicy4);
        burnBtn->setMinimumSize(QSize(40, 23));
        burnBtn->setMaximumSize(QSize(100, 23));

        gridLayout_5->addWidget(burnBtn, 1, 0, 1, 1);


        gridLayout_12->addWidget(maBurnGroup, 0, 0, 1, 1);

        maMeasureGroup = new QGroupBox(MATab);
        maMeasureGroup->setObjectName(QStringLiteral("maMeasureGroup"));
        sizePolicy.setHeightForWidth(maMeasureGroup->sizePolicy().hasHeightForWidth());
        maMeasureGroup->setSizePolicy(sizePolicy);
        maMeasureGroup->setMaximumSize(QSize(250, 16777215));
        gridLayout_3 = new QGridLayout(maMeasureGroup);
        gridLayout_3->setObjectName(QStringLiteral("gridLayout_3"));
        pruneMASlider1 = new QSlider(maMeasureGroup);
        pruneMASlider1->setObjectName(QStringLiteral("pruneMASlider1"));
        sizePolicy3.setHeightForWidth(pruneMASlider1->sizePolicy().hasHeightForWidth());
        pruneMASlider1->setSizePolicy(sizePolicy3);
        pruneMASlider1->setMaximumSize(QSize(100, 16777215));
        pruneMASlider1->setMaximum(1000000);
        pruneMASlider1->setOrientation(Qt::Horizontal);

        gridLayout_3->addWidget(pruneMASlider1, 0, 0, 1, 1);

        pruneMADistCombo1 = new QComboBox(maMeasureGroup);
        pruneMADistCombo1->setObjectName(QStringLiteral("pruneMADistCombo1"));
        sizePolicy3.setHeightForWidth(pruneMADistCombo1->sizePolicy().hasHeightForWidth());
        pruneMADistCombo1->setSizePolicy(sizePolicy3);

        gridLayout_3->addWidget(pruneMADistCombo1, 0, 1, 1, 1);

        pruneMASpin1 = new QDoubleSpinBox(maMeasureGroup);
        pruneMASpin1->setObjectName(QStringLiteral("pruneMASpin1"));
        pruneMASpin1->setMaximumSize(QSize(100, 16777215));
        pruneMASpin1->setDecimals(8);
        pruneMASpin1->setMaximum(99999);

        gridLayout_3->addWidget(pruneMASpin1, 1, 0, 1, 1);

        enableFinePruneMA = new QCheckBox(maMeasureGroup);
        enableFinePruneMA->setObjectName(QStringLiteral("enableFinePruneMA"));

        gridLayout_3->addWidget(enableFinePruneMA, 1, 1, 1, 1);


        gridLayout_12->addWidget(maMeasureGroup, 1, 0, 1, 1);

        exportBTGroup = new QGroupBox(MATab);
        exportBTGroup->setObjectName(QStringLiteral("exportBTGroup"));
        gridLayout_32 = new QGridLayout(exportBTGroup);
        gridLayout_32->setObjectName(QStringLiteral("gridLayout_32"));
        exportBTBtn = new QPushButton(exportBTGroup);
        exportBTBtn->setObjectName(QStringLiteral("exportBTBtn"));

        gridLayout_32->addWidget(exportBTBtn, 2, 1, 1, 1);

        exportPerSectorBTBox = new QCheckBox(exportBTGroup);
        exportPerSectorBTBox->setObjectName(QStringLiteral("exportPerSectorBTBox"));

        gridLayout_32->addWidget(exportPerSectorBTBox, 2, 0, 1, 1);


        gridLayout_12->addWidget(exportBTGroup, 2, 0, 1, 1);

        verticalSpacer = new QSpacerItem(20, 40, QSizePolicy::Minimum, QSizePolicy::Expanding);

        gridLayout_12->addItem(verticalSpacer, 3, 0, 1, 1);

        toolTabs->addTab(MATab, QString());
        HSTab = new QWidget();
        HSTab->setObjectName(QStringLiteral("HSTab"));
        gridLayout_34 = new QGridLayout(HSTab);
        gridLayout_34->setObjectName(QStringLiteral("gridLayout_34"));
        mcPruneGroup = new QGroupBox(HSTab);
        mcPruneGroup->setObjectName(QStringLiteral("mcPruneGroup"));
        gridLayout_15 = new QGridLayout(mcPruneGroup);
        gridLayout_15->setObjectName(QStringLiteral("gridLayout_15"));
        pruneMCSlider1 = new QSlider(mcPruneGroup);
        pruneMCSlider1->setObjectName(QStringLiteral("pruneMCSlider1"));
        sizePolicy3.setHeightForWidth(pruneMCSlider1->sizePolicy().hasHeightForWidth());
        pruneMCSlider1->setSizePolicy(sizePolicy3);
        pruneMCSlider1->setOrientation(Qt::Horizontal);

        gridLayout_15->addWidget(pruneMCSlider1, 2, 0, 1, 1);

        pruneMCDistSpin1 = new QDoubleSpinBox(mcPruneGroup);
        pruneMCDistSpin1->setObjectName(QStringLiteral("pruneMCDistSpin1"));
        sizePolicy3.setHeightForWidth(pruneMCDistSpin1->sizePolicy().hasHeightForWidth());
        pruneMCDistSpin1->setSizePolicy(sizePolicy3);
        pruneMCDistSpin1->setDecimals(8);

        gridLayout_15->addWidget(pruneMCDistSpin1, 4, 0, 1, 1);

        pruneMCDistCombo1 = new QComboBox(mcPruneGroup);
        pruneMCDistCombo1->setObjectName(QStringLiteral("pruneMCDistCombo1"));
        sizePolicy3.setHeightForWidth(pruneMCDistCombo1->sizePolicy().hasHeightForWidth());
        pruneMCDistCombo1->setSizePolicy(sizePolicy3);
        pruneMCDistCombo1->setMinimumSize(QSize(0, 0));

        gridLayout_15->addWidget(pruneMCDistCombo1, 2, 1, 1, 1);


        gridLayout_34->addWidget(mcPruneGroup, 1, 0, 1, 2);

        exportHSGroup = new QGroupBox(HSTab);
        exportHSGroup->setObjectName(QStringLiteral("exportHSGroup"));
        gridLayout_35 = new QGridLayout(exportHSGroup);
        gridLayout_35->setObjectName(QStringLiteral("gridLayout_35"));
        exportHSBtn = new QPushButton(exportHSGroup);
        exportHSBtn->setObjectName(QStringLiteral("exportHSBtn"));

        gridLayout_35->addWidget(exportHSBtn, 0, 0, 1, 1);


        gridLayout_34->addWidget(exportHSGroup, 3, 0, 1, 1);

        hsPruneGroup = new QGroupBox(HSTab);
        hsPruneGroup->setObjectName(QStringLiteral("hsPruneGroup"));
        gridLayout_18 = new QGridLayout(hsPruneGroup);
        gridLayout_18->setObjectName(QStringLiteral("gridLayout_18"));
        label_5 = new QLabel(hsPruneGroup);
        label_5->setObjectName(QStringLiteral("label_5"));
        QSizePolicy sizePolicy7(QSizePolicy::Minimum, QSizePolicy::Preferred);
        sizePolicy7.setHorizontalStretch(0);
        sizePolicy7.setVerticalStretch(0);
        sizePolicy7.setHeightForWidth(label_5->sizePolicy().hasHeightForWidth());
        label_5->setSizePolicy(sizePolicy7);
        label_5->setMaximumSize(QSize(200, 16777215));

        gridLayout_18->addWidget(label_5, 0, 0, 1, 1);

        removeSmallCmpnts = new QCheckBox(hsPruneGroup);
        removeSmallCmpnts->setObjectName(QStringLiteral("removeSmallCmpnts"));

        gridLayout_18->addWidget(removeSmallCmpnts, 7, 0, 1, 1);

        hsBt2Bt3Slider = new QSlider(hsPruneGroup);
        hsBt2Bt3Slider->setObjectName(QStringLiteral("hsBt2Bt3Slider"));
        sizePolicy3.setHeightForWidth(hsBt2Bt3Slider->sizePolicy().hasHeightForWidth());
        hsBt2Bt3Slider->setSizePolicy(sizePolicy3);
        hsBt2Bt3Slider->setMaximumSize(QSize(200, 16777215));
        hsBt2Bt3Slider->setMaximum(1000);
        hsBt2Bt3Slider->setOrientation(Qt::Horizontal);

        gridLayout_18->addWidget(hsBt2Bt3Slider, 1, 0, 1, 1);

        label_7 = new QLabel(hsPruneGroup);
        label_7->setObjectName(QStringLiteral("label_7"));
        sizePolicy7.setHeightForWidth(label_7->sizePolicy().hasHeightForWidth());
        label_7->setSizePolicy(sizePolicy7);
        label_7->setMaximumSize(QSize(200, 16777215));

        gridLayout_18->addWidget(label_7, 0, 1, 1, 1);

        hsFaceDiffSpin = new QDoubleSpinBox(hsPruneGroup);
        hsFaceDiffSpin->setObjectName(QStringLiteral("hsFaceDiffSpin"));
        hsFaceDiffSpin->setDecimals(8);

        gridLayout_18->addWidget(hsFaceDiffSpin, 3, 0, 1, 1);

        hsBt1Bt2Slider = new QSlider(hsPruneGroup);
        hsBt1Bt2Slider->setObjectName(QStringLiteral("hsBt1Bt2Slider"));
        sizePolicy3.setHeightForWidth(hsBt1Bt2Slider->sizePolicy().hasHeightForWidth());
        hsBt1Bt2Slider->setSizePolicy(sizePolicy3);
        hsBt1Bt2Slider->setMaximumSize(QSize(200, 16777215));
        hsBt1Bt2Slider->setMaximum(1000);
        hsBt1Bt2Slider->setOrientation(Qt::Horizontal);

        gridLayout_18->addWidget(hsBt1Bt2Slider, 1, 1, 1, 1);

        hsCurveDiffSpin = new QDoubleSpinBox(hsPruneGroup);
        hsCurveDiffSpin->setObjectName(QStringLiteral("hsCurveDiffSpin"));
        hsCurveDiffSpin->setDecimals(8);

        gridLayout_18->addWidget(hsCurveDiffSpin, 3, 1, 1, 1);

        hsFaceDegenThreshold = new QDoubleSpinBox(hsPruneGroup);
        hsFaceDegenThreshold->setObjectName(QStringLiteral("hsFaceDegenThreshold"));
        sizePolicy5.setHeightForWidth(hsFaceDegenThreshold->sizePolicy().hasHeightForWidth());
        hsFaceDegenThreshold->setSizePolicy(sizePolicy5);
        hsFaceDegenThreshold->setDecimals(8);

        gridLayout_18->addWidget(hsFaceDegenThreshold, 10, 0, 1, 1);

        visHSBtn = new QPushButton(hsPruneGroup);
        visHSBtn->setObjectName(QStringLiteral("visHSBtn"));
        sizePolicy5.setHeightForWidth(visHSBtn->sizePolicy().hasHeightForWidth());
        visHSBtn->setSizePolicy(sizePolicy5);
        visHSBtn->setMaximumSize(QSize(150, 16777215));

        gridLayout_18->addWidget(visHSBtn, 7, 1, 1, 1);


        gridLayout_34->addWidget(hsPruneGroup, 2, 0, 1, 2);

        hsCreateGroup = new QGroupBox(HSTab);
        hsCreateGroup->setObjectName(QStringLiteral("hsCreateGroup"));
        QSizePolicy sizePolicy8(QSizePolicy::Expanding, QSizePolicy::Preferred);
        sizePolicy8.setHorizontalStretch(0);
        sizePolicy8.setVerticalStretch(0);
        sizePolicy8.setHeightForWidth(hsCreateGroup->sizePolicy().hasHeightForWidth());
        hsCreateGroup->setSizePolicy(sizePolicy8);
        gridLayout_33 = new QGridLayout(hsCreateGroup);
        gridLayout_33->setObjectName(QStringLiteral("gridLayout_33"));
        createHSBtn = new QPushButton(hsCreateGroup);
        createHSBtn->setObjectName(QStringLiteral("createHSBtn"));
        sizePolicy3.setHeightForWidth(createHSBtn->sizePolicy().hasHeightForWidth());
        createHSBtn->setSizePolicy(sizePolicy3);

        gridLayout_33->addWidget(createHSBtn, 0, 0, 1, 1);


        gridLayout_34->addWidget(hsCreateGroup, 0, 0, 1, 1);

        verticalSpacer_8 = new QSpacerItem(20, 365, QSizePolicy::Minimum, QSizePolicy::Expanding);

        gridLayout_34->addItem(verticalSpacer_8, 4, 1, 1, 1);

        toolTabs->addTab(HSTab, QString());
        MCTab = new QWidget();
        MCTab->setObjectName(QStringLiteral("MCTab"));
        gridLayout_36 = new QGridLayout(MCTab);
        gridLayout_36->setObjectName(QStringLiteral("gridLayout_36"));
        verticalSpacer_4 = new QSpacerItem(20, 40, QSizePolicy::Minimum, QSizePolicy::Expanding);

        gridLayout_36->addItem(verticalSpacer_4, 2, 0, 1, 1);

        mcDualizeGroup = new QGroupBox(MCTab);
        mcDualizeGroup->setObjectName(QStringLiteral("mcDualizeGroup"));
        gridLayout_7 = new QGridLayout(mcDualizeGroup);
        gridLayout_7->setObjectName(QStringLiteral("gridLayout_7"));
        dualizeBtn = new QPushButton(mcDualizeGroup);
        dualizeBtn->setObjectName(QStringLiteral("dualizeBtn"));

        gridLayout_7->addWidget(dualizeBtn, 0, 0, 1, 1);


        gridLayout_36->addWidget(mcDualizeGroup, 0, 0, 1, 1);

        mpGroup = new QGroupBox(MCTab);
        mpGroup->setObjectName(QStringLiteral("mpGroup"));
        gridLayout_21 = new QGridLayout(mpGroup);
        gridLayout_21->setObjectName(QStringLiteral("gridLayout_21"));
        showMPBox = new QCheckBox(mpGroup);
        showMPBox->setObjectName(QStringLiteral("showMPBox"));

        gridLayout_21->addWidget(showMPBox, 0, 0, 1, 1);


        gridLayout_36->addWidget(mpGroup, 1, 0, 1, 1);

        toolTabs->addTab(MCTab, QString());
        ExportTab = new QWidget();
        ExportTab->setObjectName(QStringLiteral("ExportTab"));
        gridLayout_29 = new QGridLayout(ExportTab);
        gridLayout_29->setObjectName(QStringLiteral("gridLayout_29"));
        exportGroup = new QGroupBox(ExportTab);
        exportGroup->setObjectName(QStringLiteral("exportGroup"));
        gridLayout_13 = new QGridLayout(exportGroup);
        gridLayout_13->setObjectName(QStringLiteral("gridLayout_13"));
        ExportSkelBox = new QCheckBox(exportGroup);
        ExportSkelBox->setObjectName(QStringLiteral("ExportSkelBox"));

        gridLayout_13->addWidget(ExportSkelBox, 1, 0, 1, 1);

        exportBtn = new QPushButton(exportGroup);
        exportBtn->setObjectName(QStringLiteral("exportBtn"));

        gridLayout_13->addWidget(exportBtn, 0, 0, 1, 1);


        gridLayout_29->addWidget(exportGroup, 0, 0, 1, 1);

        verticalSpacer_7 = new QSpacerItem(20, 40, QSizePolicy::Minimum, QSizePolicy::Expanding);

        gridLayout_29->addItem(verticalSpacer_7, 1, 0, 1, 1);

        toolTabs->addTab(ExportTab, QString());
        SMetricTab = new QWidget();
        SMetricTab->setObjectName(QStringLiteral("SMetricTab"));
        gridLayout_20 = new QGridLayout(SMetricTab);
        gridLayout_20->setObjectName(QStringLiteral("gridLayout_20"));
        surfFuncGroup = new QGroupBox(SMetricTab);
        surfFuncGroup->setObjectName(QStringLiteral("surfFuncGroup"));
        gridLayout_19 = new QGridLayout(surfFuncGroup);
        gridLayout_19->setObjectName(QStringLiteral("gridLayout_19"));
        pruneSMetricCombo = new QComboBox(surfFuncGroup);
        pruneSMetricCombo->setObjectName(QStringLiteral("pruneSMetricCombo"));
        sizePolicy5.setHeightForWidth(pruneSMetricCombo->sizePolicy().hasHeightForWidth());
        pruneSMetricCombo->setSizePolicy(sizePolicy5);

        gridLayout_19->addWidget(pruneSMetricCombo, 5, 0, 1, 1);

        surfFuncProjectBtn = new QPushButton(surfFuncGroup);
        surfFuncProjectBtn->setObjectName(QStringLiteral("surfFuncProjectBtn"));
        sizePolicy5.setHeightForWidth(surfFuncProjectBtn->sizePolicy().hasHeightForWidth());
        surfFuncProjectBtn->setSizePolicy(sizePolicy5);

        gridLayout_19->addWidget(surfFuncProjectBtn, 6, 1, 1, 1);

        surfFuncCorrespSchemeCombo = new QComboBox(surfFuncGroup);
        surfFuncCorrespSchemeCombo->setObjectName(QStringLiteral("surfFuncCorrespSchemeCombo"));
        sizePolicy5.setHeightForWidth(surfFuncCorrespSchemeCombo->sizePolicy().hasHeightForWidth());
        surfFuncCorrespSchemeCombo->setSizePolicy(sizePolicy5);

        gridLayout_19->addWidget(surfFuncCorrespSchemeCombo, 2, 0, 1, 1);

        smoothSMetricSlider = new QSlider(surfFuncGroup);
        smoothSMetricSlider->setObjectName(QStringLiteral("smoothSMetricSlider"));
        smoothSMetricSlider->setMaximum(10000);
        smoothSMetricSlider->setOrientation(Qt::Horizontal);

        gridLayout_19->addWidget(smoothSMetricSlider, 5, 1, 1, 1);

        diffuseSurfFunc = new QCheckBox(surfFuncGroup);
        diffuseSurfFunc->setObjectName(QStringLiteral("diffuseSurfFunc"));

        gridLayout_19->addWidget(diffuseSurfFunc, 6, 0, 1, 1);

        surfFuncRadiusEps = new QDoubleSpinBox(surfFuncGroup);
        surfFuncRadiusEps->setObjectName(QStringLiteral("surfFuncRadiusEps"));
        surfFuncRadiusEps->setDecimals(8);
        surfFuncRadiusEps->setMinimum(-1);
        surfFuncRadiusEps->setMaximum(1);

        gridLayout_19->addWidget(surfFuncRadiusEps, 3, 0, 1, 1);

        surfFuncSetupCorrespBtn = new QPushButton(surfFuncGroup);
        surfFuncSetupCorrespBtn->setObjectName(QStringLiteral("surfFuncSetupCorrespBtn"));

        gridLayout_19->addWidget(surfFuncSetupCorrespBtn, 2, 1, 1, 1);

        surfFuncKNN = new QSpinBox(surfFuncGroup);
        surfFuncKNN->setObjectName(QStringLiteral("surfFuncKNN"));
        surfFuncKNN->setMinimum(1);

        gridLayout_19->addWidget(surfFuncKNN, 3, 1, 1, 1);


        gridLayout_20->addWidget(surfFuncGroup, 0, 0, 1, 1);

        verticalSpacer_5 = new QSpacerItem(20, 514, QSizePolicy::Minimum, QSizePolicy::Expanding);

        gridLayout_20->addItem(verticalSpacer_5, 1, 0, 1, 1);

        toolTabs->addTab(SMetricTab, QString());
        VideoTab = new QWidget();
        VideoTab->setObjectName(QStringLiteral("VideoTab"));
        gridLayout_23 = new QGridLayout(VideoTab);
        gridLayout_23->setObjectName(QStringLiteral("gridLayout_23"));
        isoSurfGroup = new QGroupBox(VideoTab);
        isoSurfGroup->setObjectName(QStringLiteral("isoSurfGroup"));
        gridLayout_24 = new QGridLayout(isoSurfGroup);
        gridLayout_24->setObjectName(QStringLiteral("gridLayout_24"));
        isoSurfPrecomputeBtn = new QPushButton(isoSurfGroup);
        isoSurfPrecomputeBtn->setObjectName(QStringLiteral("isoSurfPrecomputeBtn"));
        sizePolicy3.setHeightForWidth(isoSurfPrecomputeBtn->sizePolicy().hasHeightForWidth());
        isoSurfPrecomputeBtn->setSizePolicy(sizePolicy3);

        gridLayout_24->addWidget(isoSurfPrecomputeBtn, 0, 0, 1, 1);

        hideIsoSurf = new QCheckBox(isoSurfGroup);
        hideIsoSurf->setObjectName(QStringLiteral("hideIsoSurf"));
        sizePolicy3.setHeightForWidth(hideIsoSurf->sizePolicy().hasHeightForWidth());
        hideIsoSurf->setSizePolicy(sizePolicy3);

        gridLayout_24->addWidget(hideIsoSurf, 0, 1, 1, 1);

        recursiveRefineBox = new QCheckBox(isoSurfGroup);
        recursiveRefineBox->setObjectName(QStringLiteral("recursiveRefineBox"));

        gridLayout_24->addWidget(recursiveRefineBox, 4, 0, 1, 1);

        isoSurfSlider = new QSlider(isoSurfGroup);
        isoSurfSlider->setObjectName(QStringLiteral("isoSurfSlider"));
        sizePolicy1.setHeightForWidth(isoSurfSlider->sizePolicy().hasHeightForWidth());
        isoSurfSlider->setSizePolicy(sizePolicy1);
        isoSurfSlider->setMaximum(9999);
        isoSurfSlider->setSingleStep(5);
        isoSurfSlider->setPageStep(20);
        isoSurfSlider->setOrientation(Qt::Horizontal);

        gridLayout_24->addWidget(isoSurfSlider, 5, 1, 1, 1);

        hideSnappedFaceBox = new QCheckBox(isoSurfGroup);
        hideSnappedFaceBox->setObjectName(QStringLiteral("hideSnappedFaceBox"));

        gridLayout_24->addWidget(hideSnappedFaceBox, 5, 0, 1, 1);

        isoSurfRefineBtn = new QPushButton(isoSurfGroup);
        isoSurfRefineBtn->setObjectName(QStringLiteral("isoSurfRefineBtn"));

        gridLayout_24->addWidget(isoSurfRefineBtn, 4, 1, 1, 1);


        gridLayout_23->addWidget(isoSurfGroup, 0, 0, 1, 2);

        groupBox_7 = new QGroupBox(VideoTab);
        groupBox_7->setObjectName(QStringLiteral("groupBox_7"));
        gridLayout_25 = new QGridLayout(groupBox_7);
        gridLayout_25->setObjectName(QStringLiteral("gridLayout_25"));
        isoContPrecomputeBtn = new QPushButton(groupBox_7);
        isoContPrecomputeBtn->setObjectName(QStringLiteral("isoContPrecomputeBtn"));
        sizePolicy3.setHeightForWidth(isoContPrecomputeBtn->sizePolicy().hasHeightForWidth());
        isoContPrecomputeBtn->setSizePolicy(sizePolicy3);

        gridLayout_25->addWidget(isoContPrecomputeBtn, 1, 0, 1, 1);

        hideIsoCont = new QCheckBox(groupBox_7);
        hideIsoCont->setObjectName(QStringLiteral("hideIsoCont"));
        sizePolicy3.setHeightForWidth(hideIsoCont->sizePolicy().hasHeightForWidth());
        hideIsoCont->setSizePolicy(sizePolicy3);

        gridLayout_25->addWidget(hideIsoCont, 1, 1, 1, 1);

        isoContShowMC = new QCheckBox(groupBox_7);
        isoContShowMC->setObjectName(QStringLiteral("isoContShowMC"));

        gridLayout_25->addWidget(isoContShowMC, 2, 0, 1, 1);

        hidePastSurfBox = new QCheckBox(groupBox_7);
        hidePastSurfBox->setObjectName(QStringLiteral("hidePastSurfBox"));

        gridLayout_25->addWidget(hidePastSurfBox, 2, 1, 1, 1);

        isoContSlider = new QSlider(groupBox_7);
        isoContSlider->setObjectName(QStringLiteral("isoContSlider"));
        isoContSlider->setMaximum(9999);
        isoContSlider->setSingleStep(5);
        isoContSlider->setPageStep(20);
        isoContSlider->setOrientation(Qt::Horizontal);

        gridLayout_25->addWidget(isoContSlider, 3, 1, 1, 1);

        groupBox_9 = new QGroupBox(groupBox_7);
        groupBox_9->setObjectName(QStringLiteral("groupBox_9"));
        gridLayout_27 = new QGridLayout(groupBox_9);
        gridLayout_27->setObjectName(QStringLiteral("gridLayout_27"));
        isoContMinAlphaSlider = new QSlider(groupBox_9);
        isoContMinAlphaSlider->setObjectName(QStringLiteral("isoContMinAlphaSlider"));
        sizePolicy1.setHeightForWidth(isoContMinAlphaSlider->sizePolicy().hasHeightForWidth());
        isoContMinAlphaSlider->setSizePolicy(sizePolicy1);
        isoContMinAlphaSlider->setOrientation(Qt::Horizontal);

        gridLayout_27->addWidget(isoContMinAlphaSlider, 0, 0, 1, 1);

        isoContExpAlphaSlider = new QSlider(groupBox_9);
        isoContExpAlphaSlider->setObjectName(QStringLiteral("isoContExpAlphaSlider"));
        sizePolicy1.setHeightForWidth(isoContExpAlphaSlider->sizePolicy().hasHeightForWidth());
        isoContExpAlphaSlider->setSizePolicy(sizePolicy1);
        isoContExpAlphaSlider->setOrientation(Qt::Horizontal);

        gridLayout_27->addWidget(isoContExpAlphaSlider, 0, 1, 1, 1);


        gridLayout_25->addWidget(groupBox_9, 0, 0, 1, 2);


        gridLayout_23->addWidget(groupBox_7, 1, 0, 2, 2);

        verticalSpacer_6 = new QSpacerItem(20, 40, QSizePolicy::Minimum, QSizePolicy::Expanding);

        gridLayout_23->addItem(verticalSpacer_6, 5, 0, 1, 1);

        groupBox_8 = new QGroupBox(VideoTab);
        groupBox_8->setObjectName(QStringLiteral("groupBox_8"));
        gridLayout_26 = new QGridLayout(groupBox_8);
        gridLayout_26->setObjectName(QStringLiteral("gridLayout_26"));
        evolveAllSlider = new QSlider(groupBox_8);
        evolveAllSlider->setObjectName(QStringLiteral("evolveAllSlider"));
        evolveAllSlider->setMaximum(9999);
        evolveAllSlider->setSingleStep(5);
        evolveAllSlider->setPageStep(20);
        evolveAllSlider->setOrientation(Qt::Horizontal);

        gridLayout_26->addWidget(evolveAllSlider, 1, 1, 1, 1);

        evolveAllPrecomputeBtn = new QPushButton(groupBox_8);
        evolveAllPrecomputeBtn->setObjectName(QStringLiteral("evolveAllPrecomputeBtn"));

        gridLayout_26->addWidget(evolveAllPrecomputeBtn, 1, 0, 1, 1);

        evolveAllShowEntireMC = new QCheckBox(groupBox_8);
        evolveAllShowEntireMC->setObjectName(QStringLiteral("evolveAllShowEntireMC"));

        gridLayout_26->addWidget(evolveAllShowEntireMC, 2, 0, 1, 1);


        gridLayout_23->addWidget(groupBox_8, 4, 0, 1, 2);

        toolTabs->addTab(VideoTab, QString());
        MiscTab = new QWidget();
        MiscTab->setObjectName(QStringLiteral("MiscTab"));
        gridLayout_31 = new QGridLayout(MiscTab);
        gridLayout_31->setObjectName(QStringLiteral("gridLayout_31"));
        verticalSpacer_2 = new QSpacerItem(20, 40, QSizePolicy::Minimum, QSizePolicy::Expanding);

        gridLayout_31->addItem(verticalSpacer_2, 4, 0, 1, 1);

        preprocessGroup = new QGroupBox(MiscTab);
        preprocessGroup->setObjectName(QStringLiteral("preprocessGroup"));
        preprocessGroup->setMaximumSize(QSize(1000, 16777215));
        gridLayout = new QGridLayout(preprocessGroup);
        gridLayout->setObjectName(QStringLiteral("gridLayout"));
        label_2 = new QLabel(preprocessGroup);
        label_2->setObjectName(QStringLiteral("label_2"));

        gridLayout->addWidget(label_2, 0, 0, 1, 1);

        pertSpin = new QDoubleSpinBox(preprocessGroup);
        pertSpin->setObjectName(QStringLiteral("pertSpin"));
        pertSpin->setDecimals(15);

        gridLayout->addWidget(pertSpin, 4, 1, 1, 1);

        label = new QLabel(preprocessGroup);
        label->setObjectName(QStringLiteral("label"));

        gridLayout->addWidget(label, 4, 0, 1, 1);

        epsilonSpin = new QDoubleSpinBox(preprocessGroup);
        epsilonSpin->setObjectName(QStringLiteral("epsilonSpin"));
        epsilonSpin->setDecimals(15);

        gridLayout->addWidget(epsilonSpin, 0, 1, 1, 1);


        gridLayout_31->addWidget(preprocessGroup, 0, 0, 1, 1);

        miscBurntimeGroup = new QGroupBox(MiscTab);
        miscBurntimeGroup->setObjectName(QStringLiteral("miscBurntimeGroup"));
        gridLayout_28 = new QGridLayout(miscBurntimeGroup);
        gridLayout_28->setObjectName(QStringLiteral("gridLayout_28"));
        visSphereFromFileBtn = new QPushButton(miscBurntimeGroup);
        visSphereFromFileBtn->setObjectName(QStringLiteral("visSphereFromFileBtn"));

        gridLayout_28->addWidget(visSphereFromFileBtn, 2, 0, 1, 1);

        visSphereSpin = new QSpinBox(miscBurntimeGroup);
        visSphereSpin->setObjectName(QStringLiteral("visSphereSpin"));

        gridLayout_28->addWidget(visSphereSpin, 1, 1, 1, 1);

        randBTBtn = new QPushButton(miscBurntimeGroup);
        randBTBtn->setObjectName(QStringLiteral("randBTBtn"));

        gridLayout_28->addWidget(randBTBtn, 4, 0, 1, 1);

        printBTBtn = new QPushButton(miscBurntimeGroup);
        printBTBtn->setObjectName(QStringLiteral("printBTBtn"));

        gridLayout_28->addWidget(printBTBtn, 4, 1, 1, 1);

        visSphereSlider = new QSlider(miscBurntimeGroup);
        visSphereSlider->setObjectName(QStringLiteral("visSphereSlider"));
        sizePolicy4.setHeightForWidth(visSphereSlider->sizePolicy().hasHeightForWidth());
        visSphereSlider->setSizePolicy(sizePolicy4);
        visSphereSlider->setOrientation(Qt::Horizontal);

        gridLayout_28->addWidget(visSphereSlider, 1, 0, 1, 1);

        showSphereBox = new QCheckBox(miscBurntimeGroup);
        showSphereBox->setObjectName(QStringLiteral("showSphereBox"));

        gridLayout_28->addWidget(showSphereBox, 0, 0, 1, 1);

        printSelectVertInfoBtn = new QPushButton(miscBurntimeGroup);
        printSelectVertInfoBtn->setObjectName(QStringLiteral("printSelectVertInfoBtn"));

        gridLayout_28->addWidget(printSelectVertInfoBtn, 2, 1, 1, 1);


        gridLayout_31->addWidget(miscBurntimeGroup, 1, 0, 1, 1);

        renderParamGroup = new QGroupBox(MiscTab);
        renderParamGroup->setObjectName(QStringLiteral("renderParamGroup"));
        QSizePolicy sizePolicy9(QSizePolicy::Preferred, QSizePolicy::Minimum);
        sizePolicy9.setHorizontalStretch(0);
        sizePolicy9.setVerticalStretch(0);
        sizePolicy9.setHeightForWidth(renderParamGroup->sizePolicy().hasHeightForWidth());
        renderParamGroup->setSizePolicy(sizePolicy9);
        formLayout_3 = new QFormLayout(renderParamGroup);
        formLayout_3->setObjectName(QStringLiteral("formLayout_3"));
        label_3 = new QLabel(renderParamGroup);
        label_3->setObjectName(QStringLiteral("label_3"));

        formLayout_3->setWidget(1, QFormLayout::LabelRole, label_3);

        linesRenderOffsetSlider = new QSlider(renderParamGroup);
        linesRenderOffsetSlider->setObjectName(QStringLiteral("linesRenderOffsetSlider"));
        linesRenderOffsetSlider->setMaximum(1000000);
        linesRenderOffsetSlider->setOrientation(Qt::Horizontal);

        formLayout_3->setWidget(1, QFormLayout::FieldRole, linesRenderOffsetSlider);


        gridLayout_31->addWidget(renderParamGroup, 2, 0, 1, 1);

        qmatGroup = new QGroupBox(MiscTab);
        qmatGroup->setObjectName(QStringLiteral("qmatGroup"));
        QSizePolicy sizePolicy10(QSizePolicy::Preferred, QSizePolicy::Expanding);
        sizePolicy10.setHorizontalStretch(0);
        sizePolicy10.setVerticalStretch(0);
        sizePolicy10.setHeightForWidth(qmatGroup->sizePolicy().hasHeightForWidth());
        qmatGroup->setSizePolicy(sizePolicy10);
        formLayout_4 = new QFormLayout(qmatGroup);
        formLayout_4->setObjectName(QStringLiteral("formLayout_4"));
        label_9 = new QLabel(qmatGroup);
        label_9->setObjectName(QStringLiteral("label_9"));
        sizePolicy8.setHeightForWidth(label_9->sizePolicy().hasHeightForWidth());
        label_9->setSizePolicy(sizePolicy8);

        formLayout_4->setWidget(0, QFormLayout::LabelRole, label_9);

        outputWenpingBtn = new QPushButton(qmatGroup);
        outputWenpingBtn->setObjectName(QStringLiteral("outputWenpingBtn"));
        sizePolicy3.setHeightForWidth(outputWenpingBtn->sizePolicy().hasHeightForWidth());
        outputWenpingBtn->setSizePolicy(sizePolicy3);

        formLayout_4->setWidget(0, QFormLayout::FieldRole, outputWenpingBtn);

        readQMatMA = new QPushButton(qmatGroup);
        readQMatMA->setObjectName(QStringLiteral("readQMatMA"));
        sizePolicy3.setHeightForWidth(readQMatMA->sizePolicy().hasHeightForWidth());
        readQMatMA->setSizePolicy(sizePolicy3);

        formLayout_4->setWidget(1, QFormLayout::LabelRole, readQMatMA);

        hideQMatBox = new QCheckBox(qmatGroup);
        hideQMatBox->setObjectName(QStringLiteral("hideQMatBox"));
        sizePolicy3.setHeightForWidth(hideQMatBox->sizePolicy().hasHeightForWidth());
        hideQMatBox->setSizePolicy(sizePolicy3);

        formLayout_4->setWidget(1, QFormLayout::FieldRole, hideQMatBox);


        gridLayout_31->addWidget(qmatGroup, 3, 0, 1, 1);

        toolTabs->addTab(MiscTab, QString());
        hiddenTab = new QWidget();
        hiddenTab->setObjectName(QStringLiteral("hiddenTab"));
        stSubdivCombo = new QComboBox(hiddenTab);
        stSubdivCombo->setObjectName(QStringLiteral("stSubdivCombo"));
        stSubdivCombo->setGeometry(QRect(160, 10, 69, 20));
        QSizePolicy sizePolicy11(QSizePolicy::Minimum, QSizePolicy::Expanding);
        sizePolicy11.setHorizontalStretch(0);
        sizePolicy11.setVerticalStretch(0);
        sizePolicy11.setHeightForWidth(stSubdivCombo->sizePolicy().hasHeightForWidth());
        stSubdivCombo->setSizePolicy(sizePolicy11);
        stSubdivCombo->setMaximumSize(QSize(100, 20));
        edgeWeightCombo = new QComboBox(hiddenTab);
        edgeWeightCombo->setObjectName(QStringLiteral("edgeWeightCombo"));
        edgeWeightCombo->setEnabled(true);
        edgeWeightCombo->setGeometry(QRect(71, 6, 69, 20));
        sizePolicy11.setHeightForWidth(edgeWeightCombo->sizePolicy().hasHeightForWidth());
        edgeWeightCombo->setSizePolicy(sizePolicy11);
        edgeWeightCombo->setMaximumSize(QSize(100, 20));
        pointSizeSpin = new QDoubleSpinBox(hiddenTab);
        pointSizeSpin->setObjectName(QStringLiteral("pointSizeSpin"));
        pointSizeSpin->setGeometry(QRect(9, 6, 49, 20));
        sizePolicy4.setHeightForWidth(pointSizeSpin->sizePolicy().hasHeightForWidth());
        pointSizeSpin->setSizePolicy(sizePolicy4);
        pointSizeSpin->setMaximumSize(QSize(100, 20));
        doMAFaceScalarDiffusion = new QCheckBox(hiddenTab);
        doMAFaceScalarDiffusion->setObjectName(QStringLiteral("doMAFaceScalarDiffusion"));
        doMAFaceScalarDiffusion->setGeometry(QRect(91, 31, 103, 17));
        computeAllDistBtn = new QPushButton(hiddenTab);
        computeAllDistBtn->setObjectName(QStringLiteral("computeAllDistBtn"));
        computeAllDistBtn->setGeometry(QRect(10, 30, 75, 23));
        usePerSheetBox = new QCheckBox(hiddenTab);
        usePerSheetBox->setObjectName(QStringLiteral("usePerSheetBox"));
        usePerSheetBox->setGeometry(QRect(201, 31, 69, 17));
        maxPruneDistMASpin2 = new QDoubleSpinBox(hiddenTab);
        maxPruneDistMASpin2->setObjectName(QStringLiteral("maxPruneDistMASpin2"));
        maxPruneDistMASpin2->setGeometry(QRect(131, 131, 97, 20));
        maxPruneDistMASpin2->setDecimals(10);
        maPruneRatio = new QDoubleSpinBox(hiddenTab);
        maPruneRatio->setObjectName(QStringLiteral("maPruneRatio"));
        maPruneRatio->setGeometry(QRect(11, 61, 97, 20));
        maPruneRatio->setDecimals(8);
        minPruneDistMASpin1 = new QDoubleSpinBox(hiddenTab);
        minPruneDistMASpin1->setObjectName(QStringLiteral("minPruneDistMASpin1"));
        minPruneDistMASpin1->setGeometry(QRect(8, 83, 107, 20));
        minPruneDistMASpin1->setDecimals(10);
        pruneMASpin2 = new QDoubleSpinBox(hiddenTab);
        pruneMASpin2->setObjectName(QStringLiteral("pruneMASpin2"));
        pruneMASpin2->setGeometry(QRect(11, 111, 107, 20));
        pruneMASpin2->setDecimals(30);
        maxPruneDistMASpin1 = new QDoubleSpinBox(hiddenTab);
        maxPruneDistMASpin1->setObjectName(QStringLiteral("maxPruneDistMASpin1"));
        maxPruneDistMASpin1->setGeometry(QRect(111, 61, 97, 20));
        maxPruneDistMASpin1->setDecimals(10);
        pruneMASlider2 = new QSlider(hiddenTab);
        pruneMASlider2->setObjectName(QStringLiteral("pruneMASlider2"));
        pruneMASlider2->setGeometry(QRect(131, 81, 107, 19));
        pruneMASlider2->setOrientation(Qt::Horizontal);
        pruneMADistCombo2 = new QComboBox(hiddenTab);
        pruneMADistCombo2->setObjectName(QStringLiteral("pruneMADistCombo2"));
        pruneMADistCombo2->setGeometry(QRect(131, 111, 97, 20));
        minPruneDistMASpin2 = new QDoubleSpinBox(hiddenTab);
        minPruneDistMASpin2->setObjectName(QStringLiteral("minPruneDistMASpin2"));
        minPruneDistMASpin2->setGeometry(QRect(11, 131, 107, 20));
        minPruneDistMASpin2->setDecimals(10);
        edgeDualOptCombo = new QComboBox(hiddenTab);
        edgeDualOptCombo->setObjectName(QStringLiteral("edgeDualOptCombo"));
        edgeDualOptCombo->setGeometry(QRect(0, 161, 115, 20));
        polyDualOptCombo = new QComboBox(hiddenTab);
        polyDualOptCombo->setObjectName(QStringLiteral("polyDualOptCombo"));
        polyDualOptCombo->setGeometry(QRect(121, 161, 115, 20));
        underEstimateBox = new QCheckBox(hiddenTab);
        underEstimateBox->setObjectName(QStringLiteral("underEstimateBox"));
        underEstimateBox->setGeometry(QRect(130, 213, 115, 17));
        onlyUnburntBox = new QCheckBox(hiddenTab);
        onlyUnburntBox->setObjectName(QStringLiteral("onlyUnburntBox"));
        onlyUnburntBox->setGeometry(QRect(9, 190, 115, 17));
        stopAtJuncBox = new QCheckBox(hiddenTab);
        stopAtJuncBox->setObjectName(QStringLiteral("stopAtJuncBox"));
        stopAtJuncBox->setGeometry(QRect(9, 213, 115, 17));
        protectBT2Box = new QCheckBox(hiddenTab);
        protectBT2Box->setObjectName(QStringLiteral("protectBT2Box"));
        protectBT2Box->setGeometry(QRect(130, 190, 115, 17));
        hsBt1Bt2RelSlider = new QSlider(hiddenTab);
        hsBt1Bt2RelSlider->setObjectName(QStringLiteral("hsBt1Bt2RelSlider"));
        hsBt1Bt2RelSlider->setGeometry(QRect(130, 260, 115, 19));
        hsBt1Bt2RelSlider->setMaximum(1000);
        hsBt1Bt2RelSlider->setOrientation(Qt::Horizontal);
        hsCurveRelDiffSpin = new QDoubleSpinBox(hiddenTab);
        hsCurveRelDiffSpin->setObjectName(QStringLiteral("hsCurveRelDiffSpin"));
        hsCurveRelDiffSpin->setGeometry(QRect(120, 300, 115, 20));
        hsCurveRelDiffSpin->setDecimals(8);
        hsFaceRelDiffSpin = new QDoubleSpinBox(hiddenTab);
        hsFaceRelDiffSpin->setObjectName(QStringLiteral("hsFaceRelDiffSpin"));
        hsFaceRelDiffSpin->setGeometry(QRect(120, 278, 115, 20));
        hsFaceRelDiffSpin->setDecimals(8);
        hsBt2Bt3RelSlider = new QSlider(hiddenTab);
        hsBt2Bt3RelSlider->setObjectName(QStringLiteral("hsBt2Bt3RelSlider"));
        hsBt2Bt3RelSlider->setGeometry(QRect(131, 240, 115, 19));
        hsBt2Bt3RelSlider->setMaximum(1000);
        hsBt2Bt3RelSlider->setOrientation(Qt::Horizontal);
        label_8 = new QLabel(hiddenTab);
        label_8->setObjectName(QStringLiteral("label_8"));
        label_8->setGeometry(QRect(9, 260, 115, 19));
        label_6 = new QLabel(hiddenTab);
        label_6->setObjectName(QStringLiteral("label_6"));
        label_6->setGeometry(QRect(10, 240, 115, 19));
        hsCmpntFaceNumTSpin = new QDoubleSpinBox(hiddenTab);
        hsCmpntFaceNumTSpin->setObjectName(QStringLiteral("hsCmpntFaceNumTSpin"));
        hsCmpntFaceNumTSpin->setGeometry(QRect(0, 320, 115, 20));
        hsCmpntFaceNumTSpin->setDecimals(8);
        hsCmpntFaceNumTSpin->setMaximum(999999);
        mcBurnGroup = new QGroupBox(hiddenTab);
        mcBurnGroup->setObjectName(QStringLiteral("mcBurnGroup"));
        mcBurnGroup->setGeometry(QRect(10, 350, 71, 61));
        gridLayout_14 = new QGridLayout(mcBurnGroup);
        gridLayout_14->setObjectName(QStringLiteral("gridLayout_14"));
        burnMedialCurveBtn = new QPushButton(mcBurnGroup);
        burnMedialCurveBtn->setObjectName(QStringLiteral("burnMedialCurveBtn"));

        gridLayout_14->addWidget(burnMedialCurveBtn, 0, 0, 1, 1);

        printMCStatsBtn = new QPushButton(hiddenTab);
        printMCStatsBtn->setObjectName(QStringLiteral("printMCStatsBtn"));
        printMCStatsBtn->setGeometry(QRect(120, 330, 115, 23));
        pruneMCSlider2 = new QSlider(hiddenTab);
        pruneMCSlider2->setObjectName(QStringLiteral("pruneMCSlider2"));
        pruneMCSlider2->setGeometry(QRect(10, 424, 115, 19));
        pruneMCSlider2->setOrientation(Qt::Horizontal);
        mcPruneRatio = new QDoubleSpinBox(hiddenTab);
        mcPruneRatio->setObjectName(QStringLiteral("mcPruneRatio"));
        mcPruneRatio->setGeometry(QRect(130, 400, 115, 20));
        sizePolicy3.setHeightForWidth(mcPruneRatio->sizePolicy().hasHeightForWidth());
        mcPruneRatio->setSizePolicy(sizePolicy3);
        mcPruneRatio->setDecimals(10);
        pruneMCDistCombo2 = new QComboBox(hiddenTab);
        pruneMCDistCombo2->setObjectName(QStringLiteral("pruneMCDistCombo2"));
        pruneMCDistCombo2->setGeometry(QRect(131, 424, 115, 20));
        pruneMCDistSpin2 = new QDoubleSpinBox(hiddenTab);
        pruneMCDistSpin2->setObjectName(QStringLiteral("pruneMCDistSpin2"));
        pruneMCDistSpin2->setGeometry(QRect(10, 450, 115, 20));
        pruneMCDistSpin2->setDecimals(10);
        colorPerVert = new QRadioButton(hiddenTab);
        colorPerVert->setObjectName(QStringLiteral("colorPerVert"));
        colorPerVert->setGeometry(QRect(12, 480, 102, 17));
        colorPerFace = new QRadioButton(hiddenTab);
        colorPerFace->setObjectName(QStringLiteral("colorPerFace"));
        colorPerFace->setGeometry(QRect(120, 480, 102, 17));
        hsUseTypedInput = new QCheckBox(hiddenTab);
        hsUseTypedInput->setObjectName(QStringLiteral("hsUseTypedInput"));
        hsUseTypedInput->setGeometry(QRect(10, 500, 115, 17));
        visBallStick = new QCheckBox(hiddenTab);
        visBallStick->setObjectName(QStringLiteral("visBallStick"));
        visBallStick->setGeometry(QRect(9, 520, 115, 17));
        sizePolicy3.setHeightForWidth(visBallStick->sizePolicy().hasHeightForWidth());
        visBallStick->setSizePolicy(sizePolicy3);
        preserveMCTopo = new QCheckBox(hiddenTab);
        preserveMCTopo->setObjectName(QStringLiteral("preserveMCTopo"));
        preserveMCTopo->setGeometry(QRect(130, 520, 115, 17));
        toolTabs->addTab(hiddenTab, QString());
        splitter->addWidget(toolTabs);
        layoutWidget = new QWidget(splitter);
        layoutWidget->setObjectName(QStringLiteral("layoutWidget"));
        glLayout = new QHBoxLayout(layoutWidget);
        glLayout->setObjectName(QStringLiteral("glLayout"));
        glLayout->setContentsMargins(0, 0, 0, 0);
        horizontalSpacer = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        glLayout->addItem(horizontalSpacer);

        splitter->addWidget(layoutWidget);

        gridLayout_11->addWidget(splitter, 0, 0, 1, 1);

        MainWindow->setCentralWidget(centralwidget);
        statusbar = new QStatusBar(MainWindow);
        statusbar->setObjectName(QStringLiteral("statusbar"));
        MainWindow->setStatusBar(statusbar);

        retranslateUi(MainWindow);

        toolTabs->setCurrentIndex(2);


        QMetaObject::connectSlotsByName(MainWindow);
    } // setupUi

    void retranslateUi(QMainWindow *MainWindow)
    {
        MainWindow->setWindowTitle(QApplication::translate("MainWindow", "MainWindow", 0));
        openFile->setText(QApplication::translate("MainWindow", "&Open", 0));
        saveView->setText(QApplication::translate("MainWindow", "Save View", 0));
#ifndef QT_NO_TOOLTIP
        saveView->setToolTip(QApplication::translate("MainWindow", "Save camera & viewport to file", 0));
#endif // QT_NO_TOOLTIP
        loadView->setText(QApplication::translate("MainWindow", "Load View", 0));
#ifndef QT_NO_TOOLTIP
        loadView->setToolTip(QApplication::translate("MainWindow", "Load a saved view from file", 0));
#endif // QT_NO_TOOLTIP
        readRadii->setText(QApplication::translate("MainWindow", "Read Radii", 0));
#ifndef QT_NO_TOOLTIP
        readRadii->setToolTip(QApplication::translate("MainWindow", "read radii of MA from .r file", 0));
#endif // QT_NO_TOOLTIP
        actionHelp->setText(QApplication::translate("MainWindow", "Help", 0));
#ifndef QT_NO_TOOLTIP
        toolTabs->setToolTip(QString());
#endif // QT_NO_TOOLTIP
#ifndef QT_NO_TOOLTIP
        fileTab->setToolTip(QApplication::translate("MainWindow", "Load all required files here", 0));
#endif // QT_NO_TOOLTIP
        loadFileGroup->setTitle(QString());
        label_11->setText(QApplication::translate("MainWindow", "medial axis", 0));
#ifndef QT_NO_TOOLTIP
        maFileEdit->setToolTip(QApplication::translate("MainWindow", "the medial axis file", 0));
#endif // QT_NO_TOOLTIP
#ifndef QT_NO_TOOLTIP
        browseMAFileBtn->setToolTip(QApplication::translate("MainWindow", "select a medial axis file (.off)", 0));
#endif // QT_NO_TOOLTIP
        browseMAFileBtn->setText(QApplication::translate("MainWindow", "browse", 0));
        label_12->setText(QApplication::translate("MainWindow", "3d shape", 0));
#ifndef QT_NO_TOOLTIP
        shapeFileEdit->setToolTip(QApplication::translate("MainWindow", "the 3d shape file (optional)", 0));
#endif // QT_NO_TOOLTIP
#ifndef QT_NO_TOOLTIP
        browse3DShapeFileBtn->setToolTip(QApplication::translate("MainWindow", "select a 3d shape file (.off)", 0));
#endif // QT_NO_TOOLTIP
        browse3DShapeFileBtn->setText(QApplication::translate("MainWindow", "browse", 0));
        label_13->setText(QApplication::translate("MainWindow", "radius", 0));
#ifndef QT_NO_TOOLTIP
        radiusFileEdit->setToolTip(QApplication::translate("MainWindow", "the radius file (optional)", 0));
#endif // QT_NO_TOOLTIP
#ifndef QT_NO_TOOLTIP
        browseRadiusFileBtn->setToolTip(QApplication::translate("MainWindow", "(optionally) select a radius file (.r)", 0));
#endif // QT_NO_TOOLTIP
        browseRadiusFileBtn->setText(QApplication::translate("MainWindow", "browse", 0));
        loadFilesBtn->setText(QApplication::translate("MainWindow", "load", 0));
        toolTabs->setTabText(toolTabs->indexOf(fileTab), QApplication::translate("MainWindow", "Files", 0));
#ifndef QT_NO_TOOLTIP
        origSurfTab->setToolTip(QApplication::translate("MainWindow", "<html><head/>\n"
"<body>visualizatoin options for<br>\n"
"- 3d shape<br>\n"
"- medial axis <br>\n"
"- medial curve <br>\n"
"- skeleton <br>\n"
"</body></html>", 0));
#endif // QT_NO_TOOLTIP
        groupBox->setTitle(QApplication::translate("MainWindow", "general", 0));
#ifndef QT_NO_TOOLTIP
        DPMaxRendersSpin->setToolTip(QApplication::translate("MainWindow", "specify # layers here", 0));
#endif // QT_NO_TOOLTIP
#ifndef QT_NO_TOOLTIP
        colorBGBtn->setToolTip(QApplication::translate("MainWindow", "change background color", 0));
#endif // QT_NO_TOOLTIP
        colorBGBtn->setText(QApplication::translate("MainWindow", "BG color", 0));
#ifndef QT_NO_TOOLTIP
        drawEdge->setToolTip(QApplication::translate("MainWindow", "show wireframe on original shape and medial axis", 0));
#endif // QT_NO_TOOLTIP
        drawEdge->setText(QApplication::translate("MainWindow", "wire-frame", 0));
#ifndef QT_NO_TOOLTIP
        saveViewBtn->setToolTip(QApplication::translate("MainWindow", "Save camera & viewport to file", 0));
#endif // QT_NO_TOOLTIP
        saveViewBtn->setText(QApplication::translate("MainWindow", "save view", 0));
#ifndef QT_NO_TOOLTIP
        loadViewBtn->setToolTip(QApplication::translate("MainWindow", "Load a saved view from file", 0));
#endif // QT_NO_TOOLTIP
        loadViewBtn->setText(QApplication::translate("MainWindow", "load view", 0));
#ifndef QT_NO_TOOLTIP
        useTrueTransparencyBox->setToolTip(QApplication::translate("MainWindow", "<html><head/>\n"
"<body>\n"
"<p>use &quot;ground truth&quot; transparency (depth-peel) instead of opengl default alpha blending <br>\n"
"(need to specify # layers for which transparency will be computed)</p>\n"
"</body></html>", 0));
#endif // QT_NO_TOOLTIP
        useTrueTransparencyBox->setText(QApplication::translate("MainWindow", "true transparency", 0));
#ifndef QT_NO_TOOLTIP
        origSurfViewGroup->setToolTip(QApplication::translate("MainWindow", "view options for 3d shape", 0));
#endif // QT_NO_TOOLTIP
        origSurfViewGroup->setTitle(QApplication::translate("MainWindow", "shape", 0));
#ifndef QT_NO_TOOLTIP
        colorOrigBtn->setToolTip(QApplication::translate("MainWindow", "Give 3d surface constant color", 0));
#endif // QT_NO_TOOLTIP
        colorOrigBtn->setText(QApplication::translate("MainWindow", "color", 0));
#ifndef QT_NO_TOOLTIP
        hideOrig->setToolTip(QApplication::translate("MainWindow", "hide/show current 3d surface", 0));
#endif // QT_NO_TOOLTIP
        hideOrig->setText(QApplication::translate("MainWindow", "hide", 0));
        label_4->setText(QApplication::translate("MainWindow", "transparent", 0));
#ifndef QT_NO_TOOLTIP
        origTransparentSlider->setToolTip(QApplication::translate("MainWindow", "change transparency of the 3D surface", 0));
#endif // QT_NO_TOOLTIP
#ifndef QT_NO_TOOLTIP
        maViewGroup->setToolTip(QApplication::translate("MainWindow", "view options for medial axis", 0));
#endif // QT_NO_TOOLTIP
        maViewGroup->setTitle(QApplication::translate("MainWindow", "medial axis", 0));
#ifndef QT_NO_TOOLTIP
        hideMA->setToolTip(QApplication::translate("MainWindow", "hide/show medial axis", 0));
#endif // QT_NO_TOOLTIP
        hideMA->setText(QApplication::translate("MainWindow", "hide", 0));
#ifndef QT_NO_TOOLTIP
        turnOffLightingForMABox->setToolTip(QApplication::translate("MainWindow", "disable/enable lighting on medial axis", 0));
#endif // QT_NO_TOOLTIP
        turnOffLightingForMABox->setText(QApplication::translate("MainWindow", "no lighting", 0));
#ifndef QT_NO_TOOLTIP
        drawStPoints->setToolTip(QApplication::translate("MainWindow", "visualize all steiner points on medial axis", 0));
#endif // QT_NO_TOOLTIP
        drawStPoints->setText(QApplication::translate("MainWindow", "steiner points", 0));
#ifndef QT_NO_TOOLTIP
        colorMABtn->setToolTip(QApplication::translate("MainWindow", "Give MA a constant color", 0));
#endif // QT_NO_TOOLTIP
        colorMABtn->setText(QApplication::translate("MainWindow", "color", 0));
#ifndef QT_NO_TOOLTIP
        distVisBtnGroup->setToolTip(QApplication::translate("MainWindow", "pick a measure to visualize on the medial axis", 0));
#endif // QT_NO_TOOLTIP
        distVisBtnGroup->setTitle(QApplication::translate("MainWindow", "view measure", 0));
#ifndef QT_NO_TOOLTIP
        visMADistCombo->setToolTip(QApplication::translate("MainWindow", "choose a distance metric to visualize on medial axis", 0));
#endif // QT_NO_TOOLTIP
#ifndef QT_NO_TOOLTIP
        clampMADistBox->setToolTip(QApplication::translate("MainWindow", "map the selected measure to the right range for visualization purpose", 0));
#endif // QT_NO_TOOLTIP
        clampMADistBox->setText(QString());
#ifndef QT_NO_TOOLTIP
        minVisDistMASpin->setToolTip(QString());
#endif // QT_NO_TOOLTIP
#ifndef QT_NO_TOOLTIP
        maxVisDistMASpin->setToolTip(QString());
#endif // QT_NO_TOOLTIP
        MAAlphaGroup->setTitle(QApplication::translate("MainWindow", "transparency", 0));
#ifndef QT_NO_TOOLTIP
        MATransparentSlider->setToolTip(QApplication::translate("MainWindow", "MA min alpha", 0));
#endif // QT_NO_TOOLTIP
#ifndef QT_NO_TOOLTIP
        MATransExpSlider->setToolTip(QApplication::translate("MainWindow", "MA alpha exponent", 0));
#endif // QT_NO_TOOLTIP
#ifndef QT_NO_TOOLTIP
        mcViewGroup->setToolTip(QApplication::translate("MainWindow", "options for viewing the MC and the measures on it", 0));
#endif // QT_NO_TOOLTIP
        mcViewGroup->setTitle(QApplication::translate("MainWindow", "medial curve", 0));
#ifndef QT_NO_TOOLTIP
        visBurnt->setToolTip(QApplication::translate("MainWindow", "visualize all burn trees", 0));
#endif // QT_NO_TOOLTIP
        visBurnt->setText(QApplication::translate("MainWindow", "show burn trees", 0));
#ifndef QT_NO_TOOLTIP
        colorMCBtn->setToolTip(QApplication::translate("MainWindow", "Give MC a constant color", 0));
#endif // QT_NO_TOOLTIP
        colorMCBtn->setText(QApplication::translate("MainWindow", "color", 0));
#ifndef QT_NO_TOOLTIP
        hideMC->setToolTip(QApplication::translate("MainWindow", "hide/show medial curve", 0));
#endif // QT_NO_TOOLTIP
        hideMC->setText(QApplication::translate("MainWindow", "hide MC", 0));
        label_10->setText(QApplication::translate("MainWindow", "line width", 0));
#ifndef QT_NO_TOOLTIP
        MCWidthSpin->setToolTip(QApplication::translate("MainWindow", "medial curve line width", 0));
#endif // QT_NO_TOOLTIP
#ifndef QT_NO_TOOLTIP
        groupBox_2->setToolTip(QApplication::translate("MainWindow", "pick a measure to view on medial curve here", 0));
#endif // QT_NO_TOOLTIP
        groupBox_2->setTitle(QApplication::translate("MainWindow", "view measure", 0));
#ifndef QT_NO_TOOLTIP
        visMCDistCombo->setToolTip(QApplication::translate("MainWindow", "choose a measure on medial curve to visualize", 0));
#endif // QT_NO_TOOLTIP
#ifndef QT_NO_TOOLTIP
        clampMCDistBox->setToolTip(QApplication::translate("MainWindow", "map the the selected measure to the right range", 0));
#endif // QT_NO_TOOLTIP
        clampMCDistBox->setText(QString());
#ifndef QT_NO_TOOLTIP
        minVisDistMCSpin->setToolTip(QApplication::translate("MainWindow", "min of the clamp range", 0));
#endif // QT_NO_TOOLTIP
#ifndef QT_NO_TOOLTIP
        maxVisDistMCSpin->setToolTip(QApplication::translate("MainWindow", "max of the clamp range", 0));
#endif // QT_NO_TOOLTIP
#ifndef QT_NO_TOOLTIP
        MCAlphaGroup->setToolTip(QApplication::translate("MainWindow", "make medial curve look transparent ", 0));
#endif // QT_NO_TOOLTIP
        MCAlphaGroup->setTitle(QApplication::translate("MainWindow", "transparency", 0));
#ifndef QT_NO_TOOLTIP
        MCAlphaSlider->setToolTip(QApplication::translate("MainWindow", "min alpha", 0));
#endif // QT_NO_TOOLTIP
#ifndef QT_NO_TOOLTIP
        MCTransExpSlider->setToolTip(QApplication::translate("MainWindow", "alpha exponent", 0));
#endif // QT_NO_TOOLTIP
#ifndef QT_NO_TOOLTIP
        hsViewGroup->setToolTip(QApplication::translate("MainWindow", "view options for skeleton", 0));
#endif // QT_NO_TOOLTIP
        hsViewGroup->setTitle(QApplication::translate("MainWindow", "skeleton", 0));
        hideHS->setText(QApplication::translate("MainWindow", "hide HS", 0));
        toolTabs->setTabText(toolTabs->indexOf(origSurfTab), QApplication::translate("MainWindow", "Visualization", 0));
#ifndef QT_NO_TOOLTIP
        MATab->setToolTip(QApplication::translate("MainWindow", "creating and viewing erosion thickness on the medial axis", 0));
#endif // QT_NO_TOOLTIP
        maBurnGroup->setTitle(QApplication::translate("MainWindow", "burn", 0));
#ifndef QT_NO_TOOLTIP
        nFixedSteinerSpin->setToolTip(QApplication::translate("MainWindow", "<html><head/><body><p>parameter \317\211 (distance between two steiner points)</p></body></html>", 0));
#endif // QT_NO_TOOLTIP
#ifndef QT_NO_TOOLTIP
        cleanTopoBtn->setToolTip(QApplication::translate("MainWindow", "<html><head/><body><p>click if you want to remove all &quot;pockets&quot; (closed subcomplex) from the medial axis. </p><p>A cleaned medial axis will be written to a file.</p></body></html>", 0));
#endif // QT_NO_TOOLTIP
        cleanTopoBtn->setText(QApplication::translate("MainWindow", "clean topo", 0));
#ifndef QT_NO_TOOLTIP
        burnBtn->setToolTip(QApplication::translate("MainWindow", "compute burn time and ET", 0));
#endif // QT_NO_TOOLTIP
        burnBtn->setText(QApplication::translate("MainWindow", "burn", 0));
#ifndef QT_NO_TOOLTIP
        maMeasureGroup->setToolTip(QApplication::translate("MainWindow", "<html><head/><body>\n"
"<p>This is where you can visualize measure on the medial axis, <br>\n"
"or use measure to prune the medial axis</p>\n"
"</body></html>", 0));
#endif // QT_NO_TOOLTIP
        maMeasureGroup->setTitle(QApplication::translate("MainWindow", "thresholding", 0));
#ifndef QT_NO_TOOLTIP
        pruneMASlider1->setToolTip(QApplication::translate("MainWindow", "prune MA by the specified metric", 0));
#endif // QT_NO_TOOLTIP
#ifndef QT_NO_TOOLTIP
        pruneMADistCombo1->setToolTip(QApplication::translate("MainWindow", "choose a dist metric for medial axis thresholding", 0));
#endif // QT_NO_TOOLTIP
#ifndef QT_NO_TOOLTIP
        pruneMASpin1->setToolTip(QApplication::translate("MainWindow", "current value of the measure", 0));
#endif // QT_NO_TOOLTIP
#ifndef QT_NO_TOOLTIP
        enableFinePruneMA->setToolTip(QApplication::translate("MainWindow", "visualize the finner medial axis whose boundary is refined by the iso-contour of the measure", 0));
#endif // QT_NO_TOOLTIP
        enableFinePruneMA->setText(QApplication::translate("MainWindow", "fine medial axis", 0));
#ifndef QT_NO_TOOLTIP
        exportBTGroup->setToolTip(QApplication::translate("MainWindow", "export ET here", 0));
#endif // QT_NO_TOOLTIP
        exportBTGroup->setTitle(QApplication::translate("MainWindow", "export ET", 0));
#ifndef QT_NO_TOOLTIP
        exportBTBtn->setToolTip(QString());
#endif // QT_NO_TOOLTIP
        exportBTBtn->setText(QApplication::translate("MainWindow", "export", 0));
#ifndef QT_NO_TOOLTIP
        exportPerSectorBTBox->setToolTip(QApplication::translate("MainWindow", "<html><head/>\n"
"<body><p>Instead of writing just the ET value for each vertex into a .et file, <br>\n"
"export the ET value for each sector of each vertex into a .etps file</p>\n"
"</body></html>", 0));
#endif // QT_NO_TOOLTIP
        exportPerSectorBTBox->setText(QApplication::translate("MainWindow", "per-sector", 0));
        toolTabs->setTabText(toolTabs->indexOf(MATab), QApplication::translate("MainWindow", "ET", 0));
#ifndef QT_NO_TOOLTIP
        HSTab->setToolTip(QApplication::translate("MainWindow", "skeleton (generate, prune, etc.)", 0));
#endif // QT_NO_TOOLTIP
#ifndef QT_NO_TOOLTIP
        mcPruneGroup->setToolTip(QApplication::translate("MainWindow", "explore medial curve pruning parameter here", 0));
#endif // QT_NO_TOOLTIP
        mcPruneGroup->setTitle(QApplication::translate("MainWindow", "thresholding medial curve", 0));
#ifndef QT_NO_TOOLTIP
        pruneMCDistSpin1->setToolTip(QApplication::translate("MainWindow", "cur value of pruning threshold", 0));
#endif // QT_NO_TOOLTIP
#ifndef QT_NO_TOOLTIP
        pruneMCDistCombo1->setToolTip(QApplication::translate("MainWindow", "choose distance metric used for pruning", 0));
#endif // QT_NO_TOOLTIP
        exportHSGroup->setTitle(QApplication::translate("MainWindow", "export skeleton", 0));
#ifndef QT_NO_TOOLTIP
        exportHSBtn->setToolTip(QApplication::translate("MainWindow", "click to export the current skeleton to a .sk file", 0));
#endif // QT_NO_TOOLTIP
        exportHSBtn->setText(QApplication::translate("MainWindow", "export", 0));
#ifndef QT_NO_TOOLTIP
        hsPruneGroup->setToolTip(QApplication::translate("MainWindow", "explore pruning params (theta2 and theta1) for skeleton generation", 0));
#endif // QT_NO_TOOLTIP
        hsPruneGroup->setTitle(QApplication::translate("MainWindow", "prune skeleton", 0));
        label_5->setText(QApplication::translate("MainWindow", "<html><head/><body><p><span style=\" font-family:'Symbol'; font-size:12pt;\">q</span><span style=\" font-size:12pt; vertical-align:sub;\">2</span></p></body></html>", 0));
#ifndef QT_NO_TOOLTIP
        removeSmallCmpnts->setToolTip(QApplication::translate("MainWindow", "tick this if want to remove small components", 0));
#endif // QT_NO_TOOLTIP
        removeSmallCmpnts->setText(QApplication::translate("MainWindow", "clean debris", 0));
#ifndef QT_NO_TOOLTIP
        hsBt2Bt3Slider->setToolTip(QApplication::translate("MainWindow", "theta_2", 0));
#endif // QT_NO_TOOLTIP
        label_7->setText(QApplication::translate("MainWindow", "<html><head/><body><p><span style=\" font-family:'Symbol'; font-size:12pt;\">q</span><span style=\" font-size:12pt; vertical-align:sub;\">1</span></p></body></html>", 0));
#ifndef QT_NO_TOOLTIP
        hsFaceDiffSpin->setToolTip(QApplication::translate("MainWindow", "theta_2", 0));
#endif // QT_NO_TOOLTIP
#ifndef QT_NO_TOOLTIP
        hsBt1Bt2Slider->setToolTip(QApplication::translate("MainWindow", "theta_1", 0));
#endif // QT_NO_TOOLTIP
#ifndef QT_NO_TOOLTIP
        hsCurveDiffSpin->setToolTip(QApplication::translate("MainWindow", "theta_1", 0));
#endif // QT_NO_TOOLTIP
#ifndef QT_NO_TOOLTIP
        hsFaceDegenThreshold->setToolTip(QApplication::translate("MainWindow", "specify the area threshold. Components with smaller size will be removed", 0));
#endif // QT_NO_TOOLTIP
#ifndef QT_NO_TOOLTIP
        visHSBtn->setToolTip(QApplication::translate("MainWindow", "update and visualize skeleton with the specified parameters", 0));
#endif // QT_NO_TOOLTIP
        visHSBtn->setText(QApplication::translate("MainWindow", "create", 0));
        hsCreateGroup->setTitle(QApplication::translate("MainWindow", "precompute", 0));
#ifndef QT_NO_TOOLTIP
        createHSBtn->setToolTip(QApplication::translate("MainWindow", "<html><head/><body>precomputation: <br>- computing medial curve, <br>- dualizing subdivisions induced by medial curve, <br>- assemblying skeleton</body></html>", 0));
#endif // QT_NO_TOOLTIP
        createHSBtn->setText(QApplication::translate("MainWindow", "precompute", 0));
        toolTabs->setTabText(toolTabs->indexOf(HSTab), QApplication::translate("MainWindow", "Skeleton", 0));
#ifndef QT_NO_TOOLTIP
        MCTab->setToolTip(QApplication::translate("MainWindow", "medial curve (generate, burn, prune, etc.)", 0));
#endif // QT_NO_TOOLTIP
#ifndef QT_NO_TOOLTIP
        mcDualizeGroup->setToolTip(QApplication::translate("MainWindow", "click the button to create medial curve and compute the 1d burn time on it", 0));
#endif // QT_NO_TOOLTIP
        mcDualizeGroup->setTitle(QApplication::translate("MainWindow", "dualize", 0));
#ifndef QT_NO_TOOLTIP
        dualizeBtn->setToolTip(QString());
#endif // QT_NO_TOOLTIP
        dualizeBtn->setText(QApplication::translate("MainWindow", "create", 0));
#ifndef QT_NO_TOOLTIP
        mpGroup->setToolTip(QApplication::translate("MainWindow", "medial point (the point last burned away during burning the medial curve)", 0));
#endif // QT_NO_TOOLTIP
        mpGroup->setTitle(QApplication::translate("MainWindow", "MP", 0));
        showMPBox->setText(QApplication::translate("MainWindow", "show MP", 0));
        toolTabs->setTabText(toolTabs->indexOf(MCTab), QApplication::translate("MainWindow", "MC/MP", 0));
#ifndef QT_NO_TOOLTIP
        ExportTab->setToolTip(QApplication::translate("MainWindow", "export (burn time, skeleton)", 0));
#endif // QT_NO_TOOLTIP
        exportGroup->setTitle(QString());
        ExportSkelBox->setText(QApplication::translate("MainWindow", "Export Skeleton", 0));
        exportBtn->setText(QApplication::translate("MainWindow", "Export", 0));
        toolTabs->setTabText(toolTabs->indexOf(ExportTab), QApplication::translate("MainWindow", "Export", 0));
#ifndef QT_NO_TOOLTIP
        SMetricTab->setToolTip(QApplication::translate("MainWindow", "surface function", 0));
#endif // QT_NO_TOOLTIP
        surfFuncGroup->setTitle(QApplication::translate("MainWindow", "vis", 0));
        surfFuncProjectBtn->setText(QApplication::translate("MainWindow", "project func", 0));
#ifndef QT_NO_TOOLTIP
        surfFuncCorrespSchemeCombo->setToolTip(QApplication::translate("MainWindow", "choose which scheme to use to compute correspondence between MC and orig. surf.", 0));
#endif // QT_NO_TOOLTIP
        diffuseSurfFunc->setText(QApplication::translate("MainWindow", "fill \"holes\"", 0));
        surfFuncSetupCorrespBtn->setText(QApplication::translate("MainWindow", "setup corresp.", 0));
        toolTabs->setTabText(toolTabs->indexOf(SMetricTab), QApplication::translate("MainWindow", "S. Metric", 0));
#ifndef QT_NO_TOOLTIP
        VideoTab->setToolTip(QApplication::translate("MainWindow", "animation (for burning, iso-surf, iso-cont, etc.)", 0));
#endif // QT_NO_TOOLTIP
        isoSurfGroup->setTitle(QApplication::translate("MainWindow", "iso surface", 0));
        isoSurfPrecomputeBtn->setText(QApplication::translate("MainWindow", "pre-compute", 0));
        hideIsoSurf->setText(QApplication::translate("MainWindow", "hide", 0));
        recursiveRefineBox->setText(QApplication::translate("MainWindow", "recursive refine", 0));
        hideSnappedFaceBox->setText(QApplication::translate("MainWindow", "hide snapped faces", 0));
        isoSurfRefineBtn->setText(QApplication::translate("MainWindow", "refine", 0));
        groupBox_7->setTitle(QApplication::translate("MainWindow", "iso contour", 0));
        isoContPrecomputeBtn->setText(QApplication::translate("MainWindow", "pre-compute", 0));
        hideIsoCont->setText(QApplication::translate("MainWindow", "hide contour", 0));
#ifndef QT_NO_TOOLTIP
        isoContShowMC->setToolTip(QApplication::translate("MainWindow", "show MC exposed by iso-contour", 0));
#endif // QT_NO_TOOLTIP
        isoContShowMC->setText(QApplication::translate("MainWindow", "show MC", 0));
        hidePastSurfBox->setText(QApplication::translate("MainWindow", "hide past surface", 0));
        groupBox_9->setTitle(QApplication::translate("MainWindow", "transparency", 0));
        groupBox_8->setTitle(QApplication::translate("MainWindow", "play all", 0));
#ifndef QT_NO_TOOLTIP
        evolveAllSlider->setToolTip(QApplication::translate("MainWindow", "drag this bar to evolve all iso-objects", 0));
#endif // QT_NO_TOOLTIP
        evolveAllPrecomputeBtn->setText(QApplication::translate("MainWindow", "pre-compute", 0));
#ifndef QT_NO_TOOLTIP
        evolveAllShowEntireMC->setToolTip(QApplication::translate("MainWindow", "Only exposed MC is shown if unchecked.", 0));
#endif // QT_NO_TOOLTIP
        evolveAllShowEntireMC->setText(QApplication::translate("MainWindow", "show entire MC", 0));
        toolTabs->setTabText(toolTabs->indexOf(VideoTab), QApplication::translate("MainWindow", "Video", 0));
#ifndef QT_NO_TOOLTIP
        MiscTab->setToolTip(QApplication::translate("MainWindow", "debug related stuff", 0));
#endif // QT_NO_TOOLTIP
        preprocessGroup->setTitle(QApplication::translate("MainWindow", "preprocess", 0));
#ifndef QT_NO_TOOLTIP
        label_2->setToolTip(QApplication::translate("MainWindow", "closeness for 2 vts to be for them to be considered \"same\"", 0));
#endif // QT_NO_TOOLTIP
        label_2->setText(QApplication::translate("MainWindow", "epsilon", 0));
#ifndef QT_NO_TOOLTIP
        label->setToolTip(QApplication::translate("MainWindow", "perturbation applied to \"same\" points", 0));
#endif // QT_NO_TOOLTIP
        label->setText(QApplication::translate("MainWindow", "pert", 0));
        miscBurntimeGroup->setTitle(QApplication::translate("MainWindow", "burntime", 0));
        visSphereFromFileBtn->setText(QApplication::translate("MainWindow", "show from file", 0));
        randBTBtn->setText(QApplication::translate("MainWindow", "randomize", 0));
        printBTBtn->setText(QApplication::translate("MainWindow", "print random", 0));
        showSphereBox->setText(QApplication::translate("MainWindow", "draw spheres", 0));
        printSelectVertInfoBtn->setText(QApplication::translate("MainWindow", "print selected", 0));
        renderParamGroup->setTitle(QApplication::translate("MainWindow", "general rendering params", 0));
        label_3->setText(QApplication::translate("MainWindow", "lines offset", 0));
        qmatGroup->setTitle(QApplication::translate("MainWindow", "qmat", 0));
        label_9->setText(QApplication::translate("MainWindow", "to qmat input", 0));
        outputWenpingBtn->setText(QApplication::translate("MainWindow", "convert", 0));
        readQMatMA->setText(QApplication::translate("MainWindow", "read qmat ma", 0));
        hideQMatBox->setText(QApplication::translate("MainWindow", "hide", 0));
        toolTabs->setTabText(toolTabs->indexOf(MiscTab), QApplication::translate("MainWindow", "Misc", 0));
#ifndef QT_NO_TOOLTIP
        pointSizeSpin->setToolTip(QApplication::translate("MainWindow", "change size of points", 0));
#endif // QT_NO_TOOLTIP
        doMAFaceScalarDiffusion->setText(QApplication::translate("MainWindow", "diffuse face field", 0));
#ifndef QT_NO_TOOLTIP
        computeAllDistBtn->setToolTip(QApplication::translate("MainWindow", "compute/update all distance metrics", 0));
#endif // QT_NO_TOOLTIP
        computeAllDistBtn->setText(QApplication::translate("MainWindow", "compute all", 0));
        usePerSheetBox->setText(QApplication::translate("MainWindow", "per sheet", 0));
#ifndef QT_NO_TOOLTIP
        maPruneRatio->setToolTip(QApplication::translate("MainWindow", "ratio: cur prune value / bbox size", 0));
#endif // QT_NO_TOOLTIP
#ifndef QT_NO_TOOLTIP
        minPruneDistMASpin1->setToolTip(QApplication::translate("MainWindow", "show/set minimal allowable dist value for pruning", 0));
#endif // QT_NO_TOOLTIP
#ifndef QT_NO_TOOLTIP
        maxPruneDistMASpin1->setToolTip(QApplication::translate("MainWindow", "show/set maximal allowable dist value for pruning", 0));
#endif // QT_NO_TOOLTIP
#ifndef QT_NO_TOOLTIP
        edgeDualOptCombo->setToolTip(QApplication::translate("MainWindow", "pick a method for dualizing edges", 0));
#endif // QT_NO_TOOLTIP
#ifndef QT_NO_TOOLTIP
        polyDualOptCombo->setToolTip(QApplication::translate("MainWindow", "pick a method for dualizing poly regions", 0));
#endif // QT_NO_TOOLTIP
        underEstimateBox->setText(QApplication::translate("MainWindow", "under estimate", 0));
#ifndef QT_NO_TOOLTIP
        onlyUnburntBox->setToolTip(QApplication::translate("MainWindow", "only show unburnt part of the dual structure", 0));
#endif // QT_NO_TOOLTIP
        onlyUnburntBox->setText(QApplication::translate("MainWindow", "only unburnt", 0));
#ifndef QT_NO_TOOLTIP
        stopAtJuncBox->setToolTip(QApplication::translate("MainWindow", "stop burning the dual line structure at junction points", 0));
#endif // QT_NO_TOOLTIP
        stopAtJuncBox->setText(QApplication::translate("MainWindow", "stop at junction", 0));
        protectBT2Box->setText(QApplication::translate("MainWindow", "protect BT2", 0));
        label_8->setText(QApplication::translate("MainWindow", "line rel diff", 0));
        label_6->setText(QApplication::translate("MainWindow", "face rel diff", 0));
#ifndef QT_NO_TOOLTIP
        mcBurnGroup->setToolTip(QApplication::translate("MainWindow", "options for burning the dual structure", 0));
#endif // QT_NO_TOOLTIP
        mcBurnGroup->setTitle(QApplication::translate("MainWindow", "burn", 0));
        burnMedialCurveBtn->setText(QApplication::translate("MainWindow", "burn", 0));
        printMCStatsBtn->setText(QApplication::translate("MainWindow", "print stats", 0));
#ifndef QT_NO_TOOLTIP
        mcPruneRatio->setToolTip(QApplication::translate("MainWindow", "distance value restricted to [0, max]", 0));
#endif // QT_NO_TOOLTIP
        colorPerVert->setText(QApplication::translate("MainWindow", "color-per-vert", 0));
        colorPerFace->setText(QApplication::translate("MainWindow", "color-per-face", 0));
#ifndef QT_NO_TOOLTIP
        hsUseTypedInput->setToolTip(QApplication::translate("MainWindow", "use theta_2 and theta_1 in the right fields as params", 0));
#endif // QT_NO_TOOLTIP
        hsUseTypedInput->setText(QApplication::translate("MainWindow", "use typed input", 0));
        visBallStick->setText(QApplication::translate("MainWindow", "ball-stick", 0));
        preserveMCTopo->setText(QApplication::translate("MainWindow", "preserve topo", 0));
        toolTabs->setTabText(toolTabs->indexOf(hiddenTab), QApplication::translate("MainWindow", "hidden", 0));
    } // retranslateUi

};

namespace Ui {
    class MainWindow: public Ui_MainWindow {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_COMPACT_H
