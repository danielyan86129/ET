#include "ETMainWindow.h"
#include <QApplication>
#include <gflags/gflags.h>
#include <string>

/// Microsoft's memory leak detection
//#define _CRTDBG_MAP_ALLOC
//#include <stdlib.h>
//#include <crtdbg.h>
//#include <vld.h>

/* flags parsed from cmd line */
DEFINE_bool( nogui, false, "use cmd line instead? OPTIONAL." );
DEFINE_string( shape_file, "", "the file for the 3d shape boundary. OPTIONAL." );
DEFINE_string( ma_file, "", "the file for the medial axis. REQUIRED." );
DEFINE_string( r_file, "", "the file for the radius function over medial axis. OPTIONAL." );
DEFINE_double( omega, 0.004, "the sampling rate for steiner subdivision. OPTIONAL." );
DEFINE_int32( burn_sch, 1, "the burning scheme (0: orig-and-steiner, 1: steiner-only(default value) ). OPTIONAL." );
//DEFINE_string( et_file, "", "the output file for the computed ET value per vertex on medial axis. OPTOINAL." );

DEFINE_string( mc_msure, "null", "output specified measure on medial curves. nothing to output by default. \
							valid value: {shape-diam, shape-width, ...}. OPTIONAL." );
DEFINE_bool( run_all, false, "compute all the way to skeleton generation." );
DEFINE_bool( export_skel, false, "compute all the way to curve skeleton generation & output skel files." );
DEFINE_double( theta_2, 1.1, "skeleton threshold for face pruning. Assume unit bounding box. All faces purged when > 1. OPTIONAL." );
DEFINE_double( theta_1, 0.05, "skeleton threshold for curve pruning. Assume unit bounding box. All curves purged when > 1. OPTIONAL." );

int main(int argc, char *argv[])
{
	std::string usage("Sample usage: \n");
	(usage += argv[ 0 ]) += " ";
	google::SetUsageMessage( usage );
	google::ParseCommandLineFlags( &argc, &argv, false );
	
	QApplication a( argc, argv );
	ETMainWindow w;

	// set the program mode
	w.setDebugMode(false);
	w.setConsoleMode( FLAGS_nogui );
	w.setWindowTitle("erosion thickness");
	w.setInputs( 
		FLAGS_ma_file, FLAGS_shape_file, FLAGS_r_file, FLAGS_omega, 
		FLAGS_mc_msure, 
		FLAGS_theta_2, FLAGS_theta_1,
		FLAGS_burn_sch );

	if ( FLAGS_nogui )
	{
		w.onLoadFilesClicked(); // load in input files
		w.onBurnBtnClicked(); // burn the medial axis
		w.onCleanTopoBtnClicked(); // check if there is closed pockets
		w.onExportETBtnClicked(); // export ET per-vertex on ma to a file
		
		if ( FLAGS_export_skel || FLAGS_run_all )
		{
			w.onCreateHSClicked(); // generate skeleton
			w.onVisHSClicked(); // prune skeleton
			if ( FLAGS_export_skel )
				w.onExportHSBtnClicked(); // export skeleton
		}
		else if ( !google::GetCommandLineFlagInfoOrDie( "mc_msure" ).is_default )
			w.onDualizeClicked();
		return 0;
	}

	w.show();

	/// Microsoft's memory leak detection
	//_CrtDumpMemoryLeaks();
	//_CrtSetDbgFlag( _CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF );
	//_CrtSetReportMode( _CRT_ERROR, _CRTDBG_MODE_DEBUG );

	return a.exec();
}
