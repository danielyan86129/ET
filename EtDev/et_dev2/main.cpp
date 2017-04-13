#include "ETMainWindow.h"
#include <QApplication>
#include <gflags/gflags.h>
#include <string>

/// microsoft's memory leak detection
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
DEFINE_int32( burn_sch, 0, "the burning scheme (0: orig-and-steiner (default value), 1: steiner-only). OPTIONAL." );
//DEFINE_string( et_file, "", "the output file for the computed ET value per vertex on medial axis. OPTOINAL." );

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

	if ( FLAGS_nogui )
	{
		w.setInputs( FLAGS_ma_file, FLAGS_shape_file, FLAGS_r_file, FLAGS_omega, FLAGS_burn_sch );
		w.onLoadFilesClicked(); // load in input files
		w.onBurnBtnClicked(); // burn the medial axis
		w.onCleanTopoBtnClicked(); // check if there is closed pockets
		w.onExportETBtnClicked(); // export ET per-vertex on ma to a file
		return 0;
	}

	w.show();

	/// microsoft's memory leak detection
	//_CrtDumpMemoryLeaks();
	//_CrtSetDbgFlag( _CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF );
	//_CrtSetReportMode( _CRT_ERROR, _CRTDBG_MODE_DEBUG );

	return a.exec();
}
