#ifndef COMISODLLEXPORT
	#ifdef WIN32
		#ifdef COMISOMDLL
			#ifdef USECOMISO
				#define COMISODLLEXPORT __declspec(dllimport)
				#define COMISODLLEXPORTONLY 
			#else
				#define COMISODLLEXPORT __declspec(dllexport)
				#define COMISODLLEXPORTONLY __declspec(dllexport)
			#endif
		#else		
			#define COMISODLLEXPORT	
			#define COMISODLLEXPORTONLY
		#endif
	#else
		#define COMISODLLEXPORT
		#define COMISODLLEXPORTONLY
	#endif
#endif

#undef min
#undef max


