/*
unpackJPG v2.5h (12/07/2013)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Copyright (c) 2006...2013, HTW Aalen University and Matthias Stirner
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

____________________________________
unpackJPG by Matthias Stirner, 12/2013
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctime>

#include "bitops.h"
#include "aricoder.h"
#include "pjpgtbl.h"
#include "dct8x8.h"

#if defined BUILD_DLL // define BUILD_LIB from the compiler options if you want to compile a DLL!
	#define BUILD_LIB
#endif

#if defined BUILD_LIB // define BUILD_LIB as compiler option if you want to compile a library!
	#include "packjpglib.h"
#endif

#define INTERN static

#define INIT_MODEL_S(a,b,c) new model_s( a, b, c, 255 )
#define INIT_MODEL_B(a,b)   new model_b( a, b, 255 )

#define QUANT(cm,bp)	( cmpnfo[cm].qtable[ bp ] )
#define MAX_V(cm,bp)	( ( QUANT(cm,bp) > 0 ) ? ( ( freqmax[bp] + QUANT(cm,bp) - 1 ) /  QUANT(cm,bp) ) : 0 )
// #define QUN_V(v,cm,bp)	( ( QUANT(cm,bp) > 0 ) ? ( ( v > 0 ) ? ( v + (QUANT(cm,bp)/2) ) /  QUANT(cm,bp) : ( v - (QUANT(cm,bp)/2) ) /  QUANT(cm,bp) ) : 0 )

#define ENVLI(s,v)		( ( v > 0 ) ? v : ( v - 1 ) + ( 1 << s ) )
#define DEVLI(s,n)		( ( n >= ( 1 << (s - 1) ) ) ? n : n + 1 - ( 1 << s ) )
#define E_ENVLI(s,v)	( v - ( 1 << s ) )
#define E_DEVLI(s,n)	( n + ( 1 << s ) )

#define ABS(v1)			( (v1 < 0) ? -v1 : v1 )
#define ABSDIFF(v1,v2)	( (v1 > v2) ? (v1 - v2) : (v2 - v1) )
#define IPOS(w,v,h)		( ( v * w ) + h )
#define NPOS(n1,n2,p)	( ( ( p / n1 ) * n2 ) + ( p % n1 ) )
#define ROUND_F(v1)		( (v1 < 0) ? (int) (v1 - 0.5) : (int) (v1 + 0.5) )
#define DIV_INT(v1,v2)	( (v1 < 0) ? (v1 - (v2>>1)) / v2 : (v1 + (v2>>1)) / v2 )
#define B_SHORT(v1,v2)	( ( ((int) v1) << 8 ) + ((int) v2) )
#define BITLEN1024P(v)	( pbitlen_0_1024[ v ] )
#define BITLEN2048N(v)	( (pbitlen_n2048_2047+2048)[ v ] )
#define CLAMPED(l,h,v)	( ( v < l ) ? l : ( v > h ) ? h : v )

#define MEM_ERRMSG	"out of memory error"
#define FRD_ERRMSG	"could not read file / file not found: %s"
#define FWR_ERRMSG	"could not write file / file write-protected: %s"
#define MSG_SIZE	128
#define BARLEN		36


/* -----------------------------------------------
	struct declarations
	----------------------------------------------- */

struct componentInfo {
	unsigned short* qtable; // quantization table
	int huffdc; // no of huffman table (DC)
	int huffac; // no of huffman table (AC)
	int sfv; // sample factor vertical
	int sfh; // sample factor horizontal	
	int mbs; // blocks in mcu
	int bcv; // block count vertical (interleaved)
	int bch; // block count horizontal (interleaved)
	int bc;  // block count (all) (interleaved)
	int ncv; // block count vertical (non interleaved)
	int nch; // block count horizontal (non interleaved)
	int nc;  // block count (all) (non interleaved)
	int sid; // statistical identity
	int jid; // jpeg internal id
};

struct huffCodes {
	unsigned short cval[ 256 ];
	unsigned short clen[ 256 ];
	unsigned short max_eobrun;
};

struct huffTree {
	unsigned short l[ 256 ];
	unsigned short r[ 256 ];
};


/* -----------------------------------------------
	function declarations: main interface
	----------------------------------------------- */
#if !defined( BUILD_LIB )
INTERN void initialize_options( int argc, char** argv );
INTERN void process_ui( void );
INTERN inline const char* get_status( bool (*function)() );
INTERN void show_help( void );
#endif
INTERN void process_file( void );
INTERN void execute( bool (*function)() );


/* -----------------------------------------------
	function declarations: main functions
	----------------------------------------------- */
#if !defined( BUILD_LIB )
INTERN bool check_file( void );
#endif
INTERN bool reset_buffers( void );
INTERN bool merge_jpeg( void );
INTERN bool recode_jpeg( void );
INTERN bool adapt_icos( void );
INTERN bool unpredict_dc( void );
INTERN bool calc_zdst_lists( void );
INTERN bool unpack_pjg( void );


/* -----------------------------------------------
	function declarations: jpeg-specific
	----------------------------------------------- */

INTERN bool jpg_setup_imginfo( void );
INTERN bool jpg_parse_jfif( unsigned char type, unsigned int len, unsigned char* segment );
INTERN bool jpg_rebuild_header( void );

INTERN int jpg_encode_block_seq( abitwriter* huffw, huffCodes* dctbl, huffCodes* actbl, short* block );
INTERN int jpg_encode_dc_prg_fs( abitwriter* huffw, huffCodes* dctbl, short* block );
INTERN int jpg_encode_ac_prg_fs( abitwriter* huffw, huffCodes* actbl, short* block,
						int* eobrun, int from, int to );
INTERN int jpg_encode_dc_prg_sa( abitwriter* huffw, short* block );
INTERN int jpg_encode_ac_prg_sa( abitwriter* huffw, abytewriter* storw, huffCodes* actbl,
						short* block, int* eobrun, int from, int to );
INTERN int jpg_encode_eobrun( abitwriter* huffw, huffCodes* actbl, int* eobrun );
INTERN int jpg_encode_crbits( abitwriter* huffw, abytewriter* storw );

INTERN int jpg_next_mcupos( int* mcu, int* cmp, int* csc, int* sub, int* dpos, int* rstw );
INTERN int jpg_next_mcuposn( int* cmp, int* dpos, int* rstw );

INTERN void jpg_build_huffcodes( unsigned char *clen, unsigned char *cval,
				huffCodes *hc, huffTree *ht );

/* -----------------------------------------------
	function declarations: pjg-specific
	----------------------------------------------- */

INTERN bool pjg_decode_zstscan( aricoder* dec, int cmp );
INTERN bool pjg_decode_zdst_high( aricoder* dec, int cmp );
INTERN bool pjg_decode_zdst_low( aricoder* dec, int cmp );
INTERN bool pjg_decode_dc( aricoder* dec, int cmp );
INTERN bool pjg_decode_ac_high( aricoder* dec, int cmp );
INTERN bool pjg_decode_ac_low( aricoder* dec, int cmp );
INTERN bool pjg_decode_generic( aricoder* dec, unsigned char** data, int* len );
INTERN bool pjg_decode_bit( aricoder* dec, unsigned char* bit );

INTERN bool pjg_unoptimize_header( void );

INTERN void pjg_aavrg_prepare( unsigned short** abs_coeffs, int* weights, unsigned short* abs_store, int cmp );
INTERN int pjg_aavrg_context( unsigned short** abs_coeffs, int* weights, int pos, int p_y, int p_x, int r_x );
INTERN int pjg_lakh_context( signed short** coeffs_x, signed short** coeffs_a, int* pred_cf, int pos );
INTERN void get_context_nnb( int pos, int w, int *a, int *b );


/* -----------------------------------------------
	function declarations: DCT
	----------------------------------------------- */

INTERN int idct_2d_fst_1x8( int cmp, int dpos, int ix, int iy );
INTERN int idct_2d_fst_8x1( int cmp, int dpos, int ix, int iy );


/* -----------------------------------------------
	function declarations: prediction
	----------------------------------------------- */

INTERN int dc_1ddct_predictor( int cmp, int dpos );
INTERN inline int median_int( int* values, int size );
INTERN inline float median_float( float* values, int size );


/* -----------------------------------------------
	function declarations: miscelaneous helpers
	----------------------------------------------- */
#if !defined( BUILD_LIB )
INTERN inline void progress_bar( int current, int last );
INTERN inline char* create_filename( const char* base, const char* extension );
INTERN inline char* unique_filename( const char* base, const char* extension );
INTERN inline void set_extension( char* filename, const char* extension );
INTERN inline void add_underscore( char* filename );
#endif
INTERN inline bool file_exists( const char* filename );


/* -----------------------------------------------
	global variables: library only variables
	----------------------------------------------- */
#if defined(BUILD_LIB)
INTERN int lib_in_type  = -1;
INTERN int lib_out_type = -1;
#endif


/* -----------------------------------------------
	global variables: data storage
	----------------------------------------------- */

INTERN unsigned short qtables[4][64];				// quantization tables
INTERN huffCodes      hcodes[2][4];				// huffman codes
INTERN huffTree       htrees[2][4];				// huffman decoding trees
INTERN unsigned char  htset[2][4];					// 1 if huffman table is set

INTERN unsigned char* grbgdata		   =   NULL;	// garbage data
INTERN unsigned char* hdrdata          =   NULL;   // header data
INTERN unsigned char* huffdata         =   NULL;   // huffman coded data
INTERN int            hufs             =    0  ;   // size of huffman data
INTERN int            hdrs             =    0  ;   // size of header
INTERN int            grbs             =    0  ;   // size of garbage

INTERN unsigned int*  rstp             =   NULL;   // restart markers positions in huffdata
INTERN unsigned int*  scnp             =   NULL;   // scan start positions in huffdata
INTERN int            rstc             =    0  ;   // count of restart markers
INTERN int            scnc             =    0  ;   // count of scans
INTERN int            rsti             =    0  ;   // restart interval
INTERN char           padbit           =    -1 ;   // padbit (for huffman coding)
INTERN unsigned char* rst_err          =   NULL;   // number of wrong-set RST markers per scan

INTERN unsigned char* zdstdata[4]      = { NULL }; // zero distribution (# of non-zeroes) lists (for higher 7x7 block)
INTERN unsigned char* eobxhigh[4]      = { NULL }; // eob in x direction (for higher 7x7 block)
INTERN unsigned char* eobyhigh[4]      = { NULL }; // eob in y direction (for higher 7x7 block)
INTERN unsigned char* zdstxlow[4]		= { NULL }; // # of non zeroes for first row
INTERN unsigned char* zdstylow[4]		= { NULL }; // # of non zeroes for first collumn
INTERN signed short*  colldata[4][64]  = {{NULL}}; // collection sorted DCT coefficients

INTERN unsigned char* freqscan[4]      = { NULL }; // optimized order for frequency scans (only pointers to scans)
INTERN unsigned char  zsrtscan[4][64];				// zero optimized frequency scan

INTERN int adpt_idct_8x8[ 4 ][ 8 * 8 * 8 * 8 ];	// precalculated/adapted values for idct (8x8)
INTERN int adpt_idct_1x8[ 4 ][ 1 * 1 * 8 * 8 ];	// precalculated/adapted values for idct (1x8)
INTERN int adpt_idct_8x1[ 4 ][ 8 * 8 * 1 * 1 ];	// precalculated/adapted values for idct (8x1)


/* -----------------------------------------------
	global variables: info about image
	----------------------------------------------- */

// seperate info for each color component
INTERN componentInfo cmpnfo[ 4 ];

INTERN int cmpc        = 0; // component count
INTERN int imgwidth    = 0; // width of image
INTERN int imgheight   = 0; // height of image

INTERN int sfhm        = 0; // max horizontal sample factor
INTERN int sfvm        = 0; // max verical sample factor
INTERN int mcuv        = 0; // mcus per line
INTERN int mcuh        = 0; // mcus per collumn
INTERN int mcuc        = 0; // count of mcus


/* -----------------------------------------------
	global variables: info about current scan
	----------------------------------------------- */

INTERN int cs_cmpc      =   0  ; // component count in current scan
INTERN int cs_cmp[ 4 ]  = { 0 }; // component numbers  in current scan
INTERN int cs_from      =   0  ; // begin - band of current scan ( inclusive )
INTERN int cs_to        =   0  ; // end - band of current scan ( inclusive )
INTERN int cs_sah       =   0  ; // successive approximation bit pos high
INTERN int cs_sal       =   0  ; // successive approximation bit pos low
	

/* -----------------------------------------------
	global variables: info about files
	----------------------------------------------- */
	
INTERN char*  jpgfilename = NULL;	// name of JPEG file
INTERN char*  pjgfilename = NULL;	// name of PJG file
INTERN int    jpgfilesize;			// size of JPEG file
INTERN int    pjgfilesize;			// size of PJG file
INTERN int    jpegtype = 0;			// type of JPEG coding: 0->unknown, 1->sequential, 2->progressive
INTERN int    filetype;				// type of current file
INTERN iostream* str_in  = NULL;	// input stream
INTERN iostream* str_out = NULL;	// output stream

#if !defined(BUILD_LIB)
INTERN iostream* str_str = NULL;	// storage stream

INTERN char** filelist = NULL;		// list of files to process 
INTERN int    file_cnt = 0;			// count of files in list
INTERN int    file_no  = 0;			// number of current file

INTERN char** err_list = NULL;		// list of error messages 
INTERN int*   err_tp   = NULL;		// list of error types
#endif

#if defined(DEV_INFOS)
INTERN int    dev_size_hdr      = 0;
INTERN int    dev_size_cmp[ 4 ] = { 0 };
INTERN int    dev_size_zsr[ 4 ] = { 0 };
INTERN int    dev_size_dc[ 4 ]  = { 0 };
INTERN int    dev_size_ach[ 4 ] = { 0 };
INTERN int    dev_size_acl[ 4 ] = { 0 };
INTERN int    dev_size_zdh[ 4 ] = { 0 };
INTERN int    dev_size_zdl[ 4 ] = { 0 };
#endif


/* -----------------------------------------------
	global variables: messages
	----------------------------------------------- */

INTERN char errormessage [ MSG_SIZE ];
INTERN bool (*errorfunction)();
INTERN int  errorlevel;
// meaning of errorlevel:
// -1 -> wrong input
// 0 -> no error
// 1 -> warning
// 2 -> fatal error


/* -----------------------------------------------
	global variables: settings
	----------------------------------------------- */

#if !defined( BUILD_LIB )
INTERN int  verbosity  = -1;	// level of verbosity
INTERN bool overwrite  = false;	// overwrite files yes / no
INTERN bool wait_exit  = true;	// pause after finished yes / no
INTERN int  verify_lv  = 0;		// verification level ( none (0), simple (1), detailed output (2) )
INTERN int  err_tol    = 1;		// error threshold ( proceed on warnings yes (2) / no (1) )
INTERN bool disc_meta  = false;	// discard meta-info yes / no

INTERN bool developer  = false;	// allow developers functions yes/no
INTERN bool auto_set   = true;	// automatic find best settings yes/no
INTERN int  action = A_COMPRESS;// what to do with JPEG/PJG files

INTERN FILE*  msgout   = stdout;// stream for output of messages
INTERN bool   pipe_on  = false;	// use stdin/stdout instead of filelist
#else
INTERN int  err_tol    = 1;		// error threshold ( proceed on warnings yes (2) / no (1) )
INTERN bool disc_meta  = false;	// discard meta-info yes / no
INTERN bool auto_set   = true;	// automatic find best settings yes/no
INTERN int  action = A_COMPRESS;// what to do with JPEG/PJG files
#endif

INTERN unsigned char nois_trs[ 4 ] = {6,6,6,6}; // bit pattern noise threshold
INTERN unsigned char segm_cnt[ 4 ] = {10,10,10,10}; // number of segments
#if !defined( BUILD_LIB )
INTERN unsigned char orig_set[ 8 ] = { 0 }; // store array for settings
#endif


/* -----------------------------------------------
	global variables: info about program
	----------------------------------------------- */

INTERN const unsigned char appversion = 25;
INTERN const char*  subversion   = "h";
INTERN const char*  apptitle     = "unpackPJG";
INTERN const char*  appname      = "unpackpjg";
INTERN const char*  versiondate  = "12/07/2013";
INTERN const char*  author       = "Matthias Stirner / Se";
#if !defined(BUILD_LIB)
INTERN const char*  website      = "http://packjpg.encode.ru/";
INTERN const char*  email        = "packjpg (at) matthiasstirner.com";
INTERN const char*	copyright    = "2006-2016 HTW Aalen University & Matthias Stirner";
INTERN const char*  pjg_ext      = "pjg";
INTERN const char*  jpg_ext      = "jpg";
#endif
INTERN const char   pjg_magic[] = { 'J', 'S' };


/* -----------------------------------------------
	main-function
	----------------------------------------------- */

#if !defined(BUILD_LIB)
int main( int argc, char** argv )
{	
	sprintf( errormessage, "no errormessage specified" );
	
	clock_t begin, end;
	
	int error_cnt = 0;
	int warn_cnt  = 0;
	
	double acc_jpgsize = 0;
	double acc_pjgsize = 0;
	
	int kbps;
	double cr;
	double total;
	
	errorlevel = 0;
	
	
	// read options from command line
	initialize_options( argc, argv );
	
	// write program info to screen
	fprintf( msgout,  "\n--> %s v%i.%i%s (%s) by %s <--\n",
			apptitle, appversion / 10, appversion % 10, subversion, versiondate, author );
	fprintf( msgout, "Copyright %s\nAll rights reserved\n\n", copyright );
	
	// check if user input is wrong, show help screen if it is
	if ( ( file_cnt == 0 ) ||
		( ( !developer ) && ( (action != A_COMPRESS) || (!auto_set) || (verify_lv > 1) ) ) ) {
		show_help();
		return -1;
	}
	
	// display warning if not using automatic settings
	if ( !auto_set ) {
		fprintf( msgout,  " custom compression settings: \n" );
		fprintf( msgout,  " -------------------------------------------------\n" );
		fprintf( msgout,  " no of segments    ->  %3i[0] %3i[1] %3i[2] %3i[3]\n",
				segm_cnt[0], segm_cnt[1], segm_cnt[2], segm_cnt[3] );
		fprintf( msgout,  " noise threshold   ->  %3i[0] %3i[1] %3i[2] %3i[3]\n",
				nois_trs[0], nois_trs[1], nois_trs[2], nois_trs[3] );
		fprintf( msgout,  " -------------------------------------------------\n\n" );
	}
	
	// (re)set program has to be done first
	reset_buffers();
	
	// process file(s) - this is the main function routine
	begin = clock();
	for ( file_no = 0; file_no < file_cnt; file_no++ ) {	
		// process current file
		process_ui();
		// store error message and type if any
		if ( errorlevel > 0 ) {
			err_list[ file_no ] = (char*) calloc( MSG_SIZE, sizeof( char ) );
			err_tp[ file_no ] = errorlevel;
			if ( err_list[ file_no ] != NULL )
				strcpy( err_list[ file_no ], errormessage );
		}
		// count errors / warnings / file sizes
		if ( errorlevel >= err_tol ) error_cnt++;
		else {
			if ( errorlevel == 1 ) warn_cnt++;
			acc_jpgsize += jpgfilesize;
			acc_pjgsize += pjgfilesize;
		}
	}
	end = clock();
	
	// errors summary: only needed for -v2 or progress bar
	if ( ( verbosity == -1 ) || ( verbosity == 2 ) ) {
		// print summary of errors to screen
		if ( error_cnt > 0 ) {
			fprintf( stderr, "\n\nfiles with errors:\n" );
			fprintf( stderr, "------------------\n" );
			for ( file_no = 0; file_no < file_cnt; file_no++ ) {
				if ( err_tp[ file_no ] >= err_tol ) {
					fprintf( stderr, "%s (%s)\n", filelist[ file_no ], err_list[ file_no ] );
				}
			}
		}
		// print summary of warnings to screen
		if ( warn_cnt > 0 ) {
			fprintf( stderr, "\n\nfiles with warnings:\n" );
			fprintf( stderr, "------------------\n" );
			for ( file_no = 0; file_no < file_cnt; file_no++ ) {
				if ( err_tp[ file_no ] == 1 ) {
					fprintf( stderr, "%s (%s)\n", filelist[ file_no ], err_list[ file_no ] );
				}
			}
		}
	}
	
	// show statistics
	fprintf( msgout,  "\n\n-> %i file(s) processed, %i error(s), %i warning(s)\n",
		file_cnt, error_cnt, warn_cnt );
	if ( ( file_cnt > error_cnt ) && ( verbosity != 0 ) &&
	 ( action == A_COMPRESS ) ) {
		acc_jpgsize /= 1024.0; acc_pjgsize /= 1024.0;
		total = (double) ( end - begin ) / CLOCKS_PER_SEC; 
		kbps  = ( total > 0 ) ? ( acc_jpgsize / total ) : acc_jpgsize;
		cr    = ( acc_jpgsize > 0 ) ? ( 100.0 * acc_pjgsize / acc_jpgsize ) : 0;
		
		fprintf( msgout,  " --------------------------------- \n" );
		if ( total >= 0 ) {
			fprintf( msgout,  " total time        : %8.2f sec\n", total );
			fprintf( msgout,  " avrg. kbyte per s : %8i byte\n", kbps );
		}
		else {
			fprintf( msgout,  " total time        : %8s sec\n", "N/A" );
			fprintf( msgout,  " avrg. kbyte per s : %8s byte\n", "N/A" );
		}
		fprintf( msgout,  " avrg. comp. ratio : %8.2f %%\n", cr );		
		fprintf( msgout,  " --------------------------------- \n" );
		#if defined(DEV_INFOS)
		if ( acc_jpgsize > 0 ) { 
			fprintf( msgout,  " header %%          : %8.2f %%\n", 100.0 * dev_size_hdr / acc_jpgsize );
			if ( dev_size_cmp[0] > 0 ) fprintf( msgout,  " component [0] %%   : %8.2f %%\n", 100.0 * dev_size_cmp[0] / acc_jpgsize );
			if ( dev_size_cmp[1] > 0 ) fprintf( msgout,  " component [1] %%   : %8.2f %%\n", 100.0 * dev_size_cmp[1] / acc_jpgsize );
			if ( dev_size_cmp[2] > 0 ) fprintf( msgout,  " component [2] %%   : %8.2f %%\n", 100.0 * dev_size_cmp[2] / acc_jpgsize );
			if ( dev_size_cmp[3] > 0 ) fprintf( msgout,  " component [3] %%   : %8.2f %%\n", 100.0 * dev_size_cmp[3] / acc_jpgsize );
			fprintf( msgout,  " --------------------------------- \n" );
			for ( int i = 0; i < 4; i++ ) {
				if ( dev_size_cmp[i] == 0 ) break;
				fprintf( msgout,  " ac coeffs h [%i] %% : %8.2f %%\n", i, 100.0 * dev_size_ach[i] / acc_jpgsize );				
				fprintf( msgout,  " ac coeffs l [%i] %% : %8.2f %%\n", i, 100.0 * dev_size_acl[i] / acc_jpgsize );
				fprintf( msgout,  " dc coeffs   [%i] %% : %8.2f %%\n", i, 100.0 * dev_size_dc[i] / acc_jpgsize );
				fprintf( msgout,  " zero dist h [%i] %% : %8.2f %%\n", i, 100.0 * dev_size_zdh[i] / acc_jpgsize );
				fprintf( msgout,  " zero dist l [%i] %% : %8.2f %%\n", i, 100.0 * dev_size_zdl[i] / acc_jpgsize );
				fprintf( msgout,  " zero sort   [%i] %% : %8.2f %%\n", i, 100.0 * dev_size_zsr[i] / acc_jpgsize );
				fprintf( msgout,  " --------------------------------- \n" );
			}
		}
		#endif
	}
	
	// pause before exit
	if ( wait_exit && ( msgout != stderr ) ) {
		fprintf( msgout, "\n\n< press ENTER >\n" );
		fgetc( stdin );
	}
	
	
	return 0;
}
#endif

/* ----------------------- Begin of library only functions -------------------------- */

/* -----------------------------------------------
	DLL export converter function
	----------------------------------------------- */
	
#if defined(BUILD_LIB)
EXPORT bool pjglib_convert_stream2stream( char* msg )
{
	// process in main function
	return pjglib_convert_stream2mem( NULL, NULL, msg ); 
}
#endif


/* -----------------------------------------------
	DLL export converter function
	----------------------------------------------- */

#if defined(BUILD_LIB)
EXPORT bool pjglib_convert_file2file( char* in, char* out, char* msg )
{
	// init streams
	pjglib_init_streams( (void*) in, 0, 0, (void*) out, 0 );
	
	// process in main function
	return pjglib_convert_stream2mem( NULL, NULL, msg ); 
}
#endif


/* -----------------------------------------------
	DLL export converter function
	----------------------------------------------- */
	
#if defined(BUILD_LIB)
EXPORT bool pjglib_convert_stream2mem( unsigned char** out_file, unsigned int* out_size, char* msg )
{
	clock_t begin, end;
	int total;
	float cr;	
	
	
	// use automatic settings
	auto_set = true;
	
	// (re)set buffers
	reset_buffers();
	action = A_COMPRESS;
	
	// main compression / decompression routines
	begin = clock();
	
	// process one file
	process_file();
	
	// fetch pointer and size of output (only for memory output)
	if ( ( errorlevel < err_tol ) && ( lib_out_type == 1 ) &&
		 ( out_file != NULL ) && ( out_size != NULL ) ) {
		*out_size = str_out->getsize();
		*out_file = str_out->getptr();
	}
	
	// close iostreams
	if ( str_in  != NULL ) delete( str_in  ); str_in  = NULL;
	if ( str_out != NULL ) delete( str_out ); str_out = NULL;
	
	end = clock();
	
	// copy errormessage / remove files if error (and output is file)
	if ( errorlevel >= err_tol ) {
		if ( lib_out_type == 0 ) {
			if ( filetype == F_PJG ) {
				if ( file_exists( jpgfilename ) ) remove( jpgfilename );
			}
		}
		if ( msg != NULL ) strcpy( msg, errormessage );
		return false;
	}
	
	// get compression info
	total = (int) ( (double) (( end - begin ) * 1000) / CLOCKS_PER_SEC );
	cr    = ( jpgfilesize > 0 ) ? ( 100.0 * pjgfilesize / jpgfilesize ) : 0;
	
	// write success message else
	if ( msg != NULL ) {
		switch( filetype )
		{
			case F_PJG:
				sprintf( msg, "Decompressed to %s (%.2f%%) in %ims",
					jpgfilename, cr, ( total >= 0 ) ? total : -1 );
				break;
			case F_UNK:
				sprintf( msg, "Unknown filetype" );
				break;	
		}
	}
	
	
	return true;
}
#endif


/* -----------------------------------------------
	DLL export init input (file/mem)
	----------------------------------------------- */
	
#if defined(BUILD_LIB)
EXPORT void pjglib_init_streams( void* in_src, int in_type, int in_size, void* out_dest, int out_type )
{
	/* a short reminder about input/output stream types:
	
	if input is file
	----------------
	in_scr -> name of input file
	in_type -> 0
	in_size -> ignore
	
	if input is memory
	------------------
	in_scr -> array containg data
	in_type -> 1
	in_size -> size of data array
	
	if input is *FILE (f.e. stdin)
	------------------------------
	in_src -> stream pointer
	in_type -> 2
	in_size -> ignore
	
	vice versa for output streams! */
	
	unsigned char buffer[ 2 ];
	
	
	// (re)set errorlevel
	errorfunction = NULL;
	errorlevel = 0;
	jpgfilesize = 0;
	pjgfilesize = 0;
	
	// open input stream, check for errors
	str_in = new iostream( in_src, in_type, in_size, 0 );
	if ( str_in->chkerr() ) {
		sprintf( errormessage, "error opening input stream" );
		errorlevel = 2;
		return;
	}	
	
	// open output stream, check for errors
	str_out = new iostream( out_dest, out_type, 0, 1 );
	if ( str_out->chkerr() ) {
		sprintf( errormessage, "error opening output stream" );
		errorlevel = 2;
		return;
	}
	
	// free memory from filenames if needed
	if ( jpgfilename != NULL ) free( jpgfilename ); jpgfilename = NULL;
	if ( pjgfilename != NULL ) free( pjgfilename ); pjgfilename = NULL;
	
	// check input stream
	str_in->read( buffer, 1, 2 );
	if ( ( buffer[0] == 0xFF ) && ( buffer[1] == 0xD8 ) ) {
		// file is JPEG
		filetype = F_JPG;
		// copy filenames
		jpgfilename = (char*) calloc( (  in_type == 0 ) ? strlen( (char*) in_src   ) + 1 : 32, sizeof( char ) );
		pjgfilename = (char*) calloc( ( out_type == 0 ) ? strlen( (char*) out_dest ) + 1 : 32, sizeof( char ) );
		strcpy( jpgfilename, (  in_type == 0 ) ? (char*) in_src   : "JPG in memory" );
		strcpy( pjgfilename, ( out_type == 0 ) ? (char*) out_dest : "PJG in memory" );
	}
	else if ( (buffer[0] == pjg_magic[0]) && (buffer[1] == pjg_magic[1]) ) {
		// file is PJG
		filetype = F_PJG;
		// copy filenames
		pjgfilename = (char*) calloc( (  in_type == 0 ) ? strlen( (char*) in_src   ) + 1 : 32, sizeof( char ) );
		jpgfilename = (char*) calloc( ( out_type == 0 ) ? strlen( (char*) out_dest ) + 1 : 32, sizeof( char ) );
		strcpy( pjgfilename, (  in_type == 0 ) ? (char*) in_src   : "PJG in memory" );
		strcpy( jpgfilename, ( out_type == 0 ) ? (char*) out_dest : "JPG in memory" );
	}
	else {
		// file is neither
		filetype = F_UNK;
		sprintf( errormessage, "filetype of input stream is unknown" );
		errorlevel = 2;
		return;
	}
	
	// store types of in-/output
	lib_in_type  = in_type;
	lib_out_type = out_type;
}
#endif


/* -----------------------------------------------
	DLL export version information
	----------------------------------------------- */
	
#if defined(BUILD_LIB)
EXPORT const char* pjglib_version_info( void )
{
	static char v_info[ 256 ];
	
	// copy version info to string
	sprintf( v_info, "--> %s library v%i.%i%s (%s) by %s <--",
			apptitle, appversion / 10, appversion % 10, subversion, versiondate, author );
			
	return (const char*) v_info;
}
#endif


/* -----------------------------------------------
	DLL export version information
	----------------------------------------------- */
	
#if defined(BUILD_LIB)
EXPORT const char* pjglib_short_name( void )
{
	static char v_name[ 256 ];
	
	// copy version info to string
	sprintf( v_name, "%s v%i.%i%s",
			apptitle, appversion / 10, appversion % 10, subversion );
			
	return (const char*) v_name;
}
#endif

/* ----------------------- End of libary only functions -------------------------- */

/* ----------------------- Begin of main interface functions -------------------------- */


/* -----------------------------------------------
	reads in commandline arguments
	----------------------------------------------- */
	
#if !defined(BUILD_LIB)	
INTERN void initialize_options( int argc, char** argv )
{	
	int tmp_val;
	char** tmp_flp;
	int i;
	
	
	// get memory for filelist & preset with NULL
	filelist = (char**) calloc( argc, sizeof( char* ) );
	for ( i = 0; i < argc; i++ )
		filelist[ i ] = NULL;
	
	// preset temporary filelist pointer
	tmp_flp = filelist;
	
	
	// read in arguments
	while ( --argc > 0 ) {
		argv++;
		// switches begin with '-'
		if ( strcmp((*argv), "-d" ) == 0 ) {
			disc_meta = true;
		}		
		else if ( sscanf( (*argv), "-v%i", &tmp_val ) == 1 ){
			verbosity = tmp_val;
			verbosity = ( verbosity < 0 ) ? 0 : verbosity;
			verbosity = ( verbosity > 2 ) ? 2 : verbosity;			
		}
		else if ( strcmp((*argv), "-vp" ) == 0 ) {
			verbosity = -1;
		}
		else if ( strcmp((*argv), "-np" ) == 0 ) {
			wait_exit = false;
		}
		else if ( strcmp((*argv), "-o" ) == 0 ) {
			overwrite = true;
		}
		else if ( strcmp((*argv), "-") == 0 ) {
			// switch standard message out stream
			msgout = stderr;
			// use "-" as placeholder for stdin
			*(tmp_flp++) = (char*) "-";
		}
		else {
			// if argument is not switch, it's a filename
			*(tmp_flp++) = *argv;
		}		
	}
	
	// count number of files (or filenames) in filelist
	for ( file_cnt = 0; filelist[ file_cnt ] != NULL; file_cnt++ );
	
	// alloc arrays for error messages and types storage
	err_list = (char**) calloc( file_cnt, sizeof( char* ) );
	err_tp   = (int*) calloc( file_cnt, sizeof( int ) );
	
	// backup settings - needed to restore original setting later
	if ( !auto_set ) {
		orig_set[ 0 ] = nois_trs[ 0 ];
		orig_set[ 1 ] = nois_trs[ 1 ];
		orig_set[ 2 ] = nois_trs[ 2 ];
		orig_set[ 3 ] = nois_trs[ 3 ];
		orig_set[ 4 ] = segm_cnt[ 0 ];
		orig_set[ 5 ] = segm_cnt[ 1 ];
		orig_set[ 6 ] = segm_cnt[ 2 ];
		orig_set[ 7 ] = segm_cnt[ 3 ];
	}
	else {
		for ( i = 0; i < 8; i++ )
			orig_set[ i ] = 0;
	}	
}
#endif


/* -----------------------------------------------
	UI for processing one file
	----------------------------------------------- */
	
#if !defined(BUILD_LIB)
INTERN void process_ui( void )
{
	clock_t begin, end;
	const char* actionmsg  = NULL;
	const char* errtypemsg = NULL;
	int total, bpms;
	float cr;	
	
	
	errorfunction = NULL;
	errorlevel = 0;
	jpgfilesize = 0;
	pjgfilesize = 0;	
	#if !defined(DEV_BUILD)
	action = A_COMPRESS;
	#endif
	
	// compare file name, set pipe if needed
	if ( ( strcmp( filelist[ file_no ], "-" ) == 0 ) && ( action == A_COMPRESS ) ) {
		pipe_on = true;
		filelist[ file_no ] = (char*) "STDIN";
	}
	else {		
		pipe_on = false;
	}
	
	if ( verbosity >= 0 ) { // standard UI
		fprintf( msgout,  "\nProcessing file %i of %i \"%s\" -> ",
					file_no + 1, file_cnt, filelist[ file_no ] );
		
		if ( verbosity > 1 )
			fprintf( msgout,  "\n----------------------------------------" );
		
		// check input file and determine filetype
		execute( check_file );
		
		// get specific action message
		if ( filetype == F_UNK ) actionmsg = "unknown filetype";
		else actionmsg = "Decompressing";
		
		if ( verbosity < 2 ) fprintf( msgout, "%s -> ", actionmsg );
	}
	else { // progress bar UI
		// update progress message
		fprintf( msgout, "Processing file %2i of %2i ", file_no + 1, file_cnt );
		progress_bar( file_no, file_cnt );
		fprintf( msgout, "\r" );
		execute( check_file );
	}
	fflush( msgout );
	
	
	// main function routine
	begin = clock();
	
	// streams are initiated, start processing file
	process_file();
	
	// close iostreams
	if ( str_in  != NULL ) delete( str_in  ); str_in  = NULL;
	if ( str_out != NULL ) delete( str_out ); str_out = NULL;
	if ( str_str != NULL ) delete( str_str ); str_str = NULL;
	// delete if broken or if output not needed
	if ( ( !pipe_on ) && ( ( errorlevel >= err_tol ) || ( action != A_COMPRESS ) ) ) {
		if ( filetype == F_PJG ) {
			if ( file_exists( jpgfilename ) ) remove( jpgfilename );
		}
	}
	
	end = clock();	
	
	// speed and compression ratio calculation
	total = (int) ( (double) (( end - begin ) * 1000) / CLOCKS_PER_SEC );
	bpms  = ( total > 0 ) ? ( jpgfilesize / total ) : jpgfilesize;
	cr    = ( jpgfilesize > 0 ) ? ( 100.0 * pjgfilesize / jpgfilesize ) : 0;

	
	if ( verbosity >= 0 ) { // standard UI
		if ( verbosity > 1 )
			fprintf( msgout,  "\n----------------------------------------" );
		
		// display success/failure message
		switch ( verbosity ) {
			case 0:			
				if ( errorlevel < err_tol ) {
					if ( action == A_COMPRESS ) fprintf( msgout,  "%.2f%%", cr );
					else fprintf( msgout, "DONE" );
				}
				else fprintf( msgout,  "ERROR" );
				if ( errorlevel > 0 ) fprintf( msgout,  "\n" );
				break;
			
			case 1:
				fprintf( msgout, "%s\n",  ( errorlevel < err_tol ) ? "DONE" : "ERROR" );
				break;
			
			case 2:
				if ( errorlevel < err_tol ) fprintf( msgout,  "\n-> %s OK\n", actionmsg );
				else  fprintf( msgout,  "\n-> %s ERROR\n", actionmsg );
				break;
		}
		
		// set type of error message
		switch ( errorlevel ) {
			case 0:	errtypemsg = "none"; break;
			case 1: errtypemsg = ( err_tol > 1 ) ?  "warning (ignored)" : "warning (skipped file)"; break;
			case 2: errtypemsg = "fatal error"; break;
		}
		
		// error/ warning message
		if ( errorlevel > 0 ) {
			fprintf( msgout, " %s -> %s:\n", get_status( errorfunction ), errtypemsg  );
			fprintf( msgout, " %s\n", errormessage );
		}
		if ( (verbosity > 0) && (errorlevel < err_tol) && (action == A_COMPRESS) ) {
			if ( total >= 0 ) {
				fprintf( msgout,  " time taken  : %7i msec\n", total );
				fprintf( msgout,  " byte per ms : %7i byte\n", bpms );
			}
			else {
				fprintf( msgout,  " time taken  : %7s msec\n", "N/A" );
				fprintf( msgout,  " byte per ms : %7s byte\n", "N/A" );
			}
			fprintf( msgout,  " comp. ratio : %7.2f %%\n", cr );		
		}	
		if ( ( verbosity > 1 ) && ( action == A_COMPRESS ) )
			fprintf( msgout,  "\n" );
	}
	else { // progress bar UI
		// if this is the last file, update progress bar one last time
		if ( file_no + 1 == file_cnt ) {
			// update progress message
			fprintf( msgout, "Processed %2i of %2i files ", file_no + 1, file_cnt );
			progress_bar( 1, 1 );
			fprintf( msgout, "\r" );
		}	
	}
}
#endif


/* -----------------------------------------------
	gets statusmessage for function
	----------------------------------------------- */
	
#if !defined(BUILD_LIB)
INTERN inline const char* get_status( bool (*function)() )
{	
	if ( function == NULL ) {
		return "unknown action";
	} else if ( function == *check_file ) {
		return "Determining filetype";
	} else if ( function == *merge_jpeg ) {
		return "Merging header & image data";
	} else if ( function == *recode_jpeg ) {
		return "Recompressing JPEG image data";
	} else if ( function == *adapt_icos ) {
		return "Adapting DCT precalc. tables";
	} else if ( function == *unpredict_dc ) {
		return "Removing prediction from DC";
	} else if ( function == *calc_zdst_lists ) {
		return "Calculating zero dist lists";
	} else if ( function == *unpack_pjg ) {
		return "Uncompressing data from PJG";
	} else if ( function == *reset_buffers ) {
		return "Resetting program";
	}
	else {
		return "Function description missing!";
	}
}
#endif


/* -----------------------------------------------
	shows help in case of wrong input
	----------------------------------------------- */
	
#if !defined(BUILD_LIB)
INTERN void show_help( void )
{	
	fprintf( msgout, "\n" );
	fprintf( msgout, "Website: %s\n", website );
	fprintf( msgout, "Email  : %s\n", email );
	fprintf( msgout, "\n" );
	fprintf( msgout, "Usage: %s [switches] [filename(s)]", appname );
	fprintf( msgout, "\n" );
	fprintf( msgout, "\n" );
	fprintf( msgout, " [-v?]    set level of verbosity (max: 2) (def: 0)\n" );
	fprintf( msgout, " [-np]    no pause after processing files\n" );
	fprintf( msgout, " [-o]     overwrite existing files\n" );
	fprintf( msgout, " [-d]     discard meta-info\n" );
	fprintf( msgout, "\n" );
	fprintf( msgout, "Examples: \"%s -v1 -o baboon.%s\"\n", appname, pjg_ext );
	fprintf( msgout, "          \"%s -p *.%s\"\n", appname, jpg_ext );	
}
#endif


/* -----------------------------------------------
	processes one file
	----------------------------------------------- */

INTERN void process_file( void )
{	
	if ( filetype == F_PJG )	{
		switch ( action )
		{
			case A_COMPRESS:
				execute( unpack_pjg );
				execute( adapt_icos );
				execute( unpredict_dc );
				execute( recode_jpeg );
				execute( merge_jpeg );
				break;
		}
	}
	// reset buffers
	reset_buffers();
}


/* -----------------------------------------------
	main-function execution routine
	----------------------------------------------- */

INTERN void execute( bool (*function)() )
{
	if ( errorlevel < err_tol ) {
		#if !defined BUILD_LIB
		clock_t begin, end;
		bool success;
		int total;
		
		// write statusmessage
		if ( verbosity == 2 ) {
			fprintf( msgout,  "\n%s ", get_status( function ) );
			for ( int i = strlen( get_status( function ) ); i <= 30; i++ )
				fprintf( msgout,  " " );			
		}
		
		// set starttime
		begin = clock();
		// call function
		success = ( *function )();
		// set endtime
		end = clock();
		
		if ( ( errorlevel > 0 ) && ( errorfunction == NULL ) )
			errorfunction = function;
		
		// write time or failure notice
		if ( success ) {
			total = (int) ( (double) (( end - begin ) * 1000) / CLOCKS_PER_SEC );
			if ( verbosity == 2 ) fprintf( msgout,  "%6ims", ( total >= 0 ) ? total : -1 );
		}
		else {
			errorfunction = function;
			if ( verbosity == 2 ) fprintf( msgout,  "%8s", "ERROR" );
		}
		#else
		// call function
		( *function )();
		
		// store errorfunction if needed
		if ( ( errorlevel > 0 ) && ( errorfunction == NULL ) )
			errorfunction = function;
		#endif
	}
}

/* ----------------------- End of main interface functions -------------------------- */

/* ----------------------- Begin of main functions -------------------------- */


/* -----------------------------------------------
	check file and determine filetype
	----------------------------------------------- */

#if !defined(BUILD_LIB)
INTERN bool check_file( void )
{	
	unsigned char fileid[ 2 ] = { 0, 0 };
	const char* filename = filelist[ file_no ];
	
	
	// open input stream, check for errors
	str_in = new iostream( (void*) filename, ( !pipe_on ) ? 0 : 2, 0, 0 );
	if ( str_in->chkerr() ) {
		sprintf( errormessage, FRD_ERRMSG, filename );
		errorlevel = 2;
		return false;
	}
	
	// free memory from filenames if needed
	if ( jpgfilename != NULL ) free( jpgfilename ); jpgfilename = NULL;
	if ( pjgfilename != NULL ) free( pjgfilename ); pjgfilename = NULL;
	
	// immediately return error if 2 bytes can't be read
	if ( str_in->read( fileid, 1, 2 ) != 2 ) { 
		filetype = F_UNK;
		sprintf( errormessage, "file doesn't contain enough data" );
		errorlevel = 2;
		return false;
	}
	
	// check file id, determine filetype
	if ( ( fileid[0] == pjg_magic[0] ) && ( fileid[1] == pjg_magic[1] ) ) {
		// file is PJG
		filetype = F_PJG;
		// create filenames
		if ( !pipe_on ) {
			pjgfilename = (char*) calloc( strlen( filename ) + 1, sizeof( char ) );
			strcpy( pjgfilename, filename );
			jpgfilename = ( overwrite ) ?
				create_filename( filename, (char*) jpg_ext ) :
				unique_filename( filename, (char*) jpg_ext );
		}
		else {
			jpgfilename = create_filename( "STDOUT", NULL );
			pjgfilename = create_filename( "STDIN", NULL );
		}
		// open output stream, check for errors
		str_out = new iostream( (void*) jpgfilename, ( !pipe_on ) ? 0 : 2, 0, 1 );
		if ( str_out->chkerr() ) {
			sprintf( errormessage, FWR_ERRMSG, jpgfilename );
			errorlevel = 2;
			return false;
		}
		// PJG specific settings - auto unless specified otherwise
		auto_set = true;
	}
	else {
		// file is neither
		filetype = F_UNK;
		sprintf( errormessage, "filetype of file \"%s\" is unknown", filename );
		errorlevel = 2;
		return false;		
	}
	
	
	return true;
}
#endif


/* -----------------------------------------------
	set each variable to its initial value
	----------------------------------------------- */

INTERN bool reset_buffers( void )
{
	int cmp, bpos;
	int i;
	
	
	// -- free buffers --
	
	// free buffers & set pointers NULL
	if ( hdrdata  != NULL ) free ( hdrdata );
	if ( huffdata != NULL ) free ( huffdata );
	if ( grbgdata != NULL ) free ( grbgdata );
	if ( rst_err  != NULL ) free ( rst_err );
	if ( rstp     != NULL ) free ( rstp );
	if ( scnp     != NULL ) free ( scnp );
	hdrdata   = NULL;
	huffdata  = NULL;
	grbgdata  = NULL;
	rst_err   = NULL;
	rstp      = NULL;
	scnp      = NULL;
	
	// free image arrays
	for ( cmp = 0; cmp < 4; cmp++ )	{
		if ( zdstdata[ cmp ] != NULL ) free( zdstdata[cmp] );
		if ( eobxhigh[ cmp ] != NULL ) free( eobxhigh[cmp] );
		if ( eobyhigh[ cmp ] != NULL ) free( eobyhigh[cmp] );
		if ( zdstxlow[ cmp ] != NULL ) free( zdstxlow[cmp] );
		if ( zdstylow[ cmp ] != NULL ) free( zdstylow[cmp] );
		zdstdata[ cmp ] = NULL;
		eobxhigh[ cmp ] = NULL;
		eobyhigh[ cmp ] = NULL;
		zdstxlow[ cmp ] = NULL;
		zdstylow[ cmp ] = NULL;
		freqscan[ cmp ] = (unsigned char*) stdscan;
		
		for ( bpos = 0; bpos < 64; bpos++ ) {
			if ( colldata[ cmp ][ bpos ] != NULL ) free( colldata[cmp][bpos] );
			colldata[ cmp ][ bpos ] = NULL;
		}		
	}
	
	
	// -- set variables --
	
	// preset componentinfo
	for ( cmp = 0; cmp < 4; cmp++ ) {
		cmpnfo[ cmp ].sfv = -1;
		cmpnfo[ cmp ].sfh = -1;
		cmpnfo[ cmp ].mbs = -1;
		cmpnfo[ cmp ].bcv = -1;
		cmpnfo[ cmp ].bch = -1;
		cmpnfo[ cmp ].bc  = -1;
		cmpnfo[ cmp ].ncv = -1;
		cmpnfo[ cmp ].nch = -1;
		cmpnfo[ cmp ].nc  = -1;
		cmpnfo[ cmp ].sid = -1;
		cmpnfo[ cmp ].jid = -1;
		cmpnfo[ cmp ].qtable = NULL;
		cmpnfo[ cmp ].huffdc = -1;
		cmpnfo[ cmp ].huffac = -1;
	}
	
	// preset imgwidth / imgheight / component count 
	imgwidth  = 0;
	imgheight = 0;
	cmpc      = 0;
	
	// preset mcu info variables / restart interval
	sfhm      = 0;
	sfvm      = 0;
	mcuc      = 0;
	mcuh      = 0;
	mcuv      = 0;
	rsti      = 0;
	
	// reset quantization / huffman tables
	for ( i = 0; i < 4; i++ ) {
		htset[ 0 ][ i ] = 0;
		htset[ 1 ][ i ] = 0;
		for ( bpos = 0; bpos < 64; bpos++ )
			qtables[ i ][ bpos ] = 0;
	}
	
	// preset jpegtype
	jpegtype  = 0;
	
	// reset padbit
	padbit = -1;
	
	
	return true;
}


/* -----------------------------------------------
	Merges header & image data to jpeg
	----------------------------------------------- */
	
INTERN bool merge_jpeg( void )
{
	unsigned char SOI[ 2 ] = { 0xFF, 0xD8 }; // SOI segment
	unsigned char EOI[ 2 ] = { 0xFF, 0xD9 }; // EOI segment
	unsigned char mrk = 0xFF; // marker start
	unsigned char stv = 0x00; // 0xFF stuff value
	unsigned char rst = 0xD0; // restart marker
	
	unsigned char  type = 0x00; // type of current marker segment
	unsigned int   len  = 0; // length of current marker segment
	unsigned int   hpos = 0; // current position in header
	unsigned int   ipos = 0; // current position in imagedata
	unsigned int   rpos = 0; // current restart marker position
	unsigned int   cpos = 0; // in scan corrected rst marker position
	unsigned int   scan = 1; // number of current scan
	unsigned int   tmp; // temporary storage variable
	
	
	// write SOI
	str_out->write( SOI, 1, 2 );
	
	// JPEG writing loop
	while ( true )
	{		
		// store current header position
		tmp = hpos;
		
		// seek till start-of-scan
		for ( type = 0x00; type != 0xDA; ) {
			if ( ( int ) hpos >= hdrs ) break;
			type = hdrdata[ hpos + 1 ];
			len = 2 + B_SHORT( hdrdata[ hpos + 2 ], hdrdata[ hpos + 3 ] );
			hpos += len;
		}
		
		// write header data to file
		str_out->write( hdrdata + tmp, 1, ( hpos - tmp ) );
		
		// get out if last marker segment type was not SOS
		if ( type != 0xDA ) break;
		
		
		// (re)set corrected rst pos
		cpos = 0;
		
		// write & expand huffman coded image data
		for ( ipos = scnp[ scan - 1 ]; ipos < scnp[ scan ]; ipos++ ) {
			// write current byte
			str_out->write( huffdata + ipos, 1, 1 );
			// check current byte, stuff if needed
			if ( huffdata[ ipos ] == 0xFF )
				str_out->write( &stv, 1, 1 );
			// insert restart markers if needed
			if ( rstp != NULL ) {
				if ( ipos == rstp[ rpos ] ) {
					rst = 0xD0 + ( cpos % 8 );
					str_out->write( &mrk, 1, 1 );
					str_out->write( &rst, 1, 1 );
					rpos++; cpos++;
				}
			}
		}
		// insert false rst markers at end if needed
		if ( rst_err != NULL ) {
			while ( rst_err[ scan - 1 ] > 0 ) {
				rst = 0xD0 + ( cpos % 8 );
				str_out->write( &mrk, 1, 1 );
				str_out->write( &rst, 1, 1 );
				cpos++;	rst_err[ scan - 1 ]--;
			}
		}

		// proceed with next scan
		scan++;
	}
	
	// write EOI
	str_out->write( EOI, 1, 2 );
	
	// write garbage if needed
	if ( grbs > 0 )
		str_out->write( grbgdata, 1, grbs );
	
	// errormessage if write error
	if ( str_out->chkerr() ) {
		sprintf( errormessage, "write error, possibly drive is full" );
		errorlevel = 2;		
		return false;
	}
	
	// get filesize
	jpgfilesize = str_out->getsize();
	
	
	return true;
}


/* -----------------------------------------------
	JPEG encoding routine
	----------------------------------------------- */

INTERN bool recode_jpeg( void )
{
	abitwriter*  huffw; // bitwise writer for image data
	abytewriter* storw; // bytewise writer for storage of correction bits
	
	unsigned char  type = 0x00; // type of current marker segment
	unsigned int   len  = 0; // length of current marker segment
	unsigned int   hpos = 0; // current position in header
		
	int lastdc[ 4 ]; // last dc for each component0
	short block[ 64 ]; // store block for coeffs
	int eobrun; // run of eobs
	int rstw; // restart wait counter
	
	int cmp, bpos, dpos;
	int mcu, sub, csc;
	int eob, sta;
	int tmp;
	
	
	// open huffman coded image data in abitwriter
	huffw = new abitwriter( 0 );
	huffw->fillbit = padbit;
	
	// init storage writer
	storw = new abytewriter( 0 );
	
	// preset count of scans and restarts
	scnc = 0;
	rstc = 0;
	
	// JPEG decompression loop
	while ( true )
	{
		// seek till start-of-scan, parse only DHT, DRI and SOS
		for ( type = 0x00; type != 0xDA; ) {
			if ( ( int ) hpos >= hdrs ) break;
			type = hdrdata[ hpos + 1 ];
			len = 2 + B_SHORT( hdrdata[ hpos + 2 ], hdrdata[ hpos + 3 ] );
			if ( ( type == 0xC4 ) || ( type == 0xDA ) || ( type == 0xDD ) ) {
				if ( !jpg_parse_jfif( type, len, &( hdrdata[ hpos ] ) ) ) {
					return false;
				}
				hpos += len;
			}
			else {
				hpos += len;
				continue;
			}			
		}
		
		// get out if last marker segment type was not SOS
		if ( type != 0xDA ) break;
		
		
		// (re)alloc scan positons array
		if ( scnp == NULL ) scnp = ( unsigned int* ) calloc( scnc + 2, sizeof( int ) );
		else scnp = ( unsigned int* ) realloc( scnp, ( scnc + 2 ) * sizeof( int ) );
		if ( scnp == NULL ) {
			sprintf( errormessage, MEM_ERRMSG );
			errorlevel = 2;
			return false;
		}
		
		// (re)alloc restart marker positons array if needed
		if ( rsti > 0 ) {
			tmp = rstc + ( ( cs_cmpc > 1 ) ?
				( mcuc / rsti ) : ( cmpnfo[ cs_cmp[ 0 ] ].bc / rsti ) );
			if ( rstp == NULL ) rstp = ( unsigned int* ) calloc( tmp + 1, sizeof( int ) );
			else rstp = ( unsigned int* ) realloc( rstp, ( tmp + 1 ) * sizeof( int ) );
			if ( rstp == NULL ) {
				sprintf( errormessage, MEM_ERRMSG );
				errorlevel = 2;
				return false;
			}
		}		
		
		// intial variables set for encoding
		cmp  = cs_cmp[ 0 ];
		csc  = 0;
		mcu  = 0;
		sub  = 0;
		dpos = 0;
		
		// store scan position
		scnp[ scnc ] = huffw->getpos();
		
		// JPEG imagedata encoding routines
		while ( true )
		{
			// (re)set last DCs for diff coding
			lastdc[ 0 ] = 0;
			lastdc[ 1 ] = 0;
			lastdc[ 2 ] = 0;
			lastdc[ 3 ] = 0;
			
			// (re)set status
			sta = 0;
			
			// (re)set eobrun
			eobrun = 0;
			
			// (re)set rst wait counter
			rstw = rsti;
			
			// encoding for interleaved data
			if ( cs_cmpc > 1 )
			{				
				if ( jpegtype == 1 ) {
					// ---> sequential interleaved encoding <---
					while ( sta == 0 ) {
						// copy from colldata
						for ( bpos = 0; bpos < 64; bpos++ )
							block[ bpos ] = colldata[ cmp ][ bpos ][ dpos ];
						
						// diff coding for dc
						block[ 0 ] -= lastdc[ cmp ];
						lastdc[ cmp ] = colldata[ cmp ][ 0 ][ dpos ];
						
						// encode block
						eob = jpg_encode_block_seq( huffw,
							&(hcodes[ 0 ][ cmpnfo[cmp].huffac ]),
							&(hcodes[ 1 ][ cmpnfo[cmp].huffac ]),
							block );
						
						// check for errors, proceed if no error encountered
						if ( eob < 0 ) sta = -1;
						else sta = jpg_next_mcupos( &mcu, &cmp, &csc, &sub, &dpos, &rstw );
					}
				}
				else if ( cs_sah == 0 ) {
					// ---> progressive interleaved DC encoding <---
					// ---> succesive approximation first stage <---
					while ( sta == 0 ) {
						// diff coding & bitshifting for dc 
						tmp = colldata[ cmp ][ 0 ][ dpos ] >> cs_sal;
						block[ 0 ] = tmp - lastdc[ cmp ];
						lastdc[ cmp ] = tmp;
						
						// encode dc
						sta = jpg_encode_dc_prg_fs( huffw,
							&(hcodes[ 0 ][ cmpnfo[cmp].huffdc ]),
							block );
						
						// next mcupos if no error happened
						if ( sta != -1 )
							sta = jpg_next_mcupos( &mcu, &cmp, &csc, &sub, &dpos, &rstw );
					}
				}
				else {
					// ---> progressive interleaved DC encoding <---
					// ---> succesive approximation later stage <---
					while ( sta == 0 ) {
						// fetch bit from current bitplane
						block[ 0 ] = BITN( colldata[ cmp ][ 0 ][ dpos ], cs_sal );
						
						// encode dc correction bit
						sta = jpg_encode_dc_prg_sa( huffw, block );
						
						// next mcupos if no error happened
						if ( sta != -1 )
							sta = jpg_next_mcupos( &mcu, &cmp, &csc, &sub, &dpos, &rstw );
					}
				}
			}
			else // encoding for non interleaved data
			{
				if ( jpegtype == 1 ) {
					// ---> sequential non interleaved encoding <---
					while ( sta == 0 ) {
						// copy from colldata
						for ( bpos = 0; bpos < 64; bpos++ )
							block[ bpos ] = colldata[ cmp ][ bpos ][ dpos ];
						
						// diff coding for dc
						block[ 0 ] -= lastdc[ cmp ];
						lastdc[ cmp ] = colldata[ cmp ][ 0 ][ dpos ];
						
						// encode block
						eob = jpg_encode_block_seq( huffw,
							&(hcodes[ 0 ][ cmpnfo[cmp].huffac ]),
							&(hcodes[ 1 ][ cmpnfo[cmp].huffac ]),
							block );
						
						// check for errors, proceed if no error encountered
						if ( eob < 0 ) sta = -1;
						else sta = jpg_next_mcuposn( &cmp, &dpos, &rstw );	
					}
				}
				else if ( cs_to == 0 ) {
					if ( cs_sah == 0 ) {
						// ---> progressive non interleaved DC encoding <---
						// ---> succesive approximation first stage <---
						while ( sta == 0 ) {
							// diff coding & bitshifting for dc 
							tmp = colldata[ cmp ][ 0 ][ dpos ] >> cs_sal;
							block[ 0 ] = tmp - lastdc[ cmp ];
							lastdc[ cmp ] = tmp;
							
							// encode dc
							sta = jpg_encode_dc_prg_fs( huffw,
								&(hcodes[ 0 ][ cmpnfo[cmp].huffdc ]),
								block );							
							
							// check for errors, increment dpos otherwise
							if ( sta != -1 )
								sta = jpg_next_mcuposn( &cmp, &dpos, &rstw );
						}
					}
					else {
						// ---> progressive non interleaved DC encoding <---
						// ---> succesive approximation later stage <---
						while ( sta == 0 ) {
							// fetch bit from current bitplane
							block[ 0 ] = BITN( colldata[ cmp ][ 0 ][ dpos ], cs_sal );
							
							// encode dc correction bit
							sta = jpg_encode_dc_prg_sa( huffw, block );
							
							// next mcupos if no error happened
							if ( sta != -1 )
								sta = jpg_next_mcuposn( &cmp, &dpos, &rstw );
						}
					}
				}
				else {
					if ( cs_sah == 0 ) {
						// ---> progressive non interleaved AC encoding <---
						// ---> succesive approximation first stage <---
						while ( sta == 0 ) {
							// copy from colldata
							for ( bpos = cs_from; bpos <= cs_to; bpos++ )
								block[ bpos ] =
									FDIV2( colldata[ cmp ][ bpos ][ dpos ], cs_sal );
							
							// encode block
							eob = jpg_encode_ac_prg_fs( huffw,
								&(hcodes[ 1 ][ cmpnfo[cmp].huffac ]),
								block, &eobrun, cs_from, cs_to );
							
							// check for errors, proceed if no error encountered
							if ( eob < 0 ) sta = -1;
							else sta = jpg_next_mcuposn( &cmp, &dpos, &rstw );
						}						
						
						// encode remaining eobrun
						jpg_encode_eobrun( huffw,
							&(hcodes[ 1 ][ cmpnfo[cmp].huffac ]),
							&eobrun );
					}
					else {
						// ---> progressive non interleaved AC encoding <---
						// ---> succesive approximation later stage <---
						while ( sta == 0 ) {
							// copy from colldata
							for ( bpos = cs_from; bpos <= cs_to; bpos++ )
								block[ bpos ] =
									FDIV2( colldata[ cmp ][ bpos ][ dpos ], cs_sal );
							
							// encode block
							eob = jpg_encode_ac_prg_sa( huffw, storw,
								&(hcodes[ 1 ][ cmpnfo[cmp].huffac ]),
								block, &eobrun, cs_from, cs_to );
							
							// check for errors, proceed if no error encountered
							if ( eob < 0 ) sta = -1;
							else sta = jpg_next_mcuposn( &cmp, &dpos, &rstw );
						}						
						
						// encode remaining eobrun
						jpg_encode_eobrun( huffw,
							&(hcodes[ 1 ][ cmpnfo[cmp].huffac ]),
							&eobrun );
							
						// encode remaining correction bits
						jpg_encode_crbits( huffw, storw );
					}
				}
			}
			
			// pad huffman writer
			huffw->pad( padbit );
			
			// evaluate status
			if ( sta == -1 ) { // status -1 means error
				sprintf( errormessage, "encode error in scan%i / mcu%i",
					scnc, ( cs_cmpc > 1 ) ? mcu : dpos );
				delete huffw;
				errorlevel = 2;
				return false;
			}
			else if ( sta == 2 ) { // status 2 means done
				scnc++; // increment scan counter
				break; // leave decoding loop, everything is done here
			}
			else if ( sta == 1 ) { // status 1 means restart
				if ( rsti > 0 ) // store rstp & stay in the loop
					rstp[ rstc++ ] = huffw->getpos() - 1;
			}
		}
	}
	
	// safety check for error in huffwriter
	if ( huffw->error ) {
		delete huffw;
		sprintf( errormessage, MEM_ERRMSG );
		errorlevel = 2;
		return false;
	}
	
	// get data into huffdata
	huffdata = huffw->getptr();
	hufs = huffw->getpos();	
	delete huffw;
	
	// remove storage writer
	delete storw;
	
	// store last scan & restart positions
	scnp[ scnc ] = hufs;
	if ( rstp != NULL )
		rstp[ rstc ] = hufs;
	
	
	return true;
}


/* -----------------------------------------------
	adapt ICOS tables for quantizer tables
	----------------------------------------------- */
	
INTERN bool adapt_icos( void )
{
	unsigned short quant[ 64 ]; // local copy of quantization
	int ipos;
	int cmp;
	
	
	for ( cmp = 0; cmp < cmpc; cmp++ ) {
		// make a local copy of the quantization values, check
		for ( ipos = 0; ipos < 64; ipos++ ) {
			quant[ ipos ] = QUANT( cmp, zigzag[ ipos ] );
			if ( quant[ ipos ] >= 2048 ) // if this is true, it can be safely assumed (for 8 bit JPEG), that all coefficients are zero
				quant[ ipos ] = 0;
		}
		// adapt idct 8x8 table
		for ( ipos = 0; ipos < 64 * 64; ipos++ )
			adpt_idct_8x8[ cmp ][ ipos ] = icos_idct_8x8[ ipos ] * quant[ ipos % 64 ];
		// adapt idct 1x8 table
		for ( ipos = 0; ipos < 8 * 8; ipos++ )
			adpt_idct_1x8[ cmp ][ ipos ] = icos_idct_1x8[ ipos ] * quant[ ( ipos % 8 ) * 8 ];
		// adapt idct 8x1 table
		for ( ipos = 0; ipos < 8 * 8; ipos++ )
			adpt_idct_8x1[ cmp ][ ipos ] = icos_idct_1x8[ ipos ] * quant[ ipos % 8 ];
	}
	
	
	return true;
}


/* -----------------------------------------------
	unpredict DC coefficients
	----------------------------------------------- */

INTERN bool unpredict_dc( void )
{	
	signed short* coef;
	int absmaxp;
	int absmaxn;
	int corr_f;
	int cmp, dpos;
	
	
	// remove prediction, store DC instead of prediction error
	for ( cmp = 0; cmp < cmpc; cmp++ ) {
		absmaxp = MAX_V( cmp, 0 );
		absmaxn = -absmaxp;
		corr_f = ( ( 2 * absmaxp ) + 1 );
		
		for ( dpos = 1; dpos < cmpnfo[cmp].bc; dpos++ ) {
			coef = &(colldata[cmp][0][dpos]);
			(*coef) += dc_1ddct_predictor( cmp, dpos ); // 1d dct predictor
			
			// fix range
			if ( (*coef) > absmaxp ) (*coef) -= corr_f;
			else if ( (*coef) < absmaxn ) (*coef) += corr_f;
		}
	}
	
	
	return true;
}


/* -----------------------------------------------
	calculate zero distribution lists
	----------------------------------------------- */
	
INTERN bool calc_zdst_lists( void )
{
	int cmp, bpos, dpos;
	int b_x, b_y;
	
	
	// this functions counts, for each DCT block, the number of non-zero coefficients
	for ( cmp = 0; cmp < cmpc; cmp++ )
	{
		// preset zdstlist
		memset( zdstdata[cmp], 0, cmpnfo[cmp].bc * sizeof( char ) );
		
		// calculate # on non-zeroes per block (separately for lower 7x7 block & first row/collumn)
		for ( bpos = 1; bpos < 64; bpos++ ) {
			b_x = unzigzag[ bpos ] % 8;
			b_y = unzigzag[ bpos ] / 8;
			if ( b_x == 0 ) {
				for ( dpos = 0; dpos < cmpnfo[cmp].bc; dpos++ )
					if ( colldata[cmp][bpos][dpos] != 0 ) zdstylow[cmp][dpos]++;
			}
			else if ( b_y == 0 ) {
				for ( dpos = 0; dpos < cmpnfo[cmp].bc; dpos++ )
					if ( colldata[cmp][bpos][dpos] != 0 ) zdstxlow[cmp][dpos]++;
			}
			else {
				for ( dpos = 0; dpos < cmpnfo[cmp].bc; dpos++ )
					if ( colldata[cmp][bpos][dpos] != 0 ) zdstdata[cmp][dpos]++;
			}
		}
	}
	
	
	return true;
}


/* -----------------------------------------------
	unpacks compressed pjg to colldata
	----------------------------------------------- */
	
INTERN bool unpack_pjg( void )
{
	aricoder* decoder;
	unsigned char hcode;
	unsigned char cb;
	int cmp;
	
	
	// check header codes ( maybe position in other function ? )
	while( true ) {
		str_in->read( &hcode, 1, 1 );
		if ( hcode == 0x00 ) {
			// retrieve compression settings from file
			str_in->read( nois_trs, 1, 4 );
			str_in->read( segm_cnt, 1, 4 );
			auto_set = false;
		}
		else if ( hcode >= 0x14 ) {
			// compare version number
			if ( hcode != appversion ) {
				sprintf( errormessage, "incompatible file, use %s v%i.%i",
					appname, hcode / 10, hcode % 10 );
				errorlevel = 2;
				return false;
			}
			else break;
		}
		else {
			sprintf( errormessage, "unknown header code, use newer version of %s", appname );
			errorlevel = 2;
			return false;
		}
	}
	
	
	// init arithmetic compression
	decoder = new aricoder( str_in, 0 );
	
	// decode JPG header
	if ( !pjg_decode_generic( decoder, &hdrdata, &hdrs ) ) return false;
	// retrieve padbit from stream
	if ( !pjg_decode_bit( decoder, &cb ) ) return false; padbit = cb;
	// decode one bit that signals false /correct use of RST markers
	if ( !pjg_decode_bit( decoder, &cb ) ) return false;
	// decode # of false set RST markers per scan only if available
	if ( cb == 1 )
		if ( !pjg_decode_generic( decoder, &rst_err, NULL ) ) return false;
	
	// undo header optimizations
	if ( !pjg_unoptimize_header() )	return false;	
	// discard meta information from header if option set
	if ( disc_meta )
		if ( !jpg_rebuild_header() ) return false;
	// parse header for image-info
	if ( !jpg_setup_imginfo() ) return false;
	
	// decode actual components data
	for ( cmp = 0; cmp < cmpc; cmp++ ) {		
		// decode frequency scan ('zero-sort-scan')
		if ( !pjg_decode_zstscan( decoder, cmp ) ) return false;		
		// decode zero-distribution-lists for higher (7x7) ACs
		if ( !pjg_decode_zdst_high( decoder, cmp ) ) return false;
		// decode coefficients for higher (7x7) ACs
		if ( !pjg_decode_ac_high( decoder, cmp ) ) return false;
		// decode zero-distribution-lists for lower ACs
		if ( !pjg_decode_zdst_low( decoder, cmp ) ) return false;
		// decode coefficients for first row / collumn ACs
		if ( !pjg_decode_ac_low( decoder, cmp ) ) return false;	
		// decode coefficients for DC
		if ( !pjg_decode_dc( decoder, cmp ) ) return false;	
	}
	
	// retrieve checkbit for garbage (0 if no garbage, 1 if garbage has to be coded)
	if ( !pjg_decode_bit( decoder, &cb ) ) return false;
	
	// decode garbage data only if available
	if ( cb == 0 ) grbs = 0;
	else if ( !pjg_decode_generic( decoder, &grbgdata, &grbs ) ) return false;
	
	// finalize arithmetic compression
	delete( decoder );
	
	
	// get filesize
	pjgfilesize = str_in->getsize();
	
	
	return true;
}

/* ----------------------- End of main functions -------------------------- */

/* ----------------------- Begin of JPEG specific functions -------------------------- */


/* -----------------------------------------------
	Parses header for imageinfo
	----------------------------------------------- */
INTERN bool jpg_setup_imginfo( void )
{
	unsigned char  type = 0x00; // type of current marker segment
	unsigned int   len  = 0; // length of current marker segment
	unsigned int   hpos = 0; // position in header
	
	int cmp, bpos;
	int i;
	
	// header parser loop
	while ( ( int ) hpos < hdrs ) {
		type = hdrdata[ hpos + 1 ];
		len = 2 + B_SHORT( hdrdata[ hpos + 2 ], hdrdata[ hpos + 3 ] );
		// do not parse DHT & DRI
		if ( ( type != 0xDA ) && ( type != 0xC4 ) && ( type != 0xDD ) ) {
			if ( !jpg_parse_jfif( type, len, &( hdrdata[ hpos ] ) ) )
				return false;
		}
		hpos += len;
	}
	
	// check if information is complete
	if ( cmpc == 0 ) {
		sprintf( errormessage, "header contains incomplete information" );
		errorlevel = 2;
		return false;
	}
	for ( cmp = 0; cmp < cmpc; cmp++ ) {
		if ( ( cmpnfo[cmp].sfv == 0 ) ||
			 ( cmpnfo[cmp].sfh == 0 ) ||
			 ( cmpnfo[cmp].qtable == NULL ) ||
			 ( cmpnfo[cmp].qtable[0] == 0 ) ||
			 ( jpegtype == 0 ) ) {
			sprintf( errormessage, "header information is incomplete" );
			errorlevel = 2;
			return false;
		}
	}
	
	// do all remaining component info calculations
	for ( cmp = 0; cmp < cmpc; cmp++ ) {
		if ( cmpnfo[ cmp ].sfh > sfhm ) sfhm = cmpnfo[ cmp ].sfh;
		if ( cmpnfo[ cmp ].sfv > sfvm ) sfvm = cmpnfo[ cmp ].sfv;
	}
	mcuv = ( int ) ceil( (float) imgheight / (float) ( 8 * sfhm ) );
	mcuh = ( int ) ceil( (float) imgwidth  / (float) ( 8 * sfvm ) );
	mcuc  = mcuv * mcuh;
	for ( cmp = 0; cmp < cmpc; cmp++ ) {
		cmpnfo[ cmp ].mbs = cmpnfo[ cmp ].sfv * cmpnfo[ cmp ].sfh;		
		cmpnfo[ cmp ].bcv = mcuv * cmpnfo[ cmp ].sfh;
		cmpnfo[ cmp ].bch = mcuh * cmpnfo[ cmp ].sfv;
		cmpnfo[ cmp ].bc  = cmpnfo[ cmp ].bcv * cmpnfo[ cmp ].bch;
		cmpnfo[ cmp ].ncv = ( int ) ceil( (float) imgheight * 
							( (float) cmpnfo[ cmp ].sfh / ( 8.0 * sfhm ) ) );
		cmpnfo[ cmp ].nch = ( int ) ceil( (float) imgwidth * 
							( (float) cmpnfo[ cmp ].sfv / ( 8.0 * sfvm ) ) );
		cmpnfo[ cmp ].nc  = cmpnfo[ cmp ].ncv * cmpnfo[ cmp ].nch;
	}
	
	// decide components' statistical ids
	if ( cmpc <= 3 ) {
		for ( cmp = 0; cmp < cmpc; cmp++ ) cmpnfo[ cmp ].sid = cmp;
	}
	else {
		for ( cmp = 0; cmp < cmpc; cmp++ ) cmpnfo[ cmp ].sid = 0;
	}
	
	// alloc memory for further operations
	for ( cmp = 0; cmp < cmpc; cmp++ )
	{
		// alloc memory for colls
		for ( bpos = 0; bpos < 64; bpos++ ) {
			colldata[cmp][bpos] = (short int*) calloc ( cmpnfo[cmp].bc, sizeof( short ) );
			if (colldata[cmp][bpos] == NULL) {
				sprintf( errormessage, MEM_ERRMSG );
				errorlevel = 2;
				return false;
			}
		}
		
		// alloc memory for zdstlist / eob x / eob y
		zdstdata[cmp] = (unsigned char*) calloc( cmpnfo[cmp].bc, sizeof( char ) );
		eobxhigh[cmp] = (unsigned char*) calloc( cmpnfo[cmp].bc, sizeof( char ) );
		eobyhigh[cmp] = (unsigned char*) calloc( cmpnfo[cmp].bc, sizeof( char ) );
		zdstxlow[cmp] = (unsigned char*) calloc( cmpnfo[cmp].bc, sizeof( char ) );
		zdstylow[cmp] = (unsigned char*) calloc( cmpnfo[cmp].bc, sizeof( char ) );
		if ( ( zdstdata[cmp] == NULL ) ||
			( eobxhigh[cmp] == NULL ) || ( eobyhigh[cmp] == NULL ) ||
			( zdstxlow[cmp] == NULL ) || ( zdstylow[cmp] == NULL ) ) {
			sprintf( errormessage, MEM_ERRMSG );
			errorlevel = 2;
			return false;
		}
	}
	
	// also decide automatic settings here
	if ( auto_set ) {
		for ( cmp = 0; cmp < cmpc; cmp++ ) {
			for ( i = 0;
				conf_sets[ i ][ cmpnfo[cmp].sid ] > (unsigned int) cmpnfo[ cmp ].bc;
				i++ );
			segm_cnt[ cmp ] = conf_segm[ i ][ cmpnfo[cmp].sid ];
			nois_trs[ cmp ] = conf_ntrs[ i ][ cmpnfo[cmp].sid ];
		}
	}
	
	
	return true;
}


/* -----------------------------------------------
	Parse routines for JFIF segments
	----------------------------------------------- */
INTERN bool jpg_parse_jfif( unsigned char type, unsigned int len, unsigned char* segment )
{
	unsigned int hpos = 4; // current position in segment, start after segment header
	int lval, rval; // temporary variables
	int skip;
	int cmp;
	int i;
	
	
	switch ( type )
	{
		case 0xC4: // DHT segment
			// build huffman trees & codes
			while ( hpos < len ) {
				lval = LBITS( segment[ hpos ], 4 );
				rval = RBITS( segment[ hpos ], 4 );
				if ( ((lval < 0) || (lval >= 2)) || ((rval < 0) || (rval >= 4)) )
					break;
					
				hpos++;
				// build huffman codes & trees
				jpg_build_huffcodes( &(segment[ hpos + 0 ]), &(segment[ hpos + 16 ]),
					&(hcodes[ lval ][ rval ]), &(htrees[ lval ][ rval ]) );
				htset[ lval ][ rval ] = 1;
				
				skip = 16;
				for ( i = 0; i < 16; i++ )		
					skip += ( int ) segment[ hpos + i ];				
				hpos += skip;
			}
			
			if ( hpos != len ) {
				// if we get here, something went wrong
				sprintf( errormessage, "size mismatch in dht marker" );
				errorlevel = 2;
				return false;
			}
			return true;
		
		case 0xDB: // DQT segment
			// copy quantization tables to internal memory
			while ( hpos < len ) {
				lval = LBITS( segment[ hpos ], 4 );
				rval = RBITS( segment[ hpos ], 4 );
				if ( (lval < 0) || (lval >= 2) ) break;
				if ( (rval < 0) || (rval >= 4) ) break;
				hpos++;				
				if ( lval == 0 ) { // 8 bit precision
					for ( i = 0; i < 64; i++ ) {
						qtables[ rval ][ i ] = ( unsigned short ) segment[ hpos + i ];
						if ( qtables[ rval ][ i ] == 0 ) break;
					}
					hpos += 64;
				}
				else { // 16 bit precision
					for ( i = 0; i < 64; i++ ) {
						qtables[ rval ][ i ] =
							B_SHORT( segment[ hpos + (2*i) ], segment[ hpos + (2*i) + 1 ] );
						if ( qtables[ rval ][ i ] == 0 ) break;
					}
					hpos += 128;
				}
			}
			
			if ( hpos != len ) {
				// if we get here, something went wrong
				sprintf( errormessage, "size mismatch in dqt marker" );
				errorlevel = 2;
				return false;
			}
			return true;
			
		case 0xDD: // DRI segment
			// define restart interval
			rsti = B_SHORT( segment[ hpos ], segment[ hpos + 1 ] );			
			return true;
			
		case 0xDA: // SOS segment
			// prepare next scan
			cs_cmpc = segment[ hpos ];
			if ( cs_cmpc > cmpc ) {
				sprintf( errormessage, "%i components in scan, only %i are allowed",
							cs_cmpc, cmpc );
				errorlevel = 2;
				return false;
			}
			hpos++;
			for ( i = 0; i < cs_cmpc; i++ ) {
				for ( cmp = 0; ( segment[ hpos ] != cmpnfo[ cmp ].jid ) && ( cmp < cmpc ); cmp++ );
				if ( cmp == cmpc ) {
					sprintf( errormessage, "component id mismatch in start-of-scan" );
					errorlevel = 2;
					return false;
				}
				cs_cmp[ i ] = cmp;
				cmpnfo[ cmp ].huffdc = LBITS( segment[ hpos + 1 ], 4 );
				cmpnfo[ cmp ].huffac = RBITS( segment[ hpos + 1 ], 4 );
				if ( ( cmpnfo[ cmp ].huffdc < 0 ) || ( cmpnfo[ cmp ].huffdc >= 4 ) ||
					 ( cmpnfo[ cmp ].huffac < 0 ) || ( cmpnfo[ cmp ].huffac >= 4 ) ) {
					sprintf( errormessage, "huffman table number mismatch" );
					errorlevel = 2;
					return false;
				}
				hpos += 2;
			}
			cs_from = segment[ hpos + 0 ];
			cs_to   = segment[ hpos + 1 ];
			cs_sah  = LBITS( segment[ hpos + 2 ], 4 );
			cs_sal  = RBITS( segment[ hpos + 2 ], 4 );
			// check for errors
			if ( ( cs_from > cs_to ) || ( cs_from > 63 ) || ( cs_to > 63 ) ) {
				sprintf( errormessage, "spectral selection parameter out of range" );
				errorlevel = 2;
				return false;
			}
			if ( ( cs_sah >= 12 ) || ( cs_sal >= 12 ) ) {
				sprintf( errormessage, "successive approximation parameter out of range" );
				errorlevel = 2;
				return false;
			}
			return true;
		
		case 0xC0: // SOF0 segment
			// coding process: baseline DCT
			
		case 0xC1: // SOF1 segment
			// coding process: extended sequential DCT
		
		case 0xC2: // SOF2 segment
			// coding process: progressive DCT
			
			// set JPEG coding type
			if ( type == 0xC2 )
				jpegtype = 2;
			else
				jpegtype = 1;
				
			// check data precision, only 8 bit is allowed
			lval = segment[ hpos ];
			if ( lval != 8 ) {
				sprintf( errormessage, "%i bit data precision is not supported", lval );
				errorlevel = 2;
				return false;
			}
			
			// image size, height & component count
			imgheight = B_SHORT( segment[ hpos + 1 ], segment[ hpos + 2 ] );
			imgwidth  = B_SHORT( segment[ hpos + 3 ], segment[ hpos + 4 ] );
			cmpc      = segment[ hpos + 5 ];
			if ( cmpc > 4 ) {
				sprintf( errormessage, "image has %i components, max 4 are supported", cmpc );
				errorlevel = 2;
				return false;
			}
			
			hpos += 6;
			// components contained in image
			for ( cmp = 0; cmp < cmpc; cmp++ ) {
				cmpnfo[ cmp ].jid = segment[ hpos ];
				cmpnfo[ cmp ].sfv = LBITS( segment[ hpos + 1 ], 4 );
				cmpnfo[ cmp ].sfh = RBITS( segment[ hpos + 1 ], 4 );				
				cmpnfo[ cmp ].qtable = qtables[ segment[ hpos + 2 ] ];
				hpos += 3;
			}
			
			return true;
		
		case 0xC3: // SOF3 segment
			// coding process: lossless sequential
			sprintf( errormessage, "sof3 marker found, image is coded lossless" );
			errorlevel = 2;
			return false;
		
		case 0xC5: // SOF5 segment
			// coding process: differential sequential DCT
			sprintf( errormessage, "sof5 marker found, image is coded diff. sequential" );
			errorlevel = 2;
			return false;
		
		case 0xC6: // SOF6 segment
			// coding process: differential progressive DCT
			sprintf( errormessage, "sof6 marker found, image is coded diff. progressive" );
			errorlevel = 2;
			return false;
		
		case 0xC7: // SOF7 segment
			// coding process: differential lossless
			sprintf( errormessage, "sof7 marker found, image is coded diff. lossless" );
			errorlevel = 2;
			return false;
			
		case 0xC9: // SOF9 segment
			// coding process: arithmetic extended sequential DCT
			sprintf( errormessage, "sof9 marker found, image is coded arithm. sequential" );
			errorlevel = 2;
			return false;
			
		case 0xCA: // SOF10 segment
			// coding process: arithmetic extended sequential DCT
			sprintf( errormessage, "sof10 marker found, image is coded arithm. progressive" );
			errorlevel = 2;
			return false;
			
		case 0xCB: // SOF11 segment
			// coding process: arithmetic extended sequential DCT
			sprintf( errormessage, "sof11 marker found, image is coded arithm. lossless" );
			errorlevel = 2;
			return false;
			
		case 0xCD: // SOF13 segment
			// coding process: arithmetic differntial sequential DCT
			sprintf( errormessage, "sof13 marker found, image is coded arithm. diff. sequential" );
			errorlevel = 2;
			return false;
			
		case 0xCE: // SOF14 segment
			// coding process: arithmetic differential progressive DCT
			sprintf( errormessage, "sof14 marker found, image is coded arithm. diff. progressive" );
			errorlevel = 2;
			return false;
		
		case 0xCF: // SOF15 segment
			// coding process: arithmetic differntial lossless
			sprintf( errormessage, "sof15 marker found, image is coded arithm. diff. lossless" );
			errorlevel = 2;
			return false;
			
		case 0xE0: // APP0 segment	
		case 0xE1: // APP1 segment
		case 0xE2: // APP2 segment
		case 0xE3: // APP3 segment
		case 0xE4: // APP4 segment
		case 0xE5: // APP5 segment
		case 0xE6: // APP6 segment
		case 0xE7: // APP7 segment
		case 0xE8: // APP8 segment
		case 0xE9: // APP9 segment
		case 0xEA: // APP10 segment
		case 0xEB: // APP11 segment
		case 0xEC: // APP12 segment
		case 0xED: // APP13 segment
		case 0xEE: // APP14 segment
		case 0xEF: // APP15 segment
		case 0xFE: // COM segment
			// do nothing - return true
			return true;
			
		case 0xD0: // RST0 segment
		case 0xD1: // RST1 segment
		case 0xD2: // RST2 segment
		case 0xD3: // RST3 segment
		case 0xD4: // RST4 segment
		case 0xD5: // RST5 segment
		case 0xD6: // RST6 segment
		case 0xD7: // RST7 segment
			// return errormessage - RST is out of place here
			sprintf( errormessage, "rst marker found out of place" );
			errorlevel = 2;
			return false;
		
		case 0xD8: // SOI segment
			// return errormessage - start-of-image is out of place here
			sprintf( errormessage, "soi marker found out of place" );
			errorlevel = 2;
			return false;
		
		case 0xD9: // EOI segment
			// return errormessage - end-of-image is out of place here
			sprintf( errormessage, "eoi marker found out of place" );
			errorlevel = 2;
			return false;
			
		default: // unknown marker segment
			// return warning
			sprintf( errormessage, "unknown marker found: FF %2X", type );
			errorlevel = 1;
			return true;
	}
}


/* -----------------------------------------------
	JFIF header rebuilding routine
	----------------------------------------------- */
INTERN bool jpg_rebuild_header( void )
{	
	abytewriter* hdrw; // new header writer
	
	unsigned char  type = 0x00; // type of current marker segment
	unsigned int   len  = 0; // length of current marker segment
	unsigned int   hpos = 0; // position in header	
	
	
	// start headerwriter
	hdrw = new abytewriter( 4096 );
	
	// header parser loop
	while ( ( int ) hpos < hdrs ) {
		type = hdrdata[ hpos + 1 ];
		len = 2 + B_SHORT( hdrdata[ hpos + 2 ], hdrdata[ hpos + 3 ] );
		// discard any unneeded meta info
		if ( ( type == 0xDA ) || ( type == 0xC4 ) || ( type == 0xDB ) ||
			 ( type == 0xC0 ) || ( type == 0xC1 ) || ( type == 0xC2 ) ||
			 ( type == 0xDD ) ) {
			hdrw->write_n( &(hdrdata[ hpos ]), len );
		}
		hpos += len;
	}
	
	// replace current header with the new one
	free( hdrdata );
	hdrdata = hdrw->getptr();
	hdrs    = hdrw->getpos();
	delete( hdrw );
	
	
	return true;
}


/* -----------------------------------------------
	sequential block encoding routine
	----------------------------------------------- */
INTERN int jpg_encode_block_seq( abitwriter* huffw, huffCodes* dctbl, huffCodes* actbl, short* block )
{
	unsigned short n;
	unsigned char  s;
	unsigned char  z;
	int bpos;
	int hc;
	
	
	// encode DC
	s = BITLEN2048N( block[ 0 ] );
	n = ENVLI( s, block[ 0 ] );
	huffw->write( dctbl->cval[ s ], dctbl->clen[ s ] );
	huffw->write( n, s );
	
	// encode AC
	z = 0;
	for ( bpos = 1; bpos < 64; bpos++ )
	{
		// if nonzero is encountered
		if ( block[ bpos ] != 0 ) {
			// write remaining zeroes
			while ( z >= 16 ) {
				huffw->write( actbl->cval[ 0xF0 ], actbl->clen[ 0xF0 ] );
				z -= 16;
			}			
			// vli encode
			s = BITLEN2048N( block[ bpos ] );
			n = ENVLI( s, block[ bpos ] );
			hc = ( ( z << 4 ) + s );
			// write to huffman writer
			huffw->write( actbl->cval[ hc ], actbl->clen[ hc ] );
			huffw->write( n, s );
			// reset zeroes
			z = 0;
		}
		else { // increment zero counter
			z++;
		}
	}
	// write eob if needed
	if ( z > 0 )
		huffw->write( actbl->cval[ 0x00 ], actbl->clen[ 0x00 ] );
		
	
	return 64 - z;
}


/* -----------------------------------------------
	progressive DC encoding routine
	----------------------------------------------- */
INTERN int jpg_encode_dc_prg_fs( abitwriter* huffw, huffCodes* dctbl, short* block )
{
	unsigned short n;
	unsigned char  s;
	
	
	// encode DC	
	s = BITLEN2048N( block[ 0 ] );
	n = ENVLI( s, block[ 0 ] );
	huffw->write( dctbl->cval[ s ], dctbl->clen[ s ] );
	huffw->write( n, s );
	
	
	// return 0 if everything is ok
	return 0;
}


/* -----------------------------------------------
	progressive AC encoding routine
	----------------------------------------------- */
INTERN int jpg_encode_ac_prg_fs( abitwriter* huffw, huffCodes* actbl, short* block, int* eobrun, int from, int to )
{
	unsigned short n;
	unsigned char  s;
	unsigned char  z;
	int bpos;
	int hc;
	
	// encode AC
	z = 0;
	for ( bpos = from; bpos <= to; bpos++ )
	{
		// if nonzero is encountered
		if ( block[ bpos ] != 0 ) {
			// encode eobrun
			jpg_encode_eobrun( huffw, actbl, eobrun );
			// write remaining zeroes
			while ( z >= 16 ) {
				huffw->write( actbl->cval[ 0xF0 ], actbl->clen[ 0xF0 ] );
				z -= 16;
			}			
			// vli encode
			s = BITLEN2048N( block[ bpos ] );
			n = ENVLI( s, block[ bpos ] );
			hc = ( ( z << 4 ) + s );
			// write to huffman writer
			huffw->write( actbl->cval[ hc ], actbl->clen[ hc ] );
			huffw->write( n, s );
			// reset zeroes
			z = 0;
		}
		else { // increment zero counter
			z++;
		}
	}
	
	// check eob, increment eobrun if needed
	if ( z > 0 ) {
		(*eobrun)++;
		// check eobrun, encode if needed
		if ( (*eobrun) == actbl->max_eobrun )
			jpg_encode_eobrun( huffw, actbl, eobrun );
		return 1 + to - z;		
	}
	else {
		return 1 + to;
	}
}


/* -----------------------------------------------
	progressive DC SA encoding routine
	----------------------------------------------- */
INTERN int jpg_encode_dc_prg_sa( abitwriter* huffw, short* block )
{
	// enocode next bit of dc coefficient
	huffw->write( block[ 0 ], 1 );
	
	// return 0 if everything is ok
	return 0;
}


/* -----------------------------------------------
	progressive AC SA encoding routine
	----------------------------------------------- */
INTERN int jpg_encode_ac_prg_sa( abitwriter* huffw, abytewriter* storw, huffCodes* actbl, short* block, int* eobrun, int from, int to )
{
	unsigned short n;
	unsigned char  s;
	unsigned char  z;
	int eob = from;
	int bpos;
	int hc;
	
	// check if block contains any newly nonzero coefficients and find out position of eob
	for ( bpos = to; bpos >= from; bpos-- )	{
		if ( ( block[ bpos ] == 1 ) || ( block[ bpos ] == -1 ) ) {
			eob = bpos + 1;
			break;
		}
	}
	
	// encode eobrun if needed
	if ( ( eob > from ) && ( (*eobrun) > 0 ) ) {
		jpg_encode_eobrun( huffw, actbl, eobrun );
		jpg_encode_crbits( huffw, storw );
	}
	
	// encode AC
	z = 0;
	for ( bpos = from; bpos < eob; bpos++ )
	{
		// if zero is encountered
		if ( block[ bpos ] == 0 ) {
			z++; // increment zero counter
			if ( z == 16 ) { // write zeroes if needed
				huffw->write( actbl->cval[ 0xF0 ], actbl->clen[ 0xF0 ] );
				jpg_encode_crbits( huffw, storw );
				z = 0;
			}
		}
		// if nonzero is encountered
		else if ( ( block[ bpos ] == 1 ) || ( block[ bpos ] == -1 ) ) {
			// vli encode			
			s = BITLEN2048N( block[ bpos ] );
			n = ENVLI( s, block[ bpos ] );
			hc = ( ( z << 4 ) + s );
			// write to huffman writer
			huffw->write( actbl->cval[ hc ], actbl->clen[ hc ] );
			huffw->write( n, s );
			// write correction bits
			jpg_encode_crbits( huffw, storw );
			// reset zeroes
			z = 0;
		}
		else { // store correction bits
			n = block[ bpos ] & 0x1;
			storw->write( n );
		}
	}
	
	// fast processing after eob
	for ( ;bpos <= to; bpos++ )
	{
		if ( block[ bpos ] != 0 ) { // store correction bits
			n = block[ bpos ] & 0x1;
			storw->write( n );
		}
	}
	
	// check eob, increment eobrun if needed
	if ( eob <= to ) {
		(*eobrun)++;	
		// check eobrun, encode if needed
		if ( (*eobrun) == actbl->max_eobrun ) {
			jpg_encode_eobrun( huffw, actbl, eobrun );
			jpg_encode_crbits( huffw, storw );		
		}
	}	
	
	// return eob
	return eob;
}


/* -----------------------------------------------
	run of EOB encoding routine
	----------------------------------------------- */
INTERN int jpg_encode_eobrun( abitwriter* huffw, huffCodes* actbl, int* eobrun )
{
	unsigned short n;
	unsigned char  s;
	int hc;
	
	
	if ( (*eobrun) > 0 ) {
		while ( (*eobrun) > actbl->max_eobrun ) {
			huffw->write( actbl->cval[ 0xE0 ], actbl->clen[ 0xE0 ] );
			huffw->write( E_ENVLI( 14, 32767 ), 14 );
			(*eobrun) -= actbl->max_eobrun;
		}
		BITLEN( s, (*eobrun) );
		s--;
		n = E_ENVLI( s, (*eobrun) );
		hc = ( s << 4 );
		huffw->write( actbl->cval[ hc ], actbl->clen[ hc ] );
		huffw->write( n, s );
		(*eobrun) = 0;
	}

	
	return 0;
}


/* -----------------------------------------------
	correction bits encoding routine
	----------------------------------------------- */
INTERN int jpg_encode_crbits( abitwriter* huffw, abytewriter* storw )
{	
	unsigned char* data;
	int len;
	int i;
	
	
	// peek into data from abytewriter	
	len = storw->getpos();
	if ( len == 0 ) return 0;
	data = storw->peekptr();
	
	// write bits to huffwriter
	for ( i = 0; i < len; i++ )
		huffw->write( data[ i ], 1 );
	
	// reset abytewriter, discard data
	storw->reset();
	
	
	return 0;
}


/* -----------------------------------------------
	calculates next position for MCU
	----------------------------------------------- */
INTERN int jpg_next_mcupos( int* mcu, int* cmp, int* csc, int* sub, int* dpos, int* rstw )
{
	int sta = 0; // status
	
	
	// increment all counts where needed
	if ( ( ++(*sub) ) >= cmpnfo[(*cmp)].mbs ) {
		(*sub) = 0;
		
		if ( ( ++(*csc) ) >= cs_cmpc ) {
			(*csc) = 0;
			(*cmp) = cs_cmp[ 0 ];
			(*mcu)++;
			if ( (*mcu) >= mcuc ) sta = 2;
			else if ( rsti > 0 )
				if ( --(*rstw) == 0 ) sta = 1;
		}
		else {
			(*cmp) = cs_cmp[(*csc)];
		}
	}
	
	// get correct position in image ( x & y )
	if ( cmpnfo[(*cmp)].sfh > 1 ) { // to fix mcu order
		(*dpos)  = ( (*mcu) / mcuh ) * cmpnfo[(*cmp)].sfh + ( (*sub) / cmpnfo[(*cmp)].sfv );
		(*dpos) *= cmpnfo[(*cmp)].bch;
		(*dpos) += ( (*mcu) % mcuh ) * cmpnfo[(*cmp)].sfv + ( (*sub) % cmpnfo[(*cmp)].sfv );
	}
	else if ( cmpnfo[(*cmp)].sfv > 1 ) {
		// simple calculation to speed up things if simple fixing is enough
		(*dpos) = ( (*mcu) * cmpnfo[(*cmp)].mbs ) + (*sub);
	}
	else {
		// no calculations needed without subsampling
		(*dpos) = (*mcu);
	}
	
	
	return sta;
}


/* -----------------------------------------------
	calculates next position (non interleaved)
	----------------------------------------------- */
INTERN int jpg_next_mcuposn( int* cmp, int* dpos, int* rstw )
{
	// increment position
	(*dpos)++;
	
	// fix for non interleaved mcu - horizontal
	if ( cmpnfo[(*cmp)].bch != cmpnfo[(*cmp)].nch ) {
		if ( (*dpos) % cmpnfo[(*cmp)].bch == cmpnfo[(*cmp)].nch )
			(*dpos) += ( cmpnfo[(*cmp)].bch - cmpnfo[(*cmp)].nch );
	}
	
	// fix for non interleaved mcu - vertical
	if ( cmpnfo[(*cmp)].bcv != cmpnfo[(*cmp)].ncv ) {
		if ( (*dpos) / cmpnfo[(*cmp)].bch == cmpnfo[(*cmp)].ncv )
			(*dpos) = cmpnfo[(*cmp)].bc;
	}
	
	// check position
	if ( (*dpos) >= cmpnfo[(*cmp)].bc ) return 2;
	else if ( rsti > 0 )
		if ( --(*rstw) == 0 ) return 1;
	

	return 0;
}


/* -----------------------------------------------
	creates huffman-codes & -trees from dht-data
	----------------------------------------------- */
INTERN void jpg_build_huffcodes( unsigned char *clen, unsigned char *cval,	huffCodes *hc, huffTree *ht )
{
	int nextfree;	
	int code;
	int node;
	int i, j, k;
	
	
	// fill with zeroes
	memset( hc->clen, 0, 256 * sizeof( short ) );
	memset( hc->cval, 0, 256 * sizeof( short ) );
	memset( ht->l, 0, 256 * sizeof( short ) );
	memset( ht->r, 0, 256 * sizeof( short ) );
	
	// 1st part -> build huffman codes
	
	// creating huffman-codes	
	k = 0;
	code = 0;	
	
	// symbol-value of code is its position in the table
	for( i = 0; i < 16; i++ ) {
		for( j = 0; j < (int) clen[ i ]; j++ ) {
			hc->clen[ (int) cval[k] ] = 1 + i;
			hc->cval[ (int) cval[k] ] = code;
			
			k++;			
			code++;
		}		
		code = code << 1;
	}
	
	// find out eobrun max value
	hc->max_eobrun = 0;
	for ( i = 14; i >= 0; i-- ) {
		if ( hc->clen[ i << 4 ] > 0 ) {
			hc->max_eobrun = ( 2 << i ) - 1;
			break;
		}
	}
	
	// 2nd -> part use codes to build the coding tree
	
	// initial value for next free place
	nextfree = 1;

	// work through every code creating links between the nodes (represented through ints)
	for ( i = 0; i < 256; i++ )	{
		// (re)set current node
		node = 0;   		   		
		// go through each code & store path
		for ( j = hc->clen[ i ] - 1; j > 0; j-- ) {
			if ( BITN( hc->cval[ i ], j ) == 1 ) {
				if ( ht->r[ node ] == 0 )
					 ht->r[ node ] = nextfree++;
				node = ht->r[ node ];
			}
			else{
				if ( ht->l[ node ] == 0 )
					ht->l[ node ] = nextfree++;
				node = ht->l[ node ];
			}   					
		}
		// last link is number of targetvalue + 256
		if ( hc->clen[ i ] > 0 ) {
			if ( BITN( hc->cval[ i ], 0 ) == 1 )
				ht->r[ node ] = i + 256;
			else
				ht->l[ node ] = i + 256;
		}	   	
	}
}

/* ----------------------- End of JPEG specific functions -------------------------- */

/* ----------------------- Begin of PJG specific functions -------------------------- */


/* -----------------------------------------------
	encodes frequency scanorder to pjg
	----------------------------------------------- */
INTERN bool pjg_decode_zstscan( aricoder* dec, int cmp )
{	
	model_s* model;;
	
	unsigned char freqlist[ 64 ];
	int tpos; // true position
	int cpos; // coded position
	int i;
	
	
	// set first position in zero sort scan
	zsrtscan[ cmp ][ 0 ] = 0;
	
	// preset freqlist
	for ( i = 0; i < 64; i++ )
		freqlist[ i ] = stdscan[ i ];
		
	// init model
	model = INIT_MODEL_S( 64, 64, 1 );
	
	// encode scanorder
	for ( i = 1; i < 64; i++ )
	{			
		// reduce range of model
		model->exclude_symbols( 'a', 64 - i );
		
		// decode symbol
		cpos = decode_ari( dec, model );
		model->shift_context( cpos );
		
		if ( cpos == 0 ) {
			// remaining list is identical to scan
			// fill the scan & make a quick exit				
			for ( tpos = 0; i < 64; i++ ) {
				while ( freqlist[ ++tpos ] == 0 );
				zsrtscan[ cmp ][ i ] = freqlist[ tpos ];
			}
			break;
		}
		
		// decode position from list
		for ( tpos = 0; tpos < 64; tpos++ ) {
			if ( freqlist[ tpos ] != 0 ) cpos--;
			if ( cpos == 0 ) break;
		}
			
		// write decoded position to zero sort scan
		zsrtscan[ cmp ][ i ] = freqlist[ tpos ];
		// remove from list
		freqlist[ tpos ] = 0;
	}
	
	// delete model
	delete( model  );		
	
	// set zero sort scan as freqscan
	freqscan[ cmp ] = zsrtscan[ cmp ];
	
	
	return true;
}


/* -----------------------------------------------
	decodes # of non zeroes from pjg (high)
	----------------------------------------------- */
INTERN bool pjg_decode_zdst_high( aricoder* dec, int cmp )
{
	model_s* model;
	
	unsigned char* zdstls;
	int dpos;
	int a, b;
	int bc;
	int w;
	
	
	// init model, constants
	model = INIT_MODEL_S( 49 + 1, 25 + 1, 1 );
	zdstls = zdstdata[ cmp ];
	w = cmpnfo[cmp].bch;
	bc = cmpnfo[cmp].bc;
	
	// arithmetic decode zero-distribution-list
	for ( dpos = 0; dpos < bc; dpos++ )	{			
		// context modelling - use average of above and left as context		
		get_context_nnb( dpos, w, &a, &b );
		a = ( a >= 0 ) ? zdstls[ a ] : 0;
		b = ( b >= 0 ) ? zdstls[ b ] : 0;
		// shift context
		model->shift_context( ( a + b + 2 ) / 4 );
		// decode symbol
		zdstls[ dpos ] = decode_ari( dec, model );
	}
	
	// clean up
	delete( model );
	
	
	return true;
}


/* -----------------------------------------------
	decodes # of non zeroes from pjg (low)
	----------------------------------------------- */	
INTERN bool pjg_decode_zdst_low( aricoder* dec, int cmp )
{
	model_s* model;
	
	unsigned char* zdstls_x;
	unsigned char* zdstls_y;
	unsigned char* ctx_zdst;
	unsigned char* ctx_eobx;
	unsigned char* ctx_eoby;
	
	int dpos;
	int bc;
	
	
	// init model, constants
	model = INIT_MODEL_S( 8, 8, 2 );
	zdstls_x = zdstxlow[ cmp ];
	zdstls_y = zdstylow[ cmp ];
	ctx_eobx = eobxhigh[ cmp ];
	ctx_eoby = eobyhigh[ cmp ];
	ctx_zdst = zdstdata[ cmp ];
	bc = cmpnfo[cmp].bc;
	
	// arithmetic encode zero-distribution-list (first row)
	for ( dpos = 0; dpos < bc; dpos++ ) {
		model->shift_context( ( ctx_zdst[dpos] + 3 ) / 7 ); // shift context
		model->shift_context( ctx_eobx[dpos] ); // shift context
		zdstls_x[ dpos ] = decode_ari( dec, model ); // decode symbol
	}
	// arithmetic encode zero-distribution-list (first collumn)
	for ( dpos = 0; dpos < bc; dpos++ ) {
		model->shift_context( ( ctx_zdst[dpos] + 3 ) / 7 ); // shift context
		model->shift_context( ctx_eoby[dpos] ); // shift context
		zdstls_y[ dpos ] = decode_ari( dec, model ); // decode symbol
	}
	
	// clean up
	delete( model );
	
	
	return true;
}


/* -----------------------------------------------
	decodes DC coefficients from pjg
	----------------------------------------------- */
INTERN bool pjg_decode_dc( aricoder* dec, int cmp )
{
	unsigned char* segm_tab;
	
	model_s* mod_len;
	model_b* mod_sgn;
	model_b* mod_res;
	
	unsigned char* zdstls; // pointer to zero distribution list
	signed short* coeffs; // pointer to current coefficent data
	
	unsigned short* absv_store; // absolute coefficients values storage
	unsigned short* c_absc[ 6 ]; // quick access array for contexts
	int c_weight[ 6 ]; // weighting for contexts

	int ctx_avr; // 'average' context
	int ctx_len; // context for bit length
	
	int max_val; // max value
	int max_len; // max bitlength
	
	int dpos;
	int clen, absv, sgn;
	int snum;
	int bt, bp;
	
	int p_x, p_y;
	int r_x; //, r_y;
	int w, bc;
	
	
	// decide segmentation setting
	segm_tab = segm_tables[ segm_cnt[ cmp ] - 1 ];
	
	// get max absolute value/bit length
	max_val = MAX_V( cmp, 0 );
	max_len = BITLEN1024P( max_val );
	
	// init models for bitlenghts and -patterns
	mod_len = INIT_MODEL_S( max_len + 1, ( segm_cnt[cmp] > max_len ) ? segm_cnt[cmp] : max_len + 1, 2 );
	mod_res = INIT_MODEL_B( ( segm_cnt[cmp] < 16 ) ? 1 << 4 : segm_cnt[cmp], 2 );
	mod_sgn = INIT_MODEL_B( 1, 0 );
	
	// set width/height of each band
	bc = cmpnfo[cmp].bc;
	w = cmpnfo[cmp].bch;
	
	// allocate memory for absolute values storage
	absv_store = (unsigned short*) calloc ( bc, sizeof( short ) );
	if ( absv_store == NULL ) {
		sprintf( errormessage, MEM_ERRMSG );
		errorlevel = 2;
		return false;
	}
	
	// set up context quick access array
	pjg_aavrg_prepare( c_absc, c_weight, absv_store, cmp );
	
	// locally store pointer to coefficients and zero distribution list
	coeffs = colldata[ cmp ][ 0 ];
	zdstls = zdstdata[ cmp ];	
	
	// arithmetic compression loop
	for ( dpos = 0; dpos < bc; dpos++ )
	{		
		//calculate x/y positions in band
		p_y = dpos / w;
		// r_y = h - ( p_y + 1 );
		p_x = dpos % w;
		r_x = w - ( p_x + 1 );
		
		// get segment-number from zero distribution list and segmentation set
		snum = segm_tab[ zdstls[dpos] ];
		// calculate contexts (for bit length)
		ctx_avr = pjg_aavrg_context( c_absc, c_weight, dpos, p_y, p_x, r_x ); // AVERAGE context
		ctx_len = BITLEN1024P( ctx_avr ); // BITLENGTH context				
		// shift context / do context modelling (segmentation is done per context)
		shift_model( mod_len, ctx_len, snum );
		// decode bit length of current coefficient
		clen = decode_ari( dec, mod_len );
		
		// simple treatment if coefficient is zero
		if ( clen == 0 ) {
			// coeffs[ dpos ] = 0;
		}
		else {
			// decoding of residual
			absv = 1;
			// first set bit must be 1, so we start at clen - 2
			for ( bp = clen - 2; bp >= 0; bp-- ) {
				shift_model( mod_res, snum, bp ); // shift in 2 contexts
				// decode bit
				bt = decode_ari( dec, mod_res );
				// update absv
				absv = absv << 1;
				if ( bt ) absv |= 1; 
			}
			// decode sign
			sgn = decode_ari( dec, mod_sgn );
			// copy to colldata
			coeffs[ dpos ] = ( sgn == 0 ) ? absv : -absv;
			// store absolute value/sign
			absv_store[ dpos ] = absv;
		}
	}
	
	// free memory / clear models
	free( absv_store );
	delete ( mod_len );
	delete ( mod_res );
	delete ( mod_sgn );
	
	
	return true;
}


/* -----------------------------------------------
	decodes high (7x7) AC coefficients to pjg
	----------------------------------------------- */
INTERN bool pjg_decode_ac_high( aricoder* dec, int cmp )
{
	unsigned char* segm_tab;
	
	model_s* mod_len;
	model_b* mod_sgn;
	model_b* mod_res;
	
	unsigned char* zdstls; // pointer to zero distribution list
	unsigned char* eob_x; // pointer to x eobs
	unsigned char* eob_y; // pointer to y eobs
	signed short* coeffs; // pointer to current coefficent data
	
	unsigned short* absv_store; // absolute coefficients values storage
	unsigned short* c_absc[ 6 ]; // quick access array for contexts
	int c_weight[ 6 ]; // weighting for contexts
	
	unsigned char* sgn_store; // sign storage for context	
	unsigned char* sgn_nbh; // left signs neighbor
	unsigned char* sgn_nbv; // upper signs neighbor

	int ctx_avr; // 'average' context
	int ctx_len; // context for bit length
	int ctx_sgn; // context for sign
	
	int max_val; // max value
	int max_len; // max bitlength
	
	int bpos, dpos;
	int clen, absv, sgn;
	int snum;
	int bt, bp;
	int i;
	
	int b_x, b_y;
	int p_x, p_y;
	int r_x;
	int w, bc;
	
	
	// decide segmentation setting
	segm_tab = segm_tables[ segm_cnt[ cmp ] - 1 ];
	
	// init models for bitlenghts and -patterns
	mod_len = INIT_MODEL_S( 11, ( segm_cnt[cmp] > 11 ) ? segm_cnt[cmp] : 11, 2 );
	mod_res = INIT_MODEL_B( ( segm_cnt[cmp] < 16 ) ? 1 << 4 : segm_cnt[cmp], 2 );
	mod_sgn = INIT_MODEL_B( 9, 1 );
	
	// set width/height of each band
	bc = cmpnfo[cmp].bc;
	w = cmpnfo[cmp].bch;
	
	// allocate memory for absolute values & signs storage
	absv_store = (unsigned short*) calloc ( bc, sizeof( short ) );	
	sgn_store = (unsigned char*) calloc ( bc, sizeof( char ) );
	zdstls = (unsigned char*) calloc ( bc, sizeof( char ) );
	if ( ( absv_store == NULL ) || ( sgn_store == NULL ) ) {
		sprintf( errormessage, MEM_ERRMSG );
		errorlevel = 2;
		return false;
	}
	
	// set up quick access arrays for signs context
	sgn_nbh = sgn_store - 1;
	sgn_nbv = sgn_store - w;	
	
	// locally store pointer to eob x / eob y
	eob_x = eobxhigh[ cmp ];
	eob_y = eobyhigh[ cmp ];
	
	// preset x/y eobs
	memset( eob_x, 0x00, bc * sizeof( char ) );
	memset( eob_y, 0x00, bc * sizeof( char ) );
	
	// make a local copy of the zero distribution list
	for ( dpos = 0; dpos < bc; dpos++ )
		zdstls[ dpos ] = zdstdata[ cmp ][ dpos ];
	
	// work through lower 7x7 bands in order of freqscan
	for ( i = 1; i < 64; i++ )
	{		
		// work through blocks in order of frequency scan
		bpos = (int) freqscan[cmp][i];
		b_x = unzigzag[ bpos ] % 8;
		b_y = unzigzag[ bpos ] / 8;		
		
		if ( ( b_x == 0 ) || ( b_y == 0 ) )
				continue; // process remaining coefficients elsewhere
		
		// preset absolute values/sign storage
		memset( absv_store, 0x00, bc * sizeof( short ) );
		memset( sgn_store, 0x00, bc * sizeof( char ) );
		
		// set up average context quick access arrays
		pjg_aavrg_prepare( c_absc, c_weight, absv_store, cmp );
		
		// locally store pointer to coefficients
		coeffs = colldata[ cmp ][ bpos ];
		
		// get max bit length
		max_val = MAX_V( cmp, bpos );
		max_len = BITLEN1024P( max_val );
		
		// arithmetic compression loop
		for ( dpos = 0; dpos < bc; dpos++ )
		{
			// skip if beyound eob
			if ( zdstls[dpos] == 0 )
				continue;
			
			//calculate x/y positions in band
			p_y = dpos / w;
			// r_y = h - ( p_y + 1 );
			p_x = dpos % w;
			r_x = w - ( p_x + 1 );					
			
			// get segment-number from zero distribution list and segmentation set
			snum = segm_tab[ zdstls[dpos] ];
			// calculate contexts (for bit length)
			ctx_avr = pjg_aavrg_context( c_absc, c_weight, dpos, p_y, p_x, r_x ); // AVERAGE context
			ctx_len = BITLEN1024P( ctx_avr ); // BITLENGTH context				
			// shift context / do context modelling (segmentation is done per context)
			shift_model( mod_len, ctx_len, snum );
			mod_len->exclude_symbols( 'a', max_len );
			
			// decode bit length of current coefficient
			clen = decode_ari( dec, mod_len );			
			// simple treatment if coefficient is zero
			if ( clen == 0 ) {
				// coeffs[ dpos ] = 0;
			}
			else {
				// decoding of residual
				absv = 1;
				// first set bit must be 1, so we start at clen - 2
				for ( bp = clen - 2; bp >= 0; bp-- ) {
					shift_model( mod_res, snum, bp ); // shift in 2 contexts
					// decode bit
					bt = decode_ari( dec, mod_res );
					// update absv
					absv = absv << 1;
					if ( bt ) absv |= 1; 
				}
				// decode sign
				ctx_sgn = ( p_x > 0 ) ? sgn_nbh[ dpos ] : 0; // sign context
				if ( p_y > 0 ) ctx_sgn += 3 * sgn_nbv[ dpos ]; // IMPROVE! !!!!!!!!!!!
				mod_sgn->shift_context( ctx_sgn );
				sgn = decode_ari( dec, mod_sgn );
				// copy to colldata
				coeffs[ dpos ] = ( sgn == 0 ) ? absv : -absv;
				// store absolute value/sign, decrement zdst
				absv_store[ dpos ] = absv;
				sgn_store[ dpos ] = sgn + 1;
				zdstls[dpos]--;
				// recalculate x/y eob
				if ( b_x > eob_x[dpos] ) eob_x[dpos] = b_x;
				if ( b_y > eob_y[dpos] ) eob_y[dpos] = b_y;	
			}
		}
		// flush models
		mod_len->flush_model( 1 );
		mod_res->flush_model( 1 );
		mod_sgn->flush_model( 1 );
	}
	
	// free memory / clear models
	free( absv_store );
	free( sgn_store );
	free( zdstls );
	delete ( mod_len );
	delete ( mod_res );
	delete ( mod_sgn );
	
	
	return true;
}


/* -----------------------------------------------
	decodes high (7x7) AC coefficients to pjg
	----------------------------------------------- */
INTERN bool pjg_decode_ac_low( aricoder* dec, int cmp )
{
	model_s* mod_len;
	model_b* mod_sgn;
	model_b* mod_res;
	model_b* mod_top;
	
	unsigned char* zdstls; // pointer to row/col # of non-zeroes
	signed short* coeffs; // pointer to current coefficent data
	
	signed short* coeffs_x[ 8 ]; // prediction coeffs - current block
	signed short* coeffs_a[ 8 ]; // prediction coeffs - neighboring block
	int pred_cf[ 8 ]; // prediction multipliers

	int ctx_lak; // lakhani context
	int ctx_abs; // absolute context
	int ctx_len; // context for bit length
	int ctx_res; // bit plane context for residual
	int ctx_sgn; // context for sign
	
	int max_valp; // max value (+)
	int max_valn; // max value (-)
	int max_len; // max bitlength
	int thrs_bp; // residual threshold bitplane
	int* edge_c; // edge criteria
	
	int bpos, dpos;
	int clen, absv, sgn;
	int bt, bp;
	int i;
	
	int b_x, b_y;
	int p_x, p_y;
	int w, bc;
	
	
	// init models for bitlenghts and -patterns
	mod_len = INIT_MODEL_S( 11, ( segm_cnt[cmp] > 11 ) ? segm_cnt[cmp] : 11, 2 );
	mod_res = INIT_MODEL_B( 1 << 4, 2 );
	mod_top = INIT_MODEL_B( ( nois_trs[cmp] > 4 ) ? 1 << nois_trs[cmp] : 1 << 4, 3 );
	mod_sgn = INIT_MODEL_B( 11, 1 );
	
	// set width/height of each band
	bc = cmpnfo[cmp].bc;
	w = cmpnfo[cmp].bch;
	
	// work through each first row / first collumn band
	for ( i = 2; i < 16; i++ )
	{		
		// alternate between first row and first collumn
		b_x = ( i % 2 == 0 ) ? i / 2 : 0;
		b_y = ( i % 2 == 1 ) ? i / 2 : 0;
		bpos = (int) zigzag[ b_x + (8*b_y) ];
		
		// locally store pointer to band coefficients
		coeffs = colldata[ cmp ][ bpos ];
		// store pointers to prediction coefficients
		if ( b_x == 0 ) {
			for ( ; b_x < 8; b_x++ ) {
				coeffs_x[ b_x ] = colldata[ cmp ][ zigzag[b_x+(8*b_y)] ];
				coeffs_a[ b_x ] = colldata[ cmp ][ zigzag[b_x+(8*b_y)] ] - 1;
				pred_cf[ b_x ] = icos_base_8x8[ b_x * 8 ] * QUANT ( cmp, zigzag[b_x+(8*b_y)] );
			} b_x = 0;
			zdstls = zdstylow[ cmp ];
			edge_c = &p_x;
		}
		else { // if ( b_y == 0 )
			for ( ; b_y < 8; b_y++ ) {
				coeffs_x[ b_y ] = colldata[ cmp ][ zigzag[b_x+(8*b_y)] ];
				coeffs_a[ b_y ] = colldata[ cmp ][ zigzag[b_x+(8*b_y)] ] - w;
				pred_cf[ b_y ] = icos_base_8x8[ b_y * 8 ] * QUANT ( cmp, zigzag[b_x+(8*b_y)] );
			} b_y = 0;
			zdstls = zdstxlow[ cmp ];
			edge_c = &p_y;
		}
		
		// get max bit length / other info
		max_valp = MAX_V( cmp, bpos );
		max_valn = -max_valp;
		max_len = BITLEN1024P( max_valp );
		thrs_bp = ( max_len > nois_trs[cmp] ) ? max_len - nois_trs[cmp] : 0;
		
		// arithmetic compression loop
		for ( dpos = 0; dpos < bc; dpos++ )
		{
			// skip if beyound eob
			if ( zdstls[ dpos ] == 0 )
				continue;
			
			//calculate x/y positions in band
			p_y = dpos / w;
			p_x = dpos % w;
			
			// edge treatment / calculate LAKHANI context
			if ( (*edge_c) > 0 )
				ctx_lak = pjg_lakh_context( coeffs_x, coeffs_a, pred_cf, dpos );
			else ctx_lak = 0;
			ctx_lak = CLAMPED( max_valn, max_valp, ctx_lak );
			ctx_len = BITLEN2048N( ctx_lak ); // BITLENGTH context				
			// shift context / do context modelling (segmentation is done per context)
			shift_model( mod_len, ctx_len, zdstls[ dpos ] );
			mod_len->exclude_symbols( 'a', max_len );
			
			// decode bit length of current coefficient
			clen = decode_ari( dec, mod_len );
			// simple treatment if coefficients == 0
			if ( clen == 0 ) {
				// coeffs[ dpos ] = 0;
			}
			else {
				// decoding of residual
				bp = clen - 2; // first set bit must be 1, so we start at clen - 2
				ctx_res = ( bp >= thrs_bp ) ? 1 : 0;
				ctx_abs = ABS( ctx_lak );
				ctx_sgn = ( ctx_lak == 0 ) ? 0 : ( ctx_lak > 0 ) ? 1 : 2;
				for ( ; bp >= thrs_bp; bp-- ) {						
					shift_model( mod_top, ctx_abs >> thrs_bp, ctx_res, clen - thrs_bp ); // shift in 3 contexts
					// decode bit
					bt = decode_ari( dec, mod_top );
					// update context
					ctx_res = ctx_res << 1;
					if ( bt ) ctx_res |= 1; 
				}
				absv = ( ctx_res == 0 ) ? 1 : ctx_res; // !!!!
				for ( ; bp >= 0; bp-- ) {
					shift_model( mod_res, zdstls[ dpos ], bp ); // shift in 2 contexts
					// decode bit
					bt = decode_ari( dec, mod_res );
					// update absv
					absv = absv << 1;
					if ( bt ) absv |= 1; 
				}
				// decode sign
				shift_model( mod_sgn, zdstls[ dpos ], ctx_sgn );
				sgn = decode_ari( dec, mod_sgn );
				// copy to colldata
				coeffs[ dpos ] = ( sgn == 0 ) ? absv : -absv;
				// decrement # of non zeroes
				zdstls[ dpos ]--;
			}
		}
		// flush models
		mod_len->flush_model( 1 );
		mod_res->flush_model( 1 );
		mod_top->flush_model( 1 );
		mod_sgn->flush_model( 1 );
	}
	
	// free memory / clear models
	delete ( mod_len );
	delete ( mod_res );
	delete ( mod_top );
	delete ( mod_sgn );
	
	
	return true;
}


/* -----------------------------------------------
	deodes a stream of generic (8bit) data from pjg
	----------------------------------------------- */
INTERN bool pjg_decode_generic( aricoder* dec, unsigned char** data, int* len )
{
	abytewriter* bwrt;
	model_s* model;
	int c;
	
	
	// start byte writer
	bwrt = new abytewriter( 1024 );
	
	// decode header, ending with 256 symbol
	model = INIT_MODEL_S( 256 + 1, 256, 1 );
	while ( true ) {
		c = decode_ari( dec, model );
		if ( c == 256 ) break;
		bwrt->write( (unsigned char) c );
		model->shift_context( c );
	}
	delete( model );
	
	// check for out of memory
	if ( bwrt->error ) {
		delete bwrt;
		sprintf( errormessage, MEM_ERRMSG );
		errorlevel = 2;
		return false;
	}
	
	// get data/length and close byte writer
	(*data) = bwrt->getptr();
	if ( len != NULL ) (*len) = bwrt->getpos();
	delete bwrt;
	
	
	return true;
}


/* -----------------------------------------------
	decodes one bit from pjg
	----------------------------------------------- */
INTERN bool pjg_decode_bit( aricoder* dec, unsigned char* bit )
{
	model_b* model;
	
	
	model = INIT_MODEL_B( 1, -1 );
	(*bit) = decode_ari( dec, model );
	delete( model );
	
	
	return true;
}


/* -----------------------------------------------
	undoes the header optimizations
	----------------------------------------------- */
INTERN bool pjg_unoptimize_header( void )
{
	unsigned char  type = 0x00; // type of current marker segment
	unsigned int   len  = 0; // length of current marker segment
	unsigned int   hpos = 0; // position in header
	
	unsigned int fpos; // end of marker position
	unsigned int skip; // bytes to skip
	unsigned int spos; // sub position	
	int i;
	
	
	// search for DHT (0xFFC4) & DQT (0xFFDB) marker segments	
	// header parser loop
	while ( ( int ) hpos < hdrs ) {
		type = hdrdata[ hpos + 1 ];
		len = 2 + B_SHORT( hdrdata[ hpos + 2 ], hdrdata[ hpos + 3 ] );
		
		if ( type == 0xC4 ) { // for DHT
			fpos = hpos + len; // reassign length to end position
			hpos += 4; // skip marker & length
			while ( hpos < fpos ) {			
				hpos++;
				// table found - check if modified
				if ( hdrdata[ hpos ] > 2 ) {	
					// reinsert the standard table
					i = hdrdata[ hpos + 1 ];
					for ( spos = 0; spos < std_huff_lengths[ i ]; spos++ ) {
						hdrdata[ hpos + spos ] = std_huff_tables[ i ][ spos ];
					}
				}
								
				skip = 16;
				for ( i = 0; i < 16; i++ )		
					skip += ( int ) hdrdata[ hpos + i ];				
				hpos += skip;
			}
		}
		else if ( type == 0xDB ) { // for DQT
			fpos = hpos + len; // reassign length to end position
			hpos += 4; // skip marker & length
			while ( hpos < fpos ) {
				i = LBITS( hdrdata[ hpos ], 4 );				
				hpos++;
				// table found
				if ( i == 1 ) { // get out for 16 bit precision
					hpos += 128;
					continue;
				}
				// undo diff coding for 8 bit precision
				for ( spos = 1; spos < 64; spos++ )
					hdrdata[ hpos + spos ] += hdrdata[ hpos + spos - 1 ];
					
				hpos += 64;
			}
		}
		else { // skip segment
			hpos += len;
		}		
	}
	
	
	return true;
}


/* -----------------------------------------------
	preparations for special average context
	----------------------------------------------- */
INTERN void pjg_aavrg_prepare( unsigned short** abs_coeffs, int* weights, unsigned short* abs_store, int cmp )
{
	int w = cmpnfo[cmp].bch;
	
	// set up quick access arrays for all prediction positions
	abs_coeffs[ 0 ] = abs_store + (  0 + ((-2)*w) ); // top-top
	abs_coeffs[ 1 ] = abs_store + ( -1 + ((-1)*w) ); // top-left
	abs_coeffs[ 2 ] = abs_store + (  0 + ((-1)*w) ); // top
	abs_coeffs[ 3 ] = abs_store + (  1 + ((-1)*w) ); // top-right
	abs_coeffs[ 4 ] = abs_store + ( -2 + (( 0)*w) ); // left-left
	abs_coeffs[ 5 ] = abs_store + ( -1 + (( 0)*w) ); // left
	// copy context weighting factors
	weights[ 0 ] = abs_ctx_weights_lum[ 0 ][ 0 ][ 2 ]; // top-top
	weights[ 1 ] = abs_ctx_weights_lum[ 0 ][ 1 ][ 1 ]; // top-left
	weights[ 2 ] = abs_ctx_weights_lum[ 0 ][ 1 ][ 2 ]; // top
	weights[ 3 ] = abs_ctx_weights_lum[ 0 ][ 1 ][ 3 ]; // top-right
	weights[ 4 ] = abs_ctx_weights_lum[ 0 ][ 2 ][ 0 ]; // left-left
	weights[ 5 ] = abs_ctx_weights_lum[ 0 ][ 2 ][ 1 ]; // left
}


/* -----------------------------------------------
	special average context used in coeff encoding
	----------------------------------------------- */
INTERN int pjg_aavrg_context( unsigned short** abs_coeffs, int* weights, int pos, int p_y, int p_x, int r_x )
{
	int ctx_avr = 0; // AVERAGE context
	int w_ctx = 0; // accumulated weight of context
	int w_curr; // current weight of context
	
	
	// different cases due to edge treatment
	if ( p_y >= 2 ) {
		w_curr = weights[ 0 ]; ctx_avr += abs_coeffs[ 0 ][ pos ] * w_curr; w_ctx += w_curr;
		w_curr = weights[ 2 ]; ctx_avr += abs_coeffs[ 2 ][ pos ] * w_curr; w_ctx += w_curr;
		if ( p_x >= 2 ) {
			w_curr = weights[ 1 ]; ctx_avr += abs_coeffs[ 1 ][ pos ] * w_curr; w_ctx += w_curr;
			w_curr = weights[ 4 ]; ctx_avr += abs_coeffs[ 4 ][ pos ] * w_curr; w_ctx += w_curr;
			w_curr = weights[ 5 ]; ctx_avr += abs_coeffs[ 5 ][ pos ] * w_curr; w_ctx += w_curr;
		}
		else if ( p_x == 1 ) {
			w_curr = weights[ 1 ]; ctx_avr += abs_coeffs[ 1 ][ pos ] * w_curr; w_ctx += w_curr;
			w_curr = weights[ 5 ]; ctx_avr += abs_coeffs[ 5 ][ pos ] * w_curr; w_ctx += w_curr;
		}
		if ( r_x >= 1 ) {
			w_curr = weights[ 3 ]; ctx_avr += abs_coeffs[ 3 ][ pos ] * w_curr; w_ctx += w_curr;
		}
	}
	else if ( p_y == 1 ) {
		w_curr = weights[ 2 ]; ctx_avr += abs_coeffs[ 2 ][ pos ] * w_curr; w_ctx += w_curr;
		if ( p_x >= 2 ) {
			w_curr = weights[ 1 ]; ctx_avr += abs_coeffs[ 1 ][ pos ] * w_curr; w_ctx += w_curr;
			w_curr = weights[ 4 ]; ctx_avr += abs_coeffs[ 4 ][ pos ] * w_curr; w_ctx += w_curr;
			w_curr = weights[ 5 ]; ctx_avr += abs_coeffs[ 5 ][ pos ] * w_curr; w_ctx += w_curr;
		}
		else if ( p_x == 1 ) {
			w_curr = weights[ 1 ]; ctx_avr += abs_coeffs[ 1 ][ pos ] * w_curr; w_ctx += w_curr;
			w_curr = weights[ 5 ]; ctx_avr += abs_coeffs[ 5 ][ pos ] * w_curr; w_ctx += w_curr;
		}
		if ( r_x >= 1 ) {
			w_curr = weights[ 3 ]; ctx_avr += abs_coeffs[ 3 ][ pos ] * w_curr; w_ctx += w_curr;
		}
	}
	else {
		if ( p_x >= 2 ) {
			w_curr = weights[ 4 ]; ctx_avr += abs_coeffs[ 4 ][ pos ] * w_curr; w_ctx += w_curr;
			w_curr = weights[ 5 ]; ctx_avr += abs_coeffs[ 5 ][ pos ] * w_curr; w_ctx += w_curr;
		}
		else if ( p_x == 1 ) {
			w_curr = weights[ 5 ]; ctx_avr += abs_coeffs[ 5 ][ pos ] * w_curr; w_ctx += w_curr;
		}
	}
	
	// return average context
	return ( w_ctx != 0 ) ? ( ctx_avr + ( w_ctx / 2 ) ) / w_ctx : 0;
}


/* -----------------------------------------------
	lakhani ac context used in coeff encoding
	----------------------------------------------- */
INTERN int pjg_lakh_context( signed short** coeffs_x, signed short** coeffs_a, int* pred_cf, int pos )
{
	int pred = 0;
	
	// calculate partial prediction
	pred -= ( coeffs_x[ 1 ][ pos ] + coeffs_a[ 1 ][ pos ] ) * pred_cf[ 1 ];
	pred -= ( coeffs_x[ 2 ][ pos ] - coeffs_a[ 2 ][ pos ] ) * pred_cf[ 2 ];
	pred -= ( coeffs_x[ 3 ][ pos ] + coeffs_a[ 3 ][ pos ] ) * pred_cf[ 3 ];
	pred -= ( coeffs_x[ 4 ][ pos ] - coeffs_a[ 4 ][ pos ] ) * pred_cf[ 4 ];
	pred -= ( coeffs_x[ 5 ][ pos ] + coeffs_a[ 5 ][ pos ] ) * pred_cf[ 5 ];
	pred -= ( coeffs_x[ 6 ][ pos ] - coeffs_a[ 6 ][ pos ] ) * pred_cf[ 6 ];
	pred -= ( coeffs_x[ 7 ][ pos ] + coeffs_a[ 7 ][ pos ] ) * pred_cf[ 7 ];
	// normalize / quantize partial prediction
	pred = ( ( pred > 0 ) ? ( pred + (pred_cf[0]/2) ) : ( pred - (pred_cf[0]/2) ) ) / pred_cf[ 0 ];
	// complete prediction
	pred += coeffs_a[ 0 ][ pos ];
	
	return pred;
}


/* -----------------------------------------------
	Calculates coordinates for nearest neighbor context
	----------------------------------------------- */
INTERN void get_context_nnb( int pos, int w, int *a, int *b )
{
	// this function calculates and returns coordinates for
	// a simple 2D context
	if ( pos == 0 ) {
		*a = -1;
		*b = -1;
	}
	else if ( ( pos % w ) == 0 ) {
		*b = pos - w;
		if ( pos >= ( w << 1 ) )
			*a = pos - ( w << 1 );
		else
			*a = *b;
	}
	else if ( pos < w ) {
		*a = pos - 1;
		if ( pos >= 2 )
			*b = pos - 2;
		else
			*b = *a;
	}
	else {
		*a = pos - 1;
		*b = pos - w;
	}
}

/* ----------------------- End of PJG specific functions -------------------------- */

/* ----------------------- Begin ofDCT specific functions -------------------------- */



/* -----------------------------------------------
	inverse DCT transform using precalc tables (fast)
	----------------------------------------------- */
INTERN int idct_2d_fst_8x1( int cmp, int dpos, int ix, int iy )
{
	int idct = 0;
	int ixy;
	
	
	// calculate start index
	ixy = ix << 3;
	
	// begin transform
	idct += colldata[ cmp ][  0 ][ dpos ] * adpt_idct_8x1[ cmp ][ ixy + 0 ];
	idct += colldata[ cmp ][  1 ][ dpos ] * adpt_idct_8x1[ cmp ][ ixy + 1 ];
	idct += colldata[ cmp ][  5 ][ dpos ] * adpt_idct_8x1[ cmp ][ ixy + 2 ];
	idct += colldata[ cmp ][  6 ][ dpos ] * adpt_idct_8x1[ cmp ][ ixy + 3 ];
	idct += colldata[ cmp ][ 14 ][ dpos ] * adpt_idct_8x1[ cmp ][ ixy + 4 ];
	idct += colldata[ cmp ][ 15 ][ dpos ] * adpt_idct_8x1[ cmp ][ ixy + 5 ];
	idct += colldata[ cmp ][ 27 ][ dpos ] * adpt_idct_8x1[ cmp ][ ixy + 6 ];
	idct += colldata[ cmp ][ 28 ][ dpos ] * adpt_idct_8x1[ cmp ][ ixy + 7 ];
	
	
	return idct;
}


/* -----------------------------------------------
	inverse DCT transform using precalc tables (fast)
	----------------------------------------------- */
INTERN int idct_2d_fst_1x8( int cmp, int dpos, int ix, int iy )
{
	int idct = 0;
	int ixy;
	
	
	// calculate start index
	ixy = iy << 3;
	
	// begin transform
	idct += colldata[ cmp ][  0 ][ dpos ] * adpt_idct_1x8[ cmp ][ ixy + 0 ];
	idct += colldata[ cmp ][  2 ][ dpos ] * adpt_idct_1x8[ cmp ][ ixy + 1 ];
	idct += colldata[ cmp ][  3 ][ dpos ] * adpt_idct_1x8[ cmp ][ ixy + 2 ];
	idct += colldata[ cmp ][  9 ][ dpos ] * adpt_idct_1x8[ cmp ][ ixy + 3 ];
	idct += colldata[ cmp ][ 10 ][ dpos ] * adpt_idct_1x8[ cmp ][ ixy + 4 ];
	idct += colldata[ cmp ][ 20 ][ dpos ] * adpt_idct_1x8[ cmp ][ ixy + 5 ];
	idct += colldata[ cmp ][ 21 ][ dpos ] * adpt_idct_1x8[ cmp ][ ixy + 6 ];
	idct += colldata[ cmp ][ 35 ][ dpos ] * adpt_idct_1x8[ cmp ][ ixy + 7 ];
	
	
	return idct;
}

/* ----------------------- End of DCT specific functions -------------------------- */

/* ----------------------- Begin of prediction functions -------------------------- */


/* -----------------------------------------------
	returns predictor for collection data
	----------------------------------------------- */
#if defined(USE_PLOCOI)
INTERN int dc_coll_predictor( int cmp, int dpos )
{
	signed short* coefs = colldata[ cmp ][ 0 ];
	int w = cmpnfo[cmp].bch;
	int a, b, c;
	
	if ( dpos < w ) {
		a = coefs[ dpos - 1 ];
		b = 0;
		c = 0;
	}
	else if ( (dpos%w) == 0 ) {
		a = 0;
		b = coefs[ dpos - w ];
		c = 0;
	}
	else {
		a = coefs[ dpos - 1 ];
		b = coefs[ dpos - w ];
		c = coefs[ dpos - 1 - w ];
	}
	
	return plocoi( a, b, c );
}
#endif


/* -----------------------------------------------
	1D DCT predictor for DC coefficients
	----------------------------------------------- */
#if !defined(USE_PLOCOI)
INTERN int dc_1ddct_predictor( int cmp, int dpos )
{
	int w  = cmpnfo[cmp].bch;
	int px = ( dpos % w );
	int py = ( dpos / w );
	
	int pred;	
	int pa = 0;
	int pb = 0;
	int xa = 0;
	int xb = 0;
	int swap;
	
	
	// store current block DC coefficient
	swap = colldata[ cmp ][ 0 ][ dpos ];
	colldata[ cmp ][ 0 ][ dpos ] = 0;
	
	// calculate prediction
	if ( ( px > 0 ) && ( py > 0 ) ) {
		pa = idct_2d_fst_8x1( cmp, dpos - 1, 7, 0 );
		pb = idct_2d_fst_1x8( cmp, dpos - w, 0, 7 );
		xa = idct_2d_fst_8x1( cmp, dpos, 0, 0 );
		xb = idct_2d_fst_1x8( cmp, dpos, 0, 0 );
		pred = ( ( pa - xa ) + ( pb - xb ) ) * ( 8 / 2 );
	}
	else if ( px > 0 ) {
		pa = idct_2d_fst_8x1( cmp, dpos - 1, 7, 0 );
		xa = idct_2d_fst_8x1( cmp, dpos, 0, 0 );
		pred = ( pa - xa ) * 8;
	}
	else if ( py > 0 ) {
		pb = idct_2d_fst_1x8( cmp, dpos - w, 0, 7 );
		xb = idct_2d_fst_1x8( cmp, dpos, 0, 0 );
		pred = ( pb - xb ) * 8;
	}
	else {
		pred = 0;
	}
	
	// write back current DCT coefficient
	colldata[ cmp ][ 0 ][ dpos ] = swap;	
	
	// clamp and quantize predictor
	pred = CLAMPED( -( 1024 * DCT_RSC_FACTOR ), ( 1016 * DCT_RSC_FACTOR ), pred );
	pred = pred / QUANT( cmp, 0 );
	pred = DCT_RESCALE( pred );
	
	
	return pred;
}
#endif


/* -----------------------------------------------
	calculates median out of an integer array
	----------------------------------------------- */
INTERN inline int median_int( int* values, int size )
{
	int middle = ( size >> 1 );
	bool done;
	int swap;
	int i;
	
	
	// sort data first
	done = false;
	while ( !done ) {
		done = true;
		for ( i = 1; i < size; i++ )
		if ( values[ i ] < values[ i - 1 ] ) {
			swap = values[ i ];
			values[ i ] = values[ i - 1 ];
			values[ i - 1 ] = swap;
			done = false;
		}
	}
	
	// return median
	return ( ( size % 2 ) == 0 ) ?
		( values[ middle ] + values[ middle - 1 ] ) / 2 : values[ middle ];
}


/* -----------------------------------------------
	calculates median out of an float array
	----------------------------------------------- */
INTERN inline float median_float( float* values, int size )
{
	int middle = ( size >> 1 );
	bool done;
	float swap;
	int i;
	
	
	// sort data first
	done = false;
	while ( !done ) {
		done = true;
		for ( i = 1; i < size; i++ )
		if ( values[ i ] < values[ i - 1 ] ) {
			swap = values[ i ];
			values[ i ] = values[ i - 1 ];
			values[ i - 1 ] = swap;
			done = false;
		}
	}
	
	// return median	
	if ( ( size % 2 ) == 0 ) {
		return ( values[ middle ] + values[ middle - 1 ] ) / 2.0;
	}
	else
		return ( values[ middle ] );
}

/* ----------------------- End of prediction functions -------------------------- */

/* ----------------------- Begin of miscellaneous helper functions -------------------------- */


/* -----------------------------------------------
	displays progress bar on screen
	----------------------------------------------- */
#if !defined(BUILD_LIB)
INTERN inline void progress_bar( int current, int last )
{
	int barpos = ( ( current * BARLEN ) + ( last / 2 ) ) / last;
	int i;
	
	
	// generate progress bar
	fprintf( msgout, "[" );
	#if defined(_WIN32)
	for ( i = 0; i < barpos; i++ )
		fprintf( msgout, "\xFE" );
	#else
	for ( i = 0; i < barpos; i++ )
		fprintf( msgout, "X" );
	#endif
	for (  ; i < BARLEN; i++ )
		fprintf( msgout, " " );
	fprintf( msgout, "]" );
}
#endif

/* -----------------------------------------------
	creates filename, callocs memory for it
	----------------------------------------------- */
#if !defined(BUILD_LIB)
INTERN inline char* create_filename( const char* base, const char* extension )
{
	int len = strlen( base ) + ( ( extension == NULL ) ? 0 : strlen( extension ) + 1 ) + 1;	
	char* filename = (char*) calloc( len, sizeof( char ) );	
	
	// create a filename from base & extension
	strcpy( filename, base );
	set_extension( filename, extension );
	
	return filename;
}
#endif

/* -----------------------------------------------
	creates filename, callocs memory for it
	----------------------------------------------- */
#if !defined(BUILD_LIB)
INTERN inline char* unique_filename( const char* base, const char* extension )
{
	int len = strlen( base ) + ( ( extension == NULL ) ? 0 : strlen( extension ) + 1 ) + 1;	
	char* filename = (char*) calloc( len, sizeof( char ) );	
	
	// create a unique filename using underscores
	strcpy( filename, base );
	set_extension( filename, extension );
	while ( file_exists( filename ) ) {
		len += sizeof( char );
		filename = (char*) realloc( filename, len );
		add_underscore( filename );
	}
	
	return filename;
}
#endif

/* -----------------------------------------------
	changes extension of filename
	----------------------------------------------- */
#if !defined(BUILD_LIB)
INTERN inline void set_extension( char* filename, const char* extension )
{
	char* extstr;
	
	// find position of extension in filename	
	extstr = ( strrchr( filename, '.' ) == NULL ) ?
		strrchr( filename, '\0' ) : strrchr( filename, '.' );
	
	// set new extension
	if ( extension != NULL ) {
		(*extstr++) = '.';
		strcpy( extstr, extension );
	}
	else
		(*extstr) = '\0';
}
#endif

/* -----------------------------------------------
	adds underscore after filename
	----------------------------------------------- */
#if !defined(BUILD_LIB)
INTERN inline void add_underscore( char* filename )
{
	char* tmpname = (char*) calloc( strlen( filename ) + 1, sizeof( char ) );
	char* extstr;
	
	// copy filename to tmpname
	strcpy( tmpname, filename );
	// search extension in filename
	extstr = strrchr( filename, '.' );
	
	// add underscore before extension
	if ( extstr != NULL ) {
		(*extstr++) = '_';
		strcpy( extstr, strrchr( tmpname, '.' ) );
	}
	else
		sprintf( filename, "%s_", tmpname );
		
	// free memory
	free( tmpname );
}
#endif

/* -----------------------------------------------
	checks if a file exists
	----------------------------------------------- */
INTERN inline bool file_exists( const char* filename )
{
	// needed for both, executable and library
	FILE* fp = fopen( filename, "rb" );
	
	if ( fp == NULL ) return false;
	else {
		fclose( fp );
		return true;
	}
}

/* ----------------------- End of miscellaneous helper functions -------------------------- */

/* ----------------------- End of file -------------------------- */
