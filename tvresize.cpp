// Image enlarger with total variation minimization
// by y.fujii <y-fujii at mimosa-pudica.net>, public domain

#include <algorithm>
#include <exception>
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <boost/lexical_cast.hpp>
#include <boost/gil/gil_all.hpp>
#define int_p_NULL 0
#define png_infopp_NULL 0
#include <boost/gil/extension/io/png_io.hpp>

using namespace std;


template<class SrcView, class DstView, class Func>
Func transform_pixel_positions_parallel( SrcView const& src , DstView const& dst, Func f ) {
	assert( src.dimensions() == dst.dimensions() );

	#pragma omp parallel for
	for( int y = 0; y < src.height(); ++y ) {
		typename SrcView::xy_locator srcLoc = src.xy_at( 0, y );
		typename DstView::x_iterator dstIt = dst.row_begin( y );
		for( int x = 0; x < src.width(); ++x ) {
			*dstIt = f( srcLoc );
			++srcLoc.x();
			++dstIt;
		}
	}        
	return f;
}

struct Clip {
	boost::gil::rgb8_pixel_t operator()( boost::gil::rgb32f_pixel_t src ) {
		return boost::gil::rgb8_pixel_t(
			int( min( max( src[0] * 255.0f + 0.5f, 0.0f ), 255.0f ) ),
			int( min( max( src[1] * 255.0f + 0.5f, 0.0f ), 255.0f ) ),
			int( min( max( src[2] * 255.0f + 0.5f, 0.0f ), 255.0f ) )
		);
	}
};

struct RandomGray {
	boost::gil::rgb32f_pixel_t operator()() {
		double r = double( rand() ) / RAND_MAX;
		return boost::gil::rgb32f_pixel_t( r, r, r );
	}
};

template<class DstPixel, class SrcLocator>
struct TVMinimizer {
	DstPixel operator()( SrcLocator const& src ) { 
		using namespace boost::gil;

		DstPixel dst;
		for( int c = 0; c < num_channels<DstPixel>::value; ++c ) { 
			double dx = double(src(+1,  0)[c]) - double(src(-1,  0)[c]);
			double dy = double(src( 0, +1)[c]) - double(src( 0, -1)[c]);
			if( dx * dx + dy * dy <= numeric_limits<float>::min() ) { 
				dst[c] = src(0, 0)[c];
			}   
			else {
				double dxdy = ( 
					+ (double(src(+1, +1)[c]) + double(src(-1, -1)[c]))
					- (double(src(+1, -1)[c]) + double(src(-1, +1)[c]))
				);  
				double np = ( 
					+ (dx * dx) * (double(src( 0, +1)[c]) + double(src( 0, -1)[c]))
					+ (dy * dy) * (double(src(+1,  0)[c]) + double(src(-1,  0)[c]))
					- dx * dy * dxdy * (1.0 / 2.0)
				) / (dx * dx + dy * dy);
				dst[c] = (np + src(0, 0)[c]) * (1.0 / 3.0);
			}   
		}   
		return dst;
	}   
};

template<class SrcView, class OrgView, class DstView>
void addConstraint( SrcView const& src, OrgView const& org, DstView const& dst, double fac ) {
	using namespace boost::gil;
	assert( src.dimensions() == dst.dimensions() );
	assert( src.width()  % org.width()  == 0 );
	assert( src.height() % org.height() == 0 );

	int const mw = src.width() / org.width();
	int const mh = src.height() / org.height();

	#pragma omp parallel for
	for( int y = 0; y < org.height(); ++y ) {
		for( int x = 0; x < org.width(); ++x ) {
			for( int c = 0; c < num_channels<DstView>::value; ++c ) {
				double s = 0.0;
				for( int sy = y * mh; sy < y * mh + mh; ++sy ) {
					for( int sx = x * mw; sx < x * mw + mw; ++sx ) {
						s += src(sx, sy)[c];
					}
				}
				double d = (org(x, y)[c] - s * (1.0 / (mw * mh))) * fac;
				for( int sy = y * mh; sy < y * mh + mh; ++sy ) {
					for( int sx = x * mw; sx < x * mw + mw; ++sx ) {
						dst(sx, sy)[c] += d;
					}
				}
			}
		}
	}
}

template<class SrcView, class OrgView, class DstView>
void iterate( SrcView const& src, OrgView const& org, DstView const& dst, double fac ) {
	using namespace boost::gil;

	int const zw = src.width()  - 2;
	int const zh = src.height() - 2;

	assert( src.dimensions() == dst.dimensions() );
	assert( src.width()  == dst.width()  );
	assert( src.height() == dst.height() );
	assert( zw % org.width()  == 0 );
	assert( zh % org.height() == 0 );

	transform_pixel_positions_parallel(
		subimage_view( src, 1, 1, zw, zh ),
		subimage_view( dst, 1, 1, zw, zh ),
		TVMinimizer<rgb32f_pixel_t, rgb32f_loc_t>()
	);
	addConstraint(
		subimage_view( src, 1, 1, zw, zh ),
		org,
		subimage_view( dst, 1, 1, zw, zh ),
		fac
	);
	copy_pixels(
		subimage_view( dst, 1, 1, zw, 1 ),
		subimage_view( dst, 1, 0, zw, 1 )
	);
	copy_pixels(
		subimage_view( dst, 1, zh    , zw, 1 ),
		subimage_view( dst, 1, zh + 1, zw, 1 )
	);
	copy_pixels(
		subimage_view( dst, 1, 1, 1, zh ),
		subimage_view( dst, 0, 1, 1, zh )
	);
	copy_pixels(
		subimage_view( dst, zw    , 1, 1, zh ),
		subimage_view( dst, zw + 1, 1, 1, zh )
	);
	for( int c = 0; c < 3; ++c ) {
		dst(     0,      0)[c] = dst( 1, 0)[c] + dst(0,  1)[c] - dst(     1,      1)[c];
		dst(zw + 1,      0)[c] = dst(zw, 0)[c] + dst(0,  1)[c] - dst(zw + 1,      1)[c];
		dst(     0, zh + 1)[c] = dst( 1, 0)[c] + dst(0, zh)[c] - dst(     1, zh + 1)[c];
		dst(zw + 1, zh + 1)[c] = dst(zw, 0)[c] + dst(0, zh)[c] - dst(zw + 1, zh + 1)[c];
	}
}

int main( int argc, char* const* argv ) {
	using namespace boost::gil;

	int scale = -1;
	int maxCount = -1;
	double smooth = 1.0;
	string srcFile, dstFile;
	try {
		char ch;
		while( ch = getopt( argc, argv, "s:c:f:" ), ch >= 0 ) {
			switch( ch ) {
				case 's':
					scale = boost::lexical_cast<int>( optarg );
					break;
				case 'c':
					maxCount = boost::lexical_cast<int>( optarg );
					break;
				case 'f':
					smooth = boost::lexical_cast<double>( optarg );
					break;
				default:
					throw exception();
			}
		}
		if( scale < 0 || maxCount < 0 || smooth < 1.0 ) {
			throw exception();
		}
		if( argc - optind != 2 ) {
			throw exception();
		}
		srcFile = argv[optind + 0];
		dstFile = argv[optind + 1];
	}
	catch( exception ) {
		cout << "Usage: tvresize -s <scale> -c <iteration count> [-f <smoothness>] src.png dst.png\n";
		return 1;
	}

	try {
		rgb32f_image_t src;
		png_read_and_convert_image( srcFile, src );
		int const zw = src.width()  * scale;
		int const zh = src.height() * scale;
		rgb32f_image_t img0( zw + 2, zh + 2 );
		rgb32f_image_t img1( zw + 2, zh + 2 );
		rgb32f_view_t view0 = view( img0 );
		rgb32f_view_t view1 = view( img1 );

		generate_pixels( view0, RandomGray() );
		for( int i = 0; i < maxCount; ++i ) {
			cout << "\riter step: " << i;
			cout.flush();

			iterate( view0, view( src ), view1, 1.0 / smooth );
			iterate( view1, view( src ), view0, 1.0 / smooth );
		}
		cout << endl;

		rgb8_image_t dst( zw, zh );
		transform_pixels(
			subimage_view( view0, 1, 1, zw, zh ),
			view( dst ),
			Clip()
		);
		png_write_view( dstFile, view( dst ) );
	}
	catch( exception const& e ) {
		cout << endl;
		cout << e.what() << endl;
		return 1;
	}

	return 0;
}
