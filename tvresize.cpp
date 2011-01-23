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


inline double hanning2( double x ) {
	if( abs( x ) < 1.0 ) {
		return 1.0 + cos( M_PI * x );
	}
	else {
		return 0.0;
	}
}

inline double sinc( double x ) {
	if( abs( x ) <= numeric_limits<double>::min() ) {
		return 1.0;
	}
	else {
		return sin( M_PI * x ) / (M_PI * x);
	}
}

struct Clip {
	boost::gil::rgb8_pixel_t operator()( boost::gil::rgb32f_pixel_t src ) {
		return boost::gil::rgb8_pixel_t(
			min( max( int( src[0] * 255.0f + 0.5f ), 0 ), 255 ),
			min( max( int( src[1] * 255.0f + 0.5f ), 0 ), 255 ),
			min( max( int( src[2] * 255.0f + 0.5f ), 0 ), 255 )
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
					(
						+ (dx * dx) * (double(src( 0, +1)[c]) + double(src( 0, -1)[c]))
						+ (dy * dy) * (double(src(+1,  0)[c]) + double(src(-1,  0)[c]))
					) * 2.0 - dx * dy * dxdy
				) / (dx * dx + dy * dy);
				dst[c] = (src(0, 0)[c] + np) * (1.0 / 5.0);
			}   
		}   
		return dst;
	}   
};

template<class SrcView, class OrgView, class DstView>
void constrain( SrcView const& src, OrgView const& org, DstView& dst ) {
	using namespace boost::gil;
	assert( src.dimensions() == dst.dimensions() );
	assert( src.width() % org.width() == 0 );
	assert( src.height() % org.height() == 0 );

	int const mw = src.width() / org.width();
	int const mh = src.height() / org.height();

	double sum = 0.0;
	double kern[mw * 2][mh * 2];
	for( int dy = 0; dy < mh * 2; ++dy ) {
		for( int dx = 0; dx < mw * 2; ++dx ) {
			double rx = double( dx - mw ) / mw;
			double ry = double( dy - mh ) / mh;
			//kern[dx][dy] = hanning2( max( abs( rx ), abs( ry ) ) );
			kern[dx][dy] = hanning2( rx ) * hanning2( ry ) / (mw * mh * 4);
			//kern[dx][dy] = hanning2( sqrt( rx * rx + ry * ry ) );
			//kern[dx][dy] = sinc( rx ) * sinc( rx / 3.0 ) * sinc( ry ) * sinc( ry / 3.0 );
			sum += kern[dx][dy];
		}
	}
	/*
	for( int dy = 0; dy < mh * 2; ++dy ) {
		for( int dx = 0; dx < mw * 2; ++dx ) {
			kern[dx][dy] /= sum;
		}
	}
	*/

	//copy_pixels( src, dst );
	fill_pixels( dst, rgb32f_pixel_t( 0.0, 0.0, 0.0 ) );
	for( int y = 0; y < org.height() - 2; ++y ) {
		for( int x = 0; x < org.width() - 2; ++x ) {
			for( int c = 0; c < num_channels<DstView>::value; ++c ) {
				double s = 0.0;
				for( int dy = 0; dy < mh * 2; ++dy ) {
					for( int dx = 0; dx < mw * 2; ++dx ) {
						s += kern[dx][dy] * src(x * mw + dx, y * mh + dy)[c];
					}
				}
				double d = org(x, y)[c] - s;
				for( int dy = 0; dy < mh * 2; ++dy ) {
					for( int dx = 0; dx < mw * 2; ++dx ) {
						dst(x * mw + dx, y * mh + dy)[c] += kern[dx][dy] * d;
					}
				}
			}
		}
	}
	for( int y = 0; y < src.height(); ++y ) {
		for( int x = 0; x < src.width(); ++x ) {
			for( int c = 0; c < num_channels<DstView>::value; ++c ) {
				dst(x, y)[c] += src(x, y)[c];
			}
		}
	}
}

int main( int argc, char* const* argv ) {
	using namespace boost::gil;

	int scale = -1;
	int maxCount = -1;
	string srcFile, dstFile;
	try {
		char ch;
		while( ch = getopt( argc, argv, "s:c:" ), ch >= 0 ) {
			switch( ch ) {
				case 's':
					scale = boost::lexical_cast<int>( optarg );
					break;
				case 'c':
					maxCount = boost::lexical_cast<int>( optarg );
					break;
				default:
					throw exception();
			}
		}
		if( scale < 0 || maxCount < 0 ) {
			throw exception();
		}
		if( argc - optind != 2 ) {
			throw exception();
		}
		srcFile = argv[optind + 0];
		dstFile = argv[optind + 1];
	}
	catch( exception ) {
		cout << "Usage: tvresize -s [scale] -c [iteration count] src.png dst.png\n";
		return 1;
	}

	try {
		rgb32f_image_t src;
		png_read_and_convert_image( srcFile, src );
		rgb32f_image_t tmp0( src.dimensions() * scale );
		rgb32f_image_t tmp1( src.dimensions() * scale );
		rgb32f_view_t sub0 = subimage_view( view( tmp0 ), 1, 1, tmp0.width() - 2, tmp0.height() - 2 );
		rgb32f_view_t sub1 = subimage_view( view( tmp1 ), 1, 1, tmp1.width() - 2, tmp1.height() - 2 );

		generate_pixels( view( tmp0 ), RandomGray() );
		for( int i = 0; i < maxCount; ++i ) {
			cout << "\riter step: " << i;
			cout.flush();
			transform_pixel_positions( sub0, sub1,
				TVMinimizer<rgb32f_pixel_t, rgb32f_loc_t>()
			);
			constrain( view( tmp1 ), view( src ), view( tmp0 ) );
		}
		cout << endl;

		rgb8_image_t dst( tmp0.dimensions() );
		transform_pixels( view( tmp0 ), view( dst ), Clip() );
		png_write_view( dstFile, view( dst ) );
	}
	catch( exception const& e ) {
		cout << endl;
		cout << e.what() << endl;
		return 1;
	}

	return 0;
}
