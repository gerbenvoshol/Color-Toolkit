/* clrtk.h - v1.0 - simple color distance Toolkit -- PUBLIC DOMAIN
					no warranty is offered or implied; use this code at your own risk

	 This is a single header file with a bunch of useful HMM functions

 ============================================================================
	 You MUST

			#define CLRTK_DEFINE

	 in EXACTLY _one_ C or C++ file that includes this header, BEFORE the
	 include, like this:

			#define CLRTK_DEFINE
			#include "clrtk.h"

	 All other files should just #include "clrtk.h" without the #define.
 ============================================================================

 Version History
		1.00  Initial release

 CITATION

 If you use this color Toolkit in a publication, please reference:
 Voshol, G.P. (2020). clrtk: A simple color Toolkit (Version 1.0) [Software]. 
 Available from https://github.com/gerbenvoshol/Color-Toolkit

 LICENSE

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

 CREDITS


 REFERENCES

 - http://www.easyrgb.com/index.php?X=MATH&H=02text2

  Written by Gerben Voshol.
 */

// http://www.easyrgb.com/index.php?X=MATH&H=02text2
// Observer	2° (CIE 1931)	10° (CIE 1964)	Note
// Illuminant	X2	Y2	Z2	X10	Y10	Z10	 
// A	109.850	100.000	35.585	111.144	100.000	35.200	Incandescent/tungsten
// B	99.0927	100.000	85.313	99.178;	100.000	84.3493	Old direct sunlight at noon
// C	98.074	100.000	118.232	97.285	100.000	116.145	Old daylight
// D50	96.422	100.000	82.521	96.720	100.000	81.427	ICC profile PCS
// D55	95.682	100.000	92.149	95.799	100.000	90.926	Mid-morning daylight
// D65	95.047	100.000	108.883	94.811	100.000	107.304	Daylight, sRGB, Adobe-RGB
// D75	94.972	100.000	122.638	94.416	100.000	120.641	North sky daylight
// E	100.000	100.000	100.000	100.000	100.000	100.000	Equal energy
// F1	92.834	100.000	103.665	94.791	100.000	103.191	Daylight Fluorescent
// F2	99.187	100.000	67.395	103.280	100.000	69.026	Cool fluorescent
// F3	103.754	100.000	49.861	108.968	100.000	51.965	White Fluorescent
// F4	109.147	100.000	38.813	114.961	100.000	40.963	Warm White Fluorescent
// F5	90.872	100.000	98.723	93.369	100.000	98.636	Daylight Fluorescent
// F6	97.309	100.000	60.191	102.148	100.000	62.074	Lite White Fluorescent
// F7	95.044	100.000	108.755	95.792	100.000	107.687	Daylight fluorescent, D65 simulator
// F8	96.413	100.000	82.333	97.115	100.000	81.135	Sylvania F40, D50 simulator
// F9	100.365	100.000	67.868	102.116	100.000	67.826	Cool White Fluorescent
// F10	96.174	100.000	81.712	99.001	100.000	83.134	Ultralume 50, Philips TL85
// F11	100.966	100.000	64.370	103.866	100.000	65.627	Ultralume 40, Philips TL84
// F12	108.046	100.000	39.228	111.428	100.000	40.353	Ultralume 30, Philips TL83

#ifndef CLRTK__H
#define CLRTK__H

#include <stdio.h>
#include <math.h>

#ifndef MAX
#define MAX(x,y)        ((x) > (y) ? (x) : (y))
#endif
#ifndef MIN
#define MIN(x,y)        ((x) < (y) ? (x) : (y))
#endif

#ifdef __cplusplus
#define CLR_EXTERN   extern "C"
#else
#define CLR_EXTERN   extern
#endif

CLR_EXTERN void convertRGBtoXYZ(int inR, int inG, int inB, double * outX, double * outY, double * outZ);
CLR_EXTERN void convertXYZtoLab(double inX, double inY, double inZ, double * outL, double * outa, double * outb);
CLR_EXTERN void convertLabtoXYZ( double inL, double ina, double  inb, double * outX, double * outY, double * outZ);
CLR_EXTERN void convertXYZtoRGB(double inX, double inY, double inZ, int * outR, int * outG, int * outB);

CLR_EXTERN double Lab_color_difference_CIE76( double inL1, double ina1, double  inb1, double inL2, double ina2, double  inb2);
CLR_EXTERN double RGB_color_Lab_difference_CIE76( int R1, int G1, int B1, int R2, int G2, int B2);
CLR_EXTERN double Lab_color_difference_CIE94( double inL1, double ina1, double  inb1, double inL2, double ina2, double  inb2);
CLR_EXTERN double RGB_color_Lab_difference_CIE94( int R1, int G1, int B1, int R2, int G2, int B2);
CLR_EXTERN double Lab_color_difference_CIEDE2000( double inL1, double ina1, double  inb1, double inL2, double ina2, double  inb2);
CLR_EXTERN double RGB_color_Lab_difference_CIEDE2000( int R1, int G1, int B1, int R2, int G2, int B2);

#ifdef CLRTK_DEFINE

//
// double ref_X = 96.422;
// double ref_Y = 100.0;
// double ref_Z = 82.521;

/***************************************************
 *  Name:        convertRGBtoXYZ
 *
 *  Returns:     Nothing
 *
 *  Parameters:  RGB values and XYZ references
 *
 *  Description: Given a color in RGB values converts  
 *               to the corresponding XYZ values
 *
 ***************************************************/
void convertRGBtoXYZ(int inR, int inG, int inB, double * outX, double * outY, double * outZ) {


	double var_R = (inR / 255.0f); //R from 0 to 255
	double var_G = (inG / 255.0f); //G from 0 to 255
	double var_B = (inB / 255.0f); //B from 0 to 255

	if (var_R > 0.04045f)
		var_R = pow(( (var_R + 0.055f) / 1.055f), 2.4f);
	else 
		var_R = var_R / 12.92f;

	if (var_G > 0.04045)
		var_G = pow(( (var_G + 0.055f) / 1.055f), 2.4f);
	else
		var_G = var_G / 12.92f;

	if (var_B > 0.04045f)
		var_B = pow(( (var_B + 0.055f) / 1.055f), 2.4f);
	else
		var_B = var_B / 12.92f;

	var_R = var_R * 100;
	var_G = var_G * 100;
	var_B = var_B * 100;

	//Observer. = 2°, Illuminant = D65
	*outX = var_R * 0.4124564 + var_G * 0.3575761 + var_B * 0.1804375;
	*outY = var_R * 0.2126729 + var_G * 0.7151522 + var_B * 0.0721750;
	*outZ = var_R * 0.0193339 + var_G * 0.1191920 + var_B * 0.9503041;
}


/***************************************************
 *  Name:        convertXYZtoLab
 *
 *  Returns:     Nothing
 *
 *  Parameters:  XYZ values and Lab references
 *
 *  Description: Given a color in XYZ values converts  
 *               to the corresponding Lab values
 *
 ***************************************************/
void convertXYZtoLab(double inX, double inY, double inZ, double * outL, double * outa, double * outb) {
	// See table above
	double ref_X = 95.047;
	double ref_Y = 100.0;
	double ref_Z = 108.883;

	double var_X = (inX / ref_X); //ref_X = 95.047
	double var_Y = (inY / ref_Y); //ref_Y = 100.0
	double var_Z = (inZ / ref_Z); //ref_Z = 108.883

	const double eps = pow(6.0 / 29.0, 3); // 0.008856
	const double m = 1.0 / 3.0 * pow(eps, -2); //7.787
	const double c = 4.0 / 29.0; // 16.0/116

	if ( var_X > eps ) 
		var_X = pow(var_X , 1.0/3.0);
	else 
		var_X = m * var_X + c;

	if ( var_Y > eps )
		var_Y = pow(var_Y , 1.0/3.0); 
	else
	    var_Y = m * var_Y  +  c;

	if ( var_Z > eps )
		var_Z = pow(var_Z , 1.0/3.0);
	else 
		var_Z = m * var_Z + c ;

	*outL = ( 116.0 * var_Y ) - 16.0;
	*outa = 500.0 * ( var_X - var_Y );
	*outb = 200.0 * ( var_Y - var_Z );
}


/***************************************************
 *  Name:        convertLabtoXYZ
 *
 *  Returns:     Nothing
 *
 *  Parameters:  Lab values and XYZ references
 *
 *  Description: Given a color in Lab values converts  
 *               to the corresponding XYZ values
 *
 ***************************************************/
void convertLabtoXYZ( double inL, double ina, double  inb, double * outX, double * outY, double * outZ) {
	// See table above
	double ref_X = 95.047;
	double ref_Y = 100.0;
	double ref_Z = 108.883;

	double var_Y = ( inL + 16.0 ) / 116.0;
	double var_X = ina / 500.0 + var_Y;
	double var_Z = var_Y - inb / 200.0;

	const double eps = pow(6.0 / 29.0, 3);
	const double m = 1.0 / 3.0 * pow(eps, -2);
	const double c = 4.0 / 29.0;

	if ( pow(var_Y,3) > eps ) 
		var_Y = pow(var_Y, 3);
	else
		var_Y = ( var_Y - c ) / m;

	if ( pow(var_X,3) > eps ) 
		var_X = pow(var_X, 3);
	else 
		var_X = ( var_X - c ) / m;
	
	if ( pow(var_Z,3) > eps )
		var_Z = pow(var_Z, 3);
	else
		var_Z = ( var_Z - c ) / m;

	*outX = ref_X * var_X;     //ref_X =  95.047     Observer= 2°, Illuminant= D65
	*outY = ref_Y * var_Y;     //ref_Y = 100.000
	*outZ = ref_Z * var_Z;     //ref_Z = 108.883
}

/***************************************************
 *  Name:        convertXYZtoRGB
 *
 *  Returns:     Nothing
 *
 *  Parameters:  XYZ values and RGB references
 *
 *  Description: Given a color in XYZ values converts  
 *               to the corresponding RGB values
 *
 ***************************************************/

// http://www.brucelindbloom.com/index.html?Eqn_RGB_XYZ_Matrix.html
// 	sRGB	D65	
// 	RGB to XYZ [M]
//  0.4124564  0.3575761  0.1804375
//  0.2126729  0.7151522  0.0721750
//  0.0193339  0.1191920  0.9503041
//  	XYZ to RGB [M]-1
//  3.2404542 -1.5371385 -0.4985314
// -0.9692660  1.8760108  0.0415560
//  0.0556434 -0.2040259  1.0572252

//  sRGB	D50	
//  RGB to XYZ [M]
//  0.4360747  0.3850649  0.1430804
//  0.2225045  0.7168786  0.0606169
//  0.0139322  0.0971045  0.7141733
//  XYZ to RGB [M]-1
//  3.1338561 -1.6168667 -0.4906146
// -0.9787684  1.9161415  0.0334540
//  0.0719453 -0.2289914  1.4052427

void convertXYZtoRGB(double inX, double inY, double inZ, int * outR, int * outG, int * outB) {


	double var_X = inX / 100.0;
	double var_Y = inY / 100.0;
	double var_Z = inZ / 100.0;

	double var_R = var_X *  3.2404542 + var_Y * -1.5371385 + var_Z * -0.4985314;
	double var_G = var_X * -0.9692660 + var_Y *  1.8760108 + var_Z *  0.0415560;
	double var_B = var_X *  0.0556434 + var_Y * -0.2040259 + var_Z *  1.0572252;

	if ( var_R > 0.0031308 ) {
		var_R = 1.055 * pow(var_R, ( 1.0 / 2.4 ) )  - 0.055;
	} else {
		var_R = 12.92 * var_R;
	}

	if ( var_G > 0.0031308 ) 
		var_G = 1.055 * pow(var_G, ( 1.0 / 2.4 ) ) - 0.055;
	else 
		var_G = 12.92 * var_G;

	if ( var_B > 0.0031308 )
		var_B = 1.055 * pow(var_B, ( 1.0 / 2.4 ) ) - 0.055;
	else
		var_B = 12.92 * var_B;

	// Clip results to correct range
	*outR = (int) round(MIN(255.0, MAX(0.0, var_R * 255.0)));
	*outG = (int) round(MIN(255.0, MAX(0.0, var_G * 255.0)));
	*outB = (int) round(MIN(255.0, MAX(0.0, var_B * 255.0)));
}

/***************************************************
 *  Name:        Lab_color_difference_CIE76
 *
 *  Returns:     double
 *
 *  Parameters:  2 Lab color values
 *
 *  Description: Calculates and returns the difference 
 *				 between 2 Lab colors based on the CIE76 formula
 *
 ***************************************************/
double Lab_color_difference_CIE76( double inL1, double ina1, double  inb1, double inL2, double ina2, double  inb2){
	return( sqrt( pow(inL1 - inL2, 2.f) + pow(ina1 - ina2, 2.f) + pow(inb1 - inb2, 2.f) ) );//CIE76
}

/***************************************************
 *  Name:        RGB_color_Lab_difference_CIE76
 *
 *  Returns:     double
 *
 *  Parameters:  2 RGB color values
 *
 *  Description: Converts RGB values to Lab values and returns  
 *				 the difference between 2 Lab colors based on 
 *				 the CIE76 formula
 *
 ***************************************************/
double RGB_color_Lab_difference_CIE76( int R1, int G1, int B1, int R2, int G2, int B2){
	double x1=0,y1=0,z1=0;
	double x2=0,y2=0,z2=0;
	double l1=0,a1=0,b1=0;
	double l2=0,a2=0,b2=0;

	convertRGBtoXYZ(R1, G1, B1, &x1, &x1, &z1);
	convertRGBtoXYZ(R2, G2, B2, &x2, &x2, &z2);

	convertXYZtoLab(x1, y1, z1, &l1, &a1, &b1);
	convertXYZtoLab(x2, y2, z2, &l2, &a2, &b2); 

	return( Lab_color_difference_CIE76(l1 ,a1 ,b1 ,l2 ,a2 ,b2) );
}

/***************************************************
 *  Name:        Lab_color_difference_CIE94
 *
 *  Returns:     double
 *
 *  Parameters:  2 Lab color values
 *
 *  Description: Calculates and returns the difference 
 *				 between 2 Lab colors based on the CIE94 formula
 *
 ***************************************************/
double Lab_color_difference_CIE94( double inL1, double ina1, double  inb1, double inL2, double ina2, double  inb2){
	// case Application.GraphicArts:
		double Kl = 1.0;
		double K1 = 0.045;
		double K2 = 0.015;
	// 	break;
	// case Application.Textiles:
	// 	Kl = 2.0;
	// 	K1 = .048;
	// 	K2 = .014;
	// break;

	double deltaL = inL1 - inL2;
	double deltaA = ina1 - ina2;
	double deltaB = inb1 - inb2;

	double c1 = sqrt(pow(ina1, 2) + pow(inb1, 2));
	double c2 = sqrt(pow(ina2, 2) + pow(inb2, 2));
	double deltaC = c1 - c2;

	double deltaH = pow(deltaA,2) + pow(deltaB,2) - pow(deltaC,2);
	deltaH = deltaH < 0 ? 0 : sqrt(deltaH);

	const double sl = 1.f;
	const double kc = 1.f;
	const double kh = 1.f;

	double sc = 1.f + K1*c1;
	double sh = 1.f + K2*c1;

	double i = pow(deltaL/(Kl*sl), 2) +
	                pow(deltaC/(kc*sc), 2) +
	                pow(deltaH/(kh*sh), 2);

	double finalResult = i < 0 ? 0 : sqrt(i);
	return (finalResult);
}

#define pi 3.141592653589793238462643383279

/// Computes the CIEDE2000 color-difference between two Lab colors
/// Based on the article:
/// The CIEDE2000 Color-Difference Formula: Implementation Notes,
/// Supplementary Test Data, and Mathematical Observations,", G. Sharma,
/// W. Wu, E. N. Dalal, submitted to Color Research and Application,
/// January 2004.
/// Available at http://www.ece.rochester.edu/~/gsharma/ciede2000/
/// Based on the C++ implementation by Ofir Pele, The Hebrew University of Jerusalem 2010.
//
double Lab_color_difference_CIEDE2000( double inL1, double ina1, double  inb1, double inL2, double ina2, double  inb2)
{
    double Cainb1= sqrt(ina1*ina1+inb1*inb1);
    double Cainb2= sqrt(ina2*ina2+inb2*inb2);

    double Cabarithmean= (Cainb1 + Cainb2)/2.0;

    double G= 0.5*( 1.0 - sqrt( pow(Cabarithmean,7.0)/(pow(Cabarithmean,7.0) + pow(25.0,7.0))));

    double apstd= (1.0+G)*ina1; // aprime in paper
    double apsample= (1.0+G)*ina2; // aprime in paper
    double Cpsample= sqrt(apsample*apsample+inb2*inb2);

    double Cpstd= sqrt(apstd*apstd+inb1*inb1);
    // Compute product of chromas
    double Cpprod= (Cpsample*Cpstd);


    // Ensure hue is between 0 and 2pi
    double hpstd= atan2(inb1,apstd);
    if (hpstd<0) hpstd+= 2.0*pi;  // rollover ones that come -ve

    double hpsample= atan2(inb2,apsample);
    if (hpsample<0) hpsample+= 2.0*pi;
    if ( (fabs(apsample)+fabs(inb2))==0.0)  hpsample= 0.0;

    double dL= (inL2-inL1);
    double dC= (Cpsample-Cpstd);

    // Computation of hue difference
    double dhp= (hpsample-hpstd);
    if (dhp>pi)  dhp-= 2.0*pi;
    if (dhp<-pi) dhp+= 2.0*pi;
    // set chroma difference to zero if the product of chromas is zero
    if (Cpprod == 0.0) dhp= 0.0;

    // Note that the defining equations actually need
    // signed Hue and chroma differences which is different
    // from prior color difference formulae

    double dH= 2.0*sqrt(Cpprod)*sin(dhp/2.0);
    //%dH2 = 4*Cpprod.*(sin(dhp/2)).^2;

    // weighting functions
    double Lp= (inL2+inL1)/2.0;
    double Cp= (Cpstd+Cpsample)/2.0;

    // Average Hue Computation
    // This is equivalent to that in the paper but simpler programmatically.
    // Note average hue is computed in radians and converted to degrees only
    // where needed
    double hp= (hpstd+hpsample)/2.0;
    // Identify positions for which abs hue diff exceeds 180 degrees
    if ( fabs(hpstd-hpsample)  > pi ) hp-= pi;
    // rollover ones that come -ve
    if (hp<0) hp+= 2.0*pi;

    // Check if one of the chroma values is zero, in which case set
    // mean hue to the sum which is equivalent to other value
    if (Cpprod==0.0) hp= hpsample+hpstd;

    double Lpm502= (Lp-50.0)*(Lp-50.0);;
    double Sl= 1.0+0.015*Lpm502/sqrt(20.0+Lpm502);
    double Sc= 1.0+0.045*Cp;
    double T= 1.0 - 0.17*cos(hp - pi/6.0) + 0.24*cos(2.0*hp) + 0.32*cos(3.0*hp+pi/30.0) - 0.20*cos(4.0*hp-63.0*pi/180.0);
    double Sh= 1.0 + 0.015*Cp*T;
    double delthetarad= (30.0*pi/180.0)*exp(- pow(( (180.0/pi*hp-275.0)/25.0),2.0));
    double Rc=  2.0*sqrt(pow(Cp,7.0)/(pow(Cp,7.0) + pow(25.0,7.0)));
    double RT= -sin(2.0*delthetarad)*Rc;

    // The CIE 00 color difference
    return sqrt( pow((dL/Sl),2.0) + pow((dC/Sc),2.0) + pow((dH/Sh),2.0) + RT*(dC/Sc)*(dH/Sh) );
}

/***************************************************
 *  Name:        RGB_color_Lab_difference_CIE94
 *
 *  Returns:     double
 *
 *  Parameters:  2 RGB color values
 *
 *  Description: Converts RGB values to Lab values and returns  
 *				 the difference between 2 Lab colors based on 
 *				 the CIE94 formula
 *
 ***************************************************/
double RGB_color_Lab_difference_CIE94( int R1, int G1, int B1, int R2, int G2, int B2){
	double x1=0,y1=0,z1=0;
	double x2=0,y2=0,z2=0;
	double l1=0,a1=0,b1=0;
	double l2=0,a2=0,b2=0;

	convertRGBtoXYZ(R1, G1, B1, &x1, &y1, &z1);
	convertRGBtoXYZ(R2, G2, B2, &x2, &y2, &z2);

	convertXYZtoLab(x1, y1, z1, &l1, &a1, &b1);
	convertXYZtoLab(x2, y2, z2, &l2, &a2, &b2); 

	return( Lab_color_difference_CIE94(l1 ,a1 ,b1 ,l2 ,a2 ,b2) );
}

double RGB_color_Lab_difference_CIEDE2000( int R1, int G1, int B1, int R2, int G2, int B2){
	double x1=0,y1=0,z1=0;
	double x2=0,y2=0,z2=0;
	double l1=0,a1=0,b1=0;
	double l2=0,a2=0,b2=0;

	convertRGBtoXYZ(R1, G1, B1, &x1, &y1, &z1);
	convertRGBtoXYZ(R2, G2, B2, &x2, &y2, &z2);

	convertXYZtoLab(x1, y1, z1, &l1, &a1, &b1);
	convertXYZtoLab(x2, y2, z2, &l2, &a2, &b2); 

	return( Lab_color_difference_CIEDE2000(l1 ,a1 ,b1 ,l2 ,a2 ,b2) );
}

/***************************************************
 *  Name:        main
 *
 *  Returns:     Nothing
 *
 *  Parameters:  argc, **argv
 *
 *  Description: Running some tests
 *
 ***************************************************/

//X11 colors
// #define DISTINCT_COLORS (140)
// char *distcolors[DISTINCT_COLORS] = {
// 	"7CFC00", "808000", "BDB76B", "FDF5E6", "BA55D3", "E0FFFF", "FFEFD5", "FFFFFF", "FAFAD2", "B22222", "90EE90", "A52A2A", "6495ED", "0000CD", "F4A460", "008080", "6A5ACD", "F0F8FF", "FA8072", "8B0000", "CD5C5C", "B8860B", "20B2AA", "FFFFF0", "7B68EE", "F8F8FF", "CD853F", "7FFFD4", "FFFAFA", "9400D3", "708090", "008000", "F5DEB3", "800080", "800000", "F0FFFF", "F5FFFA", "778899", "B0E0E6", "FFD700", "8FBC8F", "E9967A", "FAEBD7", "87CEEB", "FFDEAD", "8B4513", "7FFF00", "B0C4DE", "000080", "FFE4E1", "EEE8AA", "ADD8E6", "FFA500", "483D8B", "2E8B57", "FF4500", "ADFF2F", "A9A9A9", "808080", "FF0000", "0000FF", "DAA520", "FF69B4", "D8BFD8", "8A2BE2", "00FFFF", "FF8C00", "00CED1", "00008B", "008B8B", "4682B4", "FFFF00", "D3D3D3", "98FB98", "DB7093", "F0FFF0", "BC8F8F", "8B008B", "DEB887", "F0E68C", "D2B48C", "00FF00", "FFF0F5", "00FFFF", "191970", "9932CC", "1E90FF", "FF1493", "FF00FF", "C0C0C0", "5F9EA0", "A0522D", "696969", "3CB371", "DA70D6", "00FA9A", "AFEEEE", "FF00FF", "6B8E23", "FF7F50", "FFC0CB", "9ACD32", "00FF7F", "FFFACD", "C71585", "87CEFA", "FFB6C1", "40E0D0", "48D1CC", "006400", "FFA07A", "FFFAF0", "66CDAA", "00BFFF", "4B0082", "D2691E", "000000", "556B2F", "FFE4C4", "FFF5EE", "FFDAB9", "DCDCDC", "2F4F4F", "9370DB", "FF6347", "F5F5F5", "FAF0E6", "F08080", "DC143C", "32CD32", "EE82EE", "FFFFE0", "E6E6FA", "F5F5DC", "FFEBCD", "FFF8DC", "228B22", "FFE4B5", "DDA0DD", "4169E1"
// };

// int main(int argc, char const *argv[])
// {
// 	int R1,G1,B1,R2,G2,B2;

// 	R1 = 200;
// 	G1 = 2;
// 	B1 = 50; 

// 	R2 = 200;
// 	G2 = 2;
// 	B2 = 70; 

// 	// // convert hex to rgb
// 	// char *str = "0000FF";
// 	// int r, g, b;
// 	// sscanf(str, "%02x%02x%02x", &r, &g, &b);
// 	// // or
// 	// struct RGB colorConverter(int hexValue)
// 	// {
// 	//   struct RGB rgbColor;
// 	//   rgbColor.r = ((hexValue >> 16) & 0xFF) / 255.0;  // Extract the RR byte
// 	//   rgbColor.g = ((hexValue >> 8) & 0xFF) / 255.0;   // Extract the GG byte
// 	//   rgbColor.b = ((hexValue) & 0xFF) / 255.0;        // Extract the BB byte

// 	//   return rgbColor; 
// 	// }

// 	printf("LAB DISTANCE CIE94= %lf \n", RGB_color_Lab_difference_CIE94(R1,G1,B1,R2,G2,B2) );
// 	printf("LAB DISTANCE CIEDE2000= %lf \n", RGB_color_Lab_difference_CIEDE2000(R1,G1,B1,R2,G2,B2) );
	
// 	double X, Y, Z;
// 	double l, a, b;
// 	printf("rgb(%i, %i, %i)\n", R1, G1, B1);
// 	convertRGBtoXYZ(R1, G1, B1, &X, &Y, &Z);
// 	printf("xyz(%f, %f, %f)\n", X, Y, Z);
// 	convertXYZtoLab(X, Y, Z, &l, &a, &b);
// 	printf("lab(%f, %f, %f)\n", l, a, b);

// 	convertLabtoXYZ(l, a, b, &X, &Y, &Z);
// 	convertXYZtoRGB(X, Y, Z, &R1, &G1, &B1);
// 	printf("xyz(%f, %f, %f)\n", X, Y, Z);
// 	printf("rgb(%i, %i, %i)\n", R1, G1, B1);
// 	printf("\n");

// 	l = 50.0;
// 	a = 2.6772;
// 	b = -79.7751;
// 	convertLabtoXYZ(l, a, b, &X, &Y, &Z);
// 	convertXYZtoRGB(X, Y, Z, &R1, &G1, &B1);
// 	printf("lab(%f, %f, %f)\n", l, a, b);
// 	printf("xyz(%f, %f, %f)\n", X, Y, Z);
// 	printf("rgb(%i, %i, %i)\n", R1, G1, B1);

// 	printf("rgb(%i, %i, %i)\n", R1, G1, B1);
// 	convertRGBtoXYZ(R1, G1, B1, &X, &Y, &Z);
// 	printf("xyz(%f, %f, %f)\n", X, Y, Z);
// 	convertXYZtoLab(X, Y, Z, &l, &a, &b);
// 	printf("lab(%f, %f, %f)\n", l, a, b);
// 	printf("\n");

// 	l = 50.0;
// 	a = 0.0;
// 	b = -82.7485;
// 	convertLabtoXYZ(l, a, b, &X, &Y, &Z);
// 	convertXYZtoRGB(X, Y, Z, &R2, &G2, &B2);
// 	printf("lab(%f, %f, %f)\n", l, a, b);
// 	printf("xyz(%f, %f, %f)\n", X, Y, Z);
// 	printf("rgb(%i, %i, %i)\n", R2, G2, B2);
// 	printf("\n");

// 	printf("LAB DISTANCE CIEDE2000= %lf \n", Lab_color_difference_CIEDE2000(50.0, 2.6772, -79.7751, 50.0, 0.0, -82.7485) );
// 	printf("LAB DISTANCE CIEDE2000= %lf \n", RGB_color_Lab_difference_CIEDE2000(R1,G1,B1,R2,G2,B2) );
// 	// Filter X11 colors on lightness > 5 < 90
// 	// double dark = 10.0;
// 	// double light = 95.0;
// 	// for (int i = 0; i < DISTINCT_COLORS; i++) {
// 	// 	int r, g, b;
// 	// 	double x1=0,y1=0,z1=0;
// 	// 	double l1=0,a1=0,b1=0;

// 	// 	sscanf(distcolors[i], "%02x%02x%02x", &r, &g, &b);

// 	// 	convertRGBtoXYZ(r, g, b, &x1, &y1, &z1);
// 	// 	convertXYZtoLab(x1, y1, z1, &l1, &a1, &b1); 

// 	// 	if ((l1 > dark) && ( l1 < light)) {
// 	// 		printf("\"#%s\"\n", distcolors[i]);
// 	// 	} else {
// 	// 		if (l1 > light) {
// 	// 			printf("Too light #%s (l:%lf) \n", distcolors[i], l1);
// 	// 		} else {
// 	// 			printf("Too dark #%s (l:%lf) \n", distcolors[i], l1);
// 	// 		}
// 	// 	}
// 	// }
// 	return 0;
// }

#endif
#endif
