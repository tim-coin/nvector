/*
	Implements the "n-vector" based geodetic operations from

		Gade, Kenneth, "A Non-Singular Horizontal Position Representation",
		The Journal of Navigation (2010), 63, 395-417.
		doi:10.1017/S0373463309990415

*/
package nvector

import (
	"fmt"
	"math"
)

type ILonLat interface {
	ToNVector() INVector
}

type INVector interface {
	ToLonLat() ILonLat
	Magnitude() float64
}

type IPVector interface {
	ToNVector() INVector
	Magnitude() float64
}

type Vec3 [3]float64

type Matrix3 [3][3]float64

type NVector struct {
	Vec3
}

type PVector struct {
	Vec3
}

type LonLat struct {
	Lon float64
	Lat float64
}

// Ellipsoid represents a geographical ellipsoid in terms of its major and
// minor axes
type Ellipsoid struct {
	a, b float64
}

type InvalidLatitudeError struct {
	Lat float64
}

func (e InvalidLatitudeError) Error() string {
	return fmt.Sprintf("invalid latitude: %f", e.Lat)
}

type NoIntersectionError struct {
}

func (e NoIntersectionError) Error() string {
	return fmt.Sprintf("no intersection")
}

func cross(u, v *Vec3) *Vec3 {

	return &Vec3{u[1]*v[2] - u[2]*v[1], u[2]*v[0] - u[0]*v[2], u[0]*v[1] - u[1]*v[0]}
}

func dot(u, v *Vec3) float64 {
	return u[0]*v[0] + u[1]*v[1] + u[2]*v[2]
}

func (v *Vec3) Magnitude() float64 {
	return math.Sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2])
}

func (m *Matrix3) Mult(v *Vec3) Vec3 {
	var p Vec3
	p[0] = v[0]*m[0][0] + v[1]*m[0][1] + v[2]*m[0][2]
	p[1] = v[0]*m[1][0] + v[1]*m[1][1] + v[2]*m[1][2]
	p[2] = v[0]*m[2][0] + v[1]*m[2][1] + v[2]*m[2][2]
	return p
}

func (m *Matrix3) Transpose() Matrix3 {
	var tr Matrix3
	tr[0] = [3]float64{m[0][0], m[1][0], m[2][0]}
	tr[1] = [3]float64{m[0][1], m[1][1], m[2][1]}
	tr[2] = [3]float64{m[0][2], m[1][2], m[2][2]}
	return tr
}

func NewLonLat(londeg float64, latdeg float64) (*LonLat, error) {
	lon := londeg * math.Pi / 180.0
	lat := latdeg * math.Pi / 180.0
	lon = math.Mod((lon+math.Pi), 2*math.Pi) - math.Pi
	if (lat < -0.5*math.Pi) || (lat > 0.5*math.Pi) {
		lonlat := new(LonLat)
		return lonlat, InvalidLatitudeError{latdeg}
	}
	return &LonLat{lon, lat}, nil
}

// ToNVector returns a Cartesian position vector.
func (ll *LonLat) ToNVector() NVector {
	z := math.Sin(ll.Lat)
	y := math.Sin(ll.Lon) * math.Cos(ll.Lat)
	x := math.Cos(ll.Lon) * math.Cos(ll.Lat)
	return NVector{Vec3{x, y, z}}
}

func (ll *LonLat) String() string {
	londeg := ll.Lon * 180.0 / math.Pi
	latdeg := ll.Lat * 180.0 / math.Pi
	return fmt.Sprintf("(%.6f, %.6f)", londeg, latdeg)
}

// ToLonLat returns a LonLat struct, where lon: [-pi, pi) and lat: [-pi/2, pi/2].
func (nv *NVector) ToLonLat() LonLat {
	lat := math.Atan2(nv.Vec3[2], math.Sqrt(nv.Vec3[0]*nv.Vec3[0]+nv.Vec3[1]*nv.Vec3[1]))
	lon := math.Atan2(nv.Vec3[1], nv.Vec3[0])
	lon = math.Mod((lon+0.5*math.Pi), math.Pi) - 0.5*math.Pi
	return LonLat{lon, lat}
}

// ToPVector returns a surface-normal vector, given an ellipsoid.
func (nv *NVector) ToPVector(ellps *Ellipsoid) PVector {
	absq := ellps.a * ellps.a / (ellps.b * ellps.b)
	coeff := ellps.b / math.Sqrt(nv.Vec3[2]*nv.Vec3[2]+
		absq*nv.Vec3[1]*nv.Vec3[1]+
		absq*nv.Vec3[0]*nv.Vec3[0])
	return PVector{Vec3{coeff * absq * nv.Vec3[0], coeff * absq * nv.Vec3[1], coeff * nv.Vec3[2]}}
}

// ToNVector returns a Cartesian position vector, given an ellipsoid.
func (pv *PVector) ToNVector(ellps *Ellipsoid) NVector {
	eccen := math.Sqrt(1 - ellps.b*ellps.b/(ellps.a*ellps.a))
	eccen2 := eccen * eccen
	eccen4 := eccen2 * eccen2
	a2 := ellps.a * ellps.a
	q := (1 - eccen2) / a2 * pv.Vec3[2] * pv.Vec3[2]
	p := (pv.Vec3[1]*pv.Vec3[1] + pv.Vec3[0]*pv.Vec3[0]) / a2
	r := (p + q - eccen4) / 6.0
	s := eccen4 * p * q / (4 * math.Pow(r, 3))
	t := math.Cbrt(1 + s + math.Sqrt(s*(2+s)))
	u := r * (1 + t + 1.0/t)
	v := math.Sqrt(u*u + eccen4*q)
	w := 0.5 * eccen2 * (u + v - q) / v
	k := math.Sqrt(u+v+w*w) - w
	d := k * math.Sqrt(pv.Vec3[1]*pv.Vec3[1]+pv.Vec3[0]*pv.Vec3[0]) / (k + eccen2)
	coeff := 1.0 / math.Sqrt(d*d+pv.Vec3[2]*pv.Vec3[2])
	kke2 := k / (k + eccen2)
	return NVector{
		Vec3{-coeff * kke2 * pv.Vec3[0],
			-coeff * kke2 * pv.Vec3[1],
			coeff * pv.Vec3[2]}}
}

// RotationMatrix returns the 3x3 matrix relating the Earth-centered
// non-singular coordinate frame to the North-East-Down singular coordinate
// frame.
func (nv *NVector) RotationMatrix() Matrix3 {
	east := cross(&Vec3{0, 0, 1}, &nv.Vec3)
	north := cross(&nv.Vec3, east)

	a := north[0] / north.Magnitude()
	b := east[0] / east.Magnitude()
	c := -nv.Vec3[0]
	d := north[1] / north.Magnitude()
	e := east[1] / east.Magnitude()
	f := -nv.Vec3[1]
	g := north[2] / north.Magnitude()
	h := east[2] / east.Magnitude()
	i := -nv.Vec3[2]
	return Matrix3{[3]float64{a, b, c}, [3]float64{d, e, f}, [3]float64{g, h, i}}
}

// SphericalDistance returns the distance from another NVector on a sphere with
// radius *R*
func (nv *NVector) SphericalDistance(nv2 *NVector, R float64) float64 {
	s_ab := math.Atan2(cross(&nv.Vec3, &nv2.Vec3).Magnitude(),
		dot(&nv.Vec3, &nv2.Vec3)) * R
	return s_ab
}

// Azimuth returns the azimuth and back azimuth from one NVector to another
// along an ellipse
func (nv *NVector) Azimuth(nv2 *NVector, ellps *Ellipsoid) float64 {
	pv1 := nv.ToPVector(ellps)
	pv2 := nv2.ToPVector(ellps)
	delta_E := Vec3{pv2.Vec3[0] - pv1.Vec3[0], pv2.Vec3[1] - pv1.Vec3[1], pv2.Vec3[2] - pv1.Vec3[2]}

	rotMat_EN := nv.RotationMatrix()
	rotMat_NE := rotMat_EN.Transpose()
	delta_N := rotMat_NE.Mult(&delta_E)
	return math.Atan2(delta_N[1], delta_N[0])
}

// Forward returns the NVector position arrived at by moving in an azimuthal
// direction for a given distance along an ellipse
func (nv *NVector) Forward(az, distance, radius float64) NVector {
	east := cross(&nv.Vec3, &Vec3{0, 0, -1})
	north := cross(&nv.Vec3, east)

	cos_az := math.Cos(az)
	sin_az := math.Sin(az)
	vec_az := Vec3{north[0]*cos_az + east[0]*sin_az,
		north[1]*cos_az + east[1]*sin_az,
		north[2]*cos_az + east[2]*sin_az}

	// Great circle angle travelled
	sab := distance / radius
	cos_sab := math.Cos(sab)
	sin_sab := math.Sin(sab)
	resultant := Vec3{nv.Vec3[0]*cos_sab + vec_az[0]*sin_sab,
		nv.Vec3[1]*cos_sab + vec_az[1]*sin_sab,
		nv.Vec3[2]*cos_sab + vec_az[2]*sin_sab}
	return NVector{resultant}
}

func interpLinear(x, x0, x1, y0, y1 float64) float64 {
	return (x-x0)/(x1-x0)*(y1-y0) + y0
}

// Interpolate returns the NVector representing the intermediate position
// between two other NVectors. *frac* is the fractional distance between *nv*
// and *nv2*.
func (nv *NVector) Interpolate(nv2 *NVector, frac float64) NVector {
	result := new(NVector)
	result.Vec3[0] = interpLinear(frac, 0, 1, nv.Vec3[0], nv2.Vec3[0])
	result.Vec3[1] = interpLinear(frac, 0, 1, nv.Vec3[1], nv2.Vec3[1])
	result.Vec3[2] = interpLinear(frac, 0, 1, nv.Vec3[2], nv2.Vec3[2])
	return *result
}

// Intersection returns the spheroidal intersection point between two geodesics
// defined by an NVector pair, if it exists. If no intersection exists,
// NoIntersectionError is returned
func Intersection(nv1a, nv1b, nv2a, nv2b *NVector) (NVector, error) {
	var normalA, normalB, intersection *Vec3
	var err error

	normalA = cross(&nv1a.Vec3, &nv1b.Vec3)
	normalB = cross(&nv2a.Vec3, &nv2b.Vec3)
	intersection = cross(normalA, normalB)

	// Select the intersection on the right side of the spheroid
	if dot(intersection, &nv1a.Vec3) < 0 {
		intersection[0] = -intersection[0]
		intersection[1] = -intersection[1]
		intersection[2] = -intersection[2]
	}

	result := NVector{*intersection}

	// Tests whether intersection is between segment endpoints to within ~4cm
	var dab, dai, dbi float64
	dab = nv1a.SphericalDistance(nv1b, 1.0)
	dai = nv1a.SphericalDistance(&result, 1.0)
	dbi = nv1b.SphericalDistance(&result, 1.0)

	if math.Abs(dab-dai-dbi) > 1e-9 {
		err = NoIntersectionError{}
	}

	dab = nv2a.SphericalDistance(nv2b, 1.0)
	dai = nv2a.SphericalDistance(&result, 1.0)
	dbi = nv2b.SphericalDistance(&result, 1.0)

	if math.Abs(dab-dai-dbi) > 1e-9 {
		err = NoIntersectionError{}
	}

	return NVector{*intersection}, err
}




// Intersection2 returns the spheroidal intersection point between two geodesics based on longitude range check.
//If both intersections are within the longitude, it is OK. Else...
//If one in range is shortest distance to the point, it is OK
// defined by an NVector pair, if it exists. If no intersection exists,
// NoIntersectionError is returned
func Intersection2(nv1a, nv1b, nv2a, nv2b *NVector) (NVector, error) {
	//Add a delta if both points are same for Point of Line

	var normalA, normalB, intersection *Vec3
	//fmt.Println(nv2a.ToLonLat().Lon,nv2b.ToLonLat().Lon, -1*math.Pi, delta)
	if(nv2a.ToLonLat().Lon == nv2b.ToLonLat().Lon && nv2a.ToLonLat().Lon == -1*math.Pi){
		//Fixing singularity
		/*
		delta := 1e-9
		_t , _ := NewLonLat((nv2a.ToLonLat().Lon - delta)*180/math.Pi, (nv2a.ToLonLat().Lat)*180/math.Pi)
		_t1 := _t.ToNVector()
		nv2a = &_t1
		*/
		//fmt.Println("Needs Delta", nv2a)
		normalB = &nv2a.Vec3
	}else{
	normalB = cross(&nv2a.Vec3, &nv2b.Vec3)
	}
	//nv1a is the point
	//nv1b is the pole
	//nv2a and nv2b is the line with which intersection is sought

	var err error

	normalA = cross(&nv1a.Vec3, &nv1b.Vec3)

	intersection = cross(normalA, normalB)
	intersection2 := Vec3{0,0,0}  //negative(intersection)
	intersection2[0] = -1*intersection[0]
	intersection2[1] = -1*intersection[1]
	intersection2[2] = -1*intersection[2]

	in1 := NVector{*intersection}
	in2 := NVector{intersection2}
	/*
	fmt.Println("322: ", in1, in2)
	fmt.Println("326: ", intersection, intersection2)
	fmt.Println(":", intersection[0],1*intersection2[0])
	fmt.Println(":", intersection[1],intersection2[1])
	fmt.Println(":", intersection[2],intersection2[2])
	*/

	din1 :=  nv1a.SphericalDistance(&in1, 1.0) //Distance of intersection 1
	din2 :=  nv1a.SphericalDistance(&in2, 1.0) //Distance of intersection 2

	//loin := in1.ToLonLat().Lon //Let's assume that 1st intersection is nearest to POI (point of interest)
	//lain := in1.ToLonLat().Lat //Let's assume that 1st intersection is nearest to POI (point of interest)
	fmt.Println("LOI:::",nv2a.ToLonLat().Lon*180/math.Pi,nv2a.ToLonLat().Lat*180/math.Pi,"|", nv2b.ToLonLat().Lon*180/math.Pi, nv2b.ToLonLat().Lat*180/math.Pi)
	fmt.Println("Intersects:::",in1.ToLonLat().Lon*180/math.Pi,in1.ToLonLat().Lat*180/math.Pi,"|", in2.ToLonLat().Lon*180/math.Pi, in2.ToLonLat().Lat*180/math.Pi)
	result := in1
	if(din2 < din1){
		//loin = in2.ToLonLat().Lon
		//lain = in2.ToLonLat().Lat
		result = in2
	} //Now we have the nearest intersection point. Finally check if it is in range of POL(point of Line)
	fmt.Println("Dist: ",result, ":",din1, din2)


	/* //THis is not needed as replaced by LOI check
	lorange := []float64{math.Min(nv2a.ToLonLat().Lon,nv2b.ToLonLat().Lon), math.Max(nv2a.ToLonLat().Lon,nv2b.ToLonLat().Lon)} //the line of interest
	larange := []float64{math.Min(nv2a.ToLonLat().Lat,nv2b.ToLonLat().Lat), math.Max(nv2a.ToLonLat().Lat,nv2b.ToLonLat().Lat)} //the line of interest

	//Check if it doesnt exist in range, generate error
	//if( (loin > lorange[1] || loin < lorange[0] ) || (lain > larange[1] || lain < larange[0] ) ){
	fmt.Println( (math.Cos(loin) > math.Cos(lorange[1]) || math.Cos(loin) < math.Cos(lorange[0]) ) , (math.Cos(lain) > math.Cos(larange[1]) || math.Cos(lain) < math.Cos(larange[0]) ) )
	fmt.Println( math.Cos(loin) > math.Cos(lorange[1]) , math.Cos(loin) < math.Cos(lorange[0])  , math.Cos(lain) > math.Cos(larange[1]) , math.Cos(lain) < math.Cos(larange[0])  )
	fmt.Println( math.Cos(loin) ,">", math.Cos(lorange[1]) , math.Cos(loin) ,"<", math.Cos(lorange[0])  , math.Cos(lain) ,">", math.Cos(larange[1]) , math.Cos(lain) ,"<", math.Cos(larange[0])  )
	fmt.Println( (loin)*180/math.Pi ,">", (lorange[1])*180/math.Pi , (loin)*180/math.Pi ,"<", (lorange[0])*180/math.Pi ,"|" , (lain)*180/math.Pi ,">", (larange[1])*180/math.Pi , (lain)*180/math.Pi ,"<", (larange[0])*180/math.Pi  )

	if( (math.Cos(loin) > math.Cos(lorange[1]) || math.Cos(loin) < math.Cos(lorange[0]) ) || (math.Cos(lain) > math.Cos(larange[1]) || math.Cos(lain) < math.Cos(larange[0]) ) ){
		err = NoIntersectionError{}
		fmt.Println("T353: ", loin, lorange,";",lain, larange)
		fmt.Println("T354: ", math.Cos(loin), math.Cos(lorange[0]), math.Cos(lorange[1]),";",math.Cos(lain), math.Cos(larange[0]),math.Cos(larange[1]) )
	}else{
		fmt.Println("T356: ", math.Cos(loin), math.Cos(lorange[0]), math.Cos(lorange[1]),";",math.Cos(lain), math.Cos(larange[0]),math.Cos(larange[1]) )

		fmt.Println("T358:" , loin*180/math.Pi, " is between ", lorange[0]*180/math.Pi, " & ", lorange[1]*180/math.Pi)
	}*/


	// Tests whether intersection is between segment endpoints to within ~4cm
	var dab, dai, dbi float64
	dab = nv1a.SphericalDistance(nv1b, 1.0)
	dai = nv1a.SphericalDistance(&result, 1.0)
	dbi = nv1b.SphericalDistance(&result, 1.0)
	fmt.Println("T359: ", dab, dai,dbi, (dab-dai-dbi))
	if math.Abs(dab-dai-dbi) > 1e-9 {
		err = NoIntersectionError{}
	}

	//This is  needed, as  longitude test is not correct
	dab = nv2a.SphericalDistance(nv2b, 1.0)
	dai = nv2a.SphericalDistance(&result, 1.0)
	dbi = nv2b.SphericalDistance(&result, 1.0)
	fmt.Println("T367: ", dab, dai,dbi, (dab-dai-dbi))
	if math.Abs(dab-dai-dbi) > 1e-9 {
		err = NoIntersectionError{}
	}



	fmt.Println("Int Longitude is,", result.ToLonLat().Lon*180/math.Pi, err)
	return result, err
}


// Extrapolation returns the spheroidal  point where the line will intersect
// NoIntersectionError is returned
func Extrapolation(nv1a, nv1b, nv2a, nv2b *NVector) (LonLat, error) {
	//Add a delta if both points are same for Point of Line
	delta := 1e-9
	//fmt.Println("T384: Point-1 is, ", nv2a.Vec3, nv2b.Vec3[0],nv1a.ToLonLat().Lon*180/math.Pi, nv1a.ToLonLat().Lat*180/math.Pi)
	if(nv2a.ToLonLat().Lon == nv2b.ToLonLat().Lon && nv2a.ToLonLat().Lon == -1*math.Pi){
		//Fixing singularity
		_t , _ := NewLonLat((nv2a.ToLonLat().Lon - delta)*180/math.Pi, (nv2a.ToLonLat().Lat)*180/math.Pi)
		_t1 := _t.ToNVector()
		nv2a = &_t1
		//fmt.Println("Needs Delta", nv2a)
	}
	//nv1a is the point
	//nv1b is the pole
	//nv2a and nv2b is the line with which intersection is sought
	var normalA, normalB, intersection *Vec3
	var err error

	normalA = cross(&nv1a.Vec3, &nv1b.Vec3)
	normalB = cross(&nv2a.Vec3, &nv2b.Vec3)
	intersection = cross(normalA, normalB)
	intersection2 := Vec3{0,0,0}  //negative(intersection)
	intersection2[0] = -1*intersection[0]
	intersection2[1] = -1*intersection[1]
	intersection2[2] = -1*intersection[2]

	in1 := NVector{*intersection}
	in2 := NVector{intersection2}

	loin := in1.ToLonLat().Lon //Let's assume that 1st intersection is nearest to POI (point of interest)
	lain := in1.ToLonLat().Lat //Let's assume that 1st intersection is nearest to POI (point of interest)
	//fmt.Println("lonlat of int1:::",in1.ToLonLat().Lon*180/math.Pi, in1.ToLonLat().Lat*180/math.Pi)
	//fmt.Println("lonlat of int2:::",in2.ToLonLat().Lon*180/math.Pi, in2.ToLonLat().Lat*180/math.Pi)
	result := in1
	lorange := []float64{math.Min(nv2a.ToLonLat().Lon,nv2b.ToLonLat().Lon), math.Max(nv2a.ToLonLat().Lon,nv2b.ToLonLat().Lon)} //the line of interest
	larange := []float64{math.Min(nv2a.ToLonLat().Lat,nv2b.ToLonLat().Lat), math.Max(nv2a.ToLonLat().Lat,nv2b.ToLonLat().Lat)} //the line of interest

	//if( (math.Cos(loin) > math.Cos(lorange[1]) || math.Cos(loin) < math.Cos(lorange[0]) ) || (math.Cos(lain) > math.Cos(larange[1]) || math.Cos(lain) < math.Cos(larange[0]) ) ){
	if( ((loin) > (lorange[1]) || (loin) < (lorange[0]) ) || ((lain) > (larange[1]) || (lain) < (larange[0]) ) ){

		loin = in2.ToLonLat().Lon
		lain = in2.ToLonLat().Lat
		result = in2
	} //Now we have the nearest intersection point. Finally check if it is in range of POL(point of Line)

	//Pole Problem: If the LOI is pole to meridian, for hemisphere triangles, it will intersect at equator at edes.
	//Calculate midpoint of LOI, That's it
	var  dai, dbi float64
	//dab = nv2a.SphericalDistance(nv2b, 1.0)
	dai = nv2a.SphericalDistance(&result, 1.0)
	dbi = nv2b.SphericalDistance(&result, 1.0)
	//fmt.Println("T431: ", dai, dbi)
	if dai*dbi == 0  {
		//result = (nv2a+nv2b)/2
		result = NVector{Vec3{(nv2a.Vec3[0]+nv2b.Vec3[0])/2,(nv2a.Vec3[1]+nv2b.Vec3[1])/2,(nv2a.Vec3[2]+nv2b.Vec3[2])/2}}
	}

	//fmt.Println("Degree: Lo:", (loin*180/math.Pi) ,"; LoMax: ", (lorange[1]*180/math.Pi) ,"; LoMin: ", (lorange[0]*180/math.Pi) , ", La: ", (lain*180/math.Pi) ,"; LaMax: ", (larange[1]*180/math.Pi) ,"; LaMin: ", (larange[0]*180/math.Pi) )


	result2 := result.ToLonLat()
	//fmt.Println("Point  is,", nv1a,"; Lon: ", nv1a.ToLonLat().Lon*180/math.Pi,"; Lat: ", nv1a.ToLonLat().Lat*180/math.Pi)
	/*fmt.Println("Point-1 is, ",nv1a.ToLonLat().Lon*180/math.Pi, nv1a.ToLonLat().Lat*180/math.Pi)
	fmt.Println("CG Point  is, ",nv1b.ToLonLat().Lon*180/math.Pi, nv1b.ToLonLat().Lat*180/math.Pi)
	fmt.Println("Point-3  is, ",nv2a.ToLonLat().Lon*180/math.Pi, nv2a.ToLonLat().Lat*180/math.Pi)
	fmt.Println("Point-4  is, ",nv2b.ToLonLat().Lon*180/math.Pi, nv2b.ToLonLat().Lat*180/math.Pi)
	fmt.Println("Result  is, ",result.ToLonLat().Lon*180/math.Pi, result.ToLonLat().Lat*180/math.Pi)
	*/
	return result2, err
}


//Find the graph which is mergable
func Merger(nv1a, nv1b, nv2a, nv2b *NVector) ([]float64, error) {
	//var normalA, normalB, intersection *Vec3
	var err error

	var  dab, dai, dbi float64
	dab = nv1a.SphericalDistance(nv1b, 1.0)
	dai = nv1a.SphericalDistance(nv2b, 1.0)
	dbi = nv2b.SphericalDistance(nv1b, 1.0)
	pt := nv2a // The Edge from common points
	if math.Abs(dab-dai-dbi) > 1e-9 {

		dai = nv1a.SphericalDistance(nv2a, 1.0)
		dbi = nv2a.SphericalDistance(nv1b, 1.0)
		pt = nv2b
		if math.Abs(dab-dai-dbi) > 1e-9 {
			err = NoIntersectionError{}
		}
	}
	//tri := [][]float64{{pt.ToLonLat().Lat*180/math.Pi, pt.ToLonLat().Lon*180/math.Pi},{nv2a.ToLonLat().Lat*180/math.Pi,nv2a.ToLonLat().Lon*180/math.Pi},{nv2b.ToLonLat().Lat*180/math.Pi,nv2b.ToLonLat().Lon*180/math.Pi} }
	return []float64{pt.ToLonLat().Lat*180/math.Pi, pt.ToLonLat().Lon*180/math.Pi} , err
}
