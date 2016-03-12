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
// minor exes
type Ellipsoid struct {
	a, b float64
}

type InvalidLatitudeError struct {
	Lat float64
}

func (e InvalidLatitudeError) Error() string {
	return fmt.Sprintf("invalid latitude: %f", e.Lat)
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

// ToNVector returns a cartesian position vector.
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
	coeff := ellps.b / math.Sqrt(nv.Vec3[0]*nv.Vec3[0]+
		absq*nv.Vec3[1]*nv.Vec3[1]+
		absq*nv.Vec3[2]*nv.Vec3[2])
	return PVector{Vec3{coeff * nv.Vec3[0], absq * nv.Vec3[1], absq * nv.Vec3[2]}}
}

// ToNVector returns a Cartesian position vector, given an ellipsoid.
func (pv *PVector) ToNVector(ellps *Ellipsoid) NVector {
	eccen := math.Sqrt(1 - ellps.b*ellps.b/(ellps.a*ellps.a))
	eccen2 := eccen * eccen
	eccen4 := eccen2 * eccen2
	a2 := ellps.a * ellps.a
	q := (1 - eccen2) / a2 * pv.Vec3[0] * pv.Vec3[0]
	p := (pv.Vec3[1]*pv.Vec3[1] + pv.Vec3[2]*pv.Vec3[2]) / a2
	r := (p + q - eccen4) / 6.0
	s := eccen4 * p * q / (4 * math.Pow(r, 3))
	t := math.Cbrt(1 + s + math.Sqrt(s*(2+s)))
	u := r * (1 + t + 1.0/t)
	v := math.Sqrt(u*u + eccen4*q)
	w := 0.5 * eccen2 * (u + v - q) / v
	k := math.Sqrt(u+v+w*w) - w
	d := k * math.Sqrt(pv.Vec3[1]*pv.Vec3[1]+pv.Vec3[2]*pv.Vec3[2]) / (k + eccen2)
	coeff := 1.0 / math.Sqrt(d*d+pv.Vec3[0]*pv.Vec3[0])
	kke2 := k / (k + eccen2)
	return NVector{
		Vec3{coeff * pv.Vec3[0],
			coeff * kke2 * pv.Vec3[1],
			coeff * kke2 * pv.Vec3[2]}}
}

// RotationMatrix returns the 3x3 matrix relating the Earth-centered
// non-singular coordinate frame to the North-up singular coordinate frame.
func (nv *NVector) RotationMatrix() Matrix3 {
	k_east := cross(&Vec3{0, 0, -1}, &nv.Vec3)
	k_north := cross(&nv.Vec3, k_east)
	a := k_north[0] / k_north.Magnitude()
	b := k_east[0] / k_east.Magnitude()
	c := nv.Vec3[0]
	d := k_north[1] / k_north.Magnitude()
	e := k_east[1] / k_east.Magnitude()
	f := nv.Vec3[1]
	g := k_north[2] / k_north.Magnitude()
	h := k_east[2] / k_east.Magnitude()
	i := nv.Vec3[2]
	return Matrix3{[3]float64{a, b, c}, [3]float64{d, e, f}, [3]float64{g, h, i}}
}

func (nv *NVector) SphericalDistance(nv2 *NVector, R float64) float64 {
	s_ab := math.Atan2(cross(&nv.Vec3, &nv2.Vec3).Magnitude(),
		dot(&nv.Vec3, &nv2.Vec3)) * R
	return s_ab
}

// Azimuth returns the azimuth and back azimuth from one NVector to another
// along an ellipse
func (nv *NVector) Azimuth(nv2 *NVector, ellps *Ellipsoid) (float64, float64, error) {
	return 0.0, 0.0, nil
}

// Forward returns the NVector position arrived at by moving in an azimuthal
// direction for a given distance along an ellipse
func (nv *NVector) Forward(az, distance float64, ellps *Ellipsoid) NVector {
	//north := nv.North()
	//east := nv.East()
	east := cross(&nv.Vec3, &Vec3{0, 0, -1})
	north := cross(&nv.Vec3, east)

	cos_az := math.Cos(az)
	sin_az := math.Sin(az)
	vec_az := Vec3{north[0]*cos_az + east[0]*sin_az,
		north[1]*cos_az + east[1]*sin_az,
		north[2]*cos_az + east[2]*sin_az}

	// Great circle angle travelled
	sab := distance / ellps.a
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
