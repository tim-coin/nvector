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
	lon float64
	lat float64
}

type Ellipsoid struct {
	a, b float64
}

type InvalidLatitudeError struct {
	lat float64
}

func (e InvalidLatitudeError) Error() string {
	return fmt.Sprintf("invalid latitude: %f", e.lat)
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

func (ll *LonLat) ToNVector() NVector {
	x := math.Sin(ll.lat)
	y := math.Sin(ll.lon) * math.Cos(ll.lat)
	z := -math.Cos(ll.lon) * math.Cos(ll.lat)
	return NVector{Vec3{x, y, z}}
}

func (ll *LonLat) String() string {
	londeg := ll.lon * 180.0 / math.Pi
	latdeg := ll.lat * 180.0 / math.Pi
	return fmt.Sprintf("(%.6f, %.6f)", londeg, latdeg)
}

func (nv *NVector) ToLonLat() LonLat {
	lat := math.Atan2(nv.Vec3[0], math.Sqrt(nv.Vec3[1]*nv.Vec3[1]+nv.Vec3[2]*nv.Vec3[2]))
	lon := math.Atan2(nv.Vec3[1], -nv.Vec3[2])
	lon = math.Mod((lon+0.5*math.Pi), math.Pi) - 0.5*math.Pi
	return LonLat{lon, lat}
}

func (nv *NVector) ToPVector(ellps *Ellipsoid) PVector {
	absq := ellps.a * ellps.a / (ellps.b * ellps.b)
	coeff := ellps.b / math.Sqrt(nv.Vec3[0]*nv.Vec3[0]+
		absq*nv.Vec3[1]*nv.Vec3[1]+
		absq*nv.Vec3[2]*nv.Vec3[2])
	return PVector{Vec3{coeff * nv.Vec3[0], absq * nv.Vec3[1], absq * nv.Vec3[2]}}
}

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

func (nv *NVector) RotationMatrix() Matrix3 {
	k_east := cross(&Vec3{1, 0, 0}, &nv.Vec3)
	k_north := cross(cross(&nv.Vec3, &Vec3{1, 0, 0}), &nv.Vec3)
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
