package nvector

import (
	"fmt"
	"math"
	"testing"
)

func isclose(a, b float64, places int32) bool {
	var tol float64
	tol = math.Pow(10, -float64(places))
	res := math.Abs(a-b) < tol
	if !res {
		fmt.Printf("%f != %f (to %d places)\n", a, b, places)
	}
	return res
}

func TestNewLonLatSimple(t *testing.T) {
	lonlat, err := NewLonLat(-140.0, 49.25)
	if err != nil {
		t.Error()
	}

	if math.Abs(lonlat.Lat-0.8595746566072073) > 0.0000000001 {
		t.Fail()
	}
	if math.Abs(lonlat.Lon-(-2.443460952792061)) > 0.0000000001 {
		t.Fail()
	}
}

func TestNewLonLatUnwrap(t *testing.T) {
	lonlat, err := NewLonLat(220.0, 49.25)
	if err != nil {
		t.Error()
	}

	if math.Abs(lonlat.Lat-0.8595746566072073) > 0.0000000001 {
		t.Fail()
	}
	if math.Abs(lonlat.Lon-(-2.443460952792061)) > 0.0000000001 {
		t.Fail()
	}
}

func TestNewLonLatInvalid(t *testing.T) {
	_, err := NewLonLat(-140.0, 92.0)

	if err == nil {
		t.Fail()
	}
}

func TestLonLatToNVector1(t *testing.T) {
	ll, _ := NewLonLat(0.0, 0.0)
	nv := ll.ToNVector()
	if !isclose(nv.Vec3[0], 1.0, 8) {
		t.Fail()
	}
	if !isclose(nv.Vec3[1], 0.0, 8) {
		t.Fail()
	}
	if !isclose(nv.Vec3[2], 0.0, 8) {
		t.Fail()
	}
}

func TestLonLatToNVector2(t *testing.T) {
	ll, _ := NewLonLat(-45, 0.0)
	nv := ll.ToNVector()
	if !isclose(nv.Vec3[0], math.Sin(math.Pi/4), 8) {
		t.Fail()
	}
	if !isclose(nv.Vec3[1], -math.Sin(math.Pi/4), 8) {
		t.Fail()
	}
	if !isclose(nv.Vec3[2], 0.0, 8) {
		t.Fail()
	}
}

func TestLonLatToNVector3(t *testing.T) {
	ll, _ := NewLonLat(90.0, 45.0)
	nv := ll.ToNVector()
	if !isclose(nv.Vec3[0], 0.0, 8) {
		t.Fail()
	}
	if !isclose(nv.Vec3[1], math.Sin(math.Pi/4), 8) {
		t.Fail()
	}
	if !isclose(nv.Vec3[2], math.Sin(math.Pi/4), 8) {
		t.Fail()
	}
}

func TestNVectorToLonLat(t *testing.T) {
	nv := NVector{Vec3{math.Sin(math.Pi/4) * math.Sin(math.Pi/4),
		math.Sin(math.Pi/4) * math.Sin(math.Pi/4),
		-math.Sin(math.Pi / 4)}}
	ll := nv.ToLonLat()
	if !isclose(ll.Lon*180/math.Pi, 45, 8) {
		t.Fail()
	}
	if !isclose(ll.Lat*180/math.Pi, -45, 8) {
		t.Fail()
	}
}

func Test_cross(t *testing.T) {
	v1 := Vec3{1, 0, 0}
	v2 := Vec3{0, 1, 0}
	v3 := cross(&v1, &v2)
	if (*v3 != Vec3{0, 0, 1}) {
		t.Fail()
	}

	v4 := Vec3{0, 1, 0}
	v5 := Vec3{1, 0, 0}
	v6 := cross(&v4, &v5)
	if (*v6 != Vec3{0, 0, -1}) {
		t.Fail()
	}
}

func TestTranspose(t *testing.T) {
	var m, mt, mt_ans Matrix3
	m = Matrix3{[3]float64{3, 6, -4}, [3]float64{8, -2, -1}, [3]float64{1, 1, 4}}
	mt = m.Transpose()
	mt_ans = Matrix3{[3]float64{3, 8, 1}, [3]float64{6, -2, 1}, [3]float64{-4, -1, 4}}
	if mt != mt_ans {
		t.Fail()
	}
}

func TestMatMult(t *testing.T) {
	var m Matrix3
	var v, v2 Vec3
	v = Vec3{2, 1, 3}
	m = Matrix3{[3]float64{3, 6, -4}, [3]float64{8, -2, -1}, [3]float64{1, 1, 4}}
	v2 = m.Mult(&v)
	if (v2 != Vec3{0, 11, 15}) {
		t.Fail()
	}
}

func Test_dot(t *testing.T) {
	v1 := Vec3{1, 0, 0}
	v2 := Vec3{0, 1, 0}
	if dot(&v1, &v2) != 0 {
		t.Fail()
	}

	v3 := Vec3{1, 0.5, 0.25}
	v4 := Vec3{0, 1, 4}
	if !isclose(dot(&v3, &v4), 1.5, 8) {
		t.Fail()
	}
}

func TestSphericalDistance1(t *testing.T) {
	// meridional
	pos1, _ := NewLonLat(-140, 49.25)
	pos2, _ := NewLonLat(-140, 48.25)

	nv1 := pos1.ToNVector()
	nv2 := pos2.ToNVector()
	R := 6370997.0
	d := nv1.SphericalDistance(&nv2, R)

	if !isclose(d, 111194.874, 2) {
		t.Fail()
	}
}

func TestSphericalDistance2(t *testing.T) {
	// zonal
	pos1, _ := NewLonLat(-140, 49.25)
	pos2, _ := NewLonLat(-143, 49.25)

	nv1 := pos1.ToNVector()
	nv2 := pos2.ToNVector()
	R := 6370997.0
	d := nv1.SphericalDistance(&nv2, R)

	if !isclose(d, 217736.339, 2) {
		t.Fail()
	}
}

func TestSphericalDistance3(t *testing.T) {
	// crosses dateline
	pos1, _ := NewLonLat(174, -15)
	pos2, _ := NewLonLat(-177.5, 36)

	nv1 := pos1.ToNVector()
	nv2 := pos2.ToNVector()
	R := 6370997.0
	d := nv1.SphericalDistance(&nv2, R)

	if !isclose(d, 5740995.595, 2) {
		t.Fail()
	}
}

func TestAzimuth1(t *testing.T) {
	// meridional, along prime meridian
	pos1, _ := NewLonLat(0.0, 0.0)
	pos2, _ := NewLonLat(0.0, 10.0)
	nv1 := pos1.ToNVector()
	nv2 := pos2.ToNVector()
	ellps := Ellipsoid{6378137.0, 6356752.3142}
	az := nv1.Azimuth(&nv2, &ellps)
	if !isclose(az, 0, 6) {
		t.Fail()
	}
}

func TestAzimuth2(t *testing.T) {
	// meridional, off prime meridian
	pos1, _ := NewLonLat(60, 40)
	pos2, _ := NewLonLat(60, 50)
	nv1 := pos1.ToNVector()
	nv2 := pos2.ToNVector()
	ellps := Ellipsoid{6378137.0, 6356752.3142}
	az := nv1.Azimuth(&nv2, &ellps)
	if !isclose(az, 0, 6) {
		t.Fail()
	}
}

func TestAzimuth3(t *testing.T) {
	// zonal
	pos1, _ := NewLonLat(-140, 49.25)
	pos2, _ := NewLonLat(-143, 49.25)
	nv1 := pos1.ToNVector()
	nv2 := pos2.ToNVector()
	ellps := Ellipsoid{6378137.0, 6356752.3142}
	az := nv1.Azimuth(&nv2, &ellps)
	if !isclose(az, -88.8635416/180*math.Pi, 6) {
		t.Fail()
	}
}

func TestAzimuth4(t *testing.T) {
	// crosses dateline
	pos1, _ := NewLonLat(174, -15)
	pos2, _ := NewLonLat(-177.5, 36)
	nv1 := pos1.ToNVector()
	nv2 := pos2.ToNVector()

	var ellps Ellipsoid
	var az float64

	// WGS84 - note poor agreement with geographiclib result
	ellps = Ellipsoid{6378137.0, 6356752.3142}
	az = nv1.Azimuth(&nv2, &ellps)
	if !isclose(az, 8.8262727/180*math.Pi, 3) {
		t.Fail()
	}

	// sphere
	ellps = Ellipsoid{6370997.0, 6370997.0}
	az = nv1.Azimuth(&nv2, &ellps)
	if !isclose(az, 8.7731219/180*math.Pi, 6) {
		t.Fail()
	}
}

func TestForward1(t *testing.T) {
	pos, _ := NewLonLat(0, 0)
	nv := pos.ToNVector()
	ellps := Ellipsoid{6370997.0, 6370997.0}

	nv2 := nv.Forward(0, 100000, ellps.a)
	pos2 := nv2.ToLonLat()

	if !isclose(pos2.Lat*180/math.Pi, 0.8993220, 6) {
		t.Fail()
	}

	if !isclose(pos2.Lon*180/math.Pi, 0, 6) {
		t.Fail()
	}
}

func TestForward2(t *testing.T) {
	pos, _ := NewLonLat(0, 0)
	nv := pos.ToNVector()
	ellps := Ellipsoid{6370997.0, 6370997.0}

	nv2 := nv.Forward(-45/180.0*math.Pi, 1000000, ellps.a)
	pos2 := nv2.ToLonLat()

	if !isclose(pos2.Lat*180/math.Pi, 6.3460548, 6) {
		t.Fail()
	}

	if !isclose(pos2.Lon*180/math.Pi, -6.3853428, 6) {
		t.Fail()
	}
}

func Test_interpLinear(t *testing.T) {
	y1 := interpLinear(0.5, 0, 1, 0, 1)
	if !isclose(y1, 0.5, 8) {
		t.Fail()
	}
	y2 := interpLinear(0.5, 0, 2, 0, 1)

	if !isclose(y2, 0.25, 8) {
		t.Fail()
	}

	y3 := interpLinear(-0.5, 0, 1, 0, 2)
	if !isclose(y3, -1.0, 8) {
		t.Fail()
	}

	y4 := interpLinear(1.5, 0, 2, 0, 10)
	if !isclose(y4, 7.5, 8) {
		t.Fail()
	}
}

func TestInterpolate(t *testing.T) {
	nv1 := NVector{Vec3{0, 3, 2}}
	nv2 := NVector{Vec3{-7, 5, -3}}
	nv_intermediate := nv1.Interpolate(&nv2, 0.2)

	if !isclose(nv_intermediate.Vec3[0], -1.4, 8) {
		t.Fail()
	}
	if !isclose(nv_intermediate.Vec3[1], 3.4, 8) {
		t.Fail()
	}
	if !isclose(nv_intermediate.Vec3[2], 1.0, 8) {
		t.Fail()
	}
}

func TestIntersection1(t *testing.T) {
	ll1, _ := NewLonLat(-50, 0)
	ll2, _ := NewLonLat(-30, 0)
	ll3, _ := NewLonLat(-40, -5)
	ll4, _ := NewLonLat(-40, 3)

	nv1 := ll1.ToNVector()
	nv2 := ll2.ToNVector()
	nv3 := ll3.ToNVector()
	nv4 := ll4.ToNVector()

	nv_intersection, err := Intersection(&nv1, &nv2, &nv3, &nv4)
	if err != nil {
		t.Error()
	}
	ll_intersection := nv_intersection.ToLonLat()

	fmt.Println("error", err)
	fmt.Println("lon", ll_intersection.Lon*180/math.Pi)
	fmt.Println("lat", ll_intersection.Lat*180/math.Pi)

	if (!isclose(ll_intersection.Lon*180/math.Pi, -40, 8)) ||
		(!isclose(ll_intersection.Lat*180/math.Pi, 0.0, 8)) {
		t.Fail()
	}
}

func TestIntersection2(t *testing.T) {
	// No intersection
	ll1, _ := NewLonLat(30, 20)
	ll2, _ := NewLonLat(35, 25)
	ll3, _ := NewLonLat(32, 23)
	ll4, _ := NewLonLat(37, 28)

	nv1 := ll1.ToNVector()
	nv2 := ll2.ToNVector()
	nv3 := ll3.ToNVector()
	nv4 := ll4.ToNVector()

	_, err := Intersection(&nv1, &nv2, &nv3, &nv4)
	expected_err := NoIntersectionError{}
	if err != expected_err {
		t.Fail()
	}
}
