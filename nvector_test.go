package nvector

import (
	"math"
	"testing"
)

func isclose(a, b float64, places int32) bool {
	var tol float64
	tol = 1.0 / float64(places)
	return math.Abs(a-b) < tol
}

func TestNewLonLatSimple(t *testing.T) {
	lonlat, err := NewLonLat(-140.0, 49.25)
	if err != nil {
		t.Error()
	}

	if math.Abs(lonlat.lat-0.8595746566072073) > 0.0000000001 {
		t.Fail()
	}
	if math.Abs(lonlat.lon-(-2.443460952792061)) > 0.0000000001 {
		t.Fail()
	}
}

func TestNewLonLatUnwrap(t *testing.T) {
	lonlat, err := NewLonLat(220.0, 49.25)
	if err != nil {
		t.Error()
	}

	if math.Abs(lonlat.lat-0.8595746566072073) > 0.0000000001 {
		t.Fail()
	}
	if math.Abs(lonlat.lon-(-2.443460952792061)) > 0.0000000001 {
		t.Fail()
	}
}

func TestNewLonLatInvalid(t *testing.T) {
	_, err := NewLonLat(-140.0, 92.0)

	if err == nil {
		t.Fail()
	}
}

func TestLonLatToNVector(t *testing.T) {
	ll, _ := NewLonLat(-140.0, 49.25)
	nv := ll.ToNVector()
	if !isclose(nv.Vec3[0], 0.7575649843840494, 8) {
		t.Fail()
	}
	if !isclose(nv.Vec3[1], -0.41958588098509064, 8) {
		t.Fail()
	}
	if !isclose(nv.Vec3[2], 0.5000429810657882, 8) {
		t.Fail()
	}
}

func TestNVectorToLonLat(t *testing.T) {
	nv := NVector{Vec3{0.7575649843840494,
		-0.41958588098509064,
		0.5000429810657882}}
	ll := nv.ToLonLat()
	if !isclose(ll.lon*180/math.Pi, -140, 6) {
		t.Fail()
	}
	if !isclose(ll.lat*180/math.Pi, 49.25, 6) {
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
	// meridional
	pos1, _ := NewLonLat(-140, 49.25)
	pos2, _ := NewLonLat(-140, 48.25)
	nv1 := pos1.ToNVector()
	nv2 := pos2.ToNVector()
	ellps := Ellipsoid{6378137.0, 6356752.3142}
	az, baz, err := nv1.Azimuth(&nv2, &ellps)
	if err != nil {
		t.Error()
	}
	if !isclose(az, math.Pi, 6) {
		t.Fail()
	}
	if !isclose(baz, math.Pi, 6) {
		t.Fail()
	}
}

func TestAzimuth2(t *testing.T) {
	// zonal
	pos1, _ := NewLonLat(-140, 49.25)
	pos2, _ := NewLonLat(-143, 49.25)
	nv1 := pos1.ToNVector()
	nv2 := pos2.ToNVector()
	ellps := Ellipsoid{6378137.0, 6356752.3142}
	az, baz, err := nv1.Azimuth(&nv2, &ellps)
	if err != nil {
		t.Error()
	}
	if !isclose(az, -88.8635416/180*math.Pi, 6) {
		t.Fail()
	}
	if !isclose(baz, 88.8635416/180*math.Pi, 6) {
		t.Fail()
	}
}

func TestAzimuth3(t *testing.T) {
	// crosses dateline
	pos1, _ := NewLonLat(174, -15)
	pos2, _ := NewLonLat(-177.5, 36)
	nv1 := pos1.ToNVector()
	nv2 := pos2.ToNVector()
	ellps := Ellipsoid{6378137.0, 6356752.3142}
	az, baz, err := nv1.Azimuth(&nv2, &ellps)
	if err != nil {
		t.Error()
	}
	if !isclose(az, 8.8262727/180*math.Pi, 6) {
		t.Fail()
	}
	if !isclose(baz, -169.4538460/180*math.Pi, 6) {
		t.Fail()
	}
}
