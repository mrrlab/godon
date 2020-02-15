package optimize

import (
	"encoding/json"
	"testing"
)

const (
	json1 = "{\"a\":7.2,\"b\":1.17e-22,\"c\":0,\"d \\\"!\":0.999999}"
)

func TestMarshalParameters(tst *testing.T) {
	var pars FloatParameters
	a := 7.2
	b := 1.17e-22
	c := 0.0
	d := 0.999999
	pars.Append(NewBasicFloatParameter(&a, "a"))
	pars.Append(NewBasicFloatParameter(&b, "b"))
	pars.Append(NewBasicFloatParameter(&c, "c"))
	pars.Append(NewBasicFloatParameter(&d, "d \"!"))
	j, err := json.Marshal(pars)
	if err != nil {
		tst.Error("Error: ", err)
	}
	if string(j) != json1 {
		tst.Errorf("Incorrect encoded json value. Expected:\n'%v'\n got\n'%v'", json1, string(j))
	}
}

func TestUnmarshalParameters(tst *testing.T) {
	var pars FloatParameters
	a := 1.0
	b := 1.0
	c := 1.0
	d := 1.0
	pars.Append(NewBasicFloatParameter(&a, "a"))
	pars.Append(NewBasicFloatParameter(&b, "b"))
	pars.Append(NewBasicFloatParameter(&c, "c"))
	pars.Append(NewBasicFloatParameter(&d, "d \"!"))
	err := json.Unmarshal([]byte(json1), &pars)
	if err != nil {
		tst.Error("Error: ", err)
	}
	j, err := json.Marshal(pars)
	if string(j) != json1 {
		tst.Errorf("Incorrect encoded json value. Expected:\n'%v'\n got\n'%v'", json1, string(j))
	}
}
