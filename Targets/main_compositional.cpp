#include "TPZCompositional.h"

int main() {
  TPZCompositional_H2O_CO2_NaCl compositional;
  TPZManVector<REAL,2> GradKCO2({1, 0.});
  TPZManVector<REAL,2> GradKH2O({-1, 0.05});
  TPZManVector<REAL,2> GradKNaCl({0.000002, 0.0005});
  
  TPZComponent &CO2 = compositional.GetComponent(ECO2);
  TPZComponent &H2O = compositional.GetComponent(EH2O);
  TPZComponent &NaCl = compositional.GetComponent(ENaCl);
  CO2.SetK(1.e-4, GradKCO2, true);
  H2O.SetK(1.e4, GradKH2O, false);
  NaCl.SetK(0.0, GradKNaCl, false);
  TPZManVector<REAL, 2> RefVals({1.0e6, 300.0});
  CO2.SetReferenceValues(RefVals);
  H2O.SetReferenceValues(RefVals);
  NaCl.SetReferenceValues(RefVals);

  compositional.SetZ(0.4, 0.5, 0.1);
  compositional.SetRhoM(1600., 1000., 2000.);

  TPZManVector<REAL, 5> Val({1e6, 300, 0.4, 0.4, 0.2});
  TPZManVector<REAL, 5> DeltaVal({1000.0, 10.0, 0.04, -0.03, -0.01});

  if (1)
  {
    CO2.ValidateKDerivative(Val, DeltaVal);
    H2O.ValidateKDerivative(Val, DeltaVal);
    NaCl.ValidateKDerivative(Val, DeltaVal);
  }
  if (1) 
  {
    compositional.ValidateXDerivative(Val, DeltaVal);
  }
  if (1) 
  {
    compositional.ValidateDDerivative(Val, DeltaVal);
  }
  if (1) 
  {
    compositional.ValidateSDerivative(Val, DeltaVal);
  }
  return 0;
  compositional.ComputeMolarFractions();
  compositional.PrintMolarFractions(std::cout);
  compositional.ComputeMolarDensityFractions();
  compositional.PrintMolarDensityFractions(std::cout);

  compositional.ComputeSaturationsAndTotalMolarDensity();
  compositional.PrintSaturations(std::cout);
  compositional.PrintMolarDensity(std::cout);

  TPZFNMatrix<20,REAL> DSR(4,5);
  compositional.GetDSRhoM(DSR);
  DSR.Print("DSR", std::cout);
  return 0;
}