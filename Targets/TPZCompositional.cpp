#include "TPZCompositional.h"
#include "pzfmatrix.h"

void TPZComponent::Print(const std::string& name, std::ostream &out)
{
  out << "Molar fraction of " << name << ": " << fx << std::endl;
  fDx.Print("Dx: ", out, EFormatted);
  out << "K fugacity ratio of " << name << ": " << fK << std::endl;
}

void TPZComponent::ValidateX()
{
  for (int i = 0; i < fx.size(); i++)
  {
    if (fx[i] < 0.0 || fx[i] > 1.0) DebugStop();
  }
}

void TPZComponent::ComputeK(const TPZVec<REAL>& Val)
{
  double inc = 0.0;
  if (fIsInverse)
  {
    for (int i = 0; i < 2; i++) inc += fK0 * fDK0[i] * (Val[i] - fRefVals[i]);
    fK = fK0 / (1.0 + inc);
  }
  else
  {
    for (int i = 0; i < 2; i++) inc += fDK0[i] * (Val[i] - fRefVals[i]);
    fK = fK0 + inc;
  }

  ComputeDK(Val);
}

void TPZComponent::ComputeDK(const TPZVec<REAL>& Val)
{
  if (fIsInverse)
  {
    for (int i = 0; i < 2; i++)
    {
      REAL a = 1.0 + fK0 * (Val[0] - fRefVals[0]) * fDK0[0] + fK0 * (Val[1] - fRefVals[1]) * fDK0[1];
      fDK[i] = -fK0 * fK0 * fDK0[i] / (a * a);
    }
  }
  else
  {
    for (int i = 0; i < 2; i++)
    {
      fDK[i] = fDK0[i];
    }
  }
}

void TPZComponent::ValidateKDerivative(const TPZVec<REAL>&Val, const TPZVec<REAL>& dVal)
{
  ComputeK(Val);
  ComputeDK(Val);
  double K = fK;
  TPZManVector<REAL, 2> DK = fDK;

  int nDiv = 10;
  TPZManVector<REAL> err(nDiv);
  for (int i = 0; i < nDiv; i++)
  {
    REAL alpha = (i + 1.0) / nDiv;

    TPZManVector<REAL, 2> increment = dVal;
    sscal(increment, alpha);
    ComputeK(Val + increment);
    double NewK = fK;

    double GuessK = fDK[0] * alpha * dVal[0] + fDK[1] * alpha * dVal[1]; //K + dK * alpha * dVal
    NewK -= K;
    // std::cout << "var newx " << Newx << std::endl;
    // std::cout << "var guess " << GuessX << std::endl;
    err[i] = abs(GuessK - NewK);

    // std::cout << "i = " << i << ", error = " << err[i] << std::endl;
  }
  std::cout << "K derivative errors = " << err << std::endl;
  for (int i = 0; i < nDiv - 1; i++)
  {
    std::cout << (std::log(err[i + 1]) - std::log(err[i])) / (std::log(1. * (i + 2)) - std::log(1. * (i + 1))) << " ";
  }
  std::cout << std::endl;
}

void TPZPhase::ValidateD()
{
  if (fD < 0.0 || fD > 1.0) DebugStop();
}

void TPZPhase::ValidateS()
{
  if (fS < 0.0 || fS > 1.0) DebugStop();
}

void TPZCompositional_H2O_CO2_NaCl::ApplyMethodToComponents(void(TPZComponent::*pmemfn)())
{
  (fCO2.*pmemfn)();
  (fH2O.*pmemfn)();
  (fNaCl.*pmemfn)();
}

void TPZCompositional_H2O_CO2_NaCl::ComputeMolarFractions(const int maxIter, const double &tol)
{
  // xCO2l, xCO2g, xH2Ol, xH2Og, xNaCll, xNaCls
  TPZFNMatrix<36, REAL> TangentMat{6, 6, 0.0};
  TPZFNMatrix<6, REAL> ResVec{6, 1, 0.0};

  auto calculateResVec = [this, &ResVec]() {
    REAL xCO2l = fCO2.fx[0];
    REAL xCO2g = fCO2.fx[1];
    REAL xH2Ol = fH2O.fx[0];
    REAL xH2Og = fH2O.fx[1];
    REAL xNaCll = fNaCl.fx[0];
    REAL xNaCls = fNaCl.fx[1];
    ResVec(0) = 1.0 - xCO2l - xH2Ol - xNaCll;
    ResVec(1) = 1.0 - xCO2g - xH2Og;
    ResVec(2) = 1.0 - xNaCls;
    ResVec(3) = fCO2.fK * xCO2g - xCO2l;
    ResVec(4) = fH2O.fK * xH2Og - xH2Ol;
    ResVec(5) = fNaCl.fK * xNaCls - xNaCll;
  };

  // std::cout << "component CO2\n";
  // fCO2.Print("CO2 ",std::cout);
  // fH2O.Print("H2O", std::cout);
  // fNaCl.Print("NaCl", std::cout);
  auto updateX = [this, &ResVec]() {
    TPZFNMatrix<6, REAL> xVals = GetX();
    xVals += ResVec;
    SetX(xVals);
  };

  TangentMat(0, 0) = 1.0;
  TangentMat(0, 2) = 1.0;
  TangentMat(0, 4) = 1.0;
  TangentMat(1, 1) = 1.0;
  TangentMat(1, 3) = 1.0;
  TangentMat(2, 5) = 1.0;
  TangentMat(3, 0) = 1.0;
  TangentMat(3, 1) = -fCO2.fK;
  TangentMat(4, 2) = 1.0;
  TangentMat(4, 3) = -fH2O.fK;
  TangentMat(5, 4) = 1.0;
  TangentMat(5, 5) = -fNaCl.fK;

  calculateResVec();

#ifdef PZDEBUG2
  ResVec.Print("ResVec = ", std::cout, EMathematicaInput);
  TangentMat.Print("TangentMat = ", std::cout, EMathematicaInput);
#endif

  TPZFNMatrix<6> TangentMatDec = TangentMat;
  REAL norm = Norm(ResVec);
  for (int i = 0; i < maxIter; i++)
  {
    TangentMatDec.SolveDirect(ResVec, ELU);
    // ResVec.Print("First delu ",std::cout);
    updateX();
    calculateResVec();
    norm = Norm(ResVec);

    if (norm < tol) break;
    if (i == maxIter - 1) DebugStop();
  }

  fCO2.ValidateX();
  fH2O.ValidateX();
  fNaCl.ValidateX();

  //ApplyMethodToComponents(TPZComponent::ValidateX);

  for (int i = 0; i < 2; i++) // Two derivatives (P, T)
  {
    TPZFNMatrix<36, REAL> DTanMat{6, 6, 0.0};
    DTanMat(3, 1) = fCO2.fDK[i];
    DTanMat(4, 3) = fH2O.fDK[i];
    DTanMat(5, 5) = fNaCl.fDK[i];

    TPZFNMatrix<6, REAL> xVals = GetX();
    TPZFNMatrix<6, REAL> Result;
    DTanMat.Multiply(xVals, Result);
    TangentMatDec.SolveDirect(Result, ELU);
    SetDX(i, Result);
  }
}

void TPZCompositional_H2O_CO2_NaCl::ValidateXDerivative(const TPZVec<REAL>& Val, const TPZVec<REAL>& dVal)
{
  LoadState(Val[0], Val[1]);
  ComputeMolarFractions();
  TPZFNMatrix<6, REAL> x = GetX();
  TPZFNMatrix<12, REAL> Dx = GetDX();
  // std::cout << "x " << x << std::endl;
  // std::cout << "Dx " << Dx << std::endl;

  int nDiv = 10;
  TPZManVector<REAL> err(nDiv);
  for (int i = 0; i < nDiv; i++)
  {
    REAL alpha = (i + 1.0) / nDiv;
    REAL dP = alpha * dVal[0];
    REAL dT = alpha * dVal[1];

    LoadState(Val[0] + dP, Val[1] + dT);
    // std::cout << "K CO2 " << fCO2.fK << std::endl;
    ComputeMolarFractions();
    TPZFNMatrix<6, REAL> Newx = GetX();

    TPZFNMatrix<6, REAL> GuessX{6, 1, 0.0}; //x + dx * alpha * dVal
    for (int i = 0; i < 6; i++)
    {
      GuessX(i) = Dx(i, 0) * alpha * dVal[0] + Dx(i, 1) * alpha * dVal[1];
      Newx(i) -= x(i);
    }
    // std::cout << "var newx " << Newx << std::endl;
    // std::cout << "var guess " << GuessX << std::endl;
    err[i] = Norm(GuessX - Newx);

    // std::cout << "i = " << i << ", error = " << err[i] << std::endl;
  }
  std::cout << "x derivative errors = " << err << std::endl;
  for(int i = 0; i < nDiv-1; i++) {
    std::cout << (std::log(err[i+1])-std::log(err[i]))/(std::log(1.*(i+2))-std::log(1.*(i+1))) << " ";
  }
  std::cout << std::endl;
}

void TPZCompositional_H2O_CO2_NaCl::ValidateDDerivative(const TPZVec<REAL>& Val, const TPZVec<REAL>& dVal)
{
  LoadState(Val[0], Val[1]);
  SetZ(Val[2], Val[3], Val[4]);
  ComputeMolarFractions();
  ComputeMolarDensityFractions();
  TPZFNMatrix<3, REAL> D(3,1);
  GetD(D);
  TPZFNMatrix<15, REAL> DD(3,5);
  GetDD(DD);

  int nDiv = 10;
  TPZManVector<REAL> err(nDiv);
  for (int i = 0; i < nDiv; i++)
  {
    REAL alpha = (i + 1.0) / nDiv;
    REAL dP = alpha * dVal[0];
    REAL dT = alpha * dVal[1];
    REAL dZCO2 = alpha * dVal[2];
    REAL dZH2O = alpha * dVal[3];
    REAL dZNaCl = alpha * dVal[4];

    LoadState(Val[0] + dP, Val[1] + dT);
    SetZ(Val[2] + dZCO2, Val[3] + dZH2O, Val[4] + dZNaCl);
    ComputeMolarFractions();
    ComputeMolarDensityFractions();
    
    TPZFNMatrix<3, REAL> NewD(3,1);
    GetD(NewD);

    // std::cout << "NewD obtained " << NewD << std::endl;
    TPZFNMatrix<6, REAL> GuessD{3, 1, 0.0}; //D + dD * alpha * dVal
    for (int i = 0; i < 3; i++)
    {
      GuessD(i) = DD(i, 0) * dP + DD(i, 1) * dT + DD(i, 2) * dZCO2 + DD(i, 3) * dZH2O + DD(i, 4) * dZNaCl;
      NewD(i) -= D(i);
    }
    // std::cout << "var newD " << NewD << std::endl;
    // std::cout << "var guessD " << GuessD << std::endl;
    err[i] = Norm(GuessD - NewD);

    // std::cout << "i = " << i << ", error = " << err[i] << std::endl;
  }
  std::cout << "D derivative errors = " << err << std::endl;
  for (int i = 0; i < nDiv - 1; i++) 
  {
    std::cout << (std::log(err[i+1]) - std::log(err[i])) / (std::log(1. * (i + 2)) - std::log(1. * (i + 1))) << " ";
  }
  std::cout << std::endl;
}

void TPZCompositional_H2O_CO2_NaCl::ValidateSDerivative(const TPZVec<REAL>& Val, const TPZVec<REAL>& dVal)
{
  LoadState(Val[0], Val[1]);
  SetZ(Val[2], Val[3], Val[4]);
  ComputeMolarFractions();
  ComputeMolarDensityFractions();
  ComputeSaturationsAndTotalMolarDensity();
  TPZFNMatrix<3, REAL> SR(4,1);
  GetSRhoM(SR);
  TPZFNMatrix<20, REAL> DSR(4,5);
  GetDSRhoM(DSR);
  // DSR.Print("DSR ",std::cout);

  int nDiv = 10;
  TPZManVector<REAL> err(nDiv);
  for (int i = 0; i < nDiv; i++)
  {
    REAL alpha = (i + 1.0) / nDiv;
    REAL dP = alpha * dVal[0];
    REAL dT = alpha * dVal[1];
    REAL dZCO2 = alpha * dVal[2];
    REAL dZH2O = alpha * dVal[3];
    REAL dZNaCl = alpha * dVal[4];

    LoadState(Val[0] + dP, Val[1] + dT);
    SetZ(Val[2] + dZCO2, Val[3] + dZH2O, Val[4] + dZNaCl);
    ComputeMolarFractions();
    ComputeMolarDensityFractions();
    ComputeSaturationsAndTotalMolarDensity();
    TPZFNMatrix<4, REAL> NewSR(4,1);
    GetSRhoM(NewSR);

    // std::cout << "NewD obtained " << NewSR << std::endl;
    TPZFNMatrix<6, REAL> GuessSR{4, 1, 0.0}; //D + dD * alpha * dVal
    for (int i = 0; i < 4; i++)
    {
      GuessSR(i) = DSR(i, 0) * dP + DSR(i, 1) * dT + DSR(i, 2) * dZCO2 + DSR(i, 3) * dZH2O + DSR(i, 4) * dZNaCl;
      NewSR(i) -= SR(i);
    }
    // std::cout << "var newSR " << NewSR << std::endl;
    // std::cout << "var guessSR " << GuessSR << std::endl;
    err[i] = Norm(GuessSR - NewSR);

    // std::cout << "i = " << i << ", error = " << err[i] << std::endl;
  }
  std::cout << "S derivative errors = " << err << std::endl;
  for (int i = 0; i < nDiv - 1; i++) 
  {
    std::cout << (std::log(err[i+1]) - std::log(err[i])) / (std::log(1. * (i + 2)) - std::log(1. * (i + 1))) << " ";
  }
  std::cout << std::endl;
}

void TPZCompositional_H2O_CO2_NaCl::ComputeMolarDensityFractions()
{
  // Dl, Dg, Ds
  TPZFNMatrix<9, REAL> TangentMat{3, 3, 0.0};
  TPZFNMatrix<3, REAL> ResVec{3, 1, 0.0};

  REAL xCO2l = fCO2.fx[0];
  REAL xCO2g = fCO2.fx[1];
  REAL xH2Ol = fH2O.fx[0];
  REAL xH2Og = fH2O.fx[1];
  REAL xNaCll = fNaCl.fx[0];
  REAL xNaCls = fNaCl.fx[1];

  REAL zCO2 = fCO2.fz;
  REAL zH2O = fH2O.fz;
  REAL zNaCl = fNaCl.fz;
  ResVec(0) = zCO2;
  ResVec(1) = zH2O;
  ResVec(2) = zNaCl;

  TangentMat.Identity();
  TangentMat(0, 0) = xCO2l;
  TangentMat(0, 1) = xCO2g;
  TangentMat(1, 0) = xH2Ol;
  TangentMat(1, 1) = xH2Og;
  TangentMat(2, 0) = xNaCll;
  TangentMat(2, 2) = xNaCls;

#ifdef PZDEBUG2
  ResVec.Print("ResVec = ", std::cout, EMathematicaInput);
  TangentMat.Print("TangentMat = ", std::cout, EMathematicaInput);
#endif

  TPZFNMatrix<9, REAL> TanMatDec = TangentMat;
  TanMatDec.SolveDirect(ResVec, ELU);

  fLiquidPhase.fD = ResVec(0);
  fGasPhase.fD = ResVec(1);
  fSolidPhase.fD = ResVec(2);

  fLiquidPhase.ValidateD();
  fGasPhase.ValidateD();
  fSolidPhase.ValidateD();

  for (int i = 0; i < 2; i++) // Five derivatives (P, T for now)
  {
    TPZFNMatrix<9, REAL> DTanMat{3, 3, 0.0};
  // TangentMat(0, 0) = xCO2l;
  // TangentMat(0, 1) = xCO2g;
  // TangentMat(1, 0) = xH2Ol;
  // TangentMat(1, 1) = xH2Og;
  // TangentMat(2, 0) = xNaCll;
  // TangentMat(2, 2) = xNaCls;
    DTanMat(0, 0) = fCO2.fDx(0, i);
    DTanMat(0, 1) = fCO2.fDx(1, i);
    DTanMat(1, 0) = fH2O.fDx(0, i);
    DTanMat(1, 1) = fH2O.fDx(1, i);
    DTanMat(2, 0) = fNaCl.fDx(0, i);
    DTanMat(2, 2) = fNaCl.fDx(1, i);

    TPZFNMatrix<3, REAL> DVals(3,1);
    GetD(DVals);
    TPZFNMatrix<3, REAL> Result;
    DTanMat.Multiply(DVals, Result);
    TanMatDec.SolveDirect(Result, ELU);
    Result *= -1.;
    SetDD(i, Result);
  }
  TPZFNMatrix<9, REAL> InverseMat;
  TangentMat.Inverse(InverseMat, ELU);
  for (int i = 0; i < 3; i++) // Now zCO2, zH2O and zNaCl
  {
    TPZFNMatrix<3, REAL> Result{3, 1, 0.0};
    for (int row = 0; row < TangentMat.Rows(); row++)
    {
      Result(row) = InverseMat(row, i); 
    }
    SetDD(i + 2, Result);
  }
}

void TPZCompositional_H2O_CO2_NaCl::ComputeSaturationsAndTotalMolarDensity()
{
  //Sl, Sg, Ss, RhoTm
  TPZFNMatrix<16, REAL> TangentMat{4, 4, 0.0}, TangentMatDec(4,4);
  TPZFNMatrix<4, REAL> ResVec{4, 1, 0.0};

  REAL rholM = fLiquidPhase.frhoM;
  REAL rhogM = fGasPhase.frhoM;
  REAL rhosM = fSolidPhase.frhoM;

  REAL Dl = fLiquidPhase.fD;
  REAL Dg = fGasPhase.fD;
  REAL Ds = fSolidPhase.fD;

  ResVec(3) = 1;

  TangentMat(0, 0) = -rholM;
  TangentMat(0, 1) = -rhogM;
  TangentMat(0, 2) = -rhosM;
  TangentMat(0, 3) = 1;
  TangentMat(1, 0) = -rholM;
  TangentMat(1, 3) = Dl;
  TangentMat(2, 2) = -rhosM;// should be 2,2?
  TangentMat(2, 3) = Ds;
  TangentMat(3, 0) = 1.0;
  TangentMat(3, 1) = 1.0;
  TangentMat(3, 2) = 1.0;

#ifdef PZDEBUG2
  ResVec.Print("ResVec = ", std::cout, EMathematicaInput);
  TangentMat.Print("TangentMat = ", std::cout, EMathematicaInput);
#endif

  TangentMatDec=TangentMat;
  TangentMatDec.SolveDirect(ResVec, ELU);

  fLiquidPhase.fS = ResVec(0);
  fGasPhase.fS = ResVec(1);
  fSolidPhase.fS = ResVec(2);
  frhoMT = ResVec(3);

  // ResVec.Print("Saturacoes ",std::cout);

  fLiquidPhase.ValidateS();
  fGasPhase.ValidateS();
  fSolidPhase.ValidateS();
  if(frhoMT < 0.) DebugStop();
  TPZFNMatrix<15,REAL> DD(3,5);
  GetDD(DD);
  // DD.Print(" DD = ",std::cout);
  TPZFNMatrix<16,REAL> DTangentMatNeg(4,4,0.);
  for(int i=0; i<5; i++) {
    DTangentMatNeg(1, 3) = -DD(0,i); //Dl;
    DTangentMatNeg(2, 3) = -DD(2,i); //Ds;
    TPZFNMatrix<3,REAL> temp(3,1);
    DTangentMatNeg.Multiply(ResVec,temp);
    TangentMatDec.SolveDirect(temp,ELU);
    fLiquidPhase.fDS[i] = temp(0,0);
    fGasPhase.fDS[i] = temp(1,0);
    fSolidPhase.fDS[i] = temp(2,0);
    fDrhoMT[i] = temp(3,0);
  }
}

void TPZCompositional_H2O_CO2_NaCl::GetSRhoM(TPZFMatrix<REAL> &SR) {
  if(SR.Rows() != 4) DebugStop();
  SR(0) = fLiquidPhase.fS;
  SR(1) = fGasPhase.fS;
  SR(2) = fSolidPhase.fS;
  SR(3) = frhoMT;
}

void TPZCompositional_H2O_CO2_NaCl::GetDSRhoM(TPZFMatrix<REAL> &DSR) {
  if(DSR.Rows() !=4 || DSR.Cols() != 5) DebugStop();
  for(int i=0; i<5; i++) {
    DSR(0,i) = fLiquidPhase.fDS[i];
    DSR(1,i) = fGasPhase.fDS[i];
    DSR(2,i) = fSolidPhase.fDS[i];
    DSR(3,i) = fDrhoMT[i];
  }
}



void TPZCompositional_H2O_CO2_NaCl::PrintMolarFractions(std::ostream &out)
{
  fCO2.Print("CO2", out);
  fH2O.Print("H2O", out);
  fNaCl.Print("NaCl", out);
}

void TPZCompositional_H2O_CO2_NaCl::PrintMolarDensityFractions(std::ostream &out)
{
  out << "Liquid phase D: " << fLiquidPhase.fD << std::endl;
  out << "Gas phase D: " << fGasPhase.fD << std::endl;
  out << "Solid phase D: " << fSolidPhase.fD << std::endl;
}

void TPZCompositional_H2O_CO2_NaCl::PrintSaturations(std::ostream &out)
{
  out << "Liquid phase S: " << fLiquidPhase.fS << std::endl;
  out << "Gas phase S: " << fGasPhase.fS << std::endl;
  out << "Solid phase S: " << fSolidPhase.fS << std::endl;
}

void TPZCompositional_H2O_CO2_NaCl::PrintMolarDensity(std::ostream &out)
{
  out << "Total molar density: " << frhoMT << std::endl;
}

// load the state variables and compute the derivatives of the internal variables with respect to P and T
void TPZCompositional_H2O_CO2_NaCl::LoadState(REAL P, REAL T) 
{
  fP = P;
  fT = T;

  TPZManVector<REAL, 2> PT = TPZManVector<REAL, 2>({fP, fT});

  fCO2.ComputeK(PT);
  fH2O.ComputeK(PT);
  fNaCl.ComputeK(PT);
}
