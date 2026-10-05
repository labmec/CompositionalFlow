// Three-phase compositional flow state for H2O, CO2, and NaCl.
#ifndef TPZCOMPOSITIONAL_H
#define TPZCOMPOSITIONAL_H

#include "pzmanvector.h"
#include "pzfmatrix.h"
#include "pzvec_extras.h"

enum MComponents {
  ECO2,
  EH2O,
  ENaCl
};

struct TPZComponent {
	TPZComponent() = default;
	TPZComponent(const TPZComponent&) = default;
	TPZComponent& operator=(const TPZComponent&) = default;
	~TPZComponent() = default;

  inline void SetK(const REAL &K, const TPZVec<REAL> &GradK, const bool InverseK) 
  {
    fIsInverse = InverseK;
    fK = K;
    fK0 = K;
    fDK = GradK;
    fDK0 = GradK;
  }

  inline void SetReferenceValues(const TPZVec<REAL> &RefVals)
  {
    fRefVals = RefVals;
  }

  void Print(const std::string& name, std::ostream &out);

  void ValidateX();

  void ComputeK(const TPZVec<REAL>& Val);

  void ComputeDK(const TPZVec<REAL>& Val);

  void ValidateKDerivative(const TPZVec<REAL>&Val, const TPZVec<REAL>& dVal);

  // Component K equilibrium constant between phases
	double fK = 0.0, fK0 = 0.0;
	TPZManVector<REAL, 2> fDK = TPZManVector<REAL, 2>(2, 0.), fDK0 = TPZManVector<REAL, 2>(2, 0.);;
  bool fIsInverse = false; //Formulation using k or 1/k
  
  // molar fractions based on fugacity
  TPZManVector<REAL, 2> fx = TPZManVector<REAL, 2>(2, 0.);
  TPZFNMatrix<4, REAL> fDx{2, 2, 0.0};

  // Overall molar fraction of the component in the mixture
	double fz = 0.0;
	TPZManVector<REAL, 2> fDz = TPZManVector<REAL, 2>(2, 0.);

  // Reference pressure and temperature of measured constants
  TPZManVector<REAL, 2> fRefVals = TPZManVector<REAL, 2>(2, 0.0);
  
  // Overall density of the component in the mixture
	// double frho = 0.0;
	// TPZManVector<REAL, 2> fDrho = TPZManVector<REAL, 2>(2, 0.);
};

struct TPZPhase {
  TPZPhase() = default;
  TPZPhase(const TPZPhase&) = default;
  TPZPhase& operator=(const TPZPhase&) = default;
  ~TPZPhase() = default;

  void ValidateD();

  void ValidateS();

  // Phase saturation
  double fS = 0.0;
  // Derivative of phase saturation with respect to pressure and temperature and Z values
  TPZManVector<REAL, 5> fDS = TPZManVector<REAL, 5>(5, 0.);
  // Fractional density of the phase in the mixture
  double fD = 0.0;
  // Derivative of density with respect to pressure, temperature and z_c
  TPZManVector<REAL, 5> fDD = TPZManVector<REAL, 5>(5, 0.);
  // Phase pressure
  double fP = 0.0;
  // Phase molar density
  double frhoM = 0.0;
  // Derivative of molar density with respect to pressure and temperature
  TPZManVector<REAL, 2> fDrhoM = TPZManVector<REAL, 2>(2, 0.);
};

class TPZCompositional_H2O_CO2_NaCl {
public:
  TPZCompositional_H2O_CO2_NaCl() = default;
  TPZCompositional_H2O_CO2_NaCl(const TPZCompositional_H2O_CO2_NaCl&) = default;
  TPZCompositional_H2O_CO2_NaCl& operator=(const TPZCompositional_H2O_CO2_NaCl&) = default;
  virtual ~TPZCompositional_H2O_CO2_NaCl() = default;

public:
  inline void SetZ(const REAL& zCO2, const REAL& zH2O, const REAL& zNaCl)
  {
    fCO2.fz = zCO2;
    fH2O.fz = zH2O;
    fNaCl.fz = zNaCl;
  }

  inline void SetRhoM(const REAL& RhoLM, const REAL& RhoGM, const REAL& RhoSM)
  {
    fLiquidPhase.frhoM = RhoLM;
    fGasPhase.frhoM = RhoGM;
    fSolidPhase.frhoM = RhoSM;
  }

  TPZComponent& GetComponent(const MComponents &ComponentIdx) {
    switch (ComponentIdx)
    {
    case ECO2:
      return fCO2;
    case EH2O:
      return fH2O;
    case ENaCl:
      return fNaCl;
    default:
      DebugStop();
    }
  }

  TPZFNMatrix<6, REAL> GetX()
  {
    return TPZFNMatrix<6, REAL>({fCO2.fx[0], fCO2.fx[1], fH2O.fx[0], fH2O.fx[1], fNaCl.fx[0], fNaCl.fx[1]});
  }

  void SetX(TPZFNMatrix<6, REAL>& xVals)
  {
    fCO2.fx[0] = xVals(0); 
    fCO2.fx[1] = xVals(1);
    fH2O.fx[0] = xVals(2);
    fH2O.fx[1] = xVals(3);
    fNaCl.fx[0] = xVals(4);
    fNaCl.fx[1] = xVals(5);
  }

  TPZFNMatrix<12, REAL> GetDX()
  {
    TPZFNMatrix<12, REAL> DX {6, 2, 0.0};
    for (int i = 0; i < 2; i++)
    {
      DX(0, i) = fCO2.fDx(0, i);
      DX(1, i) = fCO2.fDx(1, i);
      DX(2, i) = fH2O.fDx(0, i);
      DX(3, i) = fH2O.fDx(1, i);
      DX(4, i) = fNaCl.fDx(0, i);
      DX(5, i) = fNaCl.fDx(1, i);
    }

    return DX;
  }

  void SetDX(const int& i, TPZFNMatrix<6, REAL> DVec)
  {
    fCO2.fDx(0, i) = DVec(0);
    fCO2.fDx(1, i) = DVec(1);
    fH2O.fDx(0, i) = DVec(2);
    fH2O.fDx(1, i) = DVec(3);
    fNaCl.fDx(0, i) = DVec(4);
    fNaCl.fDx(1, i) = DVec(5);
  }

  void GetD(TPZFMatrix<REAL> &D) const
  {
    if(D.Rows() != 3 || D.Cols() != 1) DebugStop();
    D(0, 0) = fLiquidPhase.fD;
    D(1, 0) = fGasPhase.fD;
    D(2, 0) = fSolidPhase.fD;

     return;
  }

  void SetD(const TPZFMatrix<REAL> &Val)
  {
    fLiquidPhase.fD = Val(0,0);
    fGasPhase.fD = Val(1,0);
    fSolidPhase.fD = Val(2,0);
  }

  void GetDD(TPZFMatrix<REAL> &DD)
  {
    if(DD.Rows() != 3 || DD.Cols() != 5) DebugStop();

    for (int i = 0; i < 5; i++)
    {
      DD(0, i) = fLiquidPhase.fDD[i];
      DD(1, i) = fGasPhase.fDD[i];
      DD(2, i) = fSolidPhase.fDD[i];
    }

    return;
  }

  void SetDD(const int& i, TPZFNMatrix<3, REAL> DVec)
  {
    fLiquidPhase.fDD[i] = DVec(0);
    fGasPhase.fDD[i] = DVec(1);
    fSolidPhase.fDD[i] = DVec(2);
  }

  void GetSRhoM(TPZFMatrix<REAL> &SR);

  void GetDSRhoM(TPZFMatrix<REAL> &DSR);

  void ApplyMethodToComponents(void(TPZComponent::*pmemfn)());

  void ComputeMolarFractions(const int maxIter = 1000, const double &tol = 1e-10);

  void ComputeMolarDensityFractions();

  void ComputeSaturationsAndTotalMolarDensity();

  void ValidateXDerivative(const TPZVec<REAL>& Val, const TPZVec<REAL>& dVal);

  void ValidateDDerivative(const TPZVec<REAL>& Val, const TPZVec<REAL>& dVal);

  void ValidateSDerivative(const TPZVec<REAL>& Val, const TPZVec<REAL>& dVal);

  void PrintMolarFractions(std::ostream &out);

  void PrintMolarDensityFractions(std::ostream &out);

  void PrintSaturations(std::ostream &out);

  void PrintMolarDensity(std::ostream &out);

  // load the state variables and compute the derivatives of the internal variables with respect to P and T
  void LoadState(REAL P, REAL T);

protected:
  // Three components
	TPZComponent fCO2;
	TPZComponent fH2O;
	TPZComponent fNaCl;

  // Three phases
  TPZPhase fLiquidPhase;
  TPZPhase fGasPhase;
  TPZPhase fSolidPhase;

  // Molar density of the mixture
	double frhoMT = 0.0;
  // Derivative of Molar density with respect to pressure and temperature and Z values
	TPZManVector<REAL, 2> fDrhoMT = TPZManVector<REAL, 5>(5, 0.);

  // Mixture pressure and temperature
	double fP = 0.0;
	double fT = 0.0;
};

#endif // TPZCOMPOSITIONAL_H
