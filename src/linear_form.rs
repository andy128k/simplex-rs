use num_rational::Rational32;
use std::collections::BTreeMap;
use std::sync::atomic::{AtomicUsize, Ordering};

#[derive(Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Debug)]
pub struct Variable(usize);

#[derive(Clone)]
pub struct LinearForm {
    pub(crate) coeffs: BTreeMap<Variable, Rational32>,
}

impl Variable {
    pub fn new() -> Self {
        static LAST_ID: AtomicUsize = AtomicUsize::new(0);
        let n = LAST_ID.fetch_add(1, Ordering::SeqCst);
        Self(n)
    }
}

impl std::ops::Mul<i32> for Variable {
    type Output = LinearForm;

    fn mul(self, rhs: i32) -> Self::Output {
        LinearForm::mononomial(rhs.into(), self)
    }
}

impl std::ops::Mul<Rational32> for Variable {
    type Output = LinearForm;

    fn mul(self, rhs: Rational32) -> Self::Output {
        LinearForm::mononomial(rhs, self)
    }
}

impl std::ops::Mul<Variable> for i32 {
    type Output = LinearForm;

    fn mul(self, rhs: Variable) -> Self::Output {
        LinearForm::mononomial(self.into(), rhs)
    }
}

impl std::ops::Mul<Variable> for Rational32 {
    type Output = LinearForm;

    fn mul(self, rhs: Variable) -> Self::Output {
        LinearForm::mononomial(self, rhs)
    }
}

impl std::ops::Div<i32> for Variable {
    type Output = LinearForm;

    fn div(self, rhs: i32) -> Self::Output {
        LinearForm::mononomial(Rational32::new(1, rhs), self)
    }
}

impl std::ops::Div<Rational32> for Variable {
    type Output = LinearForm;

    fn div(self, rhs: Rational32) -> Self::Output {
        LinearForm::mononomial(Rational32::from(1) / rhs, self)
    }
}

impl From<Variable> for LinearForm {
    fn from(var: Variable) -> Self {
        LinearForm::mononomial(1.into(), var)
    }
}

impl std::ops::Neg for Variable {
    type Output = LinearForm;

    fn neg(self) -> Self::Output {
        LinearForm::mononomial((-1).into(), self)
    }
}

impl std::ops::Add<Variable> for Variable {
    type Output = LinearForm;

    fn add(self, rhs: Variable) -> Self::Output {
        LinearForm::mononomial(1.into(), self) + LinearForm::mononomial(1.into(), rhs)
    }
}

impl std::ops::Sub<Variable> for Variable {
    type Output = LinearForm;

    fn sub(self, rhs: Variable) -> Self::Output {
        LinearForm::mononomial(1.into(), self) - LinearForm::mononomial(1.into(), rhs)
    }
}

impl<Rhs: Into<LinearForm>> std::ops::AddAssign<Rhs> for LinearForm {
    fn add_assign(&mut self, rhs: Rhs) {
        for (var, coeff) in &rhs.into().coeffs {
            self.coeffs
                .entry(*var)
                .and_modify(|c| *c += coeff)
                .or_insert(*coeff);
        }
    }
}

impl<Rhs: Into<LinearForm>> std::ops::SubAssign<Rhs> for LinearForm {
    fn sub_assign(&mut self, rhs: Rhs) {
        for (var, coeff) in &rhs.into().coeffs {
            self.coeffs
                .entry(*var)
                .and_modify(|c| *c -= coeff)
                .or_insert(-*coeff);
        }
    }
}

impl<Rhs: Into<Rational32>> std::ops::MulAssign<Rhs> for LinearForm {
    fn mul_assign(&mut self, rhs: Rhs) {
        let rhs = rhs.into();
        for coeff in self.coeffs.values_mut() {
            *coeff *= rhs;
        }
    }
}

impl std::ops::Neg for LinearForm {
    type Output = LinearForm;

    fn neg(self) -> Self::Output {
        let mut result = self.clone();
        result *= -1;
        result
    }
}

impl<Rhs: Into<LinearForm>> std::ops::Add<Rhs> for LinearForm {
    type Output = LinearForm;

    fn add(self, rhs: Rhs) -> Self {
        let mut result = self.clone();
        result += rhs;
        result
    }
}

impl<Rhs: Into<LinearForm>> std::ops::Sub<Rhs> for LinearForm {
    type Output = LinearForm;

    fn sub(self, rhs: Rhs) -> Self {
        let mut result = self.clone();
        result -= rhs;
        result
    }
}

impl LinearForm {
    pub fn empty() -> Self {
        Self {
            coeffs: BTreeMap::default(),
        }
    }

    pub fn mononomial(coeff: Rational32, var: Variable) -> Self {
        let mut coeffs = BTreeMap::new();
        coeffs.insert(var, coeff);
        Self { coeffs }
    }

    pub fn binomial(
        coeff1: Rational32,
        var1: Variable,
        coeff2: Rational32,
        var2: Variable,
    ) -> Self {
        Self::mononomial(coeff1, var1) + Self::mononomial(coeff2, var2)
    }

    pub fn get_coeff(&self, var: Variable) -> Option<Rational32> {
        self.coeffs.get(&var).cloned()
    }
}
