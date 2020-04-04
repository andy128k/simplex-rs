use crate::linear_form::*;
use num_rational::Rational32;
use std::collections::{BTreeMap, BTreeSet};

#[derive(Clone, Copy)]
pub enum Comparation {
    LessOrEqual,
    GreaterOrEqual,
    Equal,
}

#[derive(Clone)]
pub struct Constraint {
    pub lhs: LinearForm,
    pub cmp: Comparation,
    pub rhs: Rational32,
}

impl Constraint {
    pub fn new(lhs: LinearForm, cmp: Comparation, rhs: Rational32) -> Self {
        Self { lhs, cmp, rhs }
    }
}

struct NormalizedConstraint {
    lhs: LinearForm,
    rhs: Rational32,
    slack: Option<Variable>,
    surplus: Option<Variable>,
}

fn normalize_constraints(constraints: &[Constraint]) -> Vec<NormalizedConstraint> {
    let mut normalized_constraints = Vec::new();
    for Constraint { lhs, cmp, rhs } in constraints {
        // TODO: ensure rhs is non-negative
        match cmp {
            Comparation::LessOrEqual => {
                let slack = Variable::new();
                normalized_constraints.push(NormalizedConstraint {
                    lhs: lhs.clone() + slack,
                    rhs: *rhs,
                    slack: Some(slack),
                    surplus: None,
                });
            }
            Comparation::GreaterOrEqual => {
                let slack = Variable::new();
                let surplus = Variable::new();
                normalized_constraints.push(NormalizedConstraint {
                    lhs: lhs.clone() - slack + surplus,
                    rhs: *rhs,
                    slack: Some(slack),
                    surplus: Some(surplus),
                });
            }
            Comparation::Equal => {
                let surplus = Variable::new();
                normalized_constraints.push(NormalizedConstraint {
                    lhs: lhs.clone() + surplus,
                    rhs: *rhs,
                    slack: None,
                    surplus: Some(surplus),
                });
            }
        }
    }
    normalized_constraints
}

fn extract_variables(
    goal: &LinearForm,
    constraints: &[NormalizedConstraint],
) -> (Vec<Variable>, Vec<Variable>, Vec<Variable>) {
    let mut normal_variables = BTreeSet::new();
    let mut slack_variables = BTreeSet::new();
    let mut surplus_variables = BTreeSet::new();
    for var in goal.coeffs.keys() {
        normal_variables.insert(*var);
    }
    for constraint in constraints {
        for var in constraint.lhs.coeffs.keys() {
            if Some(var) == constraint.slack.as_ref() {
                slack_variables.insert(*var);
            } else if Some(var) == constraint.surplus.as_ref() {
                surplus_variables.insert(*var);
            } else {
                normal_variables.insert(*var);
            }
        }
    }
    (
        normal_variables.into_iter().collect(),
        slack_variables.into_iter().collect(),
        surplus_variables.into_iter().collect(),
    )
}

fn phase1(variables: &[Variable], surplus_variables: &[Variable], tableau: &mut Tableau) -> bool {
    if !tableau.run(2) {
        // unreachable!();
        return false;
    }

    if tableau.goal_value() != 0.into() {
        return false;
    }

    for (col_index, variable) in variables.iter().enumerate() {
        let remove = if surplus_variables.contains(variable) {
            tableau.get_base_row(col_index).is_none()
        } else {
            tableau[(0, col_index)] > 0.into()
        };
        if remove {
            tableau.zero_column(col_index);
        }
    }
    tableau.remove_row();

    true
}

pub fn maximize(
    goal: LinearForm,
    constraints: &[Constraint],
) -> Option<(Rational32, BTreeMap<Variable, Rational32>)> {
    let normalized_constraints = normalize_constraints(constraints);

    let (normal_variables, slack_variables, surplus_variables) =
        extract_variables(&goal, &normalized_constraints);

    let mut variables = Vec::new();
    variables.extend(&normal_variables);
    variables.extend(&slack_variables);
    variables.extend(&surplus_variables);

    // initial tableau
    let mut tableau = Tableau::empty(variables.len(), normalized_constraints.len());
    let last_column = tableau.width() - 1;

    for (row_index, constraint) in normalized_constraints.iter().enumerate() {
        for (col_index, var) in variables.iter().enumerate() {
            if let Some(coeff) = constraint.lhs.get_coeff(*var) {
                tableau[(1 + row_index, col_index)] = coeff;
            }
        }
        tableau[(1 + row_index, last_column)] = constraint.rhs;
    }

    for (col_index, var) in variables.iter().enumerate() {
        if let Some(coeff) = goal.get_coeff(*var) {
            tableau[(0, col_index)] = -coeff; // Maximize goal, so negate coeffs
        }
    }

    if !surplus_variables.is_empty() {
        let mut goal_value = 0.into();
        let mut phase1_goal = surplus_variables
            .iter()
            .fold(LinearForm::empty(), |acc, v| acc + *v);
        for constraint in &normalized_constraints {
            if constraint.surplus.is_some() {
                phase1_goal -= constraint.lhs.clone();
                goal_value -= constraint.rhs;
            }
        }

        tableau.insert_row();
        for (col_index, var) in variables.iter().enumerate() {
            if let Some(coeff) = phase1_goal.get_coeff(*var) {
                tableau[(0, col_index)] = coeff;
            }
        }
        tableau.set_goal_value(goal_value);

        if !phase1(&variables, &surplus_variables, &mut tableau) {
            return None;
        }
    }

    // Second phase
    if !tableau.run(1) {
        return None;
    }

    // Get solution
    let values = tableau.get_values();
    let mut solution = BTreeMap::new();
    for (index, var) in normal_variables.iter().enumerate() {
        let value = values[index];
        solution.insert(*var, value);
    }

    Some((tableau.goal_value(), solution))
}

pub fn minimize(
    goal: LinearForm,
    constraints: &[Constraint],
) -> Option<(Rational32, BTreeMap<Variable, Rational32>)> {
    let (value, solution) = maximize(-goal, constraints)?;
    Some((-value, solution))
}

struct Tableau {
    width: usize,
    values: Vec<Vec<Rational32>>,
}

impl Tableau {
    pub fn empty(variables: usize, constraints: usize) -> Self {
        let width = variables + 1;
        let height = constraints + 1;
        Self {
            width,
            values: vec![vec![0.into(); width]; height],
        }
    }

    pub fn width(&self) -> usize {
        self.width
    }

    pub fn height(&self) -> usize {
        self.values.len()
    }

    fn insert_row(&mut self) {
        let new_row = vec![0.into(); self.width()];
        self.values.insert(0, new_row);
    }

    fn remove_row(&mut self) {
        self.values.remove(0);
    }

    fn zero_column(&mut self, col: usize) {
        for row in 0..self.height() {
            self.values[row][col] = 0.into();
        }
    }

    fn goal_value(&self) -> Rational32 {
        let last_column = self.width() - 1;
        self.values[0][last_column]
    }

    fn set_goal_value(&mut self, goal_value: Rational32) {
        let last_column = self.width() - 1;
        self.values[0][last_column] = goal_value;
    }

    fn run(&mut self, start_row: usize) -> bool {
        loop {
            if let Some(pivot_column) = self.find_pivot_column(start_row) {
                match self.find_pivot_row(pivot_column, start_row) {
                    Some(pivot_row) => {
                        let pivot = (pivot_column, pivot_row);
                        self.next(pivot);
                    }
                    None => {
                        return false;
                    }
                }
            } else {
                return true;
            }
        }
    }

    fn is_valid_pivot_column(&self, column: usize, start_row: usize) -> bool {
        let height = self.height();
        (start_row..height).any(|row| self.values[row][column] > 0.into())
    }

    fn find_valid_pivot_column(&self, start_row: usize) -> Option<usize> {
        let count = self.width() - 1;
        for i in 0..count {
            if self.values[0][i] < 0.into() && self.is_valid_pivot_column(i, start_row) {
                return Some(i);
            }
        }
        None
    }

    fn find_pivot_column(&self, start_row: usize) -> Option<usize> {
        self.find_valid_pivot_column(start_row).or_else(|| {
            let count = self.width() - 1;
            (0..count)
                .filter(|i| self.values[0][*i] < 0.into())
                .min_by_key(|i| self.values[0][*i])
        })
    }

    fn find_pivot_row(&self, entry_var: usize, start_row: usize) -> Option<usize> {
        let width = self.width();
        let height = self.height();
        let last_column = width - 1;

        let mut pairs = (start_row..height)
            .filter(|row| self.values[*row][entry_var] > 0.into())
            .map(|row| {
                (
                    row,
                    self.values[row][last_column] / self.values[row][entry_var],
                )
            })
            .collect::<Vec<_>>();

        pairs.sort_by_key(|(j, r)| (*r, *j));

        match pairs.as_slice() {
            [] => None,
            [(j0, r0), (_j1, r1), ..] if r0 == r1 => {
                // degenerate
                Some(*j0)
            }
            [(j, _r), ..] => Some(*j),
        }
    }

    fn next(&mut self, (pivot_column, pivot_row): (usize, usize)) {
        let pivot = self.values[pivot_row][pivot_column];

        let width = self.width();
        let height = self.height();

        for col in 0..width {
            self.values[pivot_row][col] /= pivot;
        }

        for row in 0..height {
            if row != pivot_row {
                let k = self.values[row][pivot_column];
                for col in 0..width {
                    let d = k * self.values[pivot_row][col];
                    self.values[row][col] -= d;
                }
            }
        }
    }

    fn get_base_row(&self, variable: usize) -> Option<usize> {
        let height = self.height();
        let mut base_row = None;
        for row in 1..height {
            let v = self.values[row][variable];
            if v == 0.into() {
                //
            } else if v == 1.into() {
                if base_row.is_some() {
                    return None;
                }
                base_row = Some(row);
            } else {
                return None;
            }
        }
        base_row
    }

    fn get_values(&self) -> Vec<Rational32> {
        let last_column = self.width() - 1;

        let mut solution = Vec::new();
        let mut used_rows = BTreeSet::new();
        for index in 0..last_column {
            let value = if let Some(row) = self.get_base_row(index) {
                if used_rows.contains(&row) {
                    0.into()
                } else {
                    used_rows.insert(row);
                    self.values[row][last_column]
                }
            } else {
                0.into()
            };
            solution.push(value);
        }
        solution
    }
}

impl std::ops::Index<(usize, usize)> for Tableau {
    type Output = Rational32;

    fn index(&self, (row, col): (usize, usize)) -> &Self::Output {
        &self.values[row][col]
    }
}

impl std::ops::IndexMut<(usize, usize)> for Tableau {
    fn index_mut(&mut self, (row, col): (usize, usize)) -> &mut Self::Output {
        &mut self.values[row][col]
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test1() {
        let x = Variable::new();
        let y = Variable::new();
        let z = Variable::new();

        let solution = maximize(
            x + y + z,
            &[
                Constraint::new(x + y, Comparation::LessOrEqual, 5.into()),
                Constraint::new(y.into(), Comparation::LessOrEqual, 3.into()),
                Constraint::new(y + z, Comparation::LessOrEqual, 7.into()),
            ],
        );
        assert!(solution.is_some(), "Solution exists");

        let (value, solution) = solution.unwrap();
        assert_eq!(value, 12.into());
        assert_eq!(solution[&x], 5.into());
        assert_eq!(solution[&y], 0.into());
        assert_eq!(solution[&z], 7.into());
    }

    #[test]
    fn test2() {
        let x = Variable::new();
        let y = Variable::new();
        let z = Variable::new();

        let solution = minimize(
            x + y + z,
            &[
                Constraint::new(x + y, Comparation::GreaterOrEqual, 5.into()),
                Constraint::new(y.into(), Comparation::GreaterOrEqual, 3.into()),
                Constraint::new(y + z, Comparation::GreaterOrEqual, 7.into()),
            ],
        );
        assert!(solution.is_some(), "Solution exists");

        let (value, solution) = solution.unwrap();
        assert_eq!(value, 7.into());
        assert_eq!(solution[&x], 0.into());
        assert_eq!(solution[&y], 5.into());
        assert_eq!(solution[&z], 2.into());
    }

    #[test]
    fn test3() {
        let x = Variable::new();
        let y = Variable::new();
        let z = Variable::new();

        let solution = maximize(
            -x - y - z,
            &[
                Constraint::new(x.into(), Comparation::GreaterOrEqual, 5.into()),
                Constraint::new(y.into(), Comparation::GreaterOrEqual, 3.into()),
                Constraint::new(z.into(), Comparation::GreaterOrEqual, 7.into()),
                Constraint::new(z.into(), Comparation::LessOrEqual, 700.into()),
            ],
        );
        assert!(solution.is_some(), "Solution exists");

        let (value, solution) = solution.unwrap();
        assert_eq!(value, (-15).into());
        assert_eq!(solution[&x], 5.into());
        assert_eq!(solution[&y], 3.into());
        assert_eq!(solution[&z], 7.into());
    }

    #[test]
    fn test_single() {
        let x = Variable::new();

        let solution = maximize(
            x.into(),
            &[Constraint::new(
                x.into(),
                Comparation::LessOrEqual,
                3.into(),
            )],
        );
        assert!(solution.is_some(), "Solution exists");
        let (value, solution) = solution.unwrap();
        assert_eq!(value, 3.into());
        assert_eq!(solution[&x], 3.into());
    }

    #[test]
    fn test_degraged() {
        let x = Variable::new();
        let y = Variable::new();

        let solution = maximize(
            x + y,
            &[Constraint::new(x + y, Comparation::LessOrEqual, 2.into())],
        );
        assert!(solution.is_some(), "Solution exists");
        let (value, solution) = solution.unwrap();
        assert_eq!(value, 2.into());
        assert_eq!(solution[&x], 2.into());
        assert_eq!(solution[&y], 0.into());
    }

    #[test]
    fn test_unbound() {
        let x = Variable::new();

        let solution = maximize(
            x.into(),
            &[Constraint::new(
                x.into(),
                Comparation::GreaterOrEqual,
                3.into(),
            )],
        );
        assert!(solution.is_none(), "Solution doesn't exists");
    }

    #[test]
    fn test_cycle() {
        let x1 = Variable::new();
        let x2 = Variable::new();
        let x3 = Variable::new();
        let x4 = Variable::new();

        let solution = maximize(
            10 * x1 - 57 * x2 - 9 * x3 - 24 * x4,
            &[
                Constraint::new(
                    x1 / 2 - 11 / 2 * x2 - 5 / 2 * x3 + 9 * x4,
                    Comparation::LessOrEqual,
                    0.into(),
                ),
                Constraint::new(
                    x1 / 2 - 3 / 2 * x2 - 1 / 2 * x3 + x4,
                    Comparation::LessOrEqual,
                    0.into(),
                ),
                Constraint::new(1 * x1, Comparation::LessOrEqual, 1.into()),
            ],
        );
        assert!(solution.is_some(), "Solution exists");
        let (value, solution) = solution.unwrap();
        assert_eq!(value, 0.into());
        assert_eq!(solution[&x1], 0.into());
        assert_eq!(solution[&x2], 0.into());
        assert_eq!(solution[&x3], 0.into());
        assert_eq!(solution[&x4], 0.into());
    }

    #[test]
    fn test_unfeasible_solution() {
        let x1 = Variable::new();
        let x2 = Variable::new();

        let solution = maximize(
            -x1,
            &[
                Constraint::new(2 * x2, Comparation::Equal, 0.into()),
                Constraint::new(10 * x2, Comparation::Equal, 10.into()),
            ],
        );
        assert!(solution.is_none(), "Solution doesn't exist");
    }

    #[test]
    fn test_ray() {
        let x1 = Variable::new();
        let x2 = Variable::new();

        let solution = maximize(
            -x1,
            &[Constraint::new(x2 / 2, Comparation::Equal, 10.into())],
        );
        assert!(solution.is_some(), "Solution exists");
        let (value, solution) = solution.unwrap();
        assert_eq!(value, 0.into());
        assert_eq!(solution[&x1], 0.into());
        assert_eq!(solution[&x2], 20.into());
    }

    #[test]
    fn test_point() {
        let x1 = Variable::new();
        let x2 = Variable::new();

        let solution = maximize(
            -x1 - x2,
            &[
                Constraint::new(x1.into(), Comparation::Equal, 4.into()),
                Constraint::new(x2.into(), Comparation::Equal, 5.into()),
            ],
        );
        assert!(solution.is_some(), "Solution exists");
        let (value, solution) = solution.unwrap();
        assert_eq!(value, (-9).into());
        assert_eq!(solution[&x1], 4.into());
        assert_eq!(solution[&x2], 5.into());
    }

    #[test]
    fn test_eq() {
        let x1 = Variable::new();
        let x2 = Variable::new();
        let x3 = Variable::new();
        let x4 = Variable::new();
        let x5 = Variable::new();
        let x6 = Variable::new();

        let solution = maximize(
            8 * x1 + 2 * x2 + 7 * x3 + 3 * x4 + 6 * x5 + 4 * x6,
            &[
                Constraint::new(x1 + x3 + x5, Comparation::Equal, 23.into()),
                Constraint::new(x2 + x4 + x6, Comparation::Equal, 23.into()),
                Constraint::new(x1.into(), Comparation::GreaterOrEqual, 10.into()),
                Constraint::new(x3.into(), Comparation::GreaterOrEqual, 8.into()),
                Constraint::new(x5.into(), Comparation::GreaterOrEqual, 5.into()),
            ],
        );
        assert!(solution.is_some(), "Solution exists");
        let (value, solution) = solution.unwrap();
        assert_eq!(value, 258.into());
        assert_eq!(solution[&x1], 10.into());
        assert_eq!(solution[&x2], 23.into());
        assert_eq!(solution[&x3], 8.into());
        assert_eq!(solution[&x4], 0.into());
        assert_eq!(solution[&x5], 5.into());
        assert_eq!(solution[&x6], 0.into());
    }

    #[test]
    fn test_degneracy() {
        let x1 = Variable::new();
        let x2 = Variable::new();

        let solution = maximize(
            -8 * x1 - 7 * x2,
            &[
                Constraint::new(x1 + x2, Comparation::LessOrEqual, 18.into()),
                Constraint::new(x1.into(), Comparation::GreaterOrEqual, 10.into()),
                Constraint::new(x2.into(), Comparation::GreaterOrEqual, 8.into()),
            ],
        );
        assert!(solution.is_some(), "Solution exists");
        let (value, solution) = solution.unwrap();
        assert_eq!(value, (-136).into());
        assert_eq!(solution[&x1], 10.into());
        assert_eq!(solution[&x2], 8.into());
    }
}
