use crate::multiset::Multiset;
use std::cmp::Ordering;

#[derive(PartialEq)]
struct MaxItem(f64);

impl Eq for MaxItem {}

// opposite of edmx
impl PartialOrd for MaxItem {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        other.0.partial_cmp(&self.0)
    }
}

impl Ord for MaxItem {
    fn cmp(&self, other: &MaxItem) -> Ordering {
        self.partial_cmp(other).unwrap()
    }
}

#[derive(PartialEq)]
struct MinItem(f64);

impl Eq for MinItem {}

// opposite of edmx
impl PartialOrd for MinItem {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        self.0.partial_cmp(&other.0)
    }
}

impl Ord for MinItem {
    fn cmp(&self, other: &MinItem) -> Ordering {
        self.partial_cmp(other).unwrap()
    }
}

fn insert_element(m: &mut Multiset<MinItem>, m2: &mut Multiset<MaxItem>, x: f64) {
    if m.is_empty() || x < m.first().unwrap().0 {
        m2.insert(MaxItem(x));
    } else {
        m.insert(MinItem(x));
    }

    if m.len() > m2.len() + 1 {
        // TODO use pop_first for performance when available
        let i = m.first().unwrap().0;
        m2.insert(MaxItem(i));
        m.remove(MinItem(i));
    } else if m2.len() > m.len() + 1 {
        // TODO use pop_first for performance when available
        let i = m2.first().unwrap().0;
        m.insert(MinItem(i));
        m2.remove(MaxItem(i));
    }
}

fn remove_element(m: &mut Multiset<MinItem>, m2: &mut Multiset<MaxItem>, x: f64) {
    if x < m.first().unwrap().0 {
        m2.remove(MaxItem(x));
    } else {
        m.remove(MinItem(x));
    }

    if m.len() > m2.len() + 1 {
        // TODO use pop_first for performance when available
        let i = m.first().unwrap().0;
        m2.insert(MaxItem(i));
        m.remove(MinItem(i));
    } else if m2.len() > m.len() + 1 {
        // TODO use pop_first for performance when available
        let i = m2.first().unwrap().0;
        m.insert(MinItem(i));
        m2.remove(MaxItem(i));
    }
}

// given a pair of trees obtain the median
fn get_median(m: &Multiset<MinItem>, m2: &Multiset<MaxItem>) -> f64 {
    if m.len() > m2.len() {
        m.first().unwrap().0
    } else if m2.len() > m.len() {
        m2.first().unwrap().0
    } else {
        (m2.first().unwrap().0 + m.first().unwrap().0) / 2.0
    }
}

fn linear(_x: f64) -> f64 {
    1.0
}

fn constant(_x: f64) -> f64 {
    0.0
}

fn quadratic(x: f64) -> f64 {
    2.0 * x + 1.0
}

pub fn edm_multi(z: &[f64], min_size: usize, beta: f64, degree: i32) -> Vec<usize> {
    // identify which type of penalization to use
    let g: fn(f64) -> f64 = match degree {
        1 => linear,
        2 => quadratic,
        _ => constant
    };

    let n = z.len();
    let mut beta = beta;
    if beta < 0.0 { // assume that beta is a positive number
        beta = -beta;
    }
    let mut prev = vec![0; n + 1];
    let mut number = vec![0; n + 1];
    let mut f = vec![-3.0; n + 1];

    // trees used to store the "upper half" of the considered observations
    let mut right_min = Multiset::new();
    let mut left_min = Multiset::new();

    // trees used to store the "lower half" of the considered observations
    let mut right_max = Multiset::new();
    let mut left_max = Multiset::new();

    // iterate over possible locations for the last change
    for s in 2 * min_size..n + 1 {
        right_max.clear();
        right_min.clear();
        left_max.clear();
        left_min.clear();

        // initialize left and right trees to account for minimum segment size
        for i in prev[min_size - 1]..min_size - 1 {
            insert_element(&mut left_min, &mut left_max, z[i]);
        }
        for i in min_size - 1..s {
            insert_element(&mut right_min, &mut right_max, z[i]);
        }

        // iterate over possible locations for the penultimate change
        for t in min_size..s - min_size + 1 { // modify limits to deal with min_size
            insert_element(&mut left_min, &mut left_max, z[t - 1]); // insert element into left tree
            remove_element(&mut right_min, &mut right_max, z[t - 1]); // remove element from right tree

            // left tree now has { Z[prev[t-1]], ..., Z[t-1] }
            // right tree now has { Z[t], ..., Z[s-1] }

            // check to see if optimal position of previous change point has changed
            // if so update the left tree
            if prev[t] > prev[t - 1] {
                for i in prev[t - 1]..prev[t] {
                    remove_element(&mut left_min, &mut left_max, z[i]);
                }
            } else if prev[t] < prev[t - 1] {
                for i in prev[t]..prev[t - 1] {
                    insert_element(&mut left_min, &mut left_max, z[i]);
                }
            }

            // calculate statistic value
            let left_median = get_median(&left_min, &left_max);
            let right_median = get_median(&right_min, &right_max);
            let normalize = ((t - prev[t]) * (s - t)) as f64 / ((s - prev[t]) as f64).powf(2.0);
            let tmp = f[t] + normalize * (left_median - right_median).powf(2.0) - beta * g(number[t] as f64);

            // check for improved optimal statistic value
            if tmp > f[s] {
                number[s] = number[t] + 1;
                f[s] = tmp;
                prev[s] = t;
            }
        }
    }

    // obtain list of optimal change point estimates
    let mut ret = Vec::new();
    let mut at = n;
    while at != 0 {
        if prev[at] != 0 { // don't insert 0 as a change point estimate
            ret.push(prev[at]);
        }
        at = prev[at];
    }
    ret.sort_unstable();

    ret
}

// Penalizes based on percent change in the statistic value.
// Linear penalty means that each new breakout must result in an at least X% increase
// Quadratic penalty means that each new breakout must result in at least an (X*k)% increase for k breakouts
pub fn edm_percent(z: &[f64], min_size: usize, percent: f64, degree: i32) -> Vec<usize> {
    // identify which type of penalization to use
    let g: fn(f64) -> f64;
    match degree {
        1 => g = linear,
        2 => g = quadratic,
        _ => g = constant
    }

    let n = z.len();
    let mut prev = vec![0; n + 1];
    let mut number = vec![0; n + 1];
    let mut f = vec![0.0; n + 1];

    // trees used to store the "upper half" of the considered observations
    let mut right_min = Multiset::new();
    let mut left_min = Multiset::new();

    // trees used to store the "lower half" of the considered observations
    let mut right_max = Multiset::new();
    let mut left_max = Multiset::new();

    // iterate over possible locations for the last change
    for s in 2 * min_size..n + 1 {
        right_max.clear();
        right_min.clear();
        left_max.clear();
        left_min.clear();

        // initialize left and right trees to account for minimum segment size
        for i in prev[min_size - 1]..min_size - 1 {
            insert_element(&mut left_min, &mut left_max, z[i]);
        }
        for i in min_size - 1..s {
            insert_element(&mut right_min, &mut right_max, z[i]);
        }

        // iterate over possible locations for the penultiamte change
        for t in min_size..s - min_size + 1 { // modify limits to deal with min_size
            insert_element(&mut left_min, &mut left_max, z[t - 1]); // insert element into left tree
            remove_element(&mut right_min, &mut right_max, z[t - 1]); // remove element from right tree

            // left tree now has { Z[prev[t-1]], ..., Z[t-1] }
            // right tree now has { Z[t], ..., Z[s-1] }

            // check to see if optimal position of previous change point has changed
            // if so update the left tree
            if prev[t] > prev[t - 1] {
                for i in prev[t - 1]..prev[t] {
                    remove_element(&mut left_min, &mut left_max, z[i]);
                }
            } else if prev[t] < prev[t - 1] {
                for i in prev[t]..prev[t - 1] {
                    insert_element(&mut left_min, &mut left_max, z[i]);
                }
            }

            // calculate statistic value
            let left_median = get_median(&left_min, &left_max);
            let right_median = get_median(&right_min, &right_max);
            let normalize = ((t - prev[t]) * (s - t)) as f64 / ((s - prev[t]) as f64).powf(2.0);
            let tmp = f[t] + normalize * (left_median - right_median).powf(2.0);

            // find best location for change point. check % condition later
            if tmp > f[s] {
                number[s] = number[t] + 1;
                f[s] = tmp;
                prev[s] = t;
            }
        }

        // check to make sure we meet the percent change requirement
        if prev[s] != 0 && f[s] - f[prev[s]] < percent * g(number[prev[s]] as f64) * f[prev[s]] {
            number[s] = number[prev[s]];
            f[s] = f[prev[s]];
            prev[s] = prev[prev[s]];
        }
    }

    // obtain list of optimal change point estimates
    let mut ret = Vec::new();
    let mut at = n;
    while at != 0 {
        if prev[at] != 0 { // don't insert 0 as a change point estimate
            ret.push(prev[at]);
        }
        at = prev[at];
    }
    ret.sort_unstable();

    ret
}
