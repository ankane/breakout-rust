// Robust estimation of 2[mean(X)-mean(Y)]^2 time normalization factor
// This is the E-Divisive E-statistic when alpha = 2
// Instead of calculating mean(X), we calculate median(X), and similarly for Y

use std::cmp::Ordering;
use std::collections::BinaryHeap;

#[derive(PartialEq)]
struct MaxItem(f64);

impl Eq for MaxItem {}

impl PartialOrd for MaxItem {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        self.0.partial_cmp(&other.0)
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

impl PartialOrd for MinItem {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        other.0.partial_cmp(&self.0)
    }
}

impl Ord for MinItem {
    fn cmp(&self, other: &MinItem) -> Ordering {
        self.partial_cmp(other).unwrap()
    }
}

fn get_median(m: &BinaryHeap<MinItem>, m2: &BinaryHeap<MaxItem>) -> f64 {
    match m.len().cmp(&m2.len()) {
        Ordering::Greater => m.peek().unwrap().0,
        Ordering::Less => m2.peek().unwrap().0,
        Ordering::Equal => (m.peek().unwrap().0 + m2.peek().unwrap().0) / 2.0,
    }
}

fn add_to_heaps(m: &mut BinaryHeap<MinItem>, m2: &mut BinaryHeap<MaxItem>, x: f64) {
    // decide on initial heap to place Item into
    if m.is_empty() || x < m.peek().unwrap().0 {
        m2.push(MaxItem(x));
    } else {
        m.push(MinItem(x));
    }

    // make sure that heaps are balanced
    if m.len() > m2.len() + 1 {
        m2.push(MaxItem(m.pop().unwrap().0));
    } else if m2.len() > m.len() + 1 {
        m.push(MinItem(m2.pop().unwrap().0));
    }
}

pub fn edmx(z: &[f64], min_size: usize, _alpha: f64) -> (usize, f64) {
    let mut left_min = BinaryHeap::new();
    let mut left_max = BinaryHeap::new();

    let mut stat_best = -3.0;
    let mut t1 = 0;

    let n = z.len();
    for i in 0..min_size - 1 {
        add_to_heaps(&mut left_min, &mut left_max, z[i]);
    }

    // iterate over breakout locations
    for tau1 in min_size..n - min_size + 1 {
        add_to_heaps(&mut left_min, &mut left_max, z[tau1 - 1]);

        let mut right_min = BinaryHeap::new();
        let mut right_max = BinaryHeap::new();

        let medl = get_median(&left_min, &left_max);

        // add first set of Items to the heaps for the right segment
        for i in tau1..tau1 + min_size - 1 {
            add_to_heaps(&mut right_min, &mut right_max, z[i]);
        }

        for tau2 in tau1 + min_size..n + 1 {
            add_to_heaps(&mut right_min, &mut right_max, z[tau2 - 1]);
            let medr = get_median(&right_min, &right_max);

            let mut stat = (medl - medr).powi(2);
            stat *= (tau1 * (tau2 - tau1)) as f64 / tau2 as f64;

            if stat > stat_best {
                t1 = tau1;
                stat_best = stat;
            }
        }
    }

    (t1, stat_best)
}
