// Calculates the between distance using the delta points around the change point estimate

// Class used to hold all the information about the
// breakout location and the interval trees
struct Information {
    a: Vec<f64>,
    bv: Vec<f64>,
    ab: Vec<f64>,
    best_stat: f64,
    best_loc: i32,
    best_t2: i32,
    min_size: usize,
    b: i32
}

impl Information {
    fn new(bb: i32, m: usize) -> Self {
        Self {
            a: vec![0.0; 1 << (bb + 1)],
            bv: vec![0.0; 1 << (bb + 1)],
            ab: vec![0.0; 1 << (bb + 1)],
            b: bb,
            best_stat: -3.0,
            best_loc: -3,
            best_t2: -3,
            min_size: m
        }
    }
}

// Get index of leaf node interval containing x
fn get_index(b: i32, x: f64) -> usize {
    ((x.abs() * (1 << b) as f64).ceil() + (1 << b) as f64 - 1.0) as usize
}

// Return approximate quantile based on the interval tree
fn get_quantile(x: &[f64], quant: f64) -> f64 {
    let n = x.len();
    let mut k = (x[1] * quant).ceil();
    let mut l = 0.0;
    let mut u = 1.0;
    let mut i = 1;
    let mut j;
    while i < n { // Make sure that we do not go beyond the array bounds
        j = i << 1;
        if j >= n {
            break;
        }
        if x[i] == k { // Exactly k elements in this node's subtree. So can terminate early
            // Return a weighted combination of the child node medians
            let l_weight = x[j] / (x[j] + x[j + 1]);
            let r_weight = 1.0 - l_weight;
            let lu = (u + l) / 2.0;
            let rl = (u + lu) / 2.0;
            return l_weight * (quant * (lu - l) + l) + r_weight * (quant * (u - rl) + rl);
        }
        else if x[j] >= k { // More than k elements in node's left child's subtree, move to left child
            i = j;
            u = (l + u) / 2.0;
        }
        else if x[j] < k { // Not enough elements in node's left child's subtree, move to right child
            k -= x[j];
            i = j + 1;
            l = (l + u) / 2.0;
        }
    }
    quant * (u - l) + l
}

pub fn edm_tail(z: &[f64], min_size: usize, alpha: f64) -> (usize, f64) {
    let quant = 0.5;
    let n = z.len();
    let mut eps = (n as f64).ln().ceil() as i32;
    eps = eps.max(10);

    let mut info = Information::new(eps, min_size);

    let mut tau1 = info.min_size;
    let mut tau2 = tau1 * 2;

    // Populate trees and calculate statistic value for starting configuration of
    // 2 min_size segments
    for i in 0..tau1 {
        for j in i + 1..tau1 {
            let mut index = get_index(info.b, z[i] - z[j]);
            while index != 0 {
                info.a[index] += 1.0;
                index /= 2;
            }
        }
    }

    // Populate trees and calculate statistic value for starting configuration of
    // 2 min_size segments
    for i in tau1..tau2 {
        for j in i + 1..tau2 {
            let mut index = get_index(info.b, z[i] - z[j]);
            while index != 0 {
                info.bv[index] += 1.0;
                index /= 2;
            }
        }
    }

    // Populate trees and calculate statistic value for starting configuration of
    // 2 min_size segments
    for i in 0..tau1 {
        for j in tau1..tau2 {
            let mut index = get_index(info.b, z[i] - z[j]);
            while index != 0 {
                info.ab[index] += 1.0;
                index /= 2;
            }
        }
    }

    let qa = get_quantile(&info.a, quant).powf(alpha);
    let mut qb = get_quantile(&info.bv, quant).powf(alpha);
    let qc = get_quantile(&info.ab, quant).powf(alpha);

    let mut stat = 2.0 * qc - qa - qb;
    stat *= (tau1 * (tau2 - tau1) / tau2) as f64;

    info.best_stat = stat;
    info.best_loc = tau1 as i32;
    info.best_t2 = tau2 as i32;

    // Increment tau2 and update trees and statistic
    tau2 += 1;
    while tau2 < n + 1 {
        let mut index = get_index(info.b, z[tau2 - 1] - z[tau2 - 2]);
        while index != 0 { // array position 0 is not used, so we exit once we reach this location
            info.bv[index] += 1.0;
            index /= 2;
        }
        qb = get_quantile(&info.bv, quant).powf(alpha);
        stat = 2.0 * qc - qa - qb;
        stat *= ((tau2 - tau1) * tau1 / tau2) as f64;

        if stat > info.best_stat {
            info.best_stat = stat;
            info.best_loc = tau1 as i32;
            info.best_t2 = tau2 as i32;
        }
        tau2 += 1;
    }

    let mut forward_move = false;
    // Initial consideration of other possible locations for tau1
    while tau1 < n - min_size {
        // "warm start" to update tree and statistic value for other prefix series
        if forward_move {
            tau1 = forward_update(z, &mut info, tau1, quant, alpha);
        } else {
            tau1 = backward_update(z, &mut info, tau1, quant, alpha);
        }
        forward_move = !forward_move;
    }

    (info.best_loc as usize, info.best_stat)
}

fn forward_update(z: &[f64], info: &mut Information, tau1: usize, quant: f64, alpha: f64) -> usize {
    let min_size = info.min_size;
    let mut tau2 = tau1 + min_size;
    let mut tau1 = tau1;
    tau1 += 1;
    let n = z.len();
    let mut index;
    // Update A tree
    for i in tau1 - min_size..tau1 - 1 {
        index = get_index(info.b, z[i] - z[tau1 - 1]);
        while index != 0 {
            info.a[index] += 1.0;
            index /= 2;
        }
    }
    for i in tau1 - min_size..tau1 {
        index = get_index(info.b, z[i] - z[tau1 - min_size - 1]);
        while index != 0 {
            info.a[index] -= 1.0;
            index /= 2;
        }
    }
    index = get_index(info.b, z[tau1 - min_size - 1] - z[tau1 - min_size]);
    while index != 0 {
        info.a[index] += 1.0;
        index /= 2;
    }
    let qa = get_quantile(&info.a, quant).powf(alpha);

    // Update AB tree
    index = get_index(info.b, z[tau1 - 1] - z[tau1 - min_size - 1]);
    while index != 0 {
        info.ab[index] -= 1.0;
        index /= 2;
    }
    for i in tau1..tau2 {
        index = get_index(info.b, z[i] - z[tau1 - min_size - 1]);
        while index != 0 {
            info.ab[index] -= 1.0;
            index /= 2;
        }
        index = get_index(info.b, z[i] - z[tau1 - 1]);
        while index != 0 {
            info.ab[index] += 1.0;
            index /= 2;
        }
    }
    for i in tau1 - min_size..tau1 - 1 {
        index = get_index(info.b, z[i] - z[tau1 - 1]);
        while index != 0 {
            info.ab[index] -= 1.0;
            index /= 2;
        }
        index = get_index(info.b, z[i] - z[tau2]);
        while index != 0 {
            info.ab[index] += 1.0;
            index /= 2;
        }
    }
    index = get_index(info.b, z[tau1 - 1] - z[tau2]);
    while index != 0 {
        info.ab[index] += 1.0;
        index /= 2;
    }
    let qc = get_quantile(&info.ab, quant).powf(alpha);

    // Update B tree
    for i in tau1..tau2 {
        index = get_index(info.b, z[i] - z[tau1 - 1]);
        while index != 0 {
            info.bv[index] += 1.0;
            index /= 2;
        }
        index = get_index(info.b, z[i] - z[tau2]);
        while index != 0 {
            info.bv[index] += 1.0;
            index /= 2;
        }
    }

    // Increment tau2 and update statistic value as we proceed
    tau2 += 1;
    while tau2 < n + 1 {
        index = get_index(info.b, z[tau2 - 1] - z[tau2 - 2]);
        while index != 0 {
            info.bv[index] += 1.0;
            index /= 2;
        }
        let qb = get_quantile(&info.bv, quant).powf(alpha);

        let mut stat = 2.0 * qc - qa - qb;
        stat *= ((tau2 - tau1) * tau1 / tau2) as f64;

        if stat > info.best_stat {
            info.best_stat = stat;
            info.best_loc = tau1 as i32;
            info.best_t2 = tau2 as i32;
        }
        tau2 += 1;
    }

    tau1
}

fn backward_update(z: &[f64], info: &mut Information, tau1: usize, quant: f64, alpha: f64) -> usize {
    let min_size = info.min_size;
    let mut tau2 = tau1 + min_size;
    let mut tau1 = tau1;
    tau1 += 1;
    let n = z.len();
    let mut index;
    // Update A tree
    for i in tau1 - min_size..tau1 - 1 {
        index = get_index(info.b, z[i] - z[tau1 - 1]);
        while index != 0 {
            info.a[index] += 1.0;
            index /= 2;
        }
    }
    for i in tau1 - min_size..tau1 {
        index = get_index(info.b, z[i] - z[tau1 - min_size - 1]);
        while index != 0 {
            info.a[index] -= 1.0;
            index /= 2;
        }
    }
    index = get_index(info.b, z[tau1 - min_size - 1] - z[tau1 - min_size]);
    while index != 0 {
        info.a[index] += 1.0;
        index /= 2;
    }
    let qa = get_quantile(&info.a, quant).powf(alpha);

    // Update AB tree
    index = get_index(info.b, z[tau1 - 1] - z[tau1 - min_size - 1]);
    while index != 0 {
        info.ab[index] -= 1.0;
        index /= 2;
    }
    for i in tau1..tau2 {
        index = get_index(info.b, z[i] - z[tau1 - min_size - 1]);
        while index != 0 {
            info.ab[index] -= 1.0;
            index /= 2;
        }
        index = get_index(info.b, z[i] - z[tau1 - 1]);
        while index != 0 {
            info.ab[index] += 1.0;
            index /= 2;
        }
    }
    for i in tau1 - min_size..tau1 - 1 {
        index = get_index(info.b, z[i] - z[tau1 - 1]);
        while index != 0 {
            info.ab[index] -= 1.0;
            index /= 2;
        }
        index = get_index(info.b, z[i] - z[tau2]);
        while index != 0 {
            info.ab[index] += 1.0;
            index /= 2;
        }
    }
    index = get_index(info.b, z[tau1 - 1] - z[tau2]);
    while index != 0 {
        info.ab[index] += 1.0;
        index /= 2;
    }
    let qc = get_quantile(&info.ab, quant).powf(alpha);

    // Update B tree
    for i in tau1..tau1 + min_size - 1 {
        index = get_index(info.b, z[tau1 + min_size - 1] - z[i]);
        while index != 0 {
            info.bv[index] += 1.0;
            index /= 2;
        }
        index = get_index(info.b, z[i] - z[tau1 - 1]);
        while index != 0 {
            info.bv[index] -= 1.0;
            index /= 2;
        }
    }
    // Move tau2 from the end of the time series to the front.
    // Update the statistic value along the way
    tau2 = n;

    while tau2 >= tau1 + min_size {
        index = get_index(info.b, z[tau2 - 1] - z[tau2 - 2]);
        while index != 0 {
            info.bv[index] += 1.0;
            index /= 2;
        }
        let qb = get_quantile(&info.bv, quant).powf(alpha);

        let mut stat = 2.0 * qc - qa - qb;
        stat *= ((tau2 - tau1) * tau1 / tau2) as f64;

        if stat > info.best_stat {
            info.best_stat = stat;
            info.best_loc = tau1 as i32;
            info.best_t2 = tau2 as i32;
        }
        tau2 -= 1;
    }

    tau1
}
