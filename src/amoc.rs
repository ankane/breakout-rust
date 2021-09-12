pub struct AmocParams {
    min_size: usize,
    alpha: f64,
    exact: bool
}

pub fn amoc() -> AmocParams {
    AmocParams {
        min_size: 30,
        alpha: 2.0,
        exact: true
    }
}

impl AmocParams {
    pub fn min_size(&mut self, value: usize) -> &mut Self {
        self.min_size = value;
        self
    }

    pub fn exact(&mut self, value: bool) -> &mut Self {
        self.exact = value;
        self
    }

    pub fn fit(&self, z: &[f64]) -> Option<usize> {
        assert!(self.min_size >= 2, "min_size must be at least 2");
        assert!(self.alpha >= 0.0 && self.alpha <= 2.0, "alpha must be between 0 and 2");

        if z.len() < self.min_size {
            return None;
        }

        // scale observations
        let min = z.iter().min_by(|i, j| i.partial_cmp(j).unwrap()).unwrap();
        let max = z.iter().max_by(|i, j| i.partial_cmp(j).unwrap()).unwrap();
        let denom = max - min;
        if denom == 0.0 {
            return None;
        }
        let zcounts: Vec<f64> = z.iter().map(|x| (x - min) / denom).collect();

        let (loc, stat) = if self.exact {
            crate::edmx::edmx(&zcounts, self.min_size, self.alpha)
        } else {
            crate::edm_tail::edm_tail(&zcounts, self.min_size, self.alpha)
        };

        if stat > 0.0 {
            Some(loc)
        } else {
            None
        }
    }
}

#[cfg(test)]
mod tests {
    fn generate_series() -> Vec<f64> {
        vec![
            3.0, 1.0, 2.0, 3.0, 2.0, 1.0, 1.0, 2.0, 2.0, 3.0,
            6.0, 4.0, 4.0, 5.0, 6.0, 4.0, 4.0, 4.0, 6.0, 5.0,
            9.0, 8.0, 7.0, 9.0, 8.0, 9.0, 9.0, 9.0, 7.0, 9.0
        ]
    }

    #[test]
    fn test_amoc() {
        let series = generate_series();
        let breakout = crate::amoc().min_size(5).fit(&series);
        assert_eq!(breakout, Some(19));
    }

    #[test]
    fn test_tail() {
        let series = generate_series();
        let breakout = crate::amoc().min_size(5).exact(false).fit(&series);
        assert_eq!(breakout, Some(20));
    }

    #[test]
    fn test_empty() {
        let series = Vec::new();
        let breakout = crate::amoc().fit(&series);
        assert_eq!(breakout, None);
    }

    #[test]
    fn test_constant() {
        let series = vec![1.0; 100];
        let breakout = crate::amoc().fit(&series);
        assert_eq!(breakout, None);
    }

    #[test]
    fn test_almost_constant() {
        let mut series = vec![1.0; 100];
        series[50] = 2.0;
        let breakout = crate::amoc().fit(&series);
        assert_eq!(breakout, None);
    }

    #[test]
    fn test_simple() {
        let series = vec![
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0
        ];
        let breakout = crate::amoc().min_size(5).fit(&series);
        assert_eq!(breakout, Some(10));
    }
}
