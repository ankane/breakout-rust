use crate::Error;

/// Parameters for detecting a single breakout (at most one change).
pub struct AmocParams {
    min_size: usize,
    alpha: f64,
    exact: bool,
}

/// Returns parameters for detecting a single breakout (at most one change).
pub fn amoc() -> AmocParams {
    AmocParams {
        min_size: 30,
        alpha: 2.0,
        exact: true,
    }
}

impl AmocParams {
    /// Sets the minimum observations between breakouts.
    pub fn min_size(&mut self, value: usize) -> &mut Self {
        self.min_size = value;
        self
    }

    /// Sets the weight of the distance between observations.
    pub fn alpha(&mut self, value: f64) -> &mut Self {
        self.alpha = value;
        self
    }

    /// Sets whether to use the exact or approximate median.
    pub fn exact(&mut self, value: bool) -> &mut Self {
        self.exact = value;
        self
    }

    /// Detects a single breakout (at most one change).
    pub fn fit(&self, z: &[f64]) -> Result<Option<usize>, Error> {
        if self.min_size < 2 {
            return Err(Error::Parameter("min_size must be at least 2".to_string()));
        }
        if self.alpha < 0.0 || self.alpha > 2.0 {
            return Err(Error::Parameter(
                "alpha must be between 0 and 2".to_string(),
            ));
        }

        if z.len() < self.min_size {
            return Ok(None);
        }

        // scale observations
        let min = z.iter().min_by(|i, j| i.partial_cmp(j).unwrap()).unwrap();
        let max = z.iter().max_by(|i, j| i.partial_cmp(j).unwrap()).unwrap();
        let denom = max - min;
        if denom == 0.0 {
            return Ok(None);
        }
        let zcounts: Vec<f64> = z.iter().map(|x| (x - min) / denom).collect();

        let (loc, stat) = if self.exact {
            crate::edmx::edmx(&zcounts, self.min_size, self.alpha)
        } else {
            crate::edm_tail::edm_tail(&zcounts, self.min_size, self.alpha)
        };

        if stat > 0.0 {
            Ok(Some(loc))
        } else {
            Ok(None)
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::Error;

    #[rustfmt::skip]
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
        let breakout = crate::amoc().min_size(5).fit(&series).unwrap();
        assert_eq!(breakout, Some(19));
    }

    #[test]
    fn test_tail() {
        let series = generate_series();
        let breakout = crate::amoc().min_size(5).exact(false).fit(&series).unwrap();
        assert_eq!(breakout, Some(20));
    }

    #[test]
    fn test_empty() {
        let series = Vec::new();
        let breakout = crate::amoc().fit(&series).unwrap();
        assert_eq!(breakout, None);
    }

    #[test]
    fn test_constant() {
        let series = vec![1.0; 100];
        let breakout = crate::amoc().fit(&series).unwrap();
        assert_eq!(breakout, None);
    }

    #[test]
    fn test_almost_constant() {
        let mut series = vec![1.0; 100];
        series[50] = 2.0;
        let breakout = crate::amoc().fit(&series).unwrap();
        assert_eq!(breakout, None);
    }

    #[test]
    #[rustfmt::skip]
    fn test_simple() {
        let series = vec![
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0
        ];
        let breakout = crate::amoc().min_size(5).fit(&series).unwrap();
        assert_eq!(breakout, Some(10));
    }

    #[test]
    fn test_bad_min_size() {
        let series = Vec::new();
        let result = crate::amoc().min_size(1).fit(&series);
        assert_eq!(
            result.unwrap_err(),
            Error::Parameter("min_size must be at least 2".to_string())
        );
    }

    #[test]
    fn test_bad_alpha() {
        let series = Vec::new();
        let result = crate::amoc().alpha(3.0).fit(&series);
        assert_eq!(
            result.unwrap_err(),
            Error::Parameter("alpha must be between 0 and 2".to_string())
        );
    }
}
