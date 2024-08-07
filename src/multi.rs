use crate::Error;

/// Parameters for detecting multiple breakouts.
pub struct MultiParams {
    min_size: usize,
    degree: i32,
    beta: Option<f64>,
    percent: Option<f64>,
}

/// Returns parameters for detecting multiple breakouts.
pub fn multi() -> MultiParams {
    MultiParams {
        min_size: 30,
        degree: 1,
        beta: None,
        percent: None,
    }
}

impl MultiParams {
    /// Sets the minimum observations between breakouts.
    pub fn min_size(&mut self, value: usize) -> &mut Self {
        self.min_size = value;
        self
    }

    /// Sets the degree of the penalization polynomial.
    pub fn degree(&mut self, value: i32) -> &mut Self {
        self.degree = value;
        self
    }

    /// Sets the penalization term.
    pub fn beta<T>(&mut self, value: T) -> &mut Self
    where
        T: Into<Option<f64>>,
    {
        self.beta = value.into();
        self
    }

    /// Sets the minimum percent change in goodness of fit statistic.
    pub fn percent<T>(&mut self, value: T) -> &mut Self
    where
        T: Into<Option<f64>>,
    {
        self.percent = value.into();
        self
    }

    /// Detects breakouts in a series.
    pub fn fit(&self, z: &[f64]) -> Result<Vec<usize>, Error> {
        if self.min_size < 2 {
            return Err(Error::Parameter("min_size must be at least 2".to_string()));
        }
        if self.beta.is_some() && self.percent.is_some() {
            return Err(Error::Parameter(
                "beta and percent cannot be passed together".to_string(),
            ));
        }
        if self.degree < 0 || self.degree > 2 {
            return Err(Error::Parameter("degree must be 0, 1, or 2".to_string()));
        }

        if z.len() < self.min_size {
            return Ok(Vec::new());
        }

        // scale observations
        let min = z.iter().min_by(|i, j| i.partial_cmp(j).unwrap()).unwrap();
        let max = z.iter().max_by(|i, j| i.partial_cmp(j).unwrap()).unwrap();
        let denom = max - min;
        if denom == 0.0 {
            return Ok(Vec::new());
        }
        let zcounts: Vec<f64> = z.iter().map(|x| (x - min) / denom).collect();

        if self.percent.is_some() {
            Ok(crate::edm_multi::edm_percent(
                &zcounts,
                self.min_size,
                self.percent.unwrap(),
                self.degree,
            ))
        } else {
            Ok(crate::edm_multi::edm_multi(
                &zcounts,
                self.min_size,
                self.beta.unwrap_or(0.008),
                self.degree,
            ))
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
    fn test_multi() {
        let series = generate_series();
        let breakouts = crate::multi().min_size(5).fit(&series).unwrap();
        assert_eq!(vec![10, 15, 20], breakouts);
    }

    #[test]
    fn test_percent() {
        let series = generate_series();
        let breakouts = crate::multi()
            .min_size(5)
            .percent(0.5)
            .fit(&series)
            .unwrap();
        assert_eq!(vec![8, 19], breakouts);
    }

    #[test]
    fn test_empty() {
        let series = Vec::new();
        let breakouts = crate::multi().fit(&series).unwrap();
        assert!(breakouts.is_empty());
    }

    #[test]
    fn test_constant() {
        let series = vec![1.0; 100];
        let breakouts = crate::multi().fit(&series).unwrap();
        assert!(breakouts.is_empty());
    }

    #[test]
    fn test_almost_constant() {
        let mut series = vec![1.0; 100];
        series[50] = 2.0;
        let breakouts = crate::multi().fit(&series).unwrap();
        assert!(breakouts.is_empty());
    }

    #[test]
    #[rustfmt::skip]
    fn test_simple() {
        let series = vec![
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0
        ];
        let breakouts = crate::multi().min_size(5).fit(&series).unwrap();
        assert_eq!(vec![10], breakouts);
    }

    #[test]
    fn test_bad_min_size() {
        let series = Vec::new();
        let result = crate::multi().min_size(1).fit(&series);
        assert_eq!(
            result.unwrap_err(),
            Error::Parameter("min_size must be at least 2".to_string())
        );
    }

    #[test]
    fn test_beta_percent() {
        let series = Vec::new();
        let result = crate::multi().beta(0.008).percent(0.5).fit(&series);
        assert_eq!(
            result.unwrap_err(),
            Error::Parameter("beta and percent cannot be passed together".to_string())
        );
    }

    #[test]
    fn test_bad_degree() {
        let series = Vec::new();
        let result = crate::multi().degree(3).fit(&series);
        assert_eq!(
            result.unwrap_err(),
            Error::Parameter("degree must be 0, 1, or 2".to_string())
        );
    }
}
