use std::collections::btree_map::Entry;
use std::collections::BTreeMap;

pub struct Multiset<T> {
    map: BTreeMap<T, usize>,
    len: usize,
}

impl<T: std::cmp::Ord> Multiset<T> {
    pub fn new() -> Self {
        Self {
            map: BTreeMap::new(),
            len: 0,
        }
    }

    pub fn len(&self) -> usize {
        self.len
    }

    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    pub fn first(&self) -> Option<&T> {
        self.map.first_key_value().map(|v| v.0)
    }

    pub fn insert(&mut self, value: T) {
        *self.map.entry(value).or_insert(0) += 1;
        self.len += 1;
    }

    pub fn remove(&mut self, value: T) -> bool {
        if let Entry::Occupied(mut o) = self.map.entry(value) {
            *o.get_mut() -= 1;
            if *o.get() == 0 {
                o.remove_entry();
            }
            self.len -= 1;
            true
        } else {
            false
        }
    }

    pub fn clear(&mut self) {
        self.map.clear();
        self.len = 0;
    }
}

#[cfg(test)]
mod tests {
    use crate::multiset::Multiset;

    #[test]
    fn test_works() {
        let mut mset = Multiset::new();
        assert_eq!(mset.len(), 0);
        assert!(mset.is_empty());
        assert_eq!(mset.first(), None);

        mset.insert(2);
        mset.insert(1);
        mset.insert(1);
        mset.insert(3);

        assert_eq!(mset.len(), 4);
        assert!(!mset.is_empty());
        assert_eq!(mset.first(), Some(&1));

        assert!(mset.remove(2));
        assert!(!mset.remove(2));

        assert_eq!(mset.len(), 3);

        mset.clear();

        assert_eq!(mset.len(), 0);
        assert!(mset.is_empty());
        assert_eq!(mset.first(), None);
    }
}
