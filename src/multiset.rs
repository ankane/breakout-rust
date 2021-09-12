use std::collections::btree_map::Entry;
use std::collections::BTreeMap;

pub struct MultiSet<T> {
    map: BTreeMap<T, usize>,
    len: usize
}

impl<T: std::cmp::Ord> MultiSet<T> {
    pub fn new() -> Self {
        Self {
            map: BTreeMap::new(),
            len: 0
        }
    }

    pub fn len(&self) -> usize {
        self.len
    }

    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    pub fn first(&self) -> Option<&T> {
        self.map.keys().next()
    }

    pub fn insert(&mut self, value: T) {
        *self.map.entry(value).or_insert(0) += 1;
        self.len += 1;
    }

    pub fn remove(&mut self, value: T) {
        if let Entry::Occupied(mut o) = self.map.entry(value) {
            *o.get_mut() -= 1;
            if *o.get() == 0 {
                o.remove_entry();
            }
            self.len -= 1;
        }
    }

    pub fn clear(&mut self) {
        self.map.clear();
        self.len = 0;
    }
}
