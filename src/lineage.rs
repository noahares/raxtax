use itertools::Itertools;

pub const F64_OUTPUT_ACCURACY: u32 = 2;

pub struct Lineage<'a> {
    tree: &'a Tree,
    confidence_prefix_sum: Vec<f64>,
    confidence_vectors: Vec<(usize, Vec<f64>)>,
}

impl<'a> Lineage<'a> {
    fn new(tree: &'a Tree, confidence_values: &[f64]) -> Self {
        let mut confidence_prefix_sum = vec![0.0];
        confidence_prefix_sum.extend(confidence_values.iter().scan(0.0, |sum, i| { *sum += i; Some(*sum)}));
        Self {
            tree,
            confidence_prefix_sum,
            confidence_vectors: Vec::new(),
        }
    }

    pub fn evaluate(mut self) -> Vec<(&'a String, Vec<f64>)> {
        self.eval_recurse(&self.tree.root, &[]);
        self.confidence_vectors.into_iter().map(|(idx, conf_values)| (&self.tree.lineages[idx], conf_values)).collect_vec()
    }

    fn get_confidence(&self, node: &TreeNode) -> f64 {
        self.confidence_prefix_sum[node.confidence_range.1] - self.confidence_prefix_sum[node.confidence_range.0]
    }

    fn eval_recurse(&mut self, node: &TreeNode, confidence_prefix: &[f64]) {
        node.children.iter().for_each(|c| {
            let mut conf_prefix = confidence_prefix.to_vec();
            conf_prefix.push(self.get_confidence(c));
            self.eval_recurse(c, &conf_prefix);
            if c.level == self.tree.num_levels { self.confidence_vectors.push((c.confidence_range.0, conf_prefix)) };
        });
    }
}

pub struct Tree {
    root: TreeNode,
    lineages: Vec<String>,
    num_levels: usize,
    num_tips: usize,
}

impl Tree {
    fn new(mut lineages: Vec<String>) -> Self {
        let mut root = TreeNode::new(String::from("root"), 0, 0);
        lineages.sort();
        let num_levels = lineages[0].split(',').count();
        let mut confidence_idx = 0_usize;
        lineages.iter().for_each(|lineage| {
            let levels = lineage.split(',').collect_vec();
            let mut current_node = &mut root;
            for (level, &label) in levels.iter().enumerate() {
                match &current_node.get_last_child_label() {
                    Some(name) => {
                        if name.as_str() != label {
                            current_node.add_child(TreeNode::new(
                                    label.to_owned(),
                                    level + 1,
                                    confidence_idx,
                            ));
                        }
                        current_node.confidence_range.1 = confidence_idx + 1;
                    }
                    None => {
                        current_node.add_child(TreeNode::new(
                                label.to_owned(),
                                level + 1,
                                confidence_idx,
                        ));
                        current_node.confidence_range.1 = confidence_idx + 1;
                    }
                };
                if level + 1 == num_levels {
                    confidence_idx += 1;
                }
                current_node = current_node.children.last_mut().unwrap();
            }
        });
        root.confidence_range.1 = confidence_idx;
        Self {
            root,
            lineages,
            num_levels,
            num_tips: confidence_idx - 1,
        }
    }

    pub fn print(&self) {
        self.root.print();
    }
}

pub struct TreeNode {
    label: String,
    level: usize,
    confidence_range: (usize, usize),
    children: Vec<TreeNode>,
}

impl TreeNode {
    fn new(label: String, level: usize, confidence_idx: usize) -> Self {
        Self {
            label,
            level,
            confidence_range: (confidence_idx, confidence_idx + 1),
            children: vec![],
        }
    }

    fn add_child(&mut self, child: TreeNode) {
        self.children.push(child);
    }

    fn get_last_child_label(&self) -> Option<&String> {
        match &self.children.last() {
            Some(c) => Some(&c.label),
            None => None,
        }
    }

    fn print(&self) {
        println!(
            "{}{} {:?}",
            "  ".repeat(self.level),
            self.label,
            self.confidence_range
        );
        for child in &self.children {
            child.print();
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::lineage::{Lineage, Tree};

    #[test]
    fn test_tree_construction() {
        let lineages = vec![
            String::from("Animalia,Chordata,Mammalia,Primates,Hominidae,Homo"),
            "Animalia,Chordata,Mammalia,Primates,Hominidae,Pan".into(),
            "Animalia,Chordata,Mammalia,Carnivora,Canidae,Canis".into(),
            "Animalia,Chordata,Mammalia,Carnivora,Felidae,Felis".into(),
        ];
        let tree = Tree::new(lineages);
        let confidence_values = &[0.1, 0.3, 0.4, 0.2];
        tree.print();
        let lineage = Lineage::new(&tree, confidence_values);
        let result = lineage.evaluate();
        dbg!(result);
        assert_eq!(0, 1);
    }
}
