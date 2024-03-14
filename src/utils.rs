use itertools::Itertools;
use log::{debug, info};

use crate::parser::LookupTables;

pub fn sequence_to_kmers(sequence: &[u8]) -> Vec<u16> {
    let mut k_mer: u16 = 0;
    let mut k_mers = Vec::new();
    sequence[0..8]
        .iter()
        .enumerate()
        .for_each(|(j, c)| k_mer |= (*c as u16) << (14 - j * 2));
    k_mers.push(k_mer);
    sequence[8..].iter().for_each(|c| {
        k_mer = (k_mer << 2) | *c as u16;
        k_mers.push(k_mer);
    });
    k_mers
}

pub fn accumulate_and_output_results(lookup_tables: &LookupTables, hit_buffer: &[f64], cutoff: usize, query_label: &str) {
    let mut phylum_values = vec![0.0; lookup_tables.level_hierarchy_maps[0].len()];
    let mut class_values = vec![0.0; lookup_tables.level_hierarchy_maps[1].len()];
    let mut order_values = vec![0.0; lookup_tables.level_hierarchy_maps[2].len()];
    let mut family_values = vec![0.0; lookup_tables.level_hierarchy_maps[3].len()];
    let mut genus_values = vec![0.0; lookup_tables.level_hierarchy_maps[4].len()];
    // let mut species_counts = vec![
    //     vec![0.0; lookup_tables.level_hierarchy_maps[0].len()],
    //     vec![0.0; lookup_tables.level_hierarchy_maps[1].len()],
    //     vec![0.0; lookup_tables.level_hierarchy_maps[2].len()],
    //     vec![0.0; lookup_tables.level_hierarchy_maps[3].len()],
    //     vec![0.0; lookup_tables.level_hierarchy_maps[4].len()],
    // ];
    for (a, phylum) in lookup_tables.level_hierarchy_maps[0].iter().enumerate() {
        for class in phylum {
            for order in lookup_tables.level_hierarchy_maps[1][*class].iter() {
                for family in lookup_tables.level_hierarchy_maps[2][*order].iter() {
                    for genus in lookup_tables.level_hierarchy_maps[3][*family].iter() {
                        for species in lookup_tables.level_hierarchy_maps[4][*genus].iter() {
                            genus_values[*genus] += hit_buffer[*species];
                        }
                        family_values[*family] += genus_values[*genus];
                    }
                    order_values[*order] += family_values[*family];
                }
                class_values[*class] += order_values[*order];
            }
            phylum_values[a] += class_values[*class];
        }
    }

    // info!("{query_label}");
    // println!("{:?}", phylum_values);
    // println!("{:?}", class_values);
    // println!("{:?}", order_values);
    // println!("{:?}", family_values);
    // println!("{:?}", genus_values);
    // println!("{:?}", hit_buffer);
    // dbg!(&lookup_tables.level_hierarchy_maps);
    // println!("{:?}", lookup_tables.species_meta_vec);

    let mut out_count = 0_usize;
    for (a, _) in phylum_values.iter().enumerate().sorted_by(|a, b| b.1.partial_cmp(a.1).unwrap()) {
        // println!("Checking Phylum {} with index {}", lookup_tables.level_name_maps[0][a], a);
        for (_, b) in lookup_tables.level_hierarchy_maps[0][a].iter().enumerate().sorted_by(|a, b| class_values[*b.1].partial_cmp(&class_values[*a.1]).unwrap()) {
            // println!("Checking Class {} with index {}", lookup_tables.level_name_maps[1][*b], b);
            for (_, c) in lookup_tables.level_hierarchy_maps[1][*b].iter().enumerate().sorted_by(|a, b| order_values[*b.1].partial_cmp(&order_values[*a.1]).unwrap()) {
                // println!("Checking Order {} with index {}", lookup_tables.level_name_maps[2][*c], c);
                for (_, d) in lookup_tables.level_hierarchy_maps[2][*c].iter().enumerate().sorted_by(|a, b| family_values[*b.1].partial_cmp(&family_values[*a.1]).unwrap()) {
                    // println!("Checking Family {} with index {}", lookup_tables.level_name_maps[3][*d], d);
                    for (_, e) in lookup_tables.level_hierarchy_maps[3][*d].iter().enumerate().sorted_by(|a, b| genus_values[*b.1].partial_cmp(&genus_values[*a.1]).unwrap()) {
                        // println!("Checking Genus {} with index {}", lookup_tables.level_name_maps[4][*e], e);
                        for species in lookup_tables.level_hierarchy_maps[4][*e].iter().sorted_by(|a, b| hit_buffer[**b].partial_cmp(&hit_buffer[**a]).unwrap()) {
                            let species_meta = &lookup_tables.species_meta_vec[*species];
                            println!("{} {} {:.5}|{:.5}|{:.5}|{:.5}|{:.5}|{:.5}",
                                    query_label,
                                    species_meta.name,
                                    phylum_values[a],
                                    class_values[*b],
                                    order_values[*c],
                                    family_values[*d],
                                    genus_values[*e],
                                    hit_buffer[*species],
                                    );
                            out_count += 1;
                            if out_count == cutoff {
                                return;
                            }
                        }
                    }
                }
            }
        }
    }

}

#[cfg(test)]
mod tests {
    use crate::parser::parse_reference_fasta_str;

    use super::accumulate_and_output_results;

    #[test]
    fn test_accumulation() {

        let fasta_str = r">Badabing|Badabum|Phylum1|Class1|Order1|Family1|Genus1|Species1
AAACCCTTTGGGA
>Badabing|Badabum|Phylum1|Class1|Order1|Family1|Genus1|Species2
ATACGCTTTGGGA
>Badabing|Badabum|Phylum1|Class1|Order4|Family5|Genus2|Species3
ATCCGCTATGGGA
>Badabing|Badabum|Phylum1|Class2|Order2|Family3|Genus3|Species6
ATACGCTTTGCGT
>Badabing|Badabum|Phylum2|Class3|Order3|Family4|Genus4|Species5
ATACGCTTTGCGT";
        let lookup_table = parse_reference_fasta_str(fasta_str).unwrap();
        let hit_buffer = [1.0 / 8.0, 2.0 / 8.0, 0.0, 2.0 / 8.0, 3.0 / 8.0];
        accumulate_and_output_results(&lookup_table, &hit_buffer, 4, "SpeciesX");
        assert_eq!(0, 1);
    }
}
