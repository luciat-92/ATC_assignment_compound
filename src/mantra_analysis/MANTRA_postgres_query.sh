postgres > /dev/null 2>&1 &
psql -U postgres -d chembl_30 -c "COPY (
    select
        new_molecule_dictionary.molregno,
        new_molecule_dictionary.pref_name,
        new_molecule_dictionary.chembl_id,
        new_molecule_dictionary.synonyms,
        molecule_atc_classification.mol_atc_id,
        molecule_atc_classification.level5,
        atc_classification.level4,
        atc_classification.level4_description,
        atc_classification.level3,
        atc_classification.level3_description,
        atc_classification.level2,
        atc_classification.level2_description,
        atc_classification.level1,
        atc_classification.level1_description,
        atc_classification.who_name,
        profiles.name as mantra_name
    from (
        select
            molecule_dictionary.molregno,
            molecule_dictionary.pref_name,
            molecule_dictionary.chembl_id,
            molecule_dictionary.pref_name as synonyms
        from
            molecule_dictionary
        union
        select
            molecule_dictionary.molregno,
            molecule_dictionary.pref_name,
            molecule_dictionary.chembl_id,
            molecule_synonyms.synonyms
        from
            molecule_dictionary
            inner join
                molecule_synonyms 
                on 
                molecule_dictionary.molregno = molecule_synonyms.molregno
        ) AS new_molecule_dictionary
        left outer join
            molecule_atc_classification
            on
            new_molecule_dictionary.molregno = molecule_atc_classification.molregno
        left outer join
            atc_classification
            on
            molecule_atc_classification.level5 = atc_classification.level5
        left outer join
            mantra_schema.profiles as profiles
            on 
            upper(new_molecule_dictionary.synonyms) = upper(profiles.name)
    group by
        new_molecule_dictionary.molregno,
        new_molecule_dictionary.pref_name,
        new_molecule_dictionary.chembl_id,
        new_molecule_dictionary.synonyms,
        molecule_atc_classification.mol_atc_id,
        molecule_atc_classification.level5,
        atc_classification.level4,
        atc_classification.level4_description,
        atc_classification.level3,
        atc_classification.level3_description,
        atc_classification.level2,
        atc_classification.level2_description,
        atc_classification.level1,
        atc_classification.level1_description,
        atc_classification.who_name,
        profiles.name
    order by
        new_molecule_dictionary.molregno
)
    TO '/home/paolo.cremaschi/mantra/chembl_atc_mapping.tsv' 
    DELIMITER E'\t' 
    CSV 
    HEADER;"
pg_ctl stop
